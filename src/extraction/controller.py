"""
controller.py

Main orchestration script for crystallization condition extraction.

Features:
- Optional ENT → PDB conversion
- Optional FASTA extraction
- Main extraction of crystallization conditions into SQLite DB
- Optional CSV export
- Multiprocessing
"""

from pathlib import Path
from multiprocessing import Pool
import sqlite3, os, subprocess
import pandas as pd
import re
from tqdm import tqdm

from ccmodule import process_pdb_file
from convert_ent_to_pdb import ConvertToPDB
from extract_sequence import PDBPreprocessingSequences

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

# Run type values are 'testing' and 'production', sets the correct filepath
run_type = 'production'

if run_type == 'testing':
    # ---------------------------------------------------------------------
    # Testing
    # ---------------------------------------------------------------------
    # Used by convert_ent_to_pdb.py
    # For converting ENT files to PDB files (preprocessing)
    ent_file_folder = '../../data/raw/ENT_Files_Test'
    # Used by extract_sequence.py and ccmodule.py
    # Where the PDB files are stored
    pdb_file_folder = '../../data/raw/PDB_Files_Test'
    # Used by extract_sequence.py and ccmodule.py
    # Where FASTA files are written
    fasta_file_folder = '../../data/raw/Protein_Sequences/SEQRES'

elif run_type == 'production':
    # ---------------------------------------------------------------------
    # Production
    # ---------------------------------------------------------------------
    ent_file_folder = 'D:/PDB_extracted'
    pdb_file_folder = 'D:/PDB_extracted'
    fasta_file_folder = 'D:/Data/Protein_Sequences/SEQRES'


db_path = "../../data/Conditions.db"


batch_size = 100          # PDB files per worker batch
insert_batch_size = 500   # rows per SQLite insert batch
export_csv = False         # toggle CSV snapshot
run_convert_ent = False   # toggle ENT→PDB conversion
run_extract_seq = False   # toggle FASTA extraction

# Optional tuning for FASTA extraction
extract_processes = 8     # number of processes to use
extract_chunksize = 50    # number of files per worker chunk


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
from tqdm import tqdm
import multiprocessing
import re
from pathlib import Path


def process_batch(file_list):
    """
    Process a batch of PDB files into row dicts.
    This is called inside multiprocessing workers.
    """
    return [process_pdb_file(f) for f in file_list]


def preload_fasta_sequences(fasta_file_folder, fasta_db_path="../../data/fasta_sequences.db"):
    """
    Ensure FASTA sequences are available in a SQLite database.
    - If the database doesn't exist, build it using build_fasta_db.py.
    - Returns a dict {Protein_ID: FASTA_Sequence}.
    """
    if not os.path.exists(fasta_db_path):
        print(f"⚠️ FASTA database not found at {fasta_db_path}. Building it now...")
        subprocess.run(
            ["python", "build_fasta_db.py", fasta_file_folder],
            check=True,
            cwd=os.path.dirname(__file__)
        )

    # Load FASTAs from the database with progress bar
    conn = sqlite3.connect(fasta_db_path)
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) FROM fasta_sequences")
    total = cur.fetchone()[0]

    cur.execute("SELECT Protein_ID, FASTA_Sequence FROM fasta_sequences")
    fasta_map = {}
    with tqdm(total=total, desc="Loading FASTA sequences", unit="seq") as pbar:
        for pid, seq in cur:
            fasta_map[pid] = seq
            pbar.update(1)

    conn.close()

    print(f"✅ Loaded {len(fasta_map)} FASTA sequences from {fasta_db_path}")
    return fasta_map


# ---------------------------------------------------------------------
# Main extraction pipeline
# ---------------------------------------------------------------------
def run_conditions_extraction():
    """
    Orchestrates optional preprocessing and main condition extraction.
    Writes results into SQLite DB (with optional CSV export).
    """
    # --- Optional ENT → PDB conversion ---
    if run_convert_ent:
        print("🔄 Running ENT → PDB conversion...")
        converter = ConvertToPDB(ent_file_folder)
        converter.run()

    # --- Optional FASTA extraction ---
    if run_extract_seq:
        print("🧬 Extracting FASTA sequences...")
        seq_extractor = PDBPreprocessingSequences(
            pdb_file_folder,
            fasta_file_folder,
            processes=extract_processes,
            chunksize=extract_chunksize,
        )
        seq_extractor.handle_multiprocessing()

    # --- Main extraction ---
    print("Getting pdb_files list...")
    pdb_files = list(Path(pdb_file_folder).glob("*.pdb"))
    print("Getting pdb_files list: COMPLETE")
    print("Calculating total number of files...")
    total_files = len(pdb_files)
    print(f"Calculating total number of files: COMPLETE with {total_files} files")
    print("Setting up batches...")
    batches = [pdb_files[i:i + batch_size] for i in range(0, total_files, batch_size)]
    print("Setting up batches: COMPLETE")

    # Pre-load FASTA sequences
    print("Preloading FASTA sequences...")
    fasta_map = preload_fasta_sequences(fasta_file_folder)
    print("Preloading FASTA sequences: COMPLETE")

    # Setup SQLite DB
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("""
    CREATE TABLE IF NOT EXISTS conditions (
        Protein_ID TEXT PRIMARY KEY,
        FASTA_Sequence TEXT,
        Solvent_Concentration TEXT,
        Matthews_Coefficient TEXT,
        Crystallography_Type TEXT,
        Vapor_Diffusion_Type TEXT,
        pH TEXT,
        Organic_Precipitates TEXT,
        Salt_Precipitates TEXT
    )
    """)

    insert_sql = """
    INSERT OR REPLACE INTO conditions
    (Protein_ID, FASTA_Sequence, Solvent_Concentration,
     Matthews_Coefficient, Crystallography_Type, Vapor_Diffusion_Type,
     pH, Organic_Precipitates, Salt_Precipitates)
    VALUES (:Protein_ID, :FASTA_Sequence, :Solvent_Concentration,
            :Matthews_Coefficient, :Crystallography_Type, :Vapor_Diffusion_Type,
            :pH, :Organic_Precipitates, :Salt_Precipitates)
    """

    # Process files in parallel batches
    print("Extraction starting...")
    with Pool() as pool, tqdm(total=total_files, desc="Processing PDB files") as pbar:
        for batch_result in pool.imap_unordered(process_batch, batches):
            # Merge FASTA into each row
            for row in batch_result:
                row["FASTA_Sequence"] = fasta_map.get(row["Protein_ID"], None)

            # Insert in chunks
            for i in range(0, len(batch_result), insert_batch_size):
                chunk = batch_result[i:i + insert_batch_size]
                cur.executemany(insert_sql, chunk)
            conn.commit()

            # Update progress bar
            pbar.update(len(batch_result))

    conn.close()
    print(f"✅ Extraction complete. Results written to {db_path}")

    if export_csv:
        export_results_to_csv()


def export_results_to_csv():
    """
    Export full SQLite conditions table to CSV.
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT * FROM conditions", conn)
    conn.close()

    out_path = "../../data/processed/CSV_Files/test_data.csv"
    df.to_csv(out_path, index=False)
    print(f"📄 Exported results to {out_path}")


if __name__ == "__main__":
    run_conditions_extraction()
