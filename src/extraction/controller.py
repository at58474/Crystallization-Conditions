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
import sqlite3
import pandas as pd
import re
from tqdm import tqdm

from ccmodule import process_pdb_file
from convert_ent_to_pdb import ConvertToPDB
from extract_sequence import PDBPreprocessingSequences

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
pdb_file_folder = '../../data/raw/PDB_Files_Test'
fasta_file_folder = '../../data/raw/Protein_Sequences/SEQRES'
ent_file_folder = '../../data/raw/ENT_Files_Test'
db_path = "../../data/Conditions.db"

batch_size = 100          # PDB files per worker batch
insert_batch_size = 500   # rows per SQLite insert batch
export_csv = True         # toggle CSV snapshot
run_convert_ent = False   # toggle ENT→PDB conversion
run_extract_seq = False   # toggle FASTA extraction

# Optional tuning for FASTA extraction
extract_processes = 8     # number of processes to use
extract_chunksize = 50    # number of files per worker chunk


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def load_fasta_sequences(fasta_directory):
    """
    Load FASTA sequences into a dict keyed by Protein_ID.
    """
    fasta_map = {}
    fasta_files = Path(fasta_directory).glob("*.fasta")
    for file in fasta_files:
        with open(file, 'r') as f:
            lines = f.read().splitlines()
        protein_id = re.search(r'[A-Z0-9]{4}:[A-Z]', lines[0])
        if protein_id:
            protein_id = protein_id.group().split(':', 1)[0]
            fasta_seq = "".join(lines[1:])
            fasta_map[protein_id] = fasta_seq
    return fasta_map


def process_batch(file_list):
    """
    Process a batch of PDB files into row dicts.
    This is called inside multiprocessing workers.
    """
    return [process_pdb_file(f) for f in file_list]


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
    pdb_files = list(Path(pdb_file_folder).glob("*.pdb"))
    total_files = len(pdb_files)
    batches = [pdb_files[i:i + batch_size] for i in range(0, total_files, batch_size)]

    # Pre-load FASTA sequences
    fasta_map = load_fasta_sequences(fasta_file_folder)

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
