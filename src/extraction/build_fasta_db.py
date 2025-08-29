"""
One-time script to build a SQLite database of FASTA sequences.

Instead of reading hundreds of thousands of tiny .fasta files every run,
this script stores them into a single database file (fasta_sequences.db).
Later, controller.py can query sequences directly from the database for speed.
"""

import re
import sqlite3
from pathlib import Path
from tqdm import tqdm
import sys

def build_fasta_db(fasta_directory, db_path="../../data/fasta_sequences.db"):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("""
        CREATE TABLE IF NOT EXISTS fasta_sequences (
            Protein_ID TEXT PRIMARY KEY,
            FASTA_Sequence TEXT
        )
    """)

    fasta_files = list(Path(fasta_directory).glob("*.fasta"))
    rows = []

    # tqdm progress bar while reading files
    for file in tqdm(fasta_files, desc="Building FASTA DB", unit="file"):
        with open(file, "r") as f:
            lines = f.read().splitlines()
        if not lines:
            continue

        # Extract protein ID from header (e.g., >1XYZ:A → "1XYZ")
        protein_id_match = re.search(r'[A-Z0-9]{4}:[A-Z]', lines[0])
        if protein_id_match:
            protein_id = protein_id_match.group().split(":")[0]
            fasta_seq = "".join(lines[1:])  # merge lines into one sequence
            rows.append((protein_id, fasta_seq))

    # tqdm progress bar while inserting
    for i in tqdm(range(0, len(rows), 1000), desc="Inserting into FASTA DB", unit="batch"):
        cur.executemany("INSERT OR REPLACE INTO fasta_sequences VALUES (?, ?)", rows[i:i+1000])
        conn.commit()

    conn.close()
    print(f"✅ Inserted {len(rows)} sequences into {db_path}")

if __name__ == "__main__":
    # Allow path override when run via subprocess
    if len(sys.argv) > 1:
        fasta_dir = sys.argv[1]
    else:
        fasta_dir = "../data/fasta_files"
    build_fasta_db(fasta_dir)
