"""
remove_invalid_fasta.py

2nd cleaning script

Removes rows where FASTA_Sequence is entirely X's (e.g. 'XXXX', 'XXXXXXXX').
Input:  ../../data/processed/FilteredConditions.db
Output: ../../data/processed/2RemovedXs.db
"""

import sqlite3
import os

def remove_invalid_fasta(
    input_db="../../data/processed/FilteredConditions.db",
    output_db="../../data/processed/2RemovedXs.db"
):
    if not os.path.exists(input_db):
        raise FileNotFoundError(f"Input database not found: {input_db}")

    conn_in = sqlite3.connect(input_db)
    cur_in = conn_in.cursor()

    conn_out = sqlite3.connect(output_db)
    cur_out = conn_out.cursor()

    # Get schema from input database
    cur_in.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name='conditions'")
    create_sql = cur_in.fetchone()[0]
    cur_out.execute(create_sql)

    # Copy rows where FASTA_Sequence is valid (not only X's)
    cur_in.execute("""
        SELECT * FROM conditions
        WHERE NOT (FASTA_Sequence GLOB 'X*' AND FASTA_Sequence LIKE 'XX%')
    """)
    rows = cur_in.fetchall()

    placeholders = ",".join("?" * len(rows[0])) if rows else ""
    if rows:
        cur_out.executemany(f"INSERT INTO conditions VALUES ({placeholders})", rows)

    conn_out.commit()
    conn_in.close()
    conn_out.close()
    print(f"Removed invalid FASTA rows. Cleaned DB saved to {output_db}")


if __name__ == "__main__":
    remove_invalid_fasta()
