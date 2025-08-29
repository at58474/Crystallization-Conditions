"""
filter_conditions.py

1st cleaning script

Standalone utility to create a filtered SQLite database from Conditions.db.

Keeps only rows that:
  - Have Protein_ID
  - Have FASTA_Sequence
  - Contain at least one chemical (Organic_Precipitates OR Salt_Precipitates not NULL)

Writes result to: ../../data/processed/FilteredConditions.db
"""

import sqlite3
from pathlib import Path

def filter_conditions(
    source_db="../../data/Conditions.db",
    target_db="../../data/processed/FilteredConditions.db"
):
    # Ensure output directory exists
    Path(target_db).parent.mkdir(parents=True, exist_ok=True)

    # Connect to source and target
    src_conn = sqlite3.connect(source_db)
    tgt_conn = sqlite3.connect(target_db)

    src_cur = src_conn.cursor()
    tgt_cur = tgt_conn.cursor()

    # Recreate the same schema
    tgt_cur.execute("""
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

    # Select only the rows of interest
    src_cur.execute("""
        SELECT *
        FROM conditions
        WHERE Protein_ID IS NOT NULL
          AND FASTA_Sequence IS NOT NULL
          AND (
                Organic_Precipitates IS NOT NULL
             OR Salt_Precipitates IS NOT NULL
          )
    """)
    rows = src_cur.fetchall()

    # Insert into target DB
    tgt_cur.executemany("""
        INSERT OR REPLACE INTO conditions VALUES (?,?,?,?,?,?,?,?,?)
    """, rows)

    tgt_conn.commit()
    src_conn.close()
    tgt_conn.close()

    print(f"✅ Filtered {len(rows)} rows into {target_db}")


if __name__ == "__main__":
    filter_conditions()
