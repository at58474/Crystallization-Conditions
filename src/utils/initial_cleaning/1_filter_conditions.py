"""
filter_conditions.py

1st cleaning script

Standalone utility to create a filtered SQLite database from Conditions.db.

Keeps only rows that:
  - Have Protein_ID
  - Have FASTA_Sequence
  - Contain at least one chemical (Organic_Precipitates OR Salt_Precipitates not NULL)

Writes result to: ../../data/processed/1FilteredConditions.db
"""

import sqlite3
from pathlib import Path

def filter_conditions(
    source_db="../../data/Conditions.db",
    target_db="../../data/processed/1FilteredConditions.db"
):
    # Ensure output directory exists
    Path(target_db).parent.mkdir(parents=True, exist_ok=True)

    # Connect to source and target
    conn_in = sqlite3.connect(source_db)
    cur_in = conn_in.cursor()

    conn_out = sqlite3.connect(target_db)
    cur_out = conn_out.cursor()

    # Recreate the same schema
    cur_out.execute("""
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
    cur_in.execute("""
        SELECT *
        FROM conditions
        WHERE Protein_ID IS NOT NULL
          AND FASTA_Sequence IS NOT NULL
          AND (
                Organic_Precipitates IS NOT NULL
             OR Salt_Precipitates IS NOT NULL
          )
    """)
    rows = cur_in.fetchall()

    # Insert into target DB
    cur_out.executemany("""
        INSERT OR REPLACE INTO conditions VALUES (?,?,?,?,?,?,?,?,?)
    """, rows)

    conn_out.commit()
    conn_in.close()
    conn_out.close()

    print(f"Filtered {len(rows)} rows into {target_db}")


if __name__ == "__main__":
    filter_conditions()
