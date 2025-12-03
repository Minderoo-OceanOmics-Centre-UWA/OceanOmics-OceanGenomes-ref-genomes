#!/usr/bin/env python3
import psycopg2
import os
import pandas as pd
from datetime import date

# =====================================
# Database connection parameters
# =====================================
db_params = {
    'dbname': 'oceanomics_genomes',
    'user': 'postgres',
    'password': 'oceanomics',
    'host': '131.217.178.144',
    'port': 5432
}

# run with singularity run $SING/psycopg2:0.1.sif python
# =====================================
# OG IDs for the samplesheet
# =====================================
og_ids = [
    'OG910'
]

# =====================================
# SQL function definition
# =====================================
drop_function_sql = "DROP FUNCTION IF EXISTS build_nfcore_samplesheet_rows(text[]);"
create_function_sql = """
CREATE OR REPLACE FUNCTION build_nfcore_samplesheet_rows(in_og_ids text[])
RETURNS TABLE (
  sample   text,
  hifi_dir text,
  hic_dir  text,
  version  text,
  date     text,
  tolid    text,
  taxid    bigint,
  species  text
)
LANGUAGE sql
AS $$
WITH p AS (
  SELECT unnest(in_og_ids) AS og_id
),
latest_seq AS (
  SELECT DISTINCT ON (seq.og_id)
         seq.og_id,
         seq.seq_date::date AS seq_date
  FROM sequencing seq
  JOIN p ON seq.og_id = p.og_id
  WHERE seq.technology = 'PacBio'
  ORDER BY seq.og_id, seq.seq_date DESC
),
smp AS (
  SELECT DISTINCT ON (s.og_id)
         s.og_id,
         s.nominal_species_id
  FROM sample s
  JOIN p ON s.og_id = p.og_id
  ORDER BY s.og_id
)
SELECT DISTINCT ON (p.og_id)
  p.og_id AS sample,
  '/scratch/pawsey0964/$USER/ref-gen/'||p.og_id||'/hifi' AS hifi_dir,
  '/scratch/pawsey0964/$USER/ref-gen/'||p.og_id||'/hic'  AS hic_dir,
  CASE WHEN rg.og_id IS NOT NULL THEN 'hic2' ELSE 'hic1' END AS version,
  CASE WHEN ls.seq_date IS NOT NULL THEN 'v'||to_char(ls.seq_date,'YYMMDD') END AS date,
  p.og_id AS tolid,
  sp.ncbi_taxon_id AS taxid,
  sp.species
FROM p
LEFT JOIN ref_genomes rg ON rg.og_id = p.og_id
LEFT JOIN latest_seq ls  ON ls.og_id = p.og_id
LEFT JOIN smp ON smp.og_id = p.og_id
LEFT JOIN species sp ON sp.species = smp.nominal_species_id
ORDER BY p.og_id;
$$;
"""

# =====================================
# Connect to the database
# =====================================
conn = psycopg2.connect(**db_params)
cur = conn.cursor()

# Drop and recreate the SQL function so schema changes are picked up
cur.execute(drop_function_sql)
cur.execute(create_function_sql)
conn.commit()

# =====================================
# Call the function with OG list
# =====================================
query = """
SELECT *
FROM build_nfcore_samplesheet_rows(%s);
"""
df = pd.read_sql_query(query, conn, params=(og_ids,))

# Remove .0 from taxid column by converting to Int64 (nullable integer)
if 'taxid' in df.columns:
    df['taxid'] = df['taxid'].astype('Int64')  # keeps NaN for missing values

# Identify and print rows with missing values
missing_rows = df[df.isnull().any(axis=1)]
if not missing_rows.empty:
    print("\nRows with missing values:\n")
    print(missing_rows)

# Close DB connection
cur.close()
conn.close()

# =====================================
# Save CSV to Pawsey scratch
# =====================================
today_str = date.today().strftime("%Y%m%d")
current_dir = os.getcwd()
output_path = os.path.join(current_dir, f"samplesheet.csv")
df.to_csv(output_path, index=False)

print(f"Samplesheet saved to: {output_path}")
