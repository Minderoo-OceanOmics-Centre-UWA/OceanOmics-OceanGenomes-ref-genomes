#!/usr/bin/env python3
"""
Create nf-core samplesheet for HiFi-ONLY assemblies (no Hi-C data).

hic_dir is always left empty. version is hifi1 (first assembly) or hifi2
(OG already has a ref_genomes entry). Date is derived from the sequencing
table; if absent, falls back to the run-date embedded in SAMPLESHEET_FILENAME_PREFIX.

Required config keys:
  POSTGRES_CFG=~/postgresql_details/oceanomics.cfg
  OG_IDS="OG38,OG39"
  STAGING_BASE_DIR=/scratch/pawsey0964/$USER/ref-gen

Optional:
  SAMPLESHEET_OUTPUT_DIR=/path/to/assets
  SAMPLESHEET_FILENAME_PREFIX=samplesheet_PACB_260408_AMD
  SAMPLESHEET_LATEST_NAME=samplesheet.csv

Run:
  singularity run $SING/psycopg2:0.1.sif python create_samplesheet_hifi_only.py pipeline.conf
"""

from __future__ import annotations

import configparser
import getpass
import os
import re
import shutil
import sys
from datetime import date
from pathlib import Path
from typing import Dict, List

import pandas as pd
import psycopg2


def load_kv_config(path: str) -> Dict[str, str]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Config file does not exist: {path}")
    cfg: Dict[str, str] = {}
    for raw in p.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"Invalid config line (expected KEY=VALUE): {raw}")
        k, v = line.split("=", 1)
        k, v = k.strip(), v.strip()
        if (v.startswith('"') and v.endswith('"')) or (v.startswith("'") and v.endswith("'")):
            v = v[1:-1]
        cfg[k] = v
    return cfg


def require(cfg: Dict[str, str], key: str) -> str:
    if key not in cfg or not cfg[key].strip():
        raise ValueError(f"Missing required config key: {key}")
    return cfg[key].strip()


def parse_og_ids(value: str) -> List[str]:
    value = value.strip().strip(",")
    return [c.strip() for c in value.replace(",", " ").split() if c.strip()]


def expand_user_placeholders(s: str, user: str) -> str:
    return s.replace("{user}", user)


def read_postgres_ini(path: str) -> Dict[str, str]:
    path = os.path.expanduser(path)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Postgres config not found: {path}")
    pg = configparser.ConfigParser()
    pg.read(path)
    if "postgres" not in pg:
        raise ValueError(f"Missing [postgres] section in {path}")
    section = pg["postgres"]
    for k in ("dbname", "user", "password", "host"):
        if k not in section or not section[k].strip():
            raise ValueError(f"Missing '{k}' in [postgres] section")
    return {
        "dbname":   section["dbname"].strip(),
        "user":     section["user"].strip(),
        "password": section["password"].strip(),
        "host":     section["host"].strip(),
        "port":     section.get("port", "5432").strip(),
    }


def build_sql(staging_base_dir: str, primary_subdir: str, hap1_subdir: str,
              hap2_subdir: str) -> str:
    base = staging_base_dir.rstrip("/").replace("'", "''")
    primary_sql = primary_subdir.replace("'", "''")
    hap1_sql    = hap1_subdir.replace("'", "''")
    hap2_sql    = hap2_subdir.replace("'", "''")
    return f"""
CREATE OR REPLACE FUNCTION build_nfcore_samplesheet_rows(in_og_ids text[])
RETURNS TABLE (
  sample            text,
  hifi_dir          text,
  hic_dir           text,
  version           text,
  date              text,
  tolid             text,
  taxid             bigint,
  species           text,
  primary_assembly  text,
  hap1_assembly     text,
  hap2_assembly     text
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
  FROM v_sequencing_browse seq  -- sequencing has no og_id since the v2 cutover (2026-10-13)
  JOIN p ON seq.og_id = p.og_id
  WHERE seq.technology IN ('PacBio HIFI', 'PacBio')
    AND seq.status IS DISTINCT FROM 'Failed'  -- a failed run is not the read set's date (migration 069)
    AND seq.seq_date ~ '^[0-9]{{6}}$'  -- a malformed date would stop the ::date cast
  ORDER BY seq.og_id, seq.seq_date DESC
),
smp AS (
  SELECT DISTINCT ON (s.og_id)
         s.og_id,
         s.nominal_species_id,
         s.tol_id
  FROM sample s
  JOIN p ON s.og_id = p.og_id
  ORDER BY s.og_id
),
genus_tax AS (
  SELECT DISTINCT ON (genus) genus, ncbi_taxon_id
  FROM species
  WHERE ncbi_taxon_id IS NOT NULL AND genus IS NOT NULL
  ORDER BY genus, ncbi_taxon_id
)
SELECT DISTINCT ON (p.og_id)
  p.og_id AS sample,
  '{base}/'||p.og_id||'/hifi' AS hifi_dir,
  ''                           AS hic_dir,
  CASE WHEN rg.og_id IS NOT NULL THEN 'hifi2' ELSE 'hifi1' END AS version,
  CASE WHEN ls.seq_date IS NOT NULL THEN 'v'||to_char(ls.seq_date,'YYMMDD') END AS date,
  COALESCE(smp.tol_id, p.og_id) AS tolid,
  COALESCE(sp.ncbi_taxon_id, gt.ncbi_taxon_id) AS taxid,
  smp.nominal_species_id AS species,
  '' AS primary_assembly,
  '' AS hap1_assembly,
  '' AS hap2_assembly
FROM p
LEFT JOIN ref_genomes rg ON rg.og_id = p.og_id
LEFT JOIN latest_seq ls  ON ls.og_id = p.og_id
LEFT JOIN smp ON smp.og_id = p.og_id
LEFT JOIN species sp ON sp.species = smp.nominal_species_id
LEFT JOIN genus_tax gt ON gt.genus = split_part(smp.nominal_species_id, ' ', 1)
ORDER BY p.og_id;
$$;
""".strip()


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: create_samplesheet_hifi_only.py <pipeline_config.conf>", file=sys.stderr)
        sys.exit(1)

    cfg = load_kv_config(sys.argv[1])
    user = (os.environ.get("USER") or getpass.getuser()).replace("'", "").replace("/", "")

    og_ids = parse_og_ids(require(cfg, "OG_IDS"))
    if not og_ids:
        raise ValueError("OG_IDS in config is empty")

    staging_base = expand_user_placeholders(require(cfg, "STAGING_BASE_DIR"), user)
    out_dir      = expand_user_placeholders(cfg.get("SAMPLESHEET_OUTPUT_DIR", "").strip(), user)
    prefix       = cfg.get("SAMPLESHEET_FILENAME_PREFIX", "samplesheet").strip() or "samplesheet"
    latest_name  = cfg.get("SAMPLESHEET_LATEST_NAME", "").strip()

    primary_subdir = cfg.get("PRIMARY_ASSEMBLY_SUBDIR", "primary_assembly").strip() or "primary_assembly"
    hap1_subdir    = cfg.get("HAP1_ASSEMBLY_SUBDIR",    "hap1").strip()             or "hap1"
    hap2_subdir    = cfg.get("HAP2_ASSEMBLY_SUBDIR",    "hap2").strip()             or "hap2"

    pg = read_postgres_ini(expand_user_placeholders(require(cfg, "POSTGRES_CFG"), user))
    sql = build_sql(staging_base, primary_subdir, hap1_subdir, hap2_subdir)

    conn = cur = None
    try:
        conn = psycopg2.connect(dbname=pg["dbname"], user=pg["user"], password=pg["password"],
                                host=pg["host"], port=int(pg["port"]))
        cur = conn.cursor()
        cur.execute("DROP FUNCTION IF EXISTS build_nfcore_samplesheet_rows(text[]);")
        conn.commit()
        cur.execute(sql)
        conn.commit()
        cur.execute("SELECT * FROM build_nfcore_samplesheet_rows(%s);", (og_ids,))
        rows = cur.fetchall()
        cols = [d[0] for d in cur.description]
        df   = pd.DataFrame(rows, columns=cols)

        if "taxid" in df.columns:
            df["taxid"] = pd.to_numeric(df["taxid"], errors="coerce").astype("Int64")

        # Warn about missing taxids
        for _, row in df.iterrows():
            if pd.isna(row.get("taxid")):
                print(f"WARNING: {row['sample']}: taxid not found at species or genus level",
                      file=sys.stderr)

        # Fall back: derive date from prefix (e.g. samplesheet_PACB_260408_AMD -> v260408)
        m = re.search(r"_(\d{6})_", prefix)
        fallback_date = f"v{m.group(1)}" if m else None
        if fallback_date and "date" in df.columns:
            mask = df["date"].isnull() | (df["date"] == "")
            if mask.any():
                df.loc[mask, "date"] = fallback_date
                print(f"INFO: used run-name fallback date {fallback_date} for "
                      f"{mask.sum()} OG(s) with no sequencing record")

        # Populate assembly columns where directory exists
        base = staging_base.rstrip("/")
        for _, row in df.iterrows():
            og = row["sample"]
            for col, subdir in [("primary_assembly", primary_subdir),
                                 ("hap1_assembly", hap1_subdir),
                                 ("hap2_assembly", hap2_subdir)]:
                p = Path(f"{base}/{og}/{subdir}")
                if p.is_dir():
                    df.loc[df["sample"] == og, col] = str(p)

        missing = df[df.isnull().any(axis=1)]
        if not missing.empty:
            print("\nRows with missing values:", file=sys.stderr)
            print(missing.to_string(index=False), file=sys.stderr)

        out_path_dir = Path(out_dir) if out_dir else Path.cwd()
        out_path_dir.mkdir(parents=True, exist_ok=True)

        dated_path = out_path_dir / f"{prefix}_{date.today().strftime('%Y%m%d')}.csv"
        df.to_csv(dated_path, index=False)
        print(f"Samplesheet (hifi_only) saved to: {dated_path}")

        if latest_name:
            latest_path = out_path_dir / latest_name
            shutil.copyfile(dated_path, latest_path)
            print(f"Latest copy saved to: {latest_path}")

    finally:
        if cur:
            cur.close()
        if conn:
            conn.close()


if __name__ == "__main__":
    main()
