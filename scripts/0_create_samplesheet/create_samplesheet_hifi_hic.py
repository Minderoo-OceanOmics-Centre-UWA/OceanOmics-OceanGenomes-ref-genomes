#!/usr/bin/env python3
"""
Create nf-core samplesheet rows by querying the OceanOmics PostgreSQL DB, using a single
bash-compatible pipeline config file (KEY=VALUE). Postgres credentials are read from an
INI file whose path is provided in the pipeline config as POSTGRES_CFG.

Intended for the ref-genomes pipeline (raw HiFi/Hi-C staging), producing columns:
  sample,hifi_dir,hic_dir,version,date,tolid,taxid,species,primary_assembly,hap1_assembly,hap2_assembly

Required config keys:
  POSTGRES_CFG=~/postgresql_details/oceanomics.cfg
  OG_IDS="OG38,OG39"
  STAGING_BASE_DIR=/scratch/pawsey0964/{user}/PIPELINE_DEV/REFGENOMES

Optional:
  SAMPLESHEET_OUTPUT_DIR=/path/to/assets
  SAMPLESHEET_FILENAME_PREFIX=samplesheet
  SAMPLESHEET_LATEST_NAME=samplesheet.csv   # writes/overwrites a stable copy in OUTPUT_DIR
  PRIMARY_ASSEMBLY_SUBDIR=primary_assembly   # subdirectory name under <OG_ID>/ (default: primary_assembly)
  HAP1_ASSEMBLY_SUBDIR=hap1                  # subdirectory name under <OG_ID>/ (default: hap1)
  HAP2_ASSEMBLY_SUBDIR=hap2                  # subdirectory name under <OG_ID>/ (default: hap2)

  Run:
  singularity run $SING/psycopg2:0.1.sif python create_samplesheet_from_config.py ../refgenomes_pipeline.conf
"""

from __future__ import annotations

import os
import sys
import shutil
import getpass
import configparser
from pathlib import Path
from datetime import date
from typing import Dict, List

import pandas as pd
import psycopg2


def load_kv_config(path: str) -> Dict[str, str]:
    """
    Load a simple KEY=VALUE config file (bash-compatible).
    - Ignores blank lines and lines starting with '#'
    - Strips surrounding quotes from values ("..." or '...')
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"❌ Config file does not exist: {path}")

    cfg: Dict[str, str] = {}
    for raw in p.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"❌ Invalid config line (expected KEY=VALUE): {raw}")
        k, v = line.split("=", 1)
        k = k.strip()
        v = v.strip()
        if (v.startswith('"') and v.endswith('"')) or (v.startswith("'") and v.endswith("'")):
            v = v[1:-1]
        cfg[k] = v
    return cfg


def require(cfg: Dict[str, str], key: str) -> str:
    if key not in cfg or cfg[key].strip() == "":
        raise ValueError(f"❌ Missing required config key: {key}")
    return cfg[key].strip()


def parse_og_ids(value: str) -> List[str]:
    # Accept "OG1,OG2" or "OG1 OG2"
    value = value.strip().strip(",")
    if not value:
        return []
    parts: List[str] = []
    for chunk in value.replace(",", " ").split():
        c = chunk.strip()
        if c:
            parts.append(c)
    return parts


def expand_user_placeholders(s: str, user: str) -> str:
    return s.replace("{user}", user)


def read_postgres_ini(postgres_cfg_path: str) -> Dict[str, str]:
    """
    Read Postgres connection details from an INI file with section [postgres].
    Required keys: dbname, user, password, host
    Optional: port
    """
    postgres_cfg_path = os.path.expanduser(postgres_cfg_path)

    if not os.path.exists(postgres_cfg_path):
        raise FileNotFoundError(f"❌ Postgres config not found: {postgres_cfg_path}")

    pg = configparser.ConfigParser()
    pg.read(postgres_cfg_path)

    if "postgres" not in pg:
        raise ValueError(f"❌ Missing [postgres] section in {postgres_cfg_path}")

    section = pg["postgres"]
    for k in ("dbname", "user", "password", "host"):
        if k not in section or section[k].strip() == "":
            raise ValueError(f"❌ Missing '{k}' in [postgres] section of {postgres_cfg_path}")

    return {
        "dbname": section["dbname"].strip(),
        "user": section["user"].strip(),
        "password": section["password"].strip(),
        "host": section["host"].strip(),
        "port": section.get("port", "5432").strip(),
    }


def build_function_sql(staging_base_dir: str, primary_subdir: str, hap1_subdir: str,
                       hap2_subdir: str, hifi_only: bool = False) -> str:
    """
    Inject staging_base_dir into the SQL function.
    When hifi_only=True, hic_dir is returned as an empty string so the pipeline
    validator does not complain about non-existent HiC directories.
    """
    base = staging_base_dir.rstrip("/")
    base_sql = base.replace("'", "''")  # SQL literal escape
    primary_sql = primary_subdir.replace("'", "''")
    hap1_sql = hap1_subdir.replace("'", "''")
    hap2_sql = hap2_subdir.replace("'", "''")

    hic_dir_expr = "''" if hifi_only else f"'{base_sql}/'||p.og_id||'/hic'"

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
-- One representative taxid per genus (for species-level fallback)
genus_tax AS (
  SELECT DISTINCT ON (genus)
         genus,
         ncbi_taxon_id
  FROM species
  WHERE ncbi_taxon_id IS NOT NULL
    AND genus IS NOT NULL
  ORDER BY genus, ncbi_taxon_id
)
SELECT DISTINCT ON (p.og_id)
  p.og_id AS sample,
  '{base_sql}/'||p.og_id||'/hifi' AS hifi_dir,
  {hic_dir_expr}                  AS hic_dir,
  CASE WHEN rg.og_id IS NOT NULL THEN 'hic2' ELSE 'hic1' END AS version,
  CASE WHEN ls.seq_date IS NOT NULL THEN 'v'||to_char(ls.seq_date,'YYMMDD') END AS date,
  COALESCE(smp.tol_id, p.og_id) AS tolid,
  COALESCE(sp.ncbi_taxon_id, gt.ncbi_taxon_id) AS taxid,
  smp.nominal_species_id AS species,
  '' AS primary_assembly,
  '' AS hap1_assembly,
  '' AS hap2_assembly
FROM p
LEFT JOIN ref_genomes rg ON rg.og_id = p.og_id AND rg.version LIKE 'hic%'
LEFT JOIN latest_seq ls  ON ls.og_id = p.og_id
LEFT JOIN smp ON smp.og_id = p.og_id
LEFT JOIN species sp ON sp.species = smp.nominal_species_id
LEFT JOIN genus_tax gt  ON gt.genus = split_part(smp.nominal_species_id, ' ', 1)
ORDER BY p.og_id;
$$;
""".strip()


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: create_samplesheet_from_config.py <pipeline_config.conf>", file=sys.stderr)
        sys.exit(1)

    conf_path = sys.argv[1]
    cfg = load_kv_config(conf_path)

    user = os.environ.get("USER") or getpass.getuser()
    user = user.replace("'", "").replace("/", "")

    og_ids = parse_og_ids(require(cfg, "OG_IDS"))
    if not og_ids:
        raise ValueError("❌ OG_IDS in config is empty")

    staging_base_dir = expand_user_placeholders(require(cfg, "STAGING_BASE_DIR"), user)

    out_dir = expand_user_placeholders(cfg.get("SAMPLESHEET_OUTPUT_DIR", "").strip(), user)
    prefix = cfg.get("SAMPLESHEET_FILENAME_PREFIX", "samplesheet").strip() or "samplesheet"
    latest_name = cfg.get("SAMPLESHEET_LATEST_NAME", "samplesheet.csv").strip() or ""

    primary_subdir = cfg.get("PRIMARY_ASSEMBLY_SUBDIR", "primary_assembly").strip() or "primary_assembly"
    hap1_subdir    = cfg.get("HAP1_ASSEMBLY_SUBDIR",    "hap1").strip()             or "hap1"
    hap2_subdir    = cfg.get("HAP2_ASSEMBLY_SUBDIR",    "hap2").strip()             or "hap2"
    assembly_mode  = cfg.get("ASSEMBLY_MODE", "").strip().lower()
    hifi_only      = assembly_mode == "hifi_only"
    if hifi_only:
        print("ℹ️  Assembly mode: hifi_only — hic_dir will be left empty")

    postgres_cfg = expand_user_placeholders(require(cfg, "POSTGRES_CFG"), user)
    pg = read_postgres_ini(postgres_cfg)

    func_sql = build_function_sql(staging_base_dir, primary_subdir, hap1_subdir, hap2_subdir,
                                  hifi_only=hifi_only)

    conn = None
    cur = None
    try:
        conn = psycopg2.connect(
            dbname=pg["dbname"],
            user=pg["user"],
            password=pg["password"],
            host=pg["host"],
            port=int(pg["port"]),
        )
        cur = conn.cursor()

        # Drop existing function first (can't change OUT columns with CREATE OR REPLACE)
        cur.execute("DROP FUNCTION IF EXISTS build_nfcore_samplesheet_rows(text[]);")
        conn.commit()

        cur.execute(func_sql)
        conn.commit()

        cur.execute("SELECT * FROM build_nfcore_samplesheet_rows(%s);", (og_ids,))
        rows = cur.fetchall()
        cols = [d[0] for d in cur.description]
        df = pd.DataFrame(rows, columns=cols)

        # taxid nullable integer
        if "taxid" in df.columns:
            df["taxid"] = pd.to_numeric(df["taxid"], errors="coerce").astype("Int64")

        # Report taxid source per sample
        cur.execute(
            "SELECT s.og_id, sp.ncbi_taxon_id "
            "FROM sample s "
            "LEFT JOIN species sp ON sp.species = s.nominal_species_id "
            "WHERE s.og_id = ANY(%s)",
            (og_ids,),
        )
        species_taxid = {r[0]: r[1] for r in cur.fetchall()}
        for _, row in df.iterrows():
            og  = row["sample"]
            tid = row.get("taxid")
            if pd.isna(tid):
                print(f"⚠️  {og}: taxid not found at species or genus level", file=sys.stderr)
            elif species_taxid.get(og) is None:
                print(f"ℹ️  {og}: taxid {tid} from genus-level fallback "
                      f"(species '{row.get('species')}' has no species-level taxid)")

        # Populate assembly columns with paths only where the directory exists on disk
        base = staging_base_dir.rstrip("/")
        for _, row in df.iterrows():
            og = row["sample"]
            for col, subdir in [
                ("primary_assembly", primary_subdir),
                ("hap1_assembly",    hap1_subdir),
                ("hap2_assembly",    hap2_subdir),
            ]:
                p = Path(f"{base}/{og}/{subdir}")
                if p.is_dir():
                    df.loc[df["sample"] == og, col] = str(p)

        missing_rows = df[df.isnull().any(axis=1)]
        if not missing_rows.empty:
            print("\nRows with missing values:\n", file=sys.stderr)
            print(missing_rows.to_string(index=False), file=sys.stderr)

        today = date.today().strftime("%Y%m%d")
        dated_filename = f"{prefix}_{today}.csv"

        if out_dir:
            out_path_dir = Path(out_dir)
            out_path_dir.mkdir(parents=True, exist_ok=True)
        else:
            out_path_dir = Path.cwd()

        dated_path = out_path_dir / dated_filename
        df.to_csv(dated_path, index=False)
        print(f"✅ Samplesheet saved to: {dated_path}")

        # Optional stable “latest” filename (e.g. samplesheet.csv) for your bash scripts
        if latest_name:
            latest_path = out_path_dir / latest_name
            shutil.copyfile(dated_path, latest_path)
            print(f"✅ Latest copy saved to: {latest_path}")

    finally:
        if cur is not None:
            cur.close()
        if conn is not None:
            conn.close()


if __name__ == "__main__":
    main()
