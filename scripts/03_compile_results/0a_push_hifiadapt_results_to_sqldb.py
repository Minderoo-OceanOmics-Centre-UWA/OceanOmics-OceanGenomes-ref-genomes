#!/usr/bin/env python3
# singularity run $SING/psycopg2:0.1.sif python 0a_push_hifiadapt_results_to_sqldb.py ~/postgresql_details/oceanomics.cfg [samplesheet.csv]
import sys
from pathlib import Path
import configparser
import pandas as pd
import numpy as np

# Phase RG: database writes go through the shared helper (one transaction, non-zero exit on error).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import refdb  # noqa: E402

def load_db_config(config_file: str) -> dict:
    p = Path(config_file)
    if not p.exists():
        raise FileNotFoundError(f"❌ Config file '{config_file}' does not exist.")
    cfg = configparser.ConfigParser()
    cfg.read(config_file)
    if 'postgres' not in cfg:
        raise ValueError("❌ Missing [postgres] section in config file.")
    required = ['dbname', 'user', 'password', 'host', 'port']
    missing = [k for k in required if not cfg.has_option('postgres', k)]
    if missing:
        raise ValueError(f"❌ Missing keys in [postgres]: {missing}")
    return {
        'dbname': cfg.get('postgres', 'dbname'),
        'user': cfg.get('postgres', 'user'),
        'password': cfg.get('postgres', 'password'),
        'host': cfg.get('postgres', 'host'),
        'port': cfg.getint('postgres', 'port')
    }

def clean_hifiadapt(df: pd.DataFrame) -> pd.DataFrame:
    if 'sample' in df.columns and 'og_id' not in df.columns:
        df = df.rename(columns={'sample': 'og_id'})
    required_cols = {'og_id', 'contam_reads'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"❌ Input file missing required columns: {sorted(missing)}")

    # Normalise types and nulls
    df = df.replace([np.inf, -np.inf], np.nan)
    df['og_id'] = df['og_id'].astype(str)
    df['contam_reads'] = pd.to_numeric(df['contam_reads'], errors='coerce')  # -> float or NaN
    # psycopg2 wants None for NULL
    df['contam_reads'] = df['contam_reads'].where(pd.notna(df['contam_reads']), None)
    return df[['og_id', 'contam_reads']]

def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: 0a_push_hifiadapt_results_to_sqldb.py /path/to/oceanomics.cfg [samplesheet.csv]")
    config_file = sys.argv[1]
    samplesheet = sys.argv[2] if len(sys.argv) > 2 else None
    tsv = "hifiadaptorfilt_stats_summary.tsv"

    print(f"Importing data from {tsv}")
    df = pd.read_csv(tsv, sep="\t")
    df = clean_hifiadapt(df)

    print("\n🔍 Dataset summary:")
    # include='all' so non-numeric columns are summarised too
    print(df.describe(include='all'))
    print(df.columns.tolist())

    # Phase RG: ref_hifi_adapter_filter, one row per HiFi read set (og_id, seq_date). seq_date comes
    # from the samplesheet's date column (vYYMMDD); without one it is NULL, "read set unknown".
    date_map = {og.split("_")[0]: d.lstrip("v")
                for og, d in refdb.load_samplesheet(samplesheet, "date").items()}
    rows = [(og_id, date_map.get(og_id), {"contam_reads": refdb.parse_int(contam_reads)})
            for og_id, contam_reads in df.itertuples(index=False, name=None)]
    if not rows:
        print("No rows to upsert. Exiting.")
        return
    refdb.write_read_set_rows(refdb.load_db_config(config_file), "ref_hifi_adapter_filter", rows)


if __name__ == "__main__":
    refdb.main_guard(main)
