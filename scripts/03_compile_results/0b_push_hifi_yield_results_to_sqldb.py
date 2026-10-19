#!/usr/bin/env python3
# singularity run $SING/psycopg2:0.1.sif python3 0b_push_hifi_yield_results_to_sqldb.py \
#     ~/postgresql_details/oceanomics.cfg /path/to/yield/*.csv
#
# Each (og_id, tissue, ext_type, lib_code, run_id) is stored as a separate row.
# Re-running the same CSV is idempotent (ON CONFLICT overwrites).
#
# CSV must include a RunName column (first column, added by 0c_backfill script or pipeline).
# Within a single run, multiple SMRT cells for the same sample are summed in Python.

import re
import sys
from pathlib import Path
import configparser
import pandas as pd
import numpy as np

# Phase RG: database writes go through the shared helper (one transaction, non-zero exit on error).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import refdb  # noqa: E402

# SampleName format: OG58W_D2_SBL  →  og_id=OG58, tissue=W, ext_type=D2, lib_code=SBL
SAMPLE_RE = re.compile(r'^(OG\d+)([A-Za-z]+)_([^_]+)_(.+)$')


def load_db_config(config_file: str) -> dict:
    p = Path(config_file)
    if not p.exists():
        raise FileNotFoundError(f"Config file '{config_file}' does not exist.")
    cfg = configparser.ConfigParser()
    cfg.read(config_file)
    if 'postgres' not in cfg:
        raise ValueError("Missing [postgres] section in config file.")
    required = ['dbname', 'user', 'password', 'host', 'port']
    missing = [k for k in required if not cfg.has_option('postgres', k)]
    if missing:
        raise ValueError(f"Missing keys in [postgres]: {missing}")
    return {
        'dbname':   cfg.get('postgres', 'dbname'),
        'user':     cfg.get('postgres', 'user'),
        'password': cfg.get('postgres', 'password'),
        'host':     cfg.get('postgres', 'host'),
        'port':     cfg.getint('postgres', 'port'),
    }


def parse_sample_name(sample_name: str) -> tuple | None:
    """Parse 'OG58W_D2_SBL' into (og_id, tissue, ext_type, lib_code). Returns None if unparseable."""
    m = SAMPLE_RE.match(str(sample_name).strip())
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3), m.group(4)


def simplify_barcode(barcode) -> str | None:
    """Convert 'bc2032--bc2032' to 'bc2032'. Pass through anything without '--'."""
    if pd.isna(barcode) or barcode is None:
        return None
    s = str(barcode)
    return s.split('--')[0] if '--' in s else s


def clean_yield(df: pd.DataFrame) -> pd.DataFrame:
    required_cols = {'RunName', 'SampleName', 'Barcode', 'BarcodeQuality', 'HiFiReads',
                     'HiFiReadLength', 'HiFiReadQuality', 'HiFiYield', 'PolymeraseReadLength'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Input file missing required columns: {sorted(missing)}")

    df = df.copy()

    parsed = df['SampleName'].apply(parse_sample_name)
    df = df[parsed.notna()].copy()
    if df.empty:
        return df

    df[['og_id', 'tissue', 'ext_type', 'lib_code']] = pd.DataFrame(
        parsed[parsed.notna()].tolist(), index=df.index
    )

    df['run_id'] = df['RunName']
    df['Barcode'] = df['Barcode'].apply(simplify_barcode)

    df = df.replace([np.inf, -np.inf], np.nan)

    for col in ['BarcodeQuality', 'HiFiReads', 'HiFiReadLength', 'HiFiYield', 'PolymeraseReadLength']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
        df[col] = df[col].where(pd.notna(df[col]), None)

    df['HiFiReadQuality'] = df['HiFiReadQuality'].astype(str).where(df['HiFiReadQuality'].notna(), None)

    return df[['og_id', 'tissue', 'ext_type', 'lib_code', 'run_id',
               'Barcode', 'BarcodeQuality', 'HiFiReads', 'HiFiReadLength',
               'HiFiReadQuality', 'HiFiYield', 'PolymeraseReadLength']]


def main():
    args = sys.argv[1:]

    if len(args) < 2:
        sys.exit("Usage: 0b_push_hifi_yield_results_to_sqldb.py "
                 "/path/to/oceanomics.cfg yield1.csv [yield2.csv ...]")

    config_file = args[0]
    csv_files   = args[1:]

    all_frames = []
    for csv_path in csv_files:
        print(f"Reading {csv_path}")
        df = pd.read_csv(csv_path)
        # If RunName column is absent, derive it from the filename (e.g. PACB_260603_AMD_yield.csv → PACB_260603_AMD)
        if 'RunName' not in df.columns:
            run_name = Path(csv_path).stem
            if run_name.endswith('_yield'):
                run_name = run_name[:-len('_yield')]
            print(f"  RunName column absent — using '{run_name}' from filename")
            df.insert(0, 'RunName', run_name)
        df = clean_yield(df)
        print(f"  -> {len(df)} barcoded rows after filtering")
        all_frames.append(df)

    if not all_frames:
        print("No data to upsert. Exiting.")
        return

    combined = pd.concat(all_frames, ignore_index=True)

    def parse_quality(q):
        if pd.isna(q) or q is None:
            return np.nan
        s = str(q).strip()
        try:
            return float(s.lstrip('Q'))
        except ValueError:
            return np.nan

    combined['_rq_num'] = combined['HiFiReadQuality'].apply(parse_quality)

    # Group by (og_id, tissue, ext_type, lib_code, run_id) — one row per run per sample.
    # Within a single run, multiple SMRT cells for the same sample are summed here.
    KEY = ['og_id', 'tissue', 'ext_type', 'lib_code', 'run_id']
    result_rows = []
    for key_vals, grp in combined.groupby(KEY, sort=True):
        og_id, tissue, ext_type, lib_code, run_id = key_vals
        reads = grp['HiFiReads'].fillna(0)
        total_reads = reads.sum()
        total_yield = int(grp['HiFiYield'].fillna(0).sum())

        def wmean(col):
            valid = grp[col].notna()
            if not valid.any() or reads[valid].sum() == 0:
                return None
            return float((grp.loc[valid, col] * reads[valid]).sum() / reads[valid].sum())

        best_idx = reads.idxmax()
        barcode = grp.loc[best_idx, 'Barcode']

        rq_num = wmean('_rq_num')
        rq_str = f"Q{round(rq_num)}" if rq_num is not None and not np.isnan(rq_num) else None

        bc_qual  = wmean('BarcodeQuality')
        rl       = wmean('HiFiReadLength')
        poly_rl  = wmean('PolymeraseReadLength')

        result_rows.append((
            og_id, tissue, ext_type, lib_code, run_id,
            barcode,
            round(bc_qual, 4) if bc_qual is not None else None,
            int(total_reads),
            round(rl, 2)      if rl      is not None else None,
            rq_str,
            total_yield,
            round(poly_rl, 2) if poly_rl is not None else None,
        ))

    n_multi = sum(1 for _, grp in combined.groupby(KEY) if len(grp) > 1)
    print(f"\n{len(result_rows)} unique run×sample rows ({n_multi} merged from multiple SMRT cells)")

    # Phase RG: hifi_read_qc, linked to the one sequencing row of this library on this run where
    # exactly one matches (scripts/lib/refdb.py). One transaction; any error exits 1.
    cols = refdb.HIFI_READ_QC_KEY + refdb.HIFI_READ_QC_COLS
    refdb.write_hifi_read_qc_rows(refdb.load_db_config(config_file),
                                  [dict(zip(cols, r)) for r in result_rows])


if __name__ == "__main__":
    refdb.main_guard(main)
