#!/usr/bin/env python3
# singularity run $SING/psycopg2:0.1.sif python 01a_push_genomescope_comp_results_to_sqldb.py ~/postgresql_details/oceanomics.cfg [samplesheet.csv]

import csv
import urllib.request
import xml.etree.ElementTree as ET
import psycopg2
import pandas as pd
import numpy as np  # Required for handling infinity values
import configparser
import sys
from pathlib import Path

# Phase RG: database writes go through the shared helper (one transaction, non-zero exit on error).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import refdb  # noqa: E402

# -------------------------------
# Load DB credentials from .cfg
# -------------------------------
def load_db_config(config_file):
    config_file = str(Path(config_file).expanduser())
    if not Path(config_file).exists():
        raise FileNotFoundError(f"❌ Config file '{config_file}' does not exist.")
    config = configparser.ConfigParser()
    config.read(config_file)
    if not config.has_section('postgres'):
        raise ValueError("❌ Missing [postgres] section in config file.")
    required_keys = ['dbname', 'user', 'password', 'host', 'port']
    for key in required_keys:
        if not config.has_option('postgres', key):
            raise ValueError(f"❌ Missing '{key}' in [postgres] section of config file.")
    return {
        'dbname': config.get('postgres', 'dbname'),
        'user': config.get('postgres', 'user'),
        'password': config.get('postgres', 'password'),
        'host': config.get('postgres', 'host'),
        'port': config.getint('postgres', 'port')
    }

FISH_GENOME_DEFAULT = 1_000_000_000  # 1 GB fallback for fish with failed GenomeScope model

# ---- load taxid → og_id mapping from samplesheet ----
def load_taxids(samplesheet_path):
    taxids = {}
    path = Path(samplesheet_path).expanduser()
    if not path.exists():
        print(f"WARNING: samplesheet not found: {samplesheet_path} — fish genome default will not apply")
        return taxids
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            og = row.get('sample', '').strip().split('_')[0]
            taxid = row.get('taxid', '').strip()
            if og and taxid:
                taxids[og] = taxid
    return taxids

# ---- check if NCBI taxid is a fish (Actinopterygii or Chondrichthyes) ----
_fish_cache = {}

def is_fish(taxid):
    if not taxid:
        return False
    taxid = str(taxid).strip()
    if taxid in _fish_cache:
        return _fish_cache[taxid]
    try:
        url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
               f"?db=taxonomy&id={taxid}&retmode=xml")
        with urllib.request.urlopen(url, timeout=15) as r:
            lineage = ET.fromstring(r.read()).findtext('.//Lineage') or ''
        result = 'Actinopterygii' in lineage or 'Chondrichthyes' in lineage
    except Exception as e:
        print(f"  WARNING: NCBI taxonomy lookup failed for taxid {taxid}: {e}")
        result = False
    _fish_cache[taxid] = result
    return result

# ---- tiny helper to parse numbers the way your code intended ----
def parse_num(val, int_like=False):
    s = str(val).strip().replace(",", "").replace("%", "")
    if s == "" or s.lower() in ("nan", "na"):
        return None
    try:
        return int(s) if int_like else float(s)
    except ValueError:
        return None

# Config file path and optional samplesheet
config_file = sys.argv[1]
samplesheet = sys.argv[2] if len(sys.argv) > 2 else None

taxid_map = load_taxids(samplesheet) if samplesheet else {}
if taxid_map:
    print(f"Loaded taxids for {len(taxid_map)} OGs from samplesheet")
else:
    print("No samplesheet provided — fish genome size default will not be applied")

genomescope_compiled_path = "genomescope_compiled_results.tsv"

print(f"Importing data from {genomescope_compiled_path}")

genomescope = pd.read_csv(genomescope_compiled_path, sep="\t")

# Require 'sample' and create og_id exactly as you did
if 'sample' in genomescope.columns:
    genomescope['og_id'] = genomescope['sample'].astype(str).str.split('_').str[0]
else:
    raise ValueError("❌ 'sample' column not found in the input file.")

print("\n🔍 Final dataset summary:")
print(genomescope.describe())
# print(genomescope.columns)  # keep if you want

# -------------------------------
# DB write (Phase RG): ref_genomescope, one row per HiFi read set (og_id, seq_date), through
# scripts/lib/refdb.py. seq_date comes from the samplesheet's date column (vYYMMDD); without a
# samplesheet it is NULL, "read set unknown", as for every row from before Phase RG.
# One transaction for the whole file; any error exits 1.
# -------------------------------
date_map = {og.split("_")[0]: d.lstrip("v")
            for og, d in refdb.load_samplesheet(samplesheet, "date").items()} if samplesheet else {}


def build_rows():
    rows = []
    for _, row in genomescope.iterrows():
        og_id = row["og_id"]
        homozygosity = parse_num(row.get("homozygosity"), int_like=False)
        genomesize   = parse_num(row.get("genomesize"), int_like=True)

        # If the model failed (homo is NULL because it was 100%) and the sample
        # is a fish, use 1 GB as the genome size estimate for coverage calculations.
        if homozygosity is None and taxid_map:
            taxid = taxid_map.get(og_id)
            if is_fish(taxid):
                print(f"  {og_id}: model failed + fish (taxid {taxid}) — applying 1 GB genome size default")
                genomesize = FISH_GENOME_DEFAULT

        params = {
            "og_id": og_id,
            "homozygosity": homozygosity,
            "heterozygosity": parse_num(row.get("heterozygosity"), int_like=False),
            "genomesize": genomesize,
            "repeatsize": parse_num(row.get("repeatsize"), int_like=True),
            "uniquesize": parse_num(row.get("uniquesize"), int_like=True),
            "modelfit": parse_num(row.get("modelfit"), int_like=False),
            "errorrate": parse_num(row.get("errorrate"), int_like=False),
        }
        rows.append((og_id, date_map.get(og_id), {k: v for k, v in params.items() if k != "og_id"}))
    return rows


refdb.main_guard(lambda: refdb.write_read_set_rows(refdb.load_db_config(config_file), "ref_genomescope", build_rows()))
