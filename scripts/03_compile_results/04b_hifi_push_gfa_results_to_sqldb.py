#!/usr/bin/env python3
# Push hifi-only (p_ctg / a_ctg) gfastats to the ref_genomes table at stage=0.
# singularity run $SING/psycopg2:0.1.sif python 04b_hifi_push_gfa_results_to_sqldb.py ~/postgresql_details/oceanomics.cfg

import psycopg2
import pandas as pd
import numpy as np
import configparser
import sys
from pathlib import Path

# Phase RG: database writes go through the shared helper (one transaction, non-zero exit on error).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import refdb  # noqa: E402


def load_db_config(config_file):
    if not Path(config_file).exists():
        raise FileNotFoundError(f"❌ Config file '{config_file}' does not exist.")
    config = configparser.ConfigParser()
    config.read(config_file)
    section = "postgres" if config.has_section("postgres") else "postgresql"
    if not config.has_section(section):
        raise ValueError("❌ Missing [postgres] section in config file.")
    params = {
        "dbname":   config.get(section, "dbname"),
        "user":     config.get(section, "user"),
        "password": config.get(section, "password"),
        "host":     config.get(section, "host"),
    }
    if config.has_option(section, "port"):
        params["port"] = config.getint(section, "port")
    return params


def parse_int(val):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None
    s = str(val).strip().replace(",", "")
    return int(float(s)) if s and s.lower() != "nan" else None


def parse_float(val):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None
    s = str(val).strip().replace(",", "").replace("%", "")
    return float(s) if s and s.lower() != "nan" else None


def parse_filename(fname):
    """
    Parse fields from a hifi gfastats filename. Handles multiple formats:
      OGxxx_vDDDDDD.hifi1.0.hifiasm.p_ctg.assembly_summary.txt       -> p_ctg
      OGxxx_vDDDDDD.hifi1.0.hifiasm.a_ctg.assembly.summary.txt       -> a_ctg
      OGxxx_vDDDDDD.hifi1.0.hifiasm.bp.hap1.p_ctg.assembly.summary.txt -> p_ctg
      OGxxx_vDDDDDD.hifi1.0.hifiasm.bp.hap2.p_ctg.assembly.summary.txt -> a_ctg
      OG855_240821.hifi1... (date without 'v' prefix)
      OG906_hifi1...       (no date — seq_date=None)
      OG90_bc2008_v240821.hifi1... (og_id=OG90_bc2008)
    """
    # Special case: OG906 — no date in filename
    if fname.startswith("OG906_hifi1"):
        haplotype = "p_ctg" if "p_ctg" in fname else ("a_ctg" if "a_ctg" in fname else None)
        return "OG906", "250212", "hifi1", 0, haplotype

    # Special case: OG90_bc2008 — barcode embedded in og_id
    if fname.startswith("OG90_bc2008"):
        haplotype = "p_ctg" if "p_ctg" in fname else ("a_ctg" if "a_ctg" in fname else None)
        return "OG90_bc2008", "240821", "hifi1", 0, haplotype

    parts = fname.split(".")
    first = parts[0].split("_")          # ["OGxxx", "vDDDDDD"] or ["OG855", "240821"]
    og_id    = first[0] if len(first) > 0 else None
    seq_date = first[1].lstrip("v") if len(first) > 1 else None
    version  = "hifi1"                   # stage 0 hifi assemblies are always hifi1
    stage    = parse_int(parts[2]) if len(parts) > 2 else None

    # Determine haplotype — handle bp.hap1/hap2 style
    if len(parts) > 4 and parts[4] == "bp":
        if len(parts) > 5 and parts[5] == "hap1":
            haplotype = "p_ctg"
        elif len(parts) > 5 and parts[5] == "hap2":
            haplotype = "a_ctg"
        else:
            haplotype = None  # bp.p_ctg — should have been filtered before reaching here
    else:
        haplotype = parts[4] if len(parts) > 4 else None

    return og_id, seq_date, version, stage, haplotype


config_file = sys.argv[1]
input_file  = sys.argv[2] if len(sys.argv) > 2 else "hifi_gfastats_report.txt"

print(f"Importing data from {input_file}")
gfa = pd.read_csv(input_file, sep="\t")

gfa[["og_id", "seq_date", "version", "stage", "haplotype"]] = gfa["filename"].apply(
    lambda x: pd.Series(parse_filename(str(x)))
)

print("\n🔍 Parsed rows:")
print(gfa[["filename", "og_id", "seq_date", "version", "stage", "haplotype"]].to_string(index=False))

# -------------------------------
# DB write (Phase RG): one row per assembly, part and tool in ref_assembly, ref_assembly_part and
# ref_gfastats, through scripts/lib/refdb.py. One transaction for the whole file; any error exits 1.
# Table names are unqualified: search_path decides (v4 during the RG window, public after RG.6).
# -------------------------------
PART_KEY = ("og_id", "seq_date", "version", "stage", "haplotype")


def build_rows():
    rows = []
    for i, row in gfa.iterrows():
        if row.get("og_id") == "OG90_bc2008":
            print("⏭️  Skipping OG90_bc2008: excluded barcode variant (Phase RG triage, RG.0 §3)")
            continue
        params = {
            "og_id":     row.get("og_id"),
            "seq_date":  row.get("seq_date"),
            "version":   row.get("version"),
            "stage":     parse_int(row.get("stage")),
            "haplotype": row.get("haplotype"),
            "num_contigs":                    parse_int(row.get("num_contigs")),
            "contig_n50":                     parse_int(row.get("contig_n50")),
            "contig_n50_size_mb":             parse_float(row.get("contig_n50_size_mb")),
            "num_scaffolds":                  parse_int(row.get("num_scaffolds")),
            "scaffold_n50":                   parse_int(row.get("scaffold_n50")),
            "scaffold_n50_size_mb":           parse_float(row.get("scaffold_n50_size_mb")),
            "largest_scaffold":               parse_int(row.get("largest_scaffold")),
            "largest_scaffold_size_mb":       parse_float(row.get("largest_scaffold_size_mb")),
            "total_scaffold_length":          parse_int(row.get("total_scaffold_length")),
            "total_scaffold_length_size_mb":  parse_float(row.get("total_scaffold_length_size_mb")),
            "gc_content_percent":             parse_float(row.get("gc_content_percent")),
        }
        rows.append(({k: params[k] for k in PART_KEY},
                     {k: v for k, v in params.items() if k not in PART_KEY}))
    return rows


refdb.main_guard(lambda: refdb.write_tool_rows(refdb.load_db_config(config_file), "ref_gfastats", build_rows()))
