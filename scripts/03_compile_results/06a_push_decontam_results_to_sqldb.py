#!/usr/bin/env python3
# singularity run $SING/psycopg2:0.1.sif python 06a_push_decontam_results_to_sqldb.py ~/postgresql_details/oceanomics.cfg [final_decontam_report.tsv]
"""Push final_decontam_report.tsv (06_decontam_compile.py) into ref_decontamination, one row per
stage-2 haplotype (genomes_db_review migration 083), through scripts/lib/refdb.py: one transaction
for the whole file, and any error exits 1. Re-pushing the same file changes nothing."""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import refdb  # noqa: E402

STAGE = 2  # decontaminated (ref_stage)
KEY = ("og_id", "seq_date", "version", "haplotype")


def build_rows(path):
    rows = []
    with open(path, newline="") as f:
        for rec in csv.DictReader(f, delimiter="\t"):
            key = {k: rec[k] for k in KEY}
            key["stage"] = STAGE
            values = {c: refdb.parse_int(rec.get(c)) for c in refdb.TOOL_TABLES["ref_decontamination"]}
            rows.append((key, values))
    return rows


def main():
    config_file = sys.argv[1]
    report = sys.argv[2] if len(sys.argv) > 2 else "final_decontam_report.tsv"
    rows = build_rows(report)
    print(f"Pushing {len(rows)} haplotype row(s) from {report}")
    refdb.write_tool_rows(refdb.load_db_config(config_file), "ref_decontamination", rows)


if __name__ == "__main__":
    refdb.main_guard(main)
