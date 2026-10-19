#!/usr/bin/env python3
"""
Query ref_genomes for QC results of a HiFi+HiC assembly batch and write a
formatted report.

Usage:
    python query_hichifi_results.py <pg_config> <og_ids_csv> <output_file>

Example:
    python query_hichifi_results.py ~/postgresql_details/oceanomics.cfg \
        OG2104,OG2164 /scratch/pawsey0964/$USER/logs/results_hichifi_20260508.txt
"""

import configparser
import sys
from datetime import datetime

import psycopg2

pg_cfg   = sys.argv[1]
og_ids   = [x.strip() for x in sys.argv[2].split(",") if x.strip()]
out_file = sys.argv[3]

cfg = configparser.ConfigParser()
cfg.read(pg_cfg)
p   = cfg["postgres"]
conn = psycopg2.connect(
    dbname=p["dbname"], user=p["user"], password=p["password"],
    host=p["host"], port=p["port"],
)
cur = conn.cursor()

# Phase RG: one row per assembly and haplotype at its latest stage, from v_ref_assembly_summary
# (unqualified: search_path decides, v4 during the RG window, public after RG.6).
cur.execute("""
    SELECT
        og_id, seq_date, version, haplotype,
        num_scaffolds,
        ROUND(scaffold_n50_size_mb::numeric, 2)          AS scaffold_n50_mb,
        ROUND(total_scaffold_length_size_mb::numeric, 2) AS total_len_mb,
        num_chromosomes,
        ROUND(pct_assigned::numeric, 1)                  AS pct_chr,
        ROUND(busco_single_copy::numeric, 1)             AS busco_s,
        ROUND(qv::numeric, 2)                            AS qv,
        ROUND(kmer_completeness::numeric, 2)             AS kmer_compl
    FROM v_ref_assembly_summary
    WHERE og_id = ANY(%s)
    ORDER BY og_id, seq_date DESC, version, haplotype
""", (og_ids,))
rows = cur.fetchall()

# Total HiC yield per OG across all runs
cur.execute("""
    SELECT
        og_id,
        array_agg(DISTINCT run_id ORDER BY run_id) AS runs,
        COUNT(DISTINCT run_id)                      AS num_runs,
        ROUND(SUM(yield_gb)::numeric, 2)            AS total_yield_gb
    FROM hic_reads_qc
    WHERE og_id = ANY(%s)
    GROUP BY og_id
    ORDER BY og_id
""", (og_ids,))
yield_rows = cur.fetchall()
conn.close()

headers = [
    "OG", "seq_date", "version", "haplotype",
    "scaffolds", "N50(Mb)", "total(Mb)", "chr", "pct_chr",
    "BUSCO_s%", "QV", "k-compl%",
]

col_w = [max(len(h), 10) for h in headers]
for row in rows:
    for i, v in enumerate(row):
        col_w[i] = max(col_w[i], len(str(v) if v is not None else "N/A"))

def fmt_row(vals):
    return "  ".join(str(v if v is not None else "N/A").ljust(col_w[i]) for i, v in enumerate(vals))

lines = [
    f"HiFi+HiC Assembly QC Results",
    f"Generated : {datetime.now().strftime('%Y-%m-%d %H:%M AWST')}",
    f"OGs       : {', '.join(og_ids)}",
    "",
    fmt_row(headers),
    "-" * (sum(col_w) + 2 * len(col_w)),
]
for row in rows:
    lines.append(fmt_row(row))

lines += ["", f"Total rows: {len(rows)}"]

# HiC yield summary
lines += ["", "HiC Yield Summary (all runs combined):"]
yield_headers = ["OG", "num_runs", "runs", "total_yield_gb"]
yield_col_w   = [max(len(h), 10) for h in yield_headers]
for r in yield_rows:
    og, runs, n, total = r
    vals = [og, str(n), ", ".join(runs), str(total)]
    for i, v in enumerate(vals):
        yield_col_w[i] = max(yield_col_w[i], len(v))

def fmt_yield_row(vals):
    return "  ".join(str(v).ljust(yield_col_w[i]) for i, v in enumerate(vals))

lines.append(fmt_yield_row(yield_headers))
lines.append("-" * (sum(yield_col_w) + 2 * len(yield_col_w)))
if yield_rows:
    for og, runs, n, total in yield_rows:
        lines.append(fmt_yield_row([og, str(n), ", ".join(runs), str(total or "N/A")]))
else:
    lines.append("  (no yield data found in hic_reads_qc)")

report = "\n".join(lines)

with open(out_file, "w") as f:
    f.write(report + "\n")

print(report)
print(f"\nSaved to: {out_file}")

# Exit 1 if no key metrics populated — signals to caller that results are NOT in DB
has_results = any(
    row[9] is not None or row[10] is not None  # busco_s or qv
    for row in rows
)
sys.exit(0 if has_results else 1)
