#!/usr/bin/env python3
"""Record which Hi-C runs each Hi-C assembly in the samplesheet was scaffolded with.

    singularity run $SING/psycopg2:0.1.sif python 05b_push_hic_runs_to_sqldb.py \\
        ~/postgresql_details/oceanomics.cfg samplesheet.csv [BASE_DIR]

Run after 05a_push_omnic_results_to_sqldb.py (the assembly row must exist). For each samplesheet row
with a hicN version it finds the assembly's Hi-C BAM under BASE_DIR/<OG>/ (default
/scratch/pawsey0964/$USER/ref-gen, as 05_omnic-compile.sh), reads the flowcells from the read names
of its first records, and gives them to the database (refdb.write_hic_flowcells), which links the
OG's Hi-C sequencing rows on those runs (genomes_db_review migration 094).

Exits 1, after committing what it could, when a BAM is missing or a flowcell does not resolve: a run
not yet in sequencing_run.flowcell, or one with no Hi-C sequencing row for the OG (the lab's
workbook). v_ref_assembly_hic_readiness lists those assemblies; once the run or row is fixed they link
without a re-push.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import refdb  # noqa: E402

# Haplotype BAMs hold the same read pairs; take the first found in this order.
BAM_SUFFIXES = (".hap1.mapped.PT.bam", ".dual.mapped.PT.bam", ".hap2.mapped.PT.bam")


def find_bam(base, og_id, seq_date, version):
    stem = f"{og_id}_v{seq_date}.{version}"
    for suffix in BAM_SUFFIXES:
        hits = sorted((Path(base) / og_id).rglob(stem + suffix))
        if hits:
            return hits[0]
    return None


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: 05b_push_hic_runs_to_sqldb.py <postgres.cfg> <samplesheet.csv> [base_dir]")
    config_file, samplesheet = sys.argv[1], sys.argv[2]
    base = sys.argv[3] if len(sys.argv) > 3 else f"/scratch/pawsey0964/{os.environ.get('USER', '')}/ref-gen"

    versions = refdb.load_samplesheet(samplesheet, "version")
    dates = refdb.load_samplesheet(samplesheet, "date")
    items, missing = [], []
    for og_id, version in sorted(versions.items()):
        if not version.startswith("hic"):
            continue
        seq_date = dates.get(og_id, "").lstrip("v")
        if not refdb.SEQ_DATE_RE.match(seq_date):
            missing.append(f"{og_id} {version}: samplesheet date {dates.get(og_id)!r} is not vYYMMDD")
            continue
        bam = find_bam(base, og_id, seq_date, version)
        if bam is None:
            missing.append(f"{og_id} {seq_date}.{version}: no Hi-C BAM ({og_id}_v{seq_date}.{version}"
                           f"{{{','.join(BAM_SUFFIXES)}}}) under {Path(base) / og_id}")
            continue
        flowcells = refdb.bam_flowcells(bam)
        print(f"{og_id} {seq_date}.{version}: {bam.name}: "
              + ", ".join(f"{fc} ({n})" for fc, n in sorted(flowcells.items(), key=lambda x: -x[1])))
        items.append(((og_id, seq_date, version), flowcells))

    if items:
        try:
            refdb.write_hic_flowcells(refdb.load_db_config(config_file), items)
        except RuntimeError as e:
            missing.append(str(e))
    else:
        print("No Hi-C assemblies with a BAM in the samplesheet.")
    if missing:
        raise RuntimeError("\n".join(missing))


refdb.main_guard(main)
