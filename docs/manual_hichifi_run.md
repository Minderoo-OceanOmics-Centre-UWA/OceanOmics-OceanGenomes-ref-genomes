# HiFi+HiC Assembly — Manual / Unattended Run Guide

This document covers the **manual entrypoint** for the HiFi+HiC assembly pipeline —
one config file + one script, self-contained inside this repo.

It is safe to run unattended: the script stages its own data, launches Nextflow
inside tmux (no need to start tmux yourself first), and on completion backs up
results to Acacia and pushes to the database.

---

## Clone the repo (do this first)

The pipeline repo must be cloned into a directory called **`ref-gen`** in your
scratch. The run script and config template both live inside the repo under
`configs/` — there is nothing to set up outside it.

```bash
cd /scratch/pawsey0964/$USER
git clone https://github.com/OceanOmics/OceanOmics-OceanGenomes-ref-genomes.git ref-gen/OceanOmics-OceanGenomes-ref-genomes
```

Your working tree will look like this:

```
/scratch/pawsey0964/$USER/
└── ref-gen/
    └── OceanOmics-OceanGenomes-ref-genomes/   ← this repo
        ├── configs/
        │   ├── hichifi_run.sh          ← the run script
        │   └── hichifi_pipeline.conf   ← the config template (edit this)
        ├── scripts/
        └── main.nf
```

> **Do not rename `ref-gen`.** The run script and config both derive output paths
> from this directory name via `/scratch/pawsey0964/$USER/ref-gen`.

---

## Quick start

```
1. Edit config   →  ref-gen/OceanOmics-OceanGenomes-ref-genomes/configs/hichifi_pipeline.conf
2. Run            →  bash ref-gen/OceanOmics-OceanGenomes-ref-genomes/configs/hichifi_run.sh
3. tmux self attached   →  tmux attach -t nf_<RUN_NAME>   (to check progress later)
```

Details for each step below.

---

## First-time setup (new user)

The BUSCO database used by **chromsyn** (post-assembly synteny step) must live in
your own scratch because `compleasm` tries to create lock/index files inside the
library directory at runtime. The shared reference copy at
`/scratch/references/busco_db/` is read-only and will cause a `PermissionError`.

Copy it once before your first run:

```bash
mkdir -p /scratch/pawsey0964/$USER/busco_db
cp -r /scratch/references/busco_db/actinopterygii_odb10 \
      /scratch/pawsey0964/$USER/busco_db/

# Prevent compleasm from trying to download placement files at runtime
# (the download code has a parsing bug; with explicit -l the files aren't needed)
mkdir -p /scratch/pawsey0964/$USER/busco_db/placement_files
touch /scratch/pawsey0964/$USER/busco_db/placement_files.done
```

This copies ~2 GB and takes a couple of minutes. You only need to do it once —
subsequent runs reuse the same directory.

The Nextflow pipeline (`--buscodb`) and `chromsyn.sh` (`LPATH`) both default to
`/scratch/pawsey0964/$USER/busco_db/actinopterygii_odb10`, so the same copy
satisfies both steps.

---

## Prerequisites

| Requirement | Location |
|---|---|
| Singularity containers | `/software/projects/pawsey0964/singularity/` |
| PostgreSQL credentials | `~/postgresql_details/oceanomics.cfg` |
| BaseSpace config | `~/.basespace/default.cfg` |
| rclone configured | `~/.config/rclone/rclone.conf` (remotes: `pawsey0964` (acacia), `s3`) |
| Asana token (optional) | `~/asana_token.config` — if absent, completion/failure notifications are silently skipped |
| Nextflow | loaded at login via `nextflow/24.10.0` module |
| BUSCO database | `/scratch/pawsey0964/$USER/busco_db/actinopterygii_odb10` — **must be in your own scratch** (see First-time setup above) |
| FCS-GX database | `/scratch/references/Foreign_Contamination_Screening/gxdb` (shared read-only, no copy needed) |
| tmux | required — the run script launches itself inside a session automatically |

---

## Step 1 — Edit the config

Everything is controlled by one file, which lives inside the repo:

```
ref-gen/OceanOmics-OceanGenomes-ref-genomes/configs/hichifi_pipeline.conf
```

The only things you need to update each run are `RUN_NAME` and `OG_IDS`:

```bash
RUN_NAME=hichifi_20260701_1
OG_IDS="OG2200,OG2205"
```

Other settings you may occasionally need to change:

| Key | Description |
|---|---|
| `STAGING_BASE_DIR` | Root dir for staged inputs + pipeline outputs (becomes Nextflow `--outdir`) |
| `SAMPLESHEET` | Path where the generated samplesheet is written |
| `HIC_BUCKET` / `HIFI_BUCKET` | rclone remotes for current-run Hi-C / HiFi reads |
| `CURRENT_HIC_RUN_IDS` | Hi-C run IDs already part of this batch (excluded from "stage previous HiC runs" so they aren't double-staged) — leave empty if this batch has no fresh Hi-C run of its own |
| `ASSEMBLY_MODE` / `SCAFFOLDER` | Nextflow pipeline mode (default `hifi_hic` / `yahs`) |
| `NOTIFY_EMAIL` | Optional — SLURM failure emails for post-assembly jobs. Leave empty to skip. |

The `{user}` placeholder in paths is replaced with `$USER`, and `{run}` is
replaced with `RUN_NAME` — both expanded automatically at runtime into a
per-run rendered config (`ref-gen/runs/<RUN_NAME>/pipeline.conf`). You don't
need to write either placeholder yourself unless you're customizing paths.

---

## Step 2 — Run the pipeline

```bash
bash ref-gen/OceanOmics-OceanGenomes-ref-genomes/configs/hichifi_run.sh
```

You do **not** need to start a tmux session first. If you're not already inside
tmux, the script detects this and relaunches itself inside a new session named
`nf_<RUN_NAME>`, then detaches — so it's safe to run from a plain login shell
and log out immediately after.

```
Launching inside new tmux session: nf_hichifi_20260701_1
Started. Attach with: tmux attach -t nf_hichifi_20260701_1
```

If a session with that name already exists (e.g. you ran it twice by mistake),
the script just prints the attach command and exits — it won't double-launch.

`hichifi_run.sh` then runs the full pipeline in order:

1. **Create samplesheet** — queries the OceanOmics PostgreSQL DB
2. **Stage previous HiC runs** — copies any prior Hi-C sequencing for these OGs from Acacia/S3 (no-op if none exist)
3. **Stage current HiFi data**
4. **Stage current HiC data**
5. **Run Nextflow** (`-resume`) — HiFiAdapterFilt → Meryl → GenomeScope2 → Hifiasm → BUSCO → Merqury → YAHS → fcs-gx → Tiara → PretextMap

On success it launches `post-hichifi-assembly.py`: chromsyn per sample, compile +
push QC results to the database, per-OG backup to Acacia, post-backup audit, and
an Asana completion notice.

On failure it launches `debug-nextflow.py`, which posts an Asana failure alert
and opens a Claude Code debug agent in a new tmux window to read the log,
diagnose the failure, and resume automatically if it's a transient error
(node failure, walltime, preemption).

### Monitoring

```bash
# Attach to watch live output
tmux attach -t nf_<RUN_NAME>

# SLURM jobs submitted by Nextflow
squeue -u $USER

# Detailed Nextflow log
tail -f ref-gen/runs/<RUN_NAME>/nextflow_*.log
```

### Extra Nextflow flags

Any arguments that aren't a `.conf` file are passed straight through to
`nextflow run`:

```bash
bash ref-gen/OceanOmics-OceanGenomes-ref-genomes/configs/hichifi_run.sh -resume -with-trace
```

---

## Outputs

- **Scratch:** `STAGING_BASE_DIR/<OG_ID>/` (Nextflow outputs)
- **Acacia:** `pawsey0964:oceanomics-refassemblies/<OG>/<OG>_<date>.<version>/`
- **Database:** on `master`, the `ref_genomes` table, upserted by the compile/push scripts. On the
  `phase-rg-v4` branch (Phase RG, not yet deployed) the scripts write the relational tables instead
  (`ref_assembly`, `ref_assembly_part` and one table per tool) through `scripts/lib/refdb.py`, one
  transaction per file, exiting non-zero on any error. Until the RG.6 swap those tables exist only in
  schema `v4`, so the branch must run with `PGOPTIONS='-c search_path=v4,public'`; with the default
  `search_path` it fails rather than writing anywhere. Tests: `python -m unittest discover tests`.

---

## Notifications

Asana notifications (completion task, failure alert) are sent only if
`~/asana_token.config` is present. If the file is absent both scripts skip
Asana silently — the pipeline still runs to completion and all other
post-assembly steps proceed normally.

SLURM failure emails for the post-assembly jobs (chromsyn, backup, audit) are
opt-in via `NOTIFY_EMAIL` in the config; leave it blank to skip them.

---

## Troubleshooting

### Nextflow failed

The Claude Code debug agent launches automatically in a `debug` tmux window —
attach with `tmux attach -t nf_<RUN_NAME>` and switch windows (`Ctrl+b` then
`w`, or `Ctrl+b` then `n`/`p`). If the agent determines the failure was
transient, it resumes the run for you. Otherwise it explains what needs manual
intervention.

To resume manually instead:

```bash
bash ref-gen/OceanOmics-OceanGenomes-ref-genomes/configs/hichifi_run.sh ref-gen/runs/<RUN_NAME>/pipeline.conf
```

(Note: pass the **rendered** per-run config here, not the original template —
it already has `{user}`/`{run}` resolved and is what the failed run actually used.)

### Staging previous HiC failed (S3 path not yet backed up)

This step warns and continues rather than failing the whole run — check the
printed warnings, but it's safe to ignore if the missing run genuinely hasn't
been backed up yet.

### Compile/push to DB failed but Nextflow succeeded

The post-assembly steps were already submitted as independent SLURM jobs
(`compile_<RUN_NAME>`, `backup_<OG>_<RUN_NAME>`, etc.) — check
`squeue -u $USER` and the corresponding `.out` log in `logs/`, fix the issue,
and re-run the failed step's script directly rather than re-running the whole
pipeline.

---

## Script reference

| Script | Purpose |
|---|---|
| `ref-gen/OceanOmics-OceanGenomes-ref-genomes/configs/hichifi_pipeline.conf` | Central config — edit `RUN_NAME`/`OG_IDS` here |
| `ref-gen/OceanOmics-OceanGenomes-ref-genomes/configs/hichifi_run.sh` | Full run: samplesheet → stage data → Nextflow → backup/DB/Asana |
| `scripts/0_create_samplesheet/create_samplesheet_hifi_hic.py` | Generate samplesheet from DB + config |
| `scripts/01_stage_data/stage_previous_hic.py` | Generate rclone commands for any prior HiC runs on these OGs |
| `scripts/01_stage_data/01_get_hic_from_config.sh` / `02_get_hifi_from_config.sh` | Stage current-run HiC / HiFi reads from Acacia |
| `scripts/06_post_assembly/post-hichifi-assembly.py` | Post-Nextflow: chromsyn, compile/push to DB, backup, Asana notify |
| `scripts/06_post_assembly/query_hichifi_results.py` | Query DB for QC results, used by the Asana completion notice |
| `scripts/06_post_assembly/debug-nextflow.py` | On failure: Asana alert + Claude Code debug agent |
