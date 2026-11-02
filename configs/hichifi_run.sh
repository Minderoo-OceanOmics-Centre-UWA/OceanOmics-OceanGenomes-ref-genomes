#!/usr/bin/env bash
# hichifi_run.sh — run the HiFi+HiC reference-genome assembly pipeline end to end.
#
# Steps:
#   [1/5] Create samplesheet from the DB (create_samplesheet_hifi_hic.py)
#   [2/5] Stage any previous HiC runs for these OGs from Acacia (stage_previous_hic.py)
#   [3/5] Stage current HiFi data from Acacia
#   [4/5] Stage current HiC data from Acacia
#   [5/5] Run Nextflow (-resume)
#   On success: post-hichifi-assembly.py — chromsyn, compile+push to DB,
#               per-OG backup to Acacia, audit, Asana completion notice
#   On failure: debug-nextflow.py — Asana failure alert + auto Claude Code debug agent
#
# Usage: bash hichifi_run.sh [config.conf] [extra nextflow args...]
#
# Edit RUN_NAME and OG_IDS in the config (default: hichifi_pipeline.conf next to
# this script) before running. No need to be inside tmux first — if you're not
# already in a tmux session, this script re-launches itself inside one
# (nf_<RUN_NAME>) so the run survives you logging out. Attach with:
#   tmux attach -t nf_<RUN_NAME>
#
# Self-contained: only needs this script + a clone of the
# OceanOmics-OceanGenomes-ref-genomes pipeline repo. This script lives inside
# the repo at configs/hichifi_run.sh; the repo must be cloned to
# /scratch/pawsey0964/$USER/ref-gen/OceanOmics-OceanGenomes-ref-genomes
# — no access to automate_hic/ or ref_gen_automation/ required.

set -uo pipefail
source /opt/cray/pe/lmod/lmod/init/bash
set +u
[[ -f ~/.bashrc ]] && source ~/.bashrc
set -eu -o pipefail

SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SELF")"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SCRATCH_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
POST_ASSEMBLY_DIR="$REPO_DIR/scripts/06_post_assembly"
SING_DIR="${SING:-/software/projects/pawsey0964/singularity}"

CONFIG="$SCRIPT_DIR/hichifi_pipeline.conf"
NF_EXTRA_ARGS=()
for arg in "$@"; do
  case "$arg" in
    *.conf) CONFIG="$arg" ;;
    *)      NF_EXTRA_ARGS+=("$arg") ;;
  esac
done

[[ -f "$CONFIG" ]] || { echo "Config not found: $CONFIG" >&2; exit 1; }

# Raw pass — just enough to validate RUN_NAME/OG_IDS and compute the run dir
source <(grep -E '^[A-Z_][A-Z0-9_]*=' "$CONFIG")
: "${RUN_NAME:?RUN_NAME not set in $CONFIG}"
: "${OG_IDS:?OG_IDS not set in $CONFIG}"
if [[ "$RUN_NAME" == "hichifi_YYYYMMDD_N" ]]; then
  echo "Edit RUN_NAME in $CONFIG to a unique batch name before running" >&2
  exit 1
fi

USER_REAL="${USER:-$(whoami)}"
RUN_DIR="$SCRATCH_ROOT/ref-gen/runs/$RUN_NAME"
mkdir -p "$RUN_DIR"

# Render {user}/{run} placeholders into a fully-resolved per-run config, so
# every downstream script (Python and bash) gets concrete paths regardless of
# who runs this or which placeholders they support natively.
RENDERED_CONF="$RUN_DIR/pipeline.conf"
sed -e "s/{user}/$USER_REAL/g" -e "s/{run}/$RUN_NAME/g" "$CONFIG" > "$RENDERED_CONF.tmp" && mv "$RENDERED_CONF.tmp" "$RENDERED_CONF"

# Resolved pass — overwrite path vars with their expanded values
source <(grep -E '^[A-Z_][A-Z0-9_]*=' "$RENDERED_CONF")

# ── Self-wrap into tmux so the run survives the runner logging out ──────────
if [[ -z "${TMUX:-}" ]]; then
  SESSION="nf_${RUN_NAME}"
  if tmux has-session -t "$SESSION" 2>/dev/null; then
    echo "tmux session '$SESSION' already exists."
    echo "Attach with: tmux attach -t $SESSION"
    exit 0
  fi
  quoted_args=()
  for a in "$CONFIG" "${NF_EXTRA_ARGS[@]}"; do
    quoted_args+=("$(printf '%q' "$a")")
  done
  echo "Launching inside new tmux session: $SESSION"
  tmux new-session -d -s "$SESSION" "bash $(printf '%q' "$SELF") ${quoted_args[*]}"
  echo "Started. Attach with: tmux attach -t $SESSION"
  exit 0
fi

EMAIL_ARGS=()
[[ -n "${NOTIFY_EMAIL:-}" ]] && EMAIL_ARGS=(--email "$NOTIFY_EMAIL")

echo "=== HiFi+HiC assembly batch: $RUN_NAME ==="
echo "OG IDs      : $OG_IDS"
echo "Samplesheet : $SAMPLESHEET"
echo "Outdir      : $STAGING_BASE_DIR"
echo "Run dir     : $RUN_DIR"
echo ""

# ────────────────────────────────────────────────────────────────────────────
# [1/5] Create samplesheet
# ────────────────────────────────────────────────────────────────────────────
echo "=== [1/5] Creating samplesheet ==="
mkdir -p "$(dirname "$SAMPLESHEET")"
singularity run "$SING_DIR/psycopg2:0.1.sif" python \
  "$REPO_DIR/scripts/0_create_samplesheet/create_samplesheet_hifi_hic.py" "$RENDERED_CONF"

# ────────────────────────────────────────────────────────────────────────────
# [2/5] Stage previous HiC runs (idempotent — no-op if none exist)
# ────────────────────────────────────────────────────────────────────────────
echo "=== [2/5] Staging previous HiC runs ==="
PREV_HIC_CMDS="$RUN_DIR/prev_hic_cmds.sh"
singularity run "$SING_DIR/psycopg2:0.1.sif" python \
  "$REPO_DIR/scripts/01_stage_data/stage_previous_hic.py" \
  "$POSTGRES_CFG" "$SAMPLESHEET" \
  --current-runs "${CURRENT_HIC_RUN_IDS:-}" \
  --cmd-output "$PREV_HIC_CMDS"
bash "$PREV_HIC_CMDS"

# ────────────────────────────────────────────────────────────────────────────
# [3/5] Stage current HiFi data
# ────────────────────────────────────────────────────────────────────────────
echo "=== [3/5] Staging HiFi data ==="
bash "$REPO_DIR/scripts/01_stage_data/02_get_hifi_from_config.sh" "$RENDERED_CONF"

# ────────────────────────────────────────────────────────────────────────────
# [4/5] Stage current HiC data
# ────────────────────────────────────────────────────────────────────────────
echo "=== [4/5] Staging HiC data ==="
bash "$REPO_DIR/scripts/01_stage_data/01_get_hic_from_config.sh" "$RENDERED_CONF"

# ────────────────────────────────────────────────────────────────────────────
# [5/5] Run Nextflow
# ────────────────────────────────────────────────────────────────────────────
echo "=== [5/5] Running Nextflow ==="
cd "$RUN_DIR"
NF_TS="$(date +%Y%m%d_%H%M%S)"
NF_LOG="$RUN_DIR/nextflow_${NF_TS}.log"
# Per-run usage reports (trace columns are set in nextflow.config). native_id in
# the trace is the SLURM job ID, for joining with sacct to calculate SU.
NF_INFO_DIR="$RUN_DIR/pipeline_info"
mkdir -p "$NF_INFO_DIR"

set +e
nextflow run "$REPO_DIR/main.nf" \
  -profile singularity \
  --input "$SAMPLESHEET" \
  --outdir "$STAGING_BASE_DIR" \
  --assembly_mode "$ASSEMBLY_MODE" \
  --scaffolder "$SCAFFOLDER" \
  --buscodb "$BUSCODB" \
  --gxdb "$GXDB" \
  --binddir "$BINDDIR" \
  -c "$REPO_DIR/pawsey_profile.config" \
  -resume \
  -with-trace "$NF_INFO_DIR/execution_trace_${NF_TS}.txt" \
  -with-report "$NF_INFO_DIR/execution_report_${NF_TS}.html" \
  -with-timeline "$NF_INFO_DIR/execution_timeline_${NF_TS}.html" \
  --tempdir "$TEMPDIR" \
  --bs_config "$BS_CONFIG" \
  --sql_config "$SQL_CONFIG" \
  "${NF_EXTRA_ARGS[@]}" \
  2>&1 | tee "$NF_LOG"
NF_EXIT=${PIPESTATUS[0]}
set -e

echo ""
if [[ $NF_EXIT -eq 0 ]]; then
  echo "Nextflow complete for $RUN_NAME — launching post-assembly steps"
  python3 "$POST_ASSEMBLY_DIR/post-hichifi-assembly.py" \
    --run "$RUN_NAME" \
    --og-ids "$OG_IDS" \
    --samplesheet "$SAMPLESHEET" \
    --staging-base "$STAGING_BASE_DIR" \
    --pipeline-dir "$REPO_DIR" \
    --log-dir "$SCRATCH_ROOT/logs" \
    "${EMAIL_ARGS[@]}"
else
  echo "Nextflow FAILED for $RUN_NAME (exit $NF_EXIT) — launching Claude Code debug agent"
  python3 "$POST_ASSEMBLY_DIR/debug-nextflow.py" \
    --run "$RUN_NAME" \
    --run-dir "$RUN_DIR" \
    --script "$SELF $RENDERED_CONF" \
    --og-ids "$OG_IDS" \
    --samplesheet "$SAMPLESHEET" \
    --log-dir "$SCRATCH_ROOT/logs" \
    --assembly-mode "$ASSEMBLY_MODE"
  exit "$NF_EXIT"
fi
