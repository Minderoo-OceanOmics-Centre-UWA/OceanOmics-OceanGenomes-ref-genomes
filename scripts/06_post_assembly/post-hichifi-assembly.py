#!/usr/bin/env python3
"""
Post-assembly orchestrator for HiFi+HiC (hifi_hic) runs.
Called from the Nextflow tmux launch script after 'nextflow run' completes successfully.

Submits SLURM jobs in order:
  1. Chromsyn (per sample, staggered 120 s apart)
  2. Compile + push to DB (final gfastats, BUSCO, Merqury QV,
     Merqury completeness, OMNIC) — runs immediately, no chromsyn dependency
  3. Per-OG backup (full_assembly_backup.sh, mode=hifi_hic) — waits for chromsyn
  4. Post-backup audit — after all backups
  5. Asana completion notification — after post-audit

This copy lives inside the pipeline repo itself (scripts/06_post_assembly/) so
it works from a plain clone of this repo alone — no sibling automate_hic/ or
ref_gen_automation/ checkout required.

Usage:
  python3 post-hichifi-assembly.py \\
      --run hichifi_20260508 \\
      --og-ids OG2104,OG2164,OG2188,OG2202 \\
      --samplesheet /path/to/samplesheet_hichifi_20260508.csv \\
      [--staging-base /scratch/pawsey0964/$USER/ref-gen] \\
      [--pipeline-dir /path/to/this/repo] \\
      [--log-dir /scratch/pawsey0964/$USER/logs] \\
      [--email you@example.com]

All path defaults are derived from $MYSCRATCH (Pawsey's per-user scratch env var)
and this script's own location, so it runs unmodified for any user — override
with the flags above if your checkout lives somewhere non-standard.
"""

import csv
import getpass
import json
import os
import subprocess
import sys
import urllib.request
from datetime import datetime
from pathlib import Path

# ── Paths (portable — derived from $MYSCRATCH / $USER / this script's own
# location in the repo, override via CLI flags) ─────────────────────────────
SCRIPT_DIR   = Path(__file__).resolve().parent
USER_REAL    = os.environ.get("USER") or getpass.getuser()
MYSCRATCH    = os.environ.get("MYSCRATCH") or f"/scratch/pawsey0964/{USER_REAL}"

STAGING_BASE  = os.environ.get("HICHIFI_STAGING_BASE", f"{MYSCRATCH}/ref-gen")
# This script lives at <repo>/scripts/06_post_assembly/ — default to the repo
# root two levels up rather than guessing a path under STAGING_BASE.
PIPELINE_DIR  = os.environ.get("HICHIFI_PIPELINE_DIR", str(SCRIPT_DIR.parent.parent))
SCRIPTS_DIR   = f"{PIPELINE_DIR}/scripts/03_compile_results"
BACKUP_SCRIPT = f"{PIPELINE_DIR}/scripts/04_backup_scripts/full_assembly_backup.sh"
AUDIT_SCRIPT  = f"{PIPELINE_DIR}/scripts/04_backup_scripts/audit_backup.sh"
CHROMSYN_SH   = f"{PIPELINE_DIR}/scripts/02_run_chromsyn/chromsyn.sh"
QUERY_SCRIPT  = str(SCRIPT_DIR / "query_hichifi_results.py")
SING          = os.environ.get("SING", "/software/projects/pawsey0964/singularity")
PSYCOPG2_SIF  = f"{SING}/psycopg2:0.1.sif"
POSTGRES_CFG  = "~/postgresql_details/oceanomics.cfg"
LOG_DIR       = os.environ.get("HICHIFI_LOG_DIR", f"{MYSCRATCH}/logs")
ACCOUNT       = "pawsey1348"
EMAIL         = os.environ.get("HICHIFI_NOTIFY_EMAIL", "")
UWA_WORKSPACE_GID = "1103015454494961"
ASANA_PROJECT_GID    = "1211147243530958"
ASANA_SECTION_RUN_QC = "1211216586258197"
CURATION_PARENT_GID  = "1212686915278515"
# Asana tasks always go to the same person regardless of who runs the
# pipeline — this is a notification routing choice, not a scratch path.
ASANA_ASSIGNEE_ME    = "1203755755031589"

# SBATCH mail directives — only included if EMAIL is set (skip otherwise, since
# whoever runs this script may not want failure emails routed to someone else)
MAIL_DIRECTIVES = f"#SBATCH --mail-type=FAIL\n#SBATCH --mail-user={EMAIL}" if EMAIL else ""


# ── Helpers ───────────────────────────────────────────────────────────────────

def get_me(token):
    req = urllib.request.Request(
        "https://app.asana.com/api/1.0/users/me",
        headers={"Authorization": f"Bearer {token}"},
    )
    with urllib.request.urlopen(req) as r:
        return json.loads(r.read())["data"]

def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def load_token():
    token = os.environ.get("ASANA_TOKEN", "")
    if not token:
        cfg = os.path.expanduser("~/asana_token.config")
        if os.path.exists(cfg):
            with open(cfg) as f:
                token = f.read().strip()
    return token


def sbatch(script_content, job_name, dependency=None, chdir=None):
    cmd = ["sbatch"]
    if dependency:
        cmd += [f"--dependency={dependency}"]
    if chdir:
        cmd += [f"--chdir={chdir}"]
    result = subprocess.run(cmd, input=script_content.encode(),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode == 0:
        jid = result.stdout.decode().strip().split()[-1]
        log(f"  Submitted {job_name}: job {jid}")
        return jid
    log(f"  ERROR submitting {job_name}: {result.stderr.decode().strip()}")
    return None


def get_project_ids(og_ids):
    """Look up project_id per OG from the sample table (run via psycopg2 singularity container)."""
    if not og_ids:
        return {}
    py = (
        "import sys, json, configparser, psycopg2\n"
        "cfg = configparser.ConfigParser(); cfg.read(sys.argv[1])\n"
        "p = cfg['postgres']\n"
        "conn = psycopg2.connect(dbname=p['dbname'], user=p['user'], password=p['password'], host=p['host'], port=p['port'])\n"
        "cur = conn.cursor()\n"
        "cur.execute('SELECT og_id, project_id FROM sample WHERE og_id = ANY(%s)', (sys.argv[2].split(','),))\n"
        "print(json.dumps(dict(cur.fetchall())))\n"
        "conn.close()\n"
    )
    cfg_path = os.path.expanduser(POSTGRES_CFG)
    result = subprocess.run(
        ["singularity", "run", PSYCOPG2_SIF, "python", "-c", py, cfg_path, ",".join(og_ids)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
    )
    if result.returncode != 0:
        log(f"  WARNING: project_id lookup failed: {result.stderr.strip()}")
        return {}
    try:
        return json.loads(result.stdout.strip())
    except Exception as e:
        log(f"  WARNING: could not parse project_id lookup output: {e}")
        return {}


def parse_samplesheet(path, og_ids):
    """Return dict: og_id → {date, version, species}."""
    og_set = set(og_ids) if og_ids else None
    rows = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            og = row["sample"].strip()
            if og_set and og not in og_set:
                continue
            rows[og] = {
                "date":    row.get("date",    "").strip(),
                "version": row.get("version", "").strip(),
                "species": row.get("species", "").strip(),
            }
    return rows


# ── SLURM job bodies ──────────────────────────────────────────────────────────

def chromsyn_job(run_name, og, date, version):
    og_dir = f"{STAGING_BASE}/{og}/09-chromsyn"
    return f"""#!/bin/bash --login
#SBATCH --job-name=chromsyn_{og}_{run_name}
#SBATCH --account={ACCOUNT}
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output={LOG_DIR}/chromsyn_{og}_{run_name}_%j.out
{MAIL_DIRECTIVES}

mkdir -p {og_dir}
cp {CHROMSYN_SH} {og_dir}/
cd {og_dir}
bash chromsyn.sh {og} {date} {version}
echo "Chromsyn complete: {og} {date}.{version}"
"""


def compile_push_job(run_name, samplesheet):
    return f"""#!/bin/bash --login
#SBATCH --job-name=compile_{run_name}
#SBATCH --account={ACCOUNT}
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output={LOG_DIR}/compile_{run_name}_%j.out
{MAIL_DIRECTIVES}

cd {SCRIPTS_DIR}
status=0  # the job exits non-zero (sacct, the compile log) if any push failed, not only the last

echo "=== [1/6] Final GFAstats ==="
bash 04_final_gfastats_compile.sh
singularity run {PSYCOPG2_SIF} python 04a_push_gfa_results_to_sqldb.py {POSTGRES_CFG} {samplesheet} || status=1

echo "=== [2/6] BUSCO ==="
bash 03_busco_compile.sh
singularity run {PSYCOPG2_SIF} python 03a_push_busco_results_to_sqldb.py {POSTGRES_CFG} {samplesheet} || status=1

echo "=== [3/6] Merqury QV ==="
bash 02a_merqury_qv_compile.sh
singularity run {PSYCOPG2_SIF} python 02b_push_merqury_qv_results_to_sqldb.py {POSTGRES_CFG} {samplesheet} || status=1

echo "=== [4/6] Merqury completeness ==="
bash 02c_merqury_completeness_compile.sh
singularity run {PSYCOPG2_SIF} python 02d_push_merqury_completeness_results_to_sqldb.py {POSTGRES_CFG} {samplesheet} || status=1

echo "=== [5/6] OMNIC ==="
bash 05_omnic-compile.sh
singularity run {PSYCOPG2_SIF} python 05a_push_omnic_results_to_sqldb.py {POSTGRES_CFG} {samplesheet} || status=1

echo "=== [6/6] Hi-C runs (flowcells in the Hi-C BAMs; genomes_db_review 094) ==="
singularity run {PSYCOPG2_SIF} python 05b_push_hic_runs_to_sqldb.py {POSTGRES_CFG} {samplesheet} || status=1

echo "Compile + push complete for {run_name} (status $status)"
exit $status
"""


def backup_job(run_name, og, date, version):
    return f"""#!/bin/bash --login
#SBATCH --job-name=backup_{og}_{run_name}
#SBATCH --account={ACCOUNT}
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output={LOG_DIR}/backup_{og}_{run_name}_%j.out
{MAIL_DIRECTIVES}

module load rclone/1.68.1

cd {STAGING_BASE}
bash {BACKUP_SCRIPT} {og} {date} {version} hifi_hic
echo "Backup complete: {og} {date}.{version} [hifi_hic]"
"""


def audit_job(run_name, samplesheet, mode, suffix):
    return f"""#!/bin/bash --login
#SBATCH --job-name=audit_{suffix}_{run_name}
#SBATCH --account={ACCOUNT}
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output={LOG_DIR}/audit_{suffix}_{run_name}_%j.out
{MAIL_DIRECTIVES}

module load rclone/1.68.1

cd {PIPELINE_DIR}/scripts/04_backup_scripts
bash {AUDIT_SCRIPT} -c {samplesheet}
"""


def notify_job(run_name, og_info, project_ids):
    samples_arg  = ",".join(sorted(og_info.keys()))
    og_list      = "\\n".join(f"  {og}: {info['species']}" for og, info in sorted(og_info.items()))
    results_file = f"{LOG_DIR}/results_{run_name}.txt"
    project_ids_json = json.dumps(project_ids)
    notify_py = f"""
import json, os, sys, urllib.request
from datetime import datetime

_cfg = os.path.expanduser('~/asana_token.config')
token = (open(_cfg).read().strip() if os.path.exists(_cfg)
         else os.environ.get('ASANA_TOKEN', ''))
if not token:
    print("No Asana token found — skipping completion notification")
    sys.exit(0)
headers = {{"Authorization": f"Bearer {{token}}", "Content-Type": "application/json"}}
project_ids = json.loads('''{project_ids_json}''')

def post(path, data):
    url  = f"https://app.asana.com/api/1.0{{path}}"
    body = json.dumps({{"data": data}}).encode()
    req  = urllib.request.Request(url, data=body, method="POST", headers=headers)
    with urllib.request.urlopen(req) as r:
        return json.loads(r.read())

assignee_gid = "{ASANA_ASSIGNEE_ME}"

samples  = "{samples_arg}".split(",")
og_lines = "\\n".join(f"  {{s}}" for s in sorted(samples))

# Read QC results and status
results_path  = "{results_file}"
status_path   = "{results_file}.status"
try:
    with open(results_path) as f:
        results_text = f.read().strip()
except Exception as e:
    results_text = f"(Results file not available: {{e}})"

try:
    db_ok = open(status_path).read().strip() == "ok"
except Exception:
    db_ok = False

task_name = (
    "HiFi+HiC assembly of {run_name} complete — results in database"
    if db_ok else
    "HiFi+HiC assembly of {run_name} complete — WARNING: compile failed, results NOT in database"
)

notes = (
    f"Run: {run_name}\\n"
    f"Assembly mode: hifi_hic\\n\\n"
    f"Samples assembled ({{len(samples)}}):\\n{{og_lines}}\\n\\n"
    f"Completed:\\n"
    f"- Nextflow hifi_hic pipeline\\n"
    f"- Chromsyn synteny plots\\n"
    f"- Backed up to Acacia (pawsey0964:oceanomics-refassemblies)\\n\\n"
    + ("- Compile + push to OceanOmics DB: SUCCESS\\n\\n" if db_ok else
       "- Compile + push to OceanOmics DB: FAILED — re-run compile scripts manually\\n\\n")
    + f"QC Results:\\n{{results_text}}\\n\\n"
    f"Full results file: {{results_path}}\\n\\n"
    f"Date: {{datetime.now().strftime('%Y-%m-%d %H:%M AWST')}}"
)

result = post("/tasks", {{
    "name":      task_name,
    "notes":     notes,
    "assignee":  assignee_gid,
    "projects":  ["{ASANA_PROJECT_GID}"],
    "memberships": [{{"project": "{ASANA_PROJECT_GID}", "section": "{ASANA_SECTION_RUN_QC}"}}],
}})
print(f"Asana task created: {{result['data']['gid']}}")

# Create curation subtasks for each OG (only if results are in DB)
CURATION_PARENT_GID = "{CURATION_PARENT_GID}"
if db_ok:
    for og in sorted(samples):
        proj = project_ids.get(og, "")
        og_task = post("/tasks", {{
            "name":   f"{{og}} Project: {{proj}}" if proj else og,
            "notes":  f"OG assignment\\nProject: {{proj}}" if proj else "OG assignment",
            "parent": CURATION_PARENT_GID,
        }})
        og_gid = og_task["data"]["gid"]
        # Create subtasks one at a time to preserve order (agp first, reviewer second)
        post("/tasks", {{"name": "agp & save state generated", "parent": og_gid}})
        import time; time.sleep(1)
        post("/tasks", {{"name": "Tag reviewer", "parent": og_gid}})
        print(f"Created curation subtask for {{og}} ({{og_gid}})")
else:
    print("Skipping curation subtasks — results not in DB")
"""
    return f"""#!/bin/bash --login
#SBATCH --job-name=notify_{run_name}
#SBATCH --account={ACCOUNT}
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output={LOG_DIR}/notify_{run_name}_%j.out
{MAIL_DIRECTIVES}

module load singularity/4.1.0-slurm

# Query DB for assembly QC results; exit code 0 = results present, 1 = empty
if singularity run {PSYCOPG2_SIF} python {QUERY_SCRIPT} \\
    {POSTGRES_CFG} "{samples_arg}" "{results_file}"; then
    echo "ok" > "{results_file}.status"
else
    echo "failed" > "{results_file}.status"
fi

python3 - <<'PYEOF'
{notify_py}
PYEOF
echo "Asana notification sent for {run_name}"
"""


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    global STAGING_BASE, PIPELINE_DIR, SCRIPTS_DIR, BACKUP_SCRIPT, AUDIT_SCRIPT
    global CHROMSYN_SH, LOG_DIR, EMAIL, MAIL_DIRECTIVES

    args = sys.argv[1:]

    run_name    = args[args.index("--run")        + 1] if "--run"        in args else None
    og_ids_arg  = args[args.index("--og-ids")     + 1] if "--og-ids"     in args else ""
    samplesheet = args[args.index("--samplesheet") + 1] if "--samplesheet" in args else None

    if not run_name or not samplesheet:
        print("Usage: post-hichifi-assembly.py --run <RUN> --og-ids <OG1,...> --samplesheet <path> "
              "[--staging-base DIR] [--pipeline-dir DIR] [--log-dir DIR] [--email ADDR]")
        sys.exit(1)

    # CLI overrides for the portable path defaults set at module load time
    if "--staging-base" in args:
        STAGING_BASE = args[args.index("--staging-base") + 1]
        PIPELINE_DIR = f"{STAGING_BASE}/OceanOmics-OceanGenomes-ref-genomes"
    if "--pipeline-dir" in args:
        PIPELINE_DIR = args[args.index("--pipeline-dir") + 1]
    SCRIPTS_DIR   = f"{PIPELINE_DIR}/scripts/03_compile_results"
    BACKUP_SCRIPT = f"{PIPELINE_DIR}/scripts/04_backup_scripts/full_assembly_backup.sh"
    AUDIT_SCRIPT  = f"{PIPELINE_DIR}/scripts/04_backup_scripts/audit_backup.sh"
    CHROMSYN_SH   = f"{PIPELINE_DIR}/scripts/02_run_chromsyn/chromsyn.sh"
    if "--log-dir" in args:
        LOG_DIR = args[args.index("--log-dir") + 1]
    if "--email" in args:
        EMAIL = args[args.index("--email") + 1]
        MAIL_DIRECTIVES = f"#SBATCH --mail-type=FAIL\n#SBATCH --mail-user={EMAIL}" if EMAIL else ""

    og_ids  = [x.strip() for x in og_ids_arg.split(",") if x.strip()]
    og_info = parse_samplesheet(samplesheet, og_ids)

    if not og_info:
        log(f"ERROR: no matching OGs in samplesheet {samplesheet}")
        sys.exit(1)

    os.makedirs(LOG_DIR, exist_ok=True)
    log(f"Post-assembly (hifi_hic) for {run_name}: {len(og_info)} OG(s): {', '.join(sorted(og_info))}")

    # 1. Chromsyn — one SLURM job per sample, staggered via array or sequential submit
    chromsyn_jids = []
    for i, (og, info) in enumerate(sorted(og_info.items())):
        if not info["date"] or not info["version"]:
            log(f"  SKIP chromsyn for {og}: missing date or version in samplesheet")
            continue
        jid = sbatch(
            chromsyn_job(run_name, og, info["date"], info["version"]),
            f"chromsyn_{og}",
            chdir=f"{STAGING_BASE}/{og}",
        )
        if jid:
            chromsyn_jids.append(jid)
        # Stagger submissions so the first job's compleasm download completes
        # before the next job starts (avoids eukaryota_odb12.tmp lock collision)
        if i < len(og_info) - 1:
            import time; time.sleep(120)

    # 2. Compile + push — runs immediately (does not depend on chromsyn)
    compile_jid = sbatch(compile_push_job(run_name, samplesheet), f"compile_{run_name}")

    # 3. Backup — wait for ALL chromsyn jobs so PDF is included
    backup_dep  = f"afterok:{':'.join(chromsyn_jids)}" if chromsyn_jids else None
    backup_jids = []
    for og, info in sorted(og_info.items()):
        if not info["date"] or not info["version"]:
            log(f"  SKIP backup for {og}: missing date or version")
            continue
        jid = sbatch(
            backup_job(run_name, og, info["date"], info["version"]),
            f"backup_{og}",
            dependency=backup_dep,
            chdir=STAGING_BASE,
        )
        if jid:
            backup_jids.append(jid)

    if not backup_jids:
        log("ERROR: no backup jobs submitted — aborting post-audit and notify")
        sys.exit(1)

    # 4. Post-backup audit
    all_backups   = ":".join(backup_jids)
    postaudit_jid = sbatch(
        audit_job(run_name, samplesheet, "remote", "post"),
        f"postaudit_{run_name}",
        dependency=f"afterok:{all_backups}",
    )

    # 5. Asana notify
    if postaudit_jid:
        token = load_token()
        if token:
            project_ids = get_project_ids(list(og_info.keys()))
            sbatch(
                notify_job(run_name, og_info, project_ids),
                f"notify_{run_name}",
                dependency=f"afterok:{postaudit_jid}",
            )
        else:
            log("No Asana token found — skipping completion notification")

    log(f"All jobs submitted for {run_name}.")
    log(f"  Chromsyn:   {', '.join(chromsyn_jids) if chromsyn_jids else 'none'}")
    log(f"  Compile:    {compile_jid}")
    log(f"  Backups:    {', '.join(backup_jids)}")
    log(f"  Post-audit: {postaudit_jid}")


if __name__ == "__main__":
    main()
