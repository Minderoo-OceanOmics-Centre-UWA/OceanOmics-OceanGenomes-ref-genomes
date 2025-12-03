#!/usr/bin/env bash
set -euo pipefail

AUDITOR="$(dirname "$0")/audit_backup.sh"
CSV_FILE="/scratch/pawsey0964/lhuet/ref-gen/OceanOmics-OceanGenomes-ref-genomes/assets/samplesheet.csv"
REMOTE="pawsey0964:oceanomics-refassemblies"

# You can override with: ./audit_loop.sh /path/to/sheet.csv
if [[ $# -ge 1 ]]; then
  CSV_FILE="$1"
fi

if [[ ! -x "${AUDITOR}" ]]; then
  echo "Error: auditor not found or not executable at ${AUDITOR}"
  exit 1
fi

"${AUDITOR}" -c "${CSV_FILE}" -r "${REMOTE}"
