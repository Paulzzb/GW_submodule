#!/bin/bash -l
#SBATCH -J run_all_shared_save
#SBATCH -t 12:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -o log/slurm-%j.out
#SBATCH -e log/slurm-%j.err
#
# SLURM job: link shared SAVE → Si_gamma/SAVE, then profiled run_all.
#
# Submit from example/:
#   mkdir -p log
#   sbatch run_all_shared_save.sh
#
# Or run interactively (no SLURM):
#   bash run_all_shared_save.sh
#   bash run_all_shared_save.sh /path/to/matlab
#
# Steps:
#   1. Soft-link Si_gamma_{ff,cohsex}_* /SAVE → ../Si_gamma/SAVE
#   2. matlab -batch run_all_profiled  (profile on + profsave)
#
# Requires: hub Si_gamma first in run_all.m. clean_case_outputs keeps SAVE
# symlinks (does not wipe hub through the link).

set -euo pipefail

# Prefer submit dir when launched via sbatch; else script location.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
cd "$SCRIPT_DIR" || exit 1

MATLAB_BIN="${1:-${RUNMATLAB:-matlab}}"
PROFILE_DIR="${SCRIPT_DIR}/profile_run_all"
LOG_DIR="${SCRIPT_DIR}/log"
LOG_FILE="${LOG_DIR}/run_all_shared_save.log"

mkdir -p "$LOG_DIR" "$PROFILE_DIR"

if ! command -v "$MATLAB_BIN" >/dev/null 2>&1; then
  echo "ERROR: MATLAB not found: ${MATLAB_BIN}" >&2
  exit 1
fi

echo "=== host $(hostname)  job=${SLURM_JOB_ID:-local} ==="
echo "SCRIPT_DIR = ${SCRIPT_DIR}"
echo "matlab     = $(command -v "$MATLAB_BIN")"
echo "cpus       = ${SLURM_CPUS_PER_TASK:-n/a}"
echo "prof       = ${PROFILE_DIR}"
echo "log        = ${LOG_FILE}"
echo

echo "=== link shared SAVE (cauchy cluster) ==="
bash "${SCRIPT_DIR}/link_shared_save.sh"

echo
echo "=== link ISDF checkpoints (16* ratio cases; vn after exact) ==="
bash "${SCRIPT_DIR}/link_isdf_checkpoints.sh"

echo
echo "=== matlab run_all_profiled ==="
# -batch exits non-zero if the MATLAB expression errors.
"$MATLAB_BIN" -nodisplay -nosplash -nodesktop -batch \
  "cd('${SCRIPT_DIR}'); run_all_profiled('${PROFILE_DIR}');" \
  2>&1 | tee "$LOG_FILE"

echo
echo "=== done ==="
echo "log:     ${LOG_FILE}"
echo "slurm:   ${LOG_DIR}/slurm-${SLURM_JOB_ID:-local}.{out,err}"
echo "profile: ${PROFILE_DIR}/file0.html"
