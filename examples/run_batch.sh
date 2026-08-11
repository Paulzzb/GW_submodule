#!/bin/bash -l
#SBATCH -J examples_run_all
#SBATCH -t 04:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -o log/slurm-%j.out
#SBATCH -e log/slurm-%j.err
#
# SLURM / batch: profiled examples/run_all (ISDF vs dense demo).
#
# Submit from examples/:
#   mkdir -p log
#   sbatch run_batch.sh
#
# Interactive:
#   bash run_batch.sh
#   bash run_batch.sh /path/to/matlab
#
# Cases share cases/qe.save and cases/SAVE (see namelist storage_dir).

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
cd "$SCRIPT_DIR" || exit 1

MATLAB_BIN="${1:-${RUNMATLAB:-matlab}}"
PROFILE_DIR="${SCRIPT_DIR}/profile_run_all"
LOG_DIR="${SCRIPT_DIR}/log"
LOG_FILE="${LOG_DIR}/run_batch.log"

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

echo "=== matlab run_all_profiled ==="
"$MATLAB_BIN" -nodisplay -nosplash -nodesktop -batch \
  "cd('${SCRIPT_DIR}'); run_all_profiled('${PROFILE_DIR}');" \
  2>&1 | tee "$LOG_FILE"

echo
echo "=== done ==="
echo "log:     ${LOG_FILE}"
echo "slurm:   ${LOG_DIR}/slurm-${SLURM_JOB_ID:-local}.{out,err}"
echo "profile: ${PROFILE_DIR}/file0.html"
