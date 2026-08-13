#!/bin/bash -l
#SBATCH -J bench_gw_small
#SBATCH -t 48:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -o log/gw_small-%j.out
#SBATCH -e log/gw_small-%j.err
#
# GW for small cells: Si8, STO3 (MATLAB): input_driver + qp.launcher.
#
# Prerequisites: QE ground state finished (see slurm_qe_gs_small.sh).
#
# Submit from benchmarks/:
#   mkdir -p log
#   sbatch slurm_gw_small.sh
#
# Overrides:
#   MATLAB_BIN=/path/to/matlab sbatch slurm_gw_small.sh
#   CASES="Si8" NAMELIST=test_dir sbatch slurm_gw_small.sh
#   REPO_ROOT=/path/to/hefeikssolvtopub sbatch slurm_gw_small.sh
#
# Namelist output_dir keeps results apart: ./isdf (test_isdf) vs ./direct (test_dir).

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  BENCH_ROOT="$SLURM_SUBMIT_DIR"
else
  BENCH_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
cd "$BENCH_ROOT" || exit 1
mkdir -p log

REPO_ROOT="${REPO_ROOT:-$(cd "$BENCH_ROOT/.." && pwd)}"
MATLAB_BIN="${MATLAB_BIN:-${RUNMATLAB:-matlab}}"
CASES="${CASES:-Si8 STO3}"
NAMELIST="${NAMELIST:-test_isdf}"

if ! command -v "$MATLAB_BIN" >/dev/null 2>&1; then
  echo "ERROR: MATLAB not found: ${MATLAB_BIN}" >&2
  exit 1
fi

echo "[$(date)] GW benchmarks (small cells)"
echo "  BENCH_ROOT=$BENCH_ROOT"
echo "  REPO_ROOT=$REPO_ROOT"
echo "  MATLAB=$MATLAB_BIN"
echo "  CASES=$CASES  NAMELIST=$NAMELIST"

for case in $CASES; do
  case_dir="$BENCH_ROOT/$case"
  namelist_path="$case_dir/$NAMELIST"
  if [[ ! -f "$namelist_path" ]]; then
    echo "ERROR: missing $namelist_path" >&2
    exit 1
  fi

  echo
  echo "===== [$case] $(date) ====="
  cd "$case_dir"

  gs_rel="$(awk -F"'" '/groundstate_dir/ {print $2; exit}' "$NAMELIST" || true)"
  if [[ -n "$gs_rel" && ! -d "$gs_rel" ]]; then
    echo "WARN: groundstate_dir='$gs_rel' not found under $case_dir" >&2
    echo "      Run slurm_qe_gs_small.sh first (expect Si8.save / STO3.save)." >&2
  fi

  logf="$BENCH_ROOT/log/gw_${case}_${NAMELIST}_$(date +%Y%m%d_%H%M%S).log"
  echo "[$case] logging to $logf"

  "$MATLAB_BIN" -batch "
    cd('${REPO_ROOT}');
    QPstartup;
    cd('${case_dir}');
    fprintf('[GW] case=%s namelist=%s\\n', '${case}', '${NAMELIST}');
    input_driver('./${NAMELIST}');
    load('./SAVE/config.mat', 'config');
    E = qp.launcher(config);
    fprintf('[GW] done case=%s  fout=%s\\n', '${case}', string(E.fout));
  " >"$logf" 2>&1

  echo "[$case] finished. qp.dat under output_dir:"
  out_rel="$(awk -F"'" '/output_dir/ {print $2; exit}' "$NAMELIST")"
  ls -lh "${out_rel}/qp.dat" 2>/dev/null || echo "  (no qp.dat — see $logf)"
done

cd "$BENCH_ROOT"
echo
echo "[$(date)] Small-cell GW finished."
