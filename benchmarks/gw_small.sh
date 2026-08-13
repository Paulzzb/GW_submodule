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
# Prerequisites: QE ground state finished (see qe_small.sh).
#
# Submit from benchmarks/:
#   mkdir -p log
#   sbatch gw_small.sh
#
# Overrides:
#   MATLAB_BIN=/path/to/matlab sbatch gw_small.sh
#   CASES="Si8" NAMELIST=test_dir sbatch gw_small.sh
#   REPO_ROOT=/path/to/repo sbatch gw_small.sh
#
# Namelist output_dir: ./isdf (test_isdf) vs ./direct (test_dir).

set -euo pipefail

# Quiet harmless "X does not support locale zh_CN.UTF-8" on old clusters.
export LANG=C
export LC_ALL=C

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

# Old MATLAB may lack -batch (R2019a+). Prefer -batch when present.
if "$MATLAB_BIN" -help 2>&1 | grep -q -- '-batch'; then
  MATLAB_MODE=batch
else
  MATLAB_MODE=r
fi

echo "[$(date)] GW benchmarks (small cells)"
echo "  BENCH_ROOT=$BENCH_ROOT"
echo "  REPO_ROOT=$REPO_ROOT"
echo "  MATLAB=$MATLAB_BIN  mode=-${MATLAB_MODE}"
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
  if [[ -z "$gs_rel" || ! -d "$gs_rel" ]]; then
    echo "ERROR: groundstate_dir='$gs_rel' missing under $case_dir" >&2
    echo "      Run qe_small.sh first (expect Si8.save / STO3.save)." >&2
    exit 1
  fi
  # pw2bgw often leaves vxc.dat in case cwd; loader needs it inside *.save/
  if [[ ! -f "$gs_rel/vxc.dat" && -f vxc.dat ]]; then
    cp -f vxc.dat "$gs_rel/vxc.dat"
    echo "[$case] copied ./vxc.dat → $gs_rel/vxc.dat"
  fi
  if [[ ! -f "$gs_rel/vxc.dat" ]]; then
    echo "ERROR: missing $gs_rel/vxc.dat (and no ./vxc.dat to copy)" >&2
    exit 1
  fi

  stamp="$(date +%Y%m%d_%H%M%S)"
  logf="$BENCH_ROOT/log/gw_${case}_${NAMELIST}_${stamp}.log"
  mlscript="$BENCH_ROOT/log/run_${case}_${NAMELIST}_${stamp}.m"
  echo "[$case] logging to $logf"
  echo "[$case] matlab script $mlscript"

  # Write a .m file so argv cannot lose the -batch/-r statement (old shells/wrappers).
  cat >"$mlscript" <<EOF
cd('${REPO_ROOT}');
QPstartup;
cd('${case_dir}');
fprintf('[GW] case=%s namelist=%s\n', '${case}', '${NAMELIST}');
input_driver('./${NAMELIST}');
load('./SAVE/config.mat', 'config');
E = qp.launcher(config);
fprintf('[GW] done case=%s  fout=%s\n', '${case}', char(E.fout));
EOF

  if [[ "$MATLAB_MODE" == batch ]]; then
    "$MATLAB_BIN" -nodisplay -nosplash -nodesktop \
      -batch "run('${mlscript}');" >"$logf" 2>&1
  else
    "$MATLAB_BIN" -nodisplay -nosplash -nodesktop \
      -r "try, run('${mlscript}'); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);" \
      >"$logf" 2>&1
  fi

  echo "[$case] finished. qp.dat under output_dir:"
  out_rel="$(awk -F"'" '/output_dir/ {print $2; exit}' "$NAMELIST")"
  ls -lh "${out_rel}/qp.dat" 2>/dev/null || echo "  (no qp.dat — see $logf)"
done

cd "$BENCH_ROOT"
echo
echo "[$(date)] Small-cell GW finished."
