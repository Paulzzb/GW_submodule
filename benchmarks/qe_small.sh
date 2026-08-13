#!/bin/bash -l
#SBATCH -J bench_qe_gs_small
#SBATCH -t 24:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -c 1
#SBATCH -o log/qe_gs_small-%j.out
#SBATCH -e log/qe_gs_small-%j.err
#
# Ground-state (QE) for small cells: Si8, STO3 — scf → nscf → pw2bgw.
#
# Submit from benchmarks/:
#   mkdir -p log
#   sbatch slurm_qe_gs_small.sh
#
# Override binaries / MPI:
#   PW=pw.x PW2BGW=pw2bgw.x NPROC=64 sbatch slurm_qe_gs_small.sh
#   CASES="Si8" sbatch slurm_qe_gs_small.sh
#   PW_FLAGS= sbatch ...          # clear default -ndiag 4 if needed

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  ROOT="$SLURM_SUBMIT_DIR"
else
  ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
cd "$ROOT" || exit 1
mkdir -p log

PW="${PW:-pw.x}"
PW2BGW="${PW2BGW:-pw2bgw.x}"
NPROC="${NPROC:-${SLURM_NTASKS:-32}}"
MPI="${MPI:-mpirun -np ${NPROC}}"
# Serial subspace diag: avoids ScaLAPACK cholesky failures on old builds.
PW_FLAGS="${PW_FLAGS--ndiag 4}"
CASES="${CASES:-Si8 STO3}"

echo "[$(date)] QE ground-state (small cells)"
echo "  ROOT=$ROOT"
echo "  PW=$PW  PW2BGW=$PW2BGW  MPI=$MPI  PW_FLAGS=$PW_FLAGS"
echo "  CASES=$CASES"
echo "  which PW: $(command -v "$PW" 2>/dev/null || echo NOT_FOUND)"
# Old QE: stdin redirect (pw.x < in), not "pw.x -in".

for case in $CASES; do
  case_dir="$ROOT/$case"
  if [[ ! -d "$case_dir" ]]; then
    echo "ERROR: missing case dir $case_dir" >&2
    exit 1
  fi
  echo
  echo "===== [$case] $(date) ====="
  cd "$case_dir" || exit 1
  echo "[$case] pwd=$(pwd)"

  for f in scf.in nscf.in pp_in; do
    if [[ ! -f "$f" ]]; then
      echo "ERROR: $case_dir/$f missing" >&2
      exit 1
    fi
  done

  echo "[$case] scf ..."
  $MPI "$PW" $PW_FLAGS < scf.in > scf.out
  echo "[$case] scf done, scf.out=$(wc -c < scf.out) bytes"

  echo "[$case] nscf ..."
  mpirun -np 4 pw.x < nscf.in > nscf.out
  echo "[$case] nscf done, nscf.out=$(wc -c < nscf.out) bytes"

  echo "[$case] pw2bgw (vxc.dat) ..."
  $MPI "$PW2BGW" < pp_in > pp.out

  # pw2bgw writes vxc.dat in cwd; MATLAB reads groundstate_dir/vxc.dat (*.save/).
  save_dir="$(ls -d ./*.save 2>/dev/null | head -n 1 || true)"
  if [[ -z "${save_dir:-}" ]]; then
    echo "ERROR: no *.save after pw2bgw" >&2
    exit 1
  fi
  if [[ ! -f vxc.dat ]]; then
    echo "ERROR: no ./vxc.dat after pw2bgw — check pp.out" >&2
    exit 1
  fi
  cp -f vxc.dat "$save_dir/vxc.dat"
  echo "[$case] copied vxc.dat → $save_dir/"

  echo "[$case] done:"
  ls -ld "$save_dir"
  ls -lh "$save_dir/vxc.dat"
done

cd "$ROOT"
echo
echo "[$(date)] Small-cell QE ground-state finished."
