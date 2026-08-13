#!/bin/bash -l
#SBATCH -J bench_qe_gs_large
#SBATCH -t 72:00:00
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -c 1
#SBATCH -o log/qe_gs_large-%j.out
#SBATCH -e log/qe_gs_large-%j.err
#
# Ground-state (QE) for large cells: Si64, STO3_8 — scf → nscf → pw2bgw.
#
# Submit from benchmarks/:
#   mkdir -p log
#   sbatch slurm_qe_gs_large.sh
#
# Override binaries / MPI:
#   PW=pw.x PW2BGW=pw2bgw.x NPROC=128 sbatch slurm_qe_gs_large.sh
#   CASES="Si64" sbatch slurm_qe_gs_large.sh

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
NPROC="${NPROC:-${SLURM_NTASKS:-64}}"
MPI="${MPI:-mpirun -np ${NPROC}}"
CASES="${CASES:-Si64 STO3_8}"

echo "[$(date)] QE ground-state (large cells)"
echo "  ROOT=$ROOT"
echo "  PW=$PW  PW2BGW=$PW2BGW  MPI=$MPI"
echo "  CASES=$CASES"

for case in $CASES; do
  case_dir="$ROOT/$case"
  if [[ ! -d "$case_dir" ]]; then
    echo "ERROR: missing case dir $case_dir" >&2
    exit 1
  fi
  echo
  echo "===== [$case] $(date) ====="
  cd "$case_dir"

  for f in scf.in nscf.in pp_in; do
    if [[ ! -f "$f" ]]; then
      echo "ERROR: $case_dir/$f missing" >&2
      exit 1
    fi
  done

  echo "[$case] scf ..."
  $MPI "$PW" -in scf.in > scf.out

  echo "[$case] nscf ..."
  $MPI "$PW" -in nscf.in > nscf.out

  echo "[$case] pw2bgw (vxc.dat) ..."
  $MPI "$PW2BGW" -in pp_in > pp.out

  echo "[$case] done. save dirs:"
  ls -d ./*.save 2>/dev/null || echo "  (no *.save found — check prefix/outdir in scf/nscf)"
  ls -lh vxc.dat 2>/dev/null || echo "  (no vxc.dat — check pp.out)"
done

cd "$ROOT"
echo
echo "[$(date)] Large-cell QE ground-state finished."
