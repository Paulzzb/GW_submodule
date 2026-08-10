#!/usr/bin/env bash
# Collect qp.dat / qp_*.dat from example/cases into $CMD/result.
#
# Usage:
#   ./collect_qp_dat.sh                  # CMD = this script's directory (example/)
#   ./collect_qp_dat.sh /path/to/example
#   CMD=/path/to/example ./collect_qp_dat.sh
#
# Output layout:
#   $CMD/result/<case_name>/qp.dat
#   $CMD/result/<case_name>/qp_*.dat

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CMD="${1:-${CMD:-$SCRIPT_DIR}}"
CASES_DIR="${CMD}/cases"
RESULT_DIR="${CMD}/result"

if [[ ! -d "$CASES_DIR" ]]; then
  echo "ERROR: cases dir not found: $CASES_DIR" >&2
  exit 1
fi

mkdir -p "$RESULT_DIR"

n_ok=0
n_skip=0

shopt -s nullglob
for case_dir in "$CASES_DIR"/*/; do
  [[ -d "$case_dir" ]] || continue
  case_name="$(basename "$case_dir")"
  hits=("$case_dir"qp.dat "$case_dir"qp_*.dat)
  found=0
  for f in "${hits[@]}"; do
    [[ -f "$f" ]] || continue
    found=1
    dest_dir="${RESULT_DIR}/${case_name}"
    mkdir -p "$dest_dir"
    cp -f "$f" "${dest_dir}/$(basename "$f")"
    echo "OK  ${case_name}/$(basename "$f")  ->  ${dest_dir}/$(basename "$f")"
    n_ok=$((n_ok + 1))
  done
  if [[ "$found" -eq 0 ]]; then
    echo "SKIP (no qp*.dat): ${case_name}"
    n_skip=$((n_skip + 1))
  fi
done
shopt -u nullglob

echo
echo "=== collected ${n_ok} file(s) into ${RESULT_DIR}  (skipped cases: ${n_skip}) ==="
ls -la "$RESULT_DIR"
