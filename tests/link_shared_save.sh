#!/usr/bin/env bash
# Link shareable Si_gamma* case SAVE dirs to the hub cases/Si_gamma/SAVE.
#
# Consumers (same GS + pseudo ISDF; differ by fdep / iscauchy / prefix):
#   Si_gamma_ff_cauchy
#   Si_gamma_ff_nocauchy
#   Si_gamma_cohsex_cauchy
#   Si_gamma_cohsex_nocauchy
#
# Usage:
#   ./link_shared_save.sh           # create / refresh links
#   ./link_shared_save.sh --unlink  # remove consumer SAVE symlinks only

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASES_DIR="${SCRIPT_DIR}/cases"
HUB_REL="../Si_gamma/SAVE"
HUB_ABS="${CASES_DIR}/Si_gamma/SAVE"

CONSUMERS=(
  Si_gamma_ff_cauchy
  Si_gamma_ff_nocauchy
  Si_gamma_cohsex_cauchy
  Si_gamma_cohsex_nocauchy
)

unlink_only=0
if [[ "${1:-}" == "--unlink" ]]; then
  unlink_only=1
fi

if [[ ! -d "${CASES_DIR}/Si_gamma" ]]; then
  echo "ERROR: hub case missing: ${CASES_DIR}/Si_gamma" >&2
  exit 1
fi

for name in "${CONSUMERS[@]}"; do
  case_dir="${CASES_DIR}/${name}"
  dest="${case_dir}/SAVE"
  if [[ ! -d "$case_dir" ]]; then
    echo "SKIP (missing case): ${name}"
    continue
  fi

  if [[ -L "$dest" ]]; then
    rm -f "$dest"
    echo "removed symlink: ${dest}"
  elif [[ -e "$dest" ]]; then
    if [[ "$unlink_only" -eq 1 ]]; then
      echo "SKIP (real SAVE, not symlink): ${dest}"
      continue
    fi
    echo "removing real SAVE before link: ${dest}"
    rm -rf "$dest"
  fi

  if [[ "$unlink_only" -eq 1 ]]; then
    continue
  fi

  ln -sfn "$HUB_REL" "$dest"
  echo "OK  ${name}/SAVE -> ${HUB_REL}"
done

if [[ "$unlink_only" -eq 1 ]]; then
  echo "=== unlinked consumer SAVE symlinks ==="
  exit 0
fi

echo
if [[ -d "$HUB_ABS" ]]; then
  echo "hub SAVE present: ${HUB_ABS}"
else
  echo "note: hub SAVE not created yet (${HUB_ABS}); run_all will build it when Si_gamma runs first."
fi
echo "=== linked ${#CONSUMERS[@]} consumer(s) to ${HUB_REL} ==="
