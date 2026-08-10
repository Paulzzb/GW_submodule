#!/usr/bin/env bash
# Link adaptive ISDF checkpoints between Si_gamma* SAVE dirs.
#
# Usage:
#   ./link_isdf_checkpoints.sh              # hub → ratio cases (pre-run)
#   ./link_isdf_checkpoints.sh --case NAME  # links needed for one case (after clean)
#   ./link_isdf_checkpoints.sh --vn         # 162416_exact vn → other 162416*
#   ./link_isdf_checkpoints.sh --unlink     # remove checkpoint symlinks in ratio cases
#
# Rules:
#   162416* (ratios 16/24/16):
#     - *vc* and *nn* from cases/Si_gamma/SAVE
#     - *vn* from Si_gamma_162416_exact/SAVE (after that case has run)
#   161616* (ratios 16/16/16):
#     - all isdf_adaptive_checkpoint_* from cases/Si_gamma/SAVE

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASES_DIR="${SCRIPT_DIR}/cases"
HUB_SAVE="${CASES_DIR}/Si_gamma/SAVE"
EXACT_162416_SAVE="${CASES_DIR}/Si_gamma_162416_exact/SAVE"

CASES_162416=(
  Si_gamma_162416_exact
  Si_gamma_162416_sum
  Si_gamma_162416_ff
)
CASES_161616=(
  Si_gamma_161616_exact
  Si_gamma_161616_sum
  Si_gamma_161616_ff
)
CASES_162416_VN_CONSUMERS=(
  Si_gamma_162416_sum
  Si_gamma_162416_ff
)

mode="all"
case_only=""
if [[ "${1:-}" == "--unlink" ]]; then
  mode="unlink"
elif [[ "${1:-}" == "--vn" ]]; then
  mode="vn"
elif [[ "${1:-}" == "--case" ]]; then
  mode="case"
  case_only="${2:-}"
  if [[ -z "$case_only" ]]; then
    echo "ERROR: --case requires a case name" >&2
    exit 1
  fi
fi

link_glob_from_hub() {
  local dest_save="$1"
  shift
  local patterns=("$@")
  mkdir -p "$dest_save"
  if [[ ! -d "$HUB_SAVE" ]]; then
    echo "WARN: hub SAVE missing: $HUB_SAVE (skip file links)"
    return 0
  fi
  local pat f base
  for pat in "${patterns[@]}"; do
    shopt -s nullglob
    for f in "${HUB_SAVE}/"${pat}; do
      [[ -e "$f" || -L "$f" ]] || continue
      base="$(basename "$f")"
      rm -f "${dest_save}/${base}"
      ln -sfn "../Si_gamma/SAVE/${base}" "${dest_save}/${base}"
      echo "OK  $(basename "$(dirname "$dest_save")")/SAVE/${base} -> ../Si_gamma/SAVE/${base}"
    done
    shopt -u nullglob
  done
}

link_vn_from_exact() {
  local dest_save="$1"
  mkdir -p "$dest_save"
  if [[ ! -d "$EXACT_162416_SAVE" ]]; then
    echo "WARN: 162416_exact SAVE missing: $EXACT_162416_SAVE"
    return 0
  fi
  local f base
  shopt -s nullglob
  for f in "${EXACT_162416_SAVE}/"isdf_adaptive_checkpoint_*vn*; do
    [[ -e "$f" || -L "$f" ]] || continue
    # do not follow if exact itself somehow pointed elsewhere incorrectly
    base="$(basename "$f")"
    rm -f "${dest_save}/${base}"
    ln -sfn "../Si_gamma_162416_exact/SAVE/${base}" "${dest_save}/${base}"
    echo "OK  $(basename "$(dirname "$dest_save")")/SAVE/${base} -> ../Si_gamma_162416_exact/SAVE/${base}"
  done
  shopt -u nullglob
}

unlink_checkpoints() {
  local dest_save="$1"
  [[ -d "$dest_save" ]] || return 0
  local f
  shopt -s nullglob
  for f in "${dest_save}/"isdf_adaptive_checkpoint_*; do
    if [[ -L "$f" ]]; then
      rm -f "$f"
      echo "removed symlink: $f"
    fi
  done
  shopt -u nullglob
}

prepare_case() {
  local name="$1"
  local dest="${CASES_DIR}/${name}/SAVE"
  case "$name" in
    Si_gamma_162416_exact)
      link_glob_from_hub "$dest" 'isdf_adaptive_checkpoint_*vc*' 'isdf_adaptive_checkpoint_*nn*'
      ;;
    Si_gamma_162416_sum|Si_gamma_162416_ff)
      link_glob_from_hub "$dest" 'isdf_adaptive_checkpoint_*vc*' 'isdf_adaptive_checkpoint_*nn*'
      link_vn_from_exact "$dest"
      ;;
    Si_gamma_161616_exact|Si_gamma_161616_sum|Si_gamma_161616_ff)
      link_glob_from_hub "$dest" 'isdf_adaptive_checkpoint_*'
      ;;
    *)
      echo "SKIP (not a ratio case): $name"
      ;;
  esac
}

if [[ "$mode" == "unlink" ]]; then
  for name in "${CASES_162416[@]}" "${CASES_161616[@]}"; do
    unlink_checkpoints "${CASES_DIR}/${name}/SAVE"
  done
  echo "=== unlinked adaptive checkpoint symlinks ==="
  exit 0
fi

if [[ "$mode" == "vn" ]]; then
  for name in "${CASES_162416_VN_CONSUMERS[@]}"; do
    link_vn_from_exact "${CASES_DIR}/${name}/SAVE"
  done
  echo "=== linked vn from Si_gamma_162416_exact ==="
  exit 0
fi

if [[ "$mode" == "case" ]]; then
  prepare_case "$case_only"
  exit 0
fi

# default: prepare all
for name in "${CASES_162416[@]}" "${CASES_161616[@]}"; do
  if [[ ! -d "${CASES_DIR}/${name}" ]]; then
    echo "SKIP (missing case): $name"
    continue
  fi
  echo "--- $name ---"
  prepare_case "$name"
done
echo "=== ISDF checkpoint links ready ==="
