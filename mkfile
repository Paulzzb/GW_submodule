#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage:"
  echo "  newfile <path>            # relative to repo root"
  echo "  newfile --here <filename> # create in current directory"
}

# ---- locate repo root (preferred) ----
repo_root=""
if repo_root="$(git rev-parse --show-toplevel 2>/dev/null)"; then
  proj_root="$repo_root/GW"
  echo "In if branch"
  : # ok
else
  # fallback: script's parent as project root
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  repo_root="$(cd "$script_dir/" && pwd)"
  proj_root="$repo_root/GW"
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "$script_dir"
echo "$proj_root"
tpl="$proj_root/licence/header"

if [ $# -lt 1 ]; then
  usage
  exit 1
fi

mode="root"
if [ "${1:-}" = "--here" ]; then
  mode="here"
  shift
fi

if [ $# -lt 1 ]; then
  usage
  exit 1
fi

arg="$1"

# ---- decide output path ----
out="$(pwd)/$arg"
# if [ "$mode" = "here" ]; then
#   out="$(pwd)/$arg"
# else
#   # treat arg as path relative to repo root (or absolute path)
#   if [[ "$arg" = /* ]]; then
#     out="$arg"
#   else
#     out="$proj_root/$arg"
#   fi
# fi

out_dir="$(dirname "$out")"
echo "$out_dir"
out_base="$(basename "$out")"

mkdir -p "$out_dir"

if [ -e "$out" ]; then
  echo "Error: file already exists: $out"
  exit 2
fi

if [ ! -f "$tpl" ]; then
  echo "Error: missing template: $tpl"
  echo "Hint: expected at $proj_root/licence/header"
  exit 3
fi

date_str="$(date -Iseconds)"
user_str="${USER:-unknown}"

# Optionally choose comment prefix by extension (simple version)
ext="${out_base##*.}"
comment="%"
case "$ext" in
  m)   comment="%" ;;
  py)  comment="#" ;;
  sh)  comment="#" ;;
  c|h|cpp|hpp) comment="//" ;;
  *)   comment="#" ;;
esac
echo "$comment"

# # Read template and prefix lines with comment if template isn't already commented
# # (Simpler: assume template is plain text and we comment it here)
{
  while IFS= read -r line; do
    echo "${comment} ${line}"
  done < "$tpl"
# 
#   echo "${comment}"
#   echo "${comment} Created: ${date_str}"
#   echo "${comment} Author : ${user_str}"
#   echo "${comment} File   : ${out_base}"
#   echo
} > "$out"

echo "Created: $out"
