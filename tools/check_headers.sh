#!/usr/bin/env bash
#
# Per-header self-containment check.
#
# Compiles every tracked src/ header on its own, once with g++ and once with
# mpic++ -DMPI_ENABLED, and reports the ones that fail.  A header that only
# compiles because another header happened to include a dependency first is an
# accident waiting to break a downstream TU.
#
# Usage:
#   tools/check_headers.sh          both sweeps
#   tools/check_headers.sh --host   g++ only
#   tools/check_headers.sh --mpi    mpic++ -DMPI_ENABLED only
#
# Exit 0 when nothing fails beyond BASELINE below, 1 when something new does,
# 2 on a usage error.  The gate is "do not add a failure", not "the baseline is
# clean".  When you fix a baseline entry the script names it; delete that line
# by hand.  Never add a line to BASELINE to silence a new failure.
#
# Why the baseline entries fail (fix them when you touch the file):
#   utils/cuda_device.h      includes cuda.h and CUDA builtins; nvcc only
#
# Not covered here (the compiler cannot see it): headers missing an include
# guard -- io/checkpoint.hh, geometry/geometry.hh, geometry/basic_geometry.hh,
# parallel/mpi_manager.hh.
#
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BASELINE='
host:src/utils/cuda_device.h
mpi:src/utils/cuda_device.h
'

WANT_HOST=1
WANT_MPI=1
for arg in "$@"; do
  case "$arg" in
    --host)    WANT_MPI=0 ;;
    --mpi)     WANT_HOST=0 ;;
    -h|--help) sed -n '2,25p' "$0"; exit 0 ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

# Include untracked headers too: a new file must be checked before it is
# committed, and `git ls-files` alone would only see it after `git add`.
mapfile -t HEADERS < <(git ls-files --cached --others --exclude-standard \
                         'src/**/*.h' 'src/*.h' | sort -u)
[ "${#HEADERS[@]}" -gt 0 ] || { echo "no headers found (run from a git checkout)"; exit 2; }

found=()

sweep() {
  local label="$1"; shift
  local bad=0
  for f in "${HEADERS[@]}"; do
    if ! "$@" "$f" 2>/dev/null; then
      found+=("$label:$f")
      bad=$((bad + 1))
    fi
  done
  echo "[$label] ${#HEADERS[@]} headers, $bad failed"
}

if [ "$WANT_HOST" = 1 ]; then
  sweep host g++ -std=c++17 -fsyntax-only -Isrc
fi
if [ "$WANT_MPI" = 1 ]; then
  if command -v mpic++ >/dev/null 2>&1; then
    sweep mpi mpic++ -std=c++17 -DMPI_ENABLED -fsyntax-only -Isrc
  else
    echo "SKIP [mpi] mpic++ not found"
  fi
fi

# Trimming matters: a stray space would make a baseline line unmatchable.
found=($(printf '%s\n' "${found[@]:-}" | sed 's/^[[:space:]]*//; /^$/d' | sort))
base=($(printf '%s\n' "$BASELINE" | sed 's/^[[:space:]]*//; /^$/d' | sort))

new=()
for e in "${found[@]:-}"; do
  [ -n "$e" ] || continue
  printf '%s\n' "${base[@]:-}" | grep -qxF -- "$e" || new+=("$e")
done
fixed=()
for e in "${base[@]:-}"; do
  [ -n "$e" ] || continue
  # Only entries from a sweep that actually ran can have been fixed.
  case "$e" in
    host:*) [ "$WANT_HOST" = 1 ] || continue ;;
    mpi:*)  [ "$WANT_MPI" = 1 ] || continue ;;
  esac
  printf '%s\n' "${found[@]:-}" | grep -qxF -- "$e" || fixed+=("$e")
done

if [ "${#fixed[@]}" -gt 0 ]; then
  echo "note: these now pass -- delete them from BASELINE: ${fixed[*]}"
fi
if [ "${#new[@]}" -gt 0 ]; then
  echo "NEW FAILURES (not in the baseline):"
  printf '  %s\n' "${new[@]}"
  exit 1
fi
exit 0
