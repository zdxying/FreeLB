#!/usr/bin/env bash
#
# Build and run every casegen golden case.
#
# Steps:
#   1. python3 -m unittest discover -s tools/casegen/tests (the summary must
#      report at least one test)
#   2. for each tests/feature/casegen*: make clean -> make -> run its TARGET
#
# Discovery is automatic, so a new golden only needs its directory to match
# tests/feature/casegen*.  Every golden must keep `--check --selftest` in its
# regen rule; otherwise the case is built without verifying the spec/output
# pair, which is exactly what regenerating the Makefile can silently remove.
#
# A directory whose Makefile sets CXXC := nvcc is skipped unless nvcc and a
# visible GPU are present (or always with --no-cuda).  A directory whose
# Makefile sets CXXC := mpic++ or FLAGS += -DMPI_ENABLED is skipped unless
# both mpic++ and mpirun are present, and runs under `mpirun -np 2`.
#
# The generated binaries assert their own fixed checksums, so a non-zero exit
# is a numeric or structural regression.
#
# Usage:
#   tools/casegen/check.sh [--no-cuda]
#
#   --no-cuda   skip CUDA goldens even when nvcc and a GPU are present
#
# Exit status: 0 if nothing failed, 1 otherwise.
#
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$ROOT"

NO_CUDA=0
while [ $# -gt 0 ]; do
  case "$1" in
    --no-cuda) NO_CUDA=1 ;;
    -h|--help) sed -n '2,/^set -u$/p' "$0" | sed '$d'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

have() { command -v "$1" >/dev/null 2>&1; }

# OpenMPI refuses to run as root unless told to; other MPI launchers do not
# understand the flag, so only add it when actually running as root.
MPIRUN="mpirun"
if [ "$(id -u)" = "0" ]; then MPIRUN="mpirun --allow-run-as-root"; fi

UNITTEST_LOG="$(mktemp)"
BUILD_LOG="$(mktemp)"
RUN_LOG="$(mktemp)"
trap 'rm -f "$UNITTEST_LOG" "$BUILD_LOG" "$RUN_LOG"' EXIT

pass=0; fail=0; skip=0
FAILED=()
SKIPPED=()

clean_case() {
  make -C "$1" clean >/dev/null 2>&1 || true
}

run_target() {
  if have timeout; then
    timeout 300 "./$1"
  else
    "./$1"
  fi
}

run_mpi_target() {
  if have timeout; then
    timeout 300 $MPIRUN -np 2 "./$1"
  else
    $MPIRUN -np 2 "./$1"
  fi
}

echo "[1/2] generator unit tests"
( cd "$ROOT" && python3 -m unittest discover -s tools/casegen/tests ) >"$UNITTEST_LOG" 2>&1
unit_rc=$?
ran_line="$(grep -m1 -E '^Ran [0-9]+ tests?' "$UNITTEST_LOG" || true)"
ran_n="$(printf '%s\n' "$ran_line" | sed -n 's/^Ran \([0-9][0-9]*\) tests\?.*/\1/p')"

if [ "$unit_rc" -ne 0 ]; then
  echo "FAIL  tools/casegen/tests"
  tail -n 15 "$UNITTEST_LOG"
  fail=$((fail+1)); FAILED+=("tools/casegen/tests")
elif [ -z "$ran_n" ]; then
  echo "FAIL  tools/casegen/tests (no 'Ran N tests' summary)"
  tail -n 15 "$UNITTEST_LOG"
  fail=$((fail+1)); FAILED+=("tools/casegen/tests:no-summary")
elif [ "$ran_n" -le 0 ]; then
  echo "FAIL  tools/casegen/tests (ran 0 tests)"
  fail=$((fail+1)); FAILED+=("tools/casegen/tests:empty")
else
  echo "PASS  tools/casegen/tests ($ran_line)"
  pass=$((pass+1))
fi

if [ "$fail" -ne 0 ]; then
  echo
  echo "-----------------------------------------"
  echo "PASS: $pass   FAIL: $fail   SKIP: $skip"
  printf '  failed: %s\n' "${FAILED[@]}"
  exit 1
fi

mapfile -t CASE_DIRS < <(find tests/feature -mindepth 1 -maxdepth 1 -type d -name 'casegen*' \
                         -printf '%p\n' | sort)

echo "[2/2] golden cases"
if [ "${#CASE_DIRS[@]}" -eq 0 ]; then
  echo "FAIL  tests/feature/casegen* (no golden directories found)"
  fail=$((fail+1)); FAILED+=("tests/feature/casegen*:empty")
fi

for rel in "${CASE_DIRS[@]}"; do
  mf="$rel/Makefile"

  if [ ! -f "$mf" ]; then
    echo "FAIL  $rel (no Makefile)"
    fail=$((fail+1)); FAILED+=("$rel:no-makefile")
    continue
  fi

  # A regenerated Makefile that drops --check/--selftest would still build,
  # so reject the drift before deciding whether this platform can run it.
  if ! grep -qE -- '--check[[:space:]]+--selftest' "$mf"; then
    echo "FAIL  $rel (regen rule lacks --check --selftest)"
    fail=$((fail+1)); FAILED+=("$rel:no-selftest")
    continue
  fi

  is_cuda=0
  grep -qE '^CXXC[[:space:]]*:?=[[:space:]]*nvcc' "$mf" && is_cuda=1
  is_mpi=0
  grep -qE '^CXXC[[:space:]]*:?=[[:space:]]*mpic\+\+' "$mf" && is_mpi=1
  grep -qE '^[[:space:]]*[^#].*-DMPI_ENABLED' "$mf" && is_mpi=1

  if [ "$is_cuda" = 1 ]; then
    reason=""
    if [ "$NO_CUDA" = 1 ]; then
      reason="--no-cuda"
    elif ! have nvcc; then
      reason="nvcc not found"
    elif ! have nvidia-smi; then
      reason="nvidia-smi not found"
    elif ! nvidia-smi -L >/dev/null 2>&1; then
      reason="no visible GPU device"
    fi
    if [ -n "$reason" ]; then
      echo "SKIP  $rel (cuda: $reason)"
      skip=$((skip+1)); SKIPPED+=("$rel: $reason")
      continue
    fi
  fi

  if [ "$is_mpi" = 1 ]; then
    reason=""
    if ! have mpic++; then
      reason="mpic++ not found"
    elif ! have mpirun; then
      reason="mpirun not found"
    fi
    if [ -n "$reason" ]; then
      echo "SKIP  $rel (mpi: $reason)"
      skip=$((skip+1)); SKIPPED+=("$rel: $reason")
      continue
    fi
  fi

  target="$(sed -n 's/^TARGET[[:space:]]*=[[:space:]]*//p' "$mf" | head -1 | tr -d '[:space:]')"
  if [ -z "$target" ]; then
    echo "FAIL  $rel (no TARGET in Makefile)"
    fail=$((fail+1)); FAILED+=("$rel:no-target")
    continue
  fi

  echo "=== $rel"
  clean_case "$rel"
  if ! make -C "$rel" >"$BUILD_LOG" 2>&1; then
    echo "FAIL  $rel (build)"
    tail -n 8 "$BUILD_LOG"
    fail=$((fail+1)); FAILED+=("$rel:build")
    clean_case "$rel"
    continue
  fi

  rc=0
  if [ "$is_mpi" = 1 ]; then
    ( cd "$rel" && run_mpi_target "$target" ) >"$RUN_LOG" 2>&1 || rc=$?
  else
    ( cd "$rel" && run_target "$target" ) >"$RUN_LOG" 2>&1 || rc=$?
  fi
  clean_case "$rel"
  if [ "$rc" -ne 0 ]; then
    echo "FAIL  $rel (run rc=$rc)"
    tail -n 10 "$RUN_LOG"
    fail=$((fail+1)); FAILED+=("$rel:run($rc)")
    continue
  fi
  echo "PASS  $rel"
  pass=$((pass+1))
done

echo
echo "-----------------------------------------"
echo "PASS: $pass   FAIL: $fail   SKIP: $skip"
if [ "${#SKIPPED[@]}" -ne 0 ]; then
  printf '  skipped: %s\n' "${SKIPPED[@]}"
fi
if [ "$fail" -ne 0 ]; then
  printf '  failed: %s\n' "${FAILED[@]}"
  exit 1
fi
