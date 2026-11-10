#!/usr/bin/env bash
#
# FreeLB test runner.
#
# Builds and runs every tests/<category>/<name> target and checks its exit
# status.  A target that fails to build or returns non-zero counts as FAIL.
#
# Usage:
#   tests/run_tests.sh [--no-cuda] [--no-mpi]
#
#   --no-cuda   skip targets that need nvcc / a GPU
#   --no-mpi    skip targets that need mpic++ (block/refblock, smoke mpi build)
#
# Exit status: 0 if no target failed, 1 otherwise.
#
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

NO_CUDA=0
NO_MPI=0
while [ $# -gt 0 ]; do
  case "$1" in
    --no-cuda) NO_CUDA=1 ;;
    --no-mpi)  NO_MPI=1 ;;
    -h|--help) sed -n '2,16p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

have() { command -v "$1" >/dev/null 2>&1; }

# OpenMPI refuses to run as root unless told to; other MPI launchers do not
# understand the flag, so only add it when actually running as root.
MPIRUN="mpirun"
if [ "$(id -u)" = "0" ]; then MPIRUN="mpirun --allow-run-as-root"; fi

BUILD_LOG="$(mktemp)"
RUN_LOG="$(mktemp)"
trap 'rm -f "$BUILD_LOG" "$RUN_LOG"' EXIT

# Discover test directories (tests/<category>/<name>/Makefile).
mapfile -t DIRS < <(find . -mindepth 3 -maxdepth 3 -name Makefile -printf '%h\n' \
                      | sed 's|^\./||' | sort)

pass=0; fail=0; skip=0
FAILED=()

run_one() {
  local dir="$1"
  local target is_cuda is_mpi
  target="$(sed -n 's/^TARGET[[:space:]]*=[[:space:]]*//p' "$dir/Makefile" | head -1 | tr -d '[:space:]')"
  is_cuda=0; grep -qE '^CXXC[[:space:]]*:?=[[:space:]]*nvcc' "$dir/Makefile" && is_cuda=1
  is_mpi=0;  grep -qE '^CXXC[[:space:]]*:?=[[:space:]]*mpic\+\+' "$dir/Makefile" && is_mpi=1

  if [ "$is_cuda" = 1 ] && { [ "$NO_CUDA" = 1 ] || ! have nvcc; }; then
    echo "SKIP  $dir (cuda)"; skip=$((skip+1)); return
  fi
  if [ "$is_mpi" = 1 ] && { [ "$NO_MPI" = 1 ] || ! have mpic++; }; then
    echo "SKIP  $dir (mpi)"; skip=$((skip+1)); return
  fi
  if [ -z "$target" ]; then
    echo "FAIL  $dir (no TARGET in Makefile)"; fail=$((fail+1)); FAILED+=("$dir:no-target"); return
  fi

  echo "=== $dir"
  make -C "$dir" clean >/dev/null 2>&1 || true
  if ! make -C "$dir" >"$BUILD_LOG" 2>&1; then
    echo "FAIL  $dir (build)"; tail -n 8 "$BUILD_LOG"
    fail=$((fail+1)); FAILED+=("$dir:build"); return
  fi

  # Compile-only MPI smoke target lives in the same directory.
  if [ "$dir" = "feature/instantiation_smoke" ] && have mpic++ && [ "$NO_MPI" = 0 ]; then
    if ! make -C "$dir" mpi >"$BUILD_LOG" 2>&1; then
      echo "FAIL  $dir (mpi build)"; tail -n 8 "$BUILD_LOG"
      fail=$((fail+1)); FAILED+=("$dir:mpi-build"); return
    fi
  fi

  local rc=0
  if [ "$dir" = "block/refblock" ] || [ "$dir" = "block/compact_mpi" ]; then
    ( cd "$dir" && timeout 300 $MPIRUN -np 2 ./"$target" ) >"$RUN_LOG" 2>&1 || rc=$?
  elif [ "$dir" = "feature/hidden_comm" ] || [ "$dir" = "block/amr_comm" ]; then
    # These spawn a two-rank child with MPI_Comm_spawn, which needs a real MPI
    # job as the parent; a singleton process is not enough everywhere.
    ( cd "$dir" && timeout 300 $MPIRUN -np 1 ./"$target" ) >"$RUN_LOG" 2>&1 || rc=$?
  else
    ( cd "$dir" && timeout 300 ./"$target" ) >"$RUN_LOG" 2>&1 || rc=$?
  fi
  if [ "$rc" -ne 0 ]; then
    echo "FAIL  $dir (run rc=$rc)"; tail -n 10 "$RUN_LOG"
    fail=$((fail+1)); FAILED+=("$dir:run($rc)"); return
  fi
  echo "PASS  $dir"
  pass=$((pass+1))
}

for d in "${DIRS[@]}"; do
  run_one "$d"
done

echo
echo "-----------------------------------------"
echo "PASS: $pass   FAIL: $fail   SKIP: $skip"
if [ "$fail" -ne 0 ]; then
  printf '  failed: %s\n' "${FAILED[@]}"
  exit 1
fi
