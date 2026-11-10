#!/usr/bin/env bash
# Block-size sweep for the plain-communication benchmark.
#
#   ./sweep.sh [2d|3d|all]        default: all
#   MPIRUN=mpiexec ./sweep.sh 3d
#
# For every block size the benchmark is run with 1 rank (shared-memory halos
# only, so the MPI path is idle) and with 2 ranks (adds the MPI exchange), and
# the better of two alternating runs is reported.
#
# Columns: block edge, blocks per rank, communicators per step, coalesced
# messages per step, exchange time with 2 ranks (us/step), MLUPS for 1 and 2
# ranks.  "comms" is what the old per-communicator implementation posted as
# messages; "msgs" is what the plan posts now.

set -u
cd "$(dirname "$0")"

DIM=${1:-all}
MPIRUN=${MPIRUN:-mpirun}
MAKE=${MAKE:-make}

N2D=${N2D:-240}
BPS2D=${BPS2D:-"4 5 6 8 10 12 15 20"}
STEPS2D=${STEPS2D:-400}
N3D=${N3D:-80}
BPS3D=${BPS3D:-"2 4 5 8 10"}
STEPS3D=${STEPS3D:-50}

$MAKE -s >/dev/null || { echo "build failed" >&2; exit 1; }

# one benchmark run -> "MLUPS comm_us comms msgs"
run_case() {  # exe N bps steps np
  local out mlups comm comms msgs
  out=$("$MPIRUN" -np "$5" "$1" "$2" "$3" "$4" 2>/dev/null)
  mlups=$(printf '%s\n' "$out" | sed -n 's/^\[result\] partial .*MLUPS= *\([0-9.]*\).*/\1/p')
  comm=$(printf '%s\n' "$out" | sed -n 's/^\[result\] partial .*comm= *\([0-9.]*\) us.*/\1/p')
  comms=$(printf '%s\n' "$out" | sed -n 's/^\[layout\] partial .*communicators\/step=\([0-9]*\).*/\1/p')
  msgs=$(printf '%s\n' "$out" | sed -n 's/^\[layout\] partial .*messages\/step=\([0-9]*\).*/\1/p')
  if [ -z "$mlups" ]; then
    echo "run failed: $1 N=$2 bps=$3 steps=$4 np=$5" >&2
    echo "0 0 0 0"
    return
  fi
  echo "${mlups:-0} ${comm:-0} ${comms:-0} ${msgs:-0}"
}

sweep() {  # dim exe N steps "bps..."
  local dim=$1 exe=$2 n=$3 steps=$4 bps=$5
  local b blocks m1 m2 c2 k2 g2 best1 best2 bc2 bk2 bg2
  echo
  echo "== $dim: domain ${n}^${n}, ${steps} steps, best of 2 alternating runs"
  printf '%-4s %-7s %10s %6s %10s %8s %8s\n' \
    "B" "blocks" "comms/step" "msgs" "comm us" "MLUPS-1" "MLUPS-2"
  for bps in $bps; do
    b=$((n / bps))
    blocks=$((bps * bps))
    [ "$dim" = 3d ] && blocks=$((bps * bps * bps))
    best1=0
    best2=0
    bc2=0
    bk2=0
    bg2=0
    for rep in 0 1; do
      if [ $((rep % 2)) -eq 0 ]; then
        read -r m1 _ _ _ <<<"$(run_case "$exe" "$n" "$bps" "$steps" 1)"
        read -r m2 c2 k2 g2 <<<"$(run_case "$exe" "$n" "$bps" "$steps" 2)"
      else
        read -r m2 c2 k2 g2 <<<"$(run_case "$exe" "$n" "$bps" "$steps" 2)"
        read -r m1 _ _ _ <<<"$(run_case "$exe" "$n" "$bps" "$steps" 1)"
      fi
      if [ "$(awk -v a="$best1" -v b="$m1" 'BEGIN { print (b>a)?1:0 }')" -eq 1 ]; then
        best1=$m1
      fi
      if [ "$(awk -v a="$best2" -v b="$m2" 'BEGIN { print (b>a)?1:0 }')" -eq 1 ]; then
        best2=$m2
        bc2=$c2
        bk2=$k2
        bg2=$g2
      fi
    done
    printf '%-4s %-7s %10s %6s %10s %8s %8s\n' \
      "$b" "$blocks" "$bk2" "$bg2" "$bc2" "$best1" "$best2"
  done
}

if [ "$DIM" = all ] || [ "$DIM" = 2d ]; then
  sweep 2d ./2d/plain_comm_bench2d.exe "$N2D" "$STEPS2D" "$BPS2D"
fi
if [ "$DIM" = all ] || [ "$DIM" = 3d ]; then
  sweep 3d ./3d/plain_comm_bench3d.exe "$N3D" "$STEPS3D" "$BPS3D"
fi
