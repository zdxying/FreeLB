#!/usr/bin/env bash
# Block-size sweep for the hidden-communication benchmark.
#
#   ./sweep.sh [2d|3d|all]        default: all
#   MPIRUN=mpiexec ./sweep.sh 3d
#
# For every block size the benchmark is run with 1 rank (collision plus the
# shared-memory halo only) and with 2 ranks (adds the MPI exchange).  Each case
# is run twice with the order alternated, and the better MLUPS is reported, so
# a thermal dip on the laptop cannot decide the winner.
#
# Columns: block edge, blocks, halo ring as a share of the lattice, MLUPS for
# 1 and 2 ranks, and -- for 2 ranks -- the communication share of the step
# (send+recv+waitall+set), the phase-1 collision share, and the requests/step.

set -u
cd "$(dirname "$0")"

DIM=${1:-all}
MPIRUN=${MPIRUN:-mpirun}
MAKE=${MAKE:-make}

# fixed domain, swept block count.  N must be divisible by every value below.
# Override with N2D=/N3D=/STEPS2D=/STEPS3D= to sweep another configuration.
N2D=${N2D:-240}
BPS2D=${BPS2D:-"4 5 6 8 10 12 15 20"}
STEPS2D=${STEPS2D:-400}
N3D=${N3D:-80}
BPS3D=${BPS3D:-"2 4 5 8 10"}
STEPS3D=${STEPS3D:-50}

$MAKE -s >/dev/null || { echo "build failed" >&2; exit 1; }

# one benchmark run -> "MLUPS requests comm% first%"
run_case() {  # exe N bps steps np
  local out mlups reqs line send recv wl setp first
  out=$("$MPIRUN" -np "$5" "$1" "$2" "$3" "$4" 2>/dev/null)
  mlups=$(printf '%s\n' "$out" | sed -n 's/.*MLUPS=\([0-9.]*\).*/\1/p')
  reqs=$(printf '%s\n' "$out" | sed -n 's/.*requests\/step=\([0-9.]*\).*/\1/p')
  line=$(printf '%s\n' "$out" | grep '^\[profile\] collideFirst=')
  send=$(printf '%s\n' "$line" | sed -n 's/.*send= *\([0-9.]*\)%.*/\1/p')
  recv=$(printf '%s\n' "$line" | sed -n 's/.*recv= *\([0-9.]*\)%.*/\1/p')
  first=$(printf '%s\n' "$line" | sed -n 's/.*collideFirst= *\([0-9.]*\)%.*/\1/p')
  wl=$(printf '%s\n' "$out" | sed -n 's/^\[profile\] waitall= *\([0-9.]*\)%.*/\1/p')
  setp=$(printf '%s\n' "$out" | sed -n 's/^\[profile\] waitall=.*set= *\([0-9.]*\)%.*/\1/p')
  if [ -z "$mlups" ]; then
    echo "run failed: $1 N=$2 bps=$3 steps=$4 np=$5" >&2
    echo "0 0 0 0"
    return
  fi
  awk -v m="$mlups" -v r="$reqs" -v a="${send:-0}" -v b="${recv:-0}" -v c="${wl:-0}" \
      -v d="${setp:-0}" -v f="${first:-0}" \
      'BEGIN { printf "%.2f %.0f %.1f %.1f\n", m, r, a+b+c+d, f }'
}

sweep() {  # dim exe N steps "bps..."
  local dim=$1 exe=$2 n=$3 steps=$4 bps=$5
  local b blocks ring m1 m2 r1 r2 c2 f2 best1 best2 br2 bc2 bf2
  echo
  echo "== $dim: domain ${n}^${n}, ${steps} steps, best of 2 alternating runs"
  printf '%-4s %-7s %-7s %8s %8s %7s %8s %6s\n' \
    "B" "blocks" "halo%" "MLUPS-1" "MLUPS-2" "comm%" "collide1%" "reqs"
  for bps in $bps; do
    b=$((n / bps))
    blocks=$((bps * bps))
    [ "$dim" = 3d ] && blocks=$((bps * bps * bps))
    if [ "$dim" = 3d ]; then
      ring=$(awk -v b="$b" 'BEGIN { t=(b+2)^3; printf "%.1f", 100*(t-b^3)/t }')
    else
      ring=$(awk -v b="$b" 'BEGIN { t=(b+2)^2; printf "%.1f", 100*(t-b^2)/t }')
    fi
    # Two repetitions, order alternated.  The profile columns come from the
    # faster of the two 2-rank runs, so they belong to the reported number.
    best1=0
    best2=0
    br2=0
    bc2=0
    bf2=0
    for rep in 0 1; do
      if [ $((rep % 2)) -eq 0 ]; then
        read -r m1 r1 _ _ <<<"$(run_case "$exe" "$n" "$bps" "$steps" 1)"
        read -r m2 r2 c2 f2 <<<"$(run_case "$exe" "$n" "$bps" "$steps" 2)"
      else
        read -r m2 r2 c2 f2 <<<"$(run_case "$exe" "$n" "$bps" "$steps" 2)"
        read -r m1 r1 _ _ <<<"$(run_case "$exe" "$n" "$bps" "$steps" 1)"
      fi
      if [ "$(awk -v a="$best1" -v b="$m1" 'BEGIN { print (b>a)?1:0 }')" -eq 1 ]; then
        best1=$m1
      fi
      if [ "$(awk -v a="$best2" -v b="$m2" 'BEGIN { print (b>a)?1:0 }')" -eq 1 ]; then
        best2=$m2
        br2=$r2
        bc2=$c2
        bf2=$f2
      fi
    done
    printf '%-4s %-7s %-7s %8s %8s %7s %8s %6s\n' \
      "$b" "$blocks" "$ring" "$best1" "$best2" "$bc2" "$bf2" "$br2"
  done
}

if [ "$DIM" = all ] || [ "$DIM" = 2d ]; then
  sweep 2d ./2d/hidden_comm_bench2d.exe "$N2D" "$STEPS2D" "$BPS2D"
fi
if [ "$DIM" = all ] || [ "$DIM" = 3d ]; then
  sweep 3d ./3d/hidden_comm_bench3d.exe "$N3D" "$STEPS3D" "$BPS3D"
fi
