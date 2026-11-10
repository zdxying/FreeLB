#!/usr/bin/env python3
"""Per-function SASS instruction census for CUDA binaries.

The quick way to verify that an optimization actually changed the generated
code (e.g. register-resident kernels: baseline LDG~438 vs RegCell LDG~78)
and that nothing spilled to local memory.

Usage: python3 tools/sass_census.py bin1 [bin2 ...]

Requires cuobjdump on PATH.  Prints, per function with >= 100 instructions:
instruction total plus the hot opcodes (LDG/STG/IMAD/FFMA/FADD/FMUL/LDS/STS).

Exits non-zero if a cell-dynamics kernel is smaller than
MIN_CELL_DYNAMICS_INSTRUCTIONS -- see the note on that constant for why.
"""
import re
import subprocess
import sys
from collections import Counter

HOT = ["LDG", "STG", "IMAD", "FFMA", "FADD", "FMUL", "LDS", "STS"]

CELL_DYNAMICS_KERNEL = "CuDevApplyCellDynamicsKernel"

# Lower bound on a real cell-dynamics kernel, used to catch a kernel that the
# optimiser deleted because it provably has no observable effect.
#
# That failure mode compiles clean, runs clean, and even produces plausible
# checksums, because a kernel that writes nothing leaves the field untouched.
# It was hit for real twice in this project:
#   * PopCache's store took the storage address by value into a local
#     (`auto r = raw(d); r = v[d];` on a T&), so the store landed in a discarded
#     temporary and the whole kernel collapsed to 16 instructions;
#   * a CSE specialisation that stopped matching the cell type would silently
#     fall back to the loop-based primary template.
# Neither shows up in a warning, a checksum comparison against a stale binary,
# or a throughput number that still looks plausible.  Only the instruction count
# betrays it.
#
# The smallest legitimate kernel measured so far is 432 (D2Q9, the rho/U-only
# task); the smallest broken one was 16.  300 leaves margin on both sides.
MIN_CELL_DYNAMICS_INSTRUCTIONS = 300


def census(path):
    sass = subprocess.run(["cuobjdump", "-sass", path],
                          capture_output=True, text=True).stdout
    funcs = {}
    cur = None
    for line in sass.splitlines():
        m = re.search(r"Function : (\S+)", line)
        if m:
            cur = m.group(1)
            funcs.setdefault(cur, Counter())
            continue
        if cur is None:
            continue
        im = re.match(r"\s*/\*[0-9a-f]+\*/\s+(.*?);", line)
        if not im:
            continue
        tokens = im.group(1).split()
        if not tokens:
            continue
        op = tokens[0].split(".")[0]
        c = funcs[cur]
        c["total"] += 1
        c[op] += 1
    return funcs


def demangle(name):
    try:
        return subprocess.run(["c++filt", name],
                              capture_output=True, text=True).stdout.strip()
    except OSError:
        return name


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    violations = 0
    for path in sys.argv[1:]:
        print(f"===== {path} =====")
        found_cell_dynamics = False
        for name, c in census(path).items():
            if CELL_DYNAMICS_KERNEL in demangle(name):
                found_cell_dynamics = True
                if c["total"] < MIN_CELL_DYNAMICS_INSTRUCTIONS:
                    print(f"  FAIL  cell-dynamics kernel has only {c['total']} "
                          f"instructions (< {MIN_CELL_DYNAMICS_INSTRUCTIONS}): "
                          f"it was almost certainly optimised away")
                    print(f"        {name}")
                    violations += 1
                    continue
                policy = "RegPop" if "RegPop" in demangle(name) else "DirectPop"
                print(f"  ok    {policy:9s} cell dynamics inst={c['total']}")
            if c["total"] < 100:
                continue
            dem = re.sub(r"<.*", "", demangle(name))[:64]
            hot = " ".join(f"{k}={c[k]:4d}" for k in HOT if c[k])
            print(f"{dem:64s} inst={c['total']:5d} {hot}")
        if not found_cell_dynamics:
            print(f"  note: no {CELL_DYNAMICS_KERNEL} found -- is this a GPU "
                  f"binary built from FreeLB?")
    if violations:
        print(f"\n{violations} cell-dynamics kernel(s) below the instruction "
              f"floor", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
