#!/usr/bin/env python3
"""Per-function SASS instruction census for CUDA binaries.

The quick way to verify that an optimization actually changed the generated
code (e.g. register-resident kernels: baseline LDG~438 vs RegCell LDG~78)
and that nothing spilled to local memory.

Usage: python3 tools/sass_census.py bin1 [bin2 ...]

Requires cuobjdump on PATH.  Prints, per function with >= 100 instructions:
instruction total plus the hot opcodes (LDG/STG/IMAD/FFMA/FADD/FMUL/LDS/STS).
"""
import re
import subprocess
import sys
from collections import Counter

HOT = ["LDG", "STG", "IMAD", "FFMA", "FADD", "FMUL", "LDS", "STS"]


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
    for path in sys.argv[1:]:
        print(f"===== {path} =====")
        for name, c in census(path).items():
            if c["total"] < 100:
                continue
            dem = re.sub(r"<.*", "", demangle(name))[:64]
            hot = " ".join(f"{k}={c[k]:4d}" for k in HOT if c[k])
            print(f"{dem:64s} inst={c['total']:5d} {hot}")


if __name__ == "__main__":
    main()
