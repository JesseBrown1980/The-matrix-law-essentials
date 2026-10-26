#!/usr/bin/env python3
"""Validate the zero-momentum C2v sectors (ymsec) against complete operators.

2x2: the union of the four sector spectra must equal the PR #14 zero-momentum
     (++) sector spectrum level for level (lowest 12, exact micro-units), and
     the four flat counts must sum to 8.
2x3: every zero-momentum sector level up to the top of the full operator's
     lowest 12 must be a level of the full 47,449-state operator, with
     multiplicity (|diff| <= 2 micro-units, the residual radius); the ground
     state must agree exactly and the flat counts must sum to 8. Full levels
     absent from the union are nonzero-momentum states.
Integer-only. usage: ymsec_validate.py 2x2|2x3 R        json=0
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = os.path.join(ROOT, "target", "release")
BUILD = os.path.join(ROOT, "build", "ymsec")
SECTORS = (("1", "1"), ("1", "-1"), ("-1", "1"), ("-1", "-1"))


def fields(line):
    return dict(kv.split("=", 1) for kv in line.strip().split("|")[1:] if "=" in kv)


def levels(matrix, r, want):
    out = subprocess.run([os.path.join(B, "cluster"), matrix, "1", str(r), str(want), "4"], capture_output=True, text=True).stdout
    lv, conv = [], None
    for line in out.splitlines():
        f = fields(line)
        if line.startswith("CLUSTER_LEVEL"):
            lv.append(int(f["theta"].split("/")[0]))
        elif line.startswith("CLUSTER|"):
            conv = f["converged"]
    if conv != "1":
        print("YMSEC_VALIDATE|error=not_converged|matrix=%s|json=0" % os.path.basename(matrix))
        sys.exit(2)
    return lv


def main():
    lat, r = sys.argv[1], int(sys.argv[2])
    lx, ly = lat.split("x")
    os.makedirs(BUILD, exist_ok=True)
    union, flats = [], 0
    for sx, sy in SECTORS:
        m = os.path.join(BUILD, "v%s_%s_%s.matrix" % (lat, sx, sy))
        out = subprocess.run([os.path.join(B, "ymsec"), lx, ly, sx, sy, m], capture_output=True, text=True).stdout
        f = fields(out)
        if f.get("w_self_adjoint") != "PASS":
            print("YMSEC_VALIDATE|error=not_self_adjoint|sector=%s,%s|json=0" % (sx, sy))
            sys.exit(2)
        flats += int(f["flat_states"])
        union += levels(m, r, 12)
    union.sort()
    if lat == "2x2":
        ref = levels(os.path.join(ROOT, "data", "matrices", "sector_pp.matrix"), r, 12)
        diffs = [u - v for u, v in zip(union[:12], ref)]
        ok = flats == 8 and all(d == 0 for d in diffs)
        print("YMSEC_VALIDATE|L=2x2|r=%d|flat_total=%d|levels_compared=12|max_abs_diff=%d|reference=sector_pp(PR14)|result=%s|json=0"
              % (r, flats, max(abs(d) for d in diffs), "PASS" if ok else "FAIL"))
    else:
        full_m = os.path.join(ROOT, "build", "torus2x3.matrix")
        if not os.path.exists(full_m):
            subprocess.run([os.path.join(B, "torusl"), "2", "3", full_m], capture_output=True, text=True)
        full = levels(full_m, r, 12)
        # Required direction: every zero-momentum level up to the full list's top
        # is a level of the full operator, WITH MULTIPLICITY (multiset inclusion,
        # |diff| <= 2 micro-units = residual radius). Full levels absent from the
        # union are nonzero-momentum states.
        top = full[-1]
        window = [u for u in union if u <= top]
        pool = list(full)
        matched = 0
        for u in window:
            hit = next((i for i, f in enumerate(pool) if abs(f - u) <= 2), None)
            if hit is not None:
                pool.pop(hit)
                matched += 1
        ok = flats == 8 and union[0] == full[0] and matched == len(window)
        print("YMSEC_VALIDATE|L=2x3|r=%d|flat_total=%d|ground_equal=%d|union_levels_in_window=%d|matched_with_multiplicity=%d|full_levels_not_in_union=%d(nonzero_momentum)|result=%s|json=0"
              % (r, flats, int(union[0] == full[0]), len(window), matched, len(pool), "PASS" if ok else "FAIL"))
    if not ok:
        sys.exit(2)


if __name__ == "__main__":
    main()
