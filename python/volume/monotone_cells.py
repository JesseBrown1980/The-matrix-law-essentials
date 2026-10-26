#!/usr/bin/env python3
"""Continuous-in-coupling band-gap certificate for the full 3x3 zero-momentum
operator from point certificates (the 2x2 ladder's monotone-cell jump).

H/JE = E + r*M with r = JB/JE. E and M are W-self-adjoint and M is diagonal
with entries 2 - chi(U_p) summed over plaquettes, all >= 0. By Courant-Fischer
in the W inner product every ordered eigenvalue E_k(r) is non-decreasing in r.
So on a cell [ra, rb], for every r in it,
    E9(r) - E8(r) >= E9(ra) - E8(rb) >= tB(ra) - t1(rb),
where tB(ra) <= E9(ra) and E8(rb) < t1(rb) come from the SHELL_GAP_CERT lines
(N(t1) = N(tB) = 8). Thresholds are micro-units of the JE = 1 energy.
A cell certifies when that bound is positive; the run certifies the union of
consecutive cells. Integer-only.
usage: monotone_cells.py CERT_FILE...        json=0
"""
import re
import sys


def main():
    pts = {}
    for path in sys.argv[1:]:
        for line in open(path):
            if not line.startswith("SHELL_GAP_CERT|") or "|result=PASS|" not in line:
                continue
            r = int(re.search(r"\|r=(\d+)\|", line).group(1))
            t1 = int(re.search(r"\|t1=(\d+)/1000000\|", line).group(1))
            tb = int(re.search(r"\|tB=(\d+)/1000000\|", line).group(1))
            pts[r] = (t1, tb)
    rs = sorted(pts)
    if len(rs) < 2:
        print("MONOTONE_CELLS|result=FAIL|reason=need_two_points|json=0")
        sys.exit(2)
    worst = None
    ok = True
    for a, b in zip(rs, rs[1:]):
        bound = pts[a][1] - pts[b][0]
        ok = ok and bound > 0
        worst = bound if worst is None else min(worst, bound)
        print("MONOTONE_CELL|L=3x3|r_in=[%d,%d]|tB(ra)=%d/1000000|t1(rb)=%d/1000000|gap_E9_minus_E8_ge=%d/1000000|certified=%d|json=0"
              % (a, b, pts[a][1], pts[b][0], bound, int(bound > 0)))
    print("MONOTONE_CELLS|L=3x3|points=%s|r_range=[%d,%d]|cells=%d|min_gap_lower_bound=%d/1000000|result=%s|json=0"
          % (",".join(map(str, rs)), rs[0], rs[-1], len(rs) - 1, worst, "PASS" if ok else "FAIL"))
    if not ok:
        sys.exit(2)


if __name__ == "__main__":
    main()
