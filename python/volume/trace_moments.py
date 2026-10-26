#!/usr/bin/env python3
"""Sector completeness by exact spectral moments (no eigensolver).

For W-self-adjoint operators the trace of A^n is similarity invariant, so if
the four C2v sectors exactly decompose the zero-momentum space,
    Tr(H_{k=0}^n) = sum over (Sx, Sy) of Tr(H_{Sx,Sy}^n)     for every n.
Reference H_{k=0}: the PR #14 zero-momentum (++) translation sector of the 2x2
torus, constructed independently. Each operator is read from its exchange file
(weights, magnetic, upper electric entries), A = JE*E + JB*M built as exact
integer sparse rows (lower entries from w_i E_ij = w_j E_ji), and Tr(A^n)
computed exactly with Python integers for n = 1..NMAX.
usage: trace_moments.py REF_MATRIX JE JB NMAX SECTOR_MATRIX...      json=0
"""
import sys


def load(path, je, jb):
    t = open(path).read().split()
    n = int(t[0])
    w = [int(x) for x in t[1:1 + n]]
    mag = [int(x) for x in t[1 + n:1 + 2 * n]]
    rows = [dict() for _ in range(n)]
    r = t[1 + 2 * n:]
    for k in range(0, len(r), 3):
        i, j, z = int(r[k]), int(r[k + 1]), int(r[k + 2])
        rows[i][j] = rows[i].get(j, 0) + je * z
        if i != j:
            num = w[i] * z
            assert num % w[j] == 0, "non-integer lower entry"
            rows[j][i] = rows[j].get(i, 0) + je * (num // w[j])
    for i in range(n):
        rows[i][i] = rows[i].get(i, 0) + jb * mag[i]
    return rows


def traces(colmap, nmax):
    """Tr(A^k), k = 1..nmax: propagate each unit vector e_i through A k times
    (A v computed from the column map) and read component i. Exact integers."""
    n = len(colmap)
    tr = [0] * (nmax + 1)
    for i in range(n):
        v = {i: 1}
        for k in range(1, nmax + 1):
            u = {}
            for j, x in v.items():
                for r_, a in colmap[j].items():
                    u[r_] = u.get(r_, 0) + a * x
            v = u
            tr[k] += v.get(i, 0)
    return tr[1:]


def main():
    ref, je, jb, nmax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    secs = sys.argv[5:]
    res = {}
    for name in [ref] + secs:
        rows = load(name, je, jb)
        colmap = [dict() for _ in range(len(rows))]
        for r_, row in enumerate(rows):
            for j, a in row.items():
                colmap[j][r_] = a
        res[name] = traces(colmap, nmax)
    tot = [sum(res[s][k] for s in secs) for k in range(nmax)]
    ok = tot == res[ref]
    for k in range(nmax):
        print("TRACE|n=%d|ref=%d|sectors_sum=%d|equal=%d|json=0" % (k + 1, res[ref][k], tot[k], int(res[ref][k] == tot[k])))
    print("TRACE_MOMENTS|JE=%d|JB=%d|nmax=%d|sectors=%d|result=%s|json=0" % (je, jb, nmax, len(secs), "PASS" if ok else "FAIL"))
    if not ok:
        sys.exit(2)


if __name__ == "__main__":
    main()
