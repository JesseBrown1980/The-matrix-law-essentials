#!/usr/bin/env python3
"""Integer-only vacuum-ball block-Krylov Ritz solver for finite S3 lattice gauge.

Takes any standard exchange .matrix file (from ball or torusl), couplings JE, JB,
and Krylov order. Uses the 8 flat vacua (mag==0) as the simultaneous 8-dimensional
block seed, avoiding single-vector Chebyshev clustering issues.

Computes exact micro-unit brackets for all 8 vacuum eigenvalues via Bareiss
fraction-free elimination.
Zero floats, Decimal, Fraction, true division, or JSON.
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path
from time import monotonic_ns


def load_matrix(path):
    with open(path, "r", encoding="ascii") as f:
        toks = f.read().split()
    n = int(toks[0])
    w = [int(x) for x in toks[1:1 + n]]
    mag = [int(x) for x in toks[1 + n:1 + 2 * n]]
    rest = [int(x) for x in toks[1 + 2 * n:]]
    assert len(rest) % 3 == 0
    return n, w, mag, rest


def build_system(n, w, mag, rest, je, jb):
    col = [[] for _ in range(n)]
    for k in range(0, len(rest), 3):
        i, j, z = rest[k], rest[k + 1], rest[k + 2]
        col[j].append((i, 24 * je * z))
        if i != j:
            num = 24 * je * w[i] * z
            assert num % w[j] == 0
            col[i].append((j, num // w[j]))
    for i in range(n):
        col[i].append((i, 24 * jb * mag[i]))
    return col


def h_apply(col, v):
    out = defaultdict(int)
    for j, vj in v.items():
        if vj == 0:
            continue
        for i, aij in col[j]:
            out[i] += aij * vj
    return {k: val for k, val in out.items() if val != 0}


def wdot(w, u, v):
    if len(u) > len(v):
        u, v = v, u
    return sum(w[i] * ui * v.get(i, 0) for i, ui in u.items())


def inertia(m):
    a = [row[:] for row in m]
    k_dim = len(a)
    prev = 1
    signs = [1]
    for k in range(k_dim):
        p = a[k][k]
        if p == 0:
            raise ValueError("zero pivot: uncertified threshold")
        signs.append(1 if p > 0 else -1)
        for i in range(k + 1, k_dim):
            x = a[i][k]
            for j in range(i, k_dim):
                top = p * a[i][j] - x * a[k][j]
                q, rem = divmod(top, prev)
                if rem:
                    raise ValueError("non-exact Bareiss division")
                a[i][j] = q
                a[j][i] = q
            a[i][k] = 0
            a[k][i] = 0
        prev = p
    return sum(s != t for s, t in zip(signs, signs[1:]))


def bracket(a, b, k, scale=1000000, upper=4800):
    k_dim = len(a)

    def count(t):
        c = [[scale * b[i][j] - t * a[i][j] for j in range(k_dim)] for i in range(k_dim)]
        return inertia(c)

    lo, hi = 0, upper * scale
    assert count(lo) <= k and count(hi) > k
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if count(mid) <= k:
            lo = mid
        else:
            hi = mid
    return lo, hi


def solve_block_ritz(matrix_path, je, jb, order=3):
    start = monotonic_ns()
    n, w, mag, rest = load_matrix(matrix_path)
    col = build_system(n, w, mag, rest, je, jb)

    flats = [i for i in range(n) if mag[i] == 0]
    assert len(flats) == 8, f"expected 8 flat vacua, found {len(flats)}"

    basis = [{f: 1} for f in flats]
    last = basis
    for _ in range(order):
        last = [h_apply(col, b) for b in last]
        basis.extend(last)

    hcols = [h_apply(col, b) for b in basis]
    dim = len(basis)
    a = [[wdot(w, u, v) for v in basis] for u in basis]
    b = [[wdot(w, u, hv) for hv in hcols] for u in basis]

    for i in range(dim):
        for j in range(dim):
            assert a[i][j] == a[j][i]
            assert b[i][j] == b[j][i]

    if inertia(a) != 0:
        raise ValueError("Gram matrix not positive definite")

    # scale = 1_000_000, 24*H operator -> divide thresholds by 24 for micro-units
    energies = []
    for k in range(8):
        lo, hi = bracket(a, b, k, scale=1000000)
        energies.append((lo // 24, hi // 24))

    e0_lo, e0_hi = energies[0]
    e1_lo, e1_hi = energies[1]
    gap_lo = e1_lo - e0_hi
    gap_hi = e1_hi - e0_lo
    elapsed_ms = (monotonic_ns() - start) // 1000000

    name = Path(matrix_path).name
    e_str = ",".join(f"{lo}" for lo, _ in energies)
    print(
        f"BLOCK_RITZ|matrix={name}|states={n}|JE={je}|JB={jb}|order={order}|dim={dim}|"
        f"E0={e0_lo}/1000000|E1={e1_lo}/1000000|gap=[{gap_lo},{gap_hi}]/1000000|"
        f"all_8_vacua_micro={e_str}|ms={elapsed_ms}|kind=VACUUM_BLOCK_RITZ_UPPER_BOUNDS|json=0"
    )
    return energies


def main():
    p = argparse.ArgumentParser(description="Vacuum-ball block-Krylov Ritz solver")
    p.add_argument("matrix", type=str, help="path to .matrix file")
    p.add_argument("je", type=int, help="JE coupling")
    p.add_argument("jb", type=int, help="JB coupling")
    p.add_argument("--order", type=int, default=3, help="Krylov block order (default 3)")
    args = p.parse_args()
    solve_block_ritz(args.matrix, args.je, args.jb, args.order)


if __name__ == "__main__":
    main()
