#!/usr/bin/env python3
"""Independent sparse-step reference for the exact long-jump control.

Reads the exchange-format matrix itself (no shared code with jumpto.rs),
builds 24H as a sparse integer operator from (weights, magnetic, i j z),
recomputes the Gershgorin shift c24 = max row |24H| + 24, regenerates the
same deterministic LCG start vector, and applies v <- (c24 - 24H) v STEPS
times one sparse step at a time, modulo the prime P. Prints the same
fingerprints as `jumpto ... exact P STEPS`:
  state_sum_mod       = sum_i v_i  (mod P)
  state_weighted_mod  = sum_i (i+1) v_i  (mod P)
  v0_component0       = v_0
Integer-only. usage: sparse_mod_reference.py MATRIX JE JB P STEPS   (json=0)
"""
import sys

MASK = (1 << 64) - 1


def load(path, je, jb):
    t = open(path).read().split()
    n = int(t[0])
    w = [int(x) for x in t[1:1 + n]]
    mag = [int(x) for x in t[1 + n:1 + 2 * n]]
    rows = [dict() for _ in range(n)]
    rest = t[1 + 2 * n:]
    for k in range(0, len(rest), 3):
        i, j, z = int(rest[k]), int(rest[k + 1]), int(rest[k + 2])
        rows[i][j] = rows[i].get(j, 0) + 24 * je * z
        if i != j:
            num = 24 * je * w[i] * z
            assert num % w[j] == 0
            rows[j][i] = rows[j].get(i, 0) + num // w[j]
    for i in range(n):
        rows[i][i] = rows[i].get(i, 0) + 24 * jb * mag[i]
    return n, rows


def lcg(n):
    seed, out = 12345, []
    for _ in range(n):
        seed = (seed * 6364136223846793005 + 1442695040888963407) & MASK
        out.append((seed >> 24) - (1 << 39))
    return out


def main():
    path, je, jb, p, steps = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    n, rows = load(path, je, jb)
    c24 = max(sum(abs(v) for v in r.values()) for r in rows) + 24
    v = [x % p for x in lcg(n)]
    for _ in range(steps):
        v = [(c24 * v[i] - sum(a * v[j] for j, a in rows[i].items())) % p for i in range(n)]
    s = sum(v) % p
    ws = sum((i + 1) * x for i, x in enumerate(v)) % p
    print("SPARSE_REF|matrix=%s|n=%d|JE=%d|JB=%d|c=(%d)/24|modulus=%d|steps=%d|state_sum_mod=%d|state_weighted_mod=%d|v0_component0=%d|json=0"
          % (path.rsplit("/", 1)[-1], n, je, jb, c24, p, steps, s, ws, v[0]))


if __name__ == "__main__":
    main()
