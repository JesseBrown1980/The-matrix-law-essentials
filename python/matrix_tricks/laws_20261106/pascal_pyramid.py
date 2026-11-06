# Pascal's pyramid and the power-of-three jump, exact integers, HBP rows (json=0).
# 1. Path counts on the 27-cube = trinomial coefficients L!/(i!j!k!) folded by digit mod 3 (Pascal's pyramid).
# 2. Mod 3, the step's 3^k-th power collapses (Lucas / Frobenius): (I + R + R^2)^(3^k) = I + R^(3^k) + R^(2*3^k) (mod 3), likewise with signs.
import itertools
from math import factorial


def counts_by_paths(L):
    c = {}
    for seq in itertools.product(range(3), repeat=L):
        s = [0, 0, 0]
        for d in seq:
            s[d] = (s[d] + 1) % 3
        c[tuple(s)] = c.get(tuple(s), 0) + 1
    return c


def counts_by_pyramid(L):
    c = {}
    for i in range(L + 1):
        for j in range(L + 1 - i):
            k = L - i - j
            key = (i % 3, j % 3, k % 3)
            c[key] = c.get(key, 0) + factorial(L) // (factorial(i) * factorial(j) * factorial(k))
    return c


for L in (6, 9):
    a, b = counts_by_paths(L), counts_by_pyramid(L)
    print("PASCAL_PYRAMID|L=%d|positions=%d|path_counts_equal_folded_trinomials=%d|counts=%s|json=0" % (
        L, len(a), int(a == b), ",".join(str(v) for v in sorted(set(a.values())))))

D = 27


def mul(A, B, p):
    return [[sum(A[i][k] * B[k][j] for k in range(D)) % p for j in range(D)] for i in range(D)]


def circ(coef, p):
    M = [[0] * D for _ in range(D)]
    for i in range(D):
        for s, c in coef.items():
            M[i][(i - s) % D] = (M[i][(i - s) % D] + c) % p
    return M


for name, coef in (("I+R+R^2", {0: 1, 1: 1, 2: 1}), ("I+R-R^2", {0: 1, 1: 1, 2: -1})):
    ok = 1
    for k in (1, 2, 3):
        n = 3 ** k
        M = circ(coef, 3)
        P = M
        for _ in range(k):
            P = mul(mul(P, P, 3), P, 3)
        Q = [[0] * D for _ in range(D)]
        for s, c in coef.items():
            for i in range(D):
                Q[i][(i - s * n) % D] = (Q[i][(i - s * n) % D] + c) % 3
        ok &= int(P == Q)
    print("LUCAS_JUMP_MOD3|step=%s|ring=27|k=1,2,3|M^(3^k)=I+R^(3^k)(+/-)R^(2*3^k)_mod_3=%d|json=0" % (name, ok))
