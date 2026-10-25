#!/usr/bin/env python3
"""Exact second-order magnetic-end law on the 8 flat vacua, integer-only.

H = JE E + JB M, M diagonal >= 0, P = span of the M = 0 (flat) orbits.
For JB >> JE:  H_eff = JE*PEP - (JE^2/JB) * K,  K_ab = sum_{i: M_i>0} E_ai E_ib / M_i.
Eigenvalues: E_k ~ JE*p - (JE^2/JB) kappa_k  (if PEP = p I), so the
magnetic-end gap ~ (kappa_0 - kappa_1) JE^2/JB.
Everything is done on L*K with L = lcm of the magnetic values met, an
integer matrix; its characteristic polynomial comes from integer
Faddeev-LeVerrier (exact divisions), its real roots are isolated by an
integer Sturm sequence and bisected to dyadic brackets. No float,
Decimal or Fraction. Output lines end json=0."""
import math, sys

def load(path):
    t = open(path).read().split()
    n = int(t[0]); w = [int(x) for x in t[1:1+n]]; mag = [int(x) for x in t[1+n:1+2*n]]
    rows = [dict() for _ in range(n)]
    r = t[1+2*n:]
    for k in range(0, len(r), 3):
        i, j, z = int(r[k]), int(r[k+1]), int(r[k+2])
        rows[i][j] = z
        if i != j:
            assert (w[i]*z) % w[j] == 0
            rows[j][i] = w[i]*z // w[j]
    return n, w, mag, rows

def charpoly(A):
    """integer Faddeev-LeVerrier: det(xI - A) = sum c[k] x^(n-k); exact."""
    n = len(A)
    c = [1]
    Mk = [[0]*n for _ in range(n)]
    for k in range(1, n+1):
        # Mk = A*M_{k-1} + c_{k-1} I
        AM = [[sum(A[i][l]*Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        Mk = [[AM[i][j] + (c[-1] if i == j else 0) for j in range(n)] for i in range(n)]
        AMk = [[sum(A[i][l]*Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        tr = sum(AMk[i][i] for i in range(n))
        assert tr % k == 0
        c.append(-tr // k)
    return c  # highest degree first

def peval(p, x_num, x_den):
    """sign-exact evaluation of p at x_num/x_den, returns integer p(x)*x_den^deg"""
    d = len(p) - 1
    return sum(a * x_num**(d-i) * x_den**i for i, a in enumerate(p))

def pdiv_rem(a, b):
    """remainder of a / b with integer pseudo-division, content removed"""
    a = a[:]
    while len(a) >= len(b) and any(a):
        if a[0] == 0:
            a.pop(0); continue
        f, g = a[0], b[0]
        a = [x*g for x in a]
        for i in range(len(b)):
            a[i] -= f * b[i]
        a.pop(0)
    while a and a[0] == 0:
        a.pop(0)
    if a:
        cont = 0
        for x in a: cont = math.gcd(cont, x)
        a = [x // cont for x in a]
    return a

def deriv(p):
    d = len(p) - 1
    return [a*(d-i) for i, a in enumerate(p[:-1])]

def sturm_chain(p):
    ch = [p, deriv(p)]
    while True:
        r = pdiv_rem(ch[-2], ch[-1])
        if not r:
            break
        ch.append([-x for x in r])
    return ch

def sign(x): return (x > 0) - (x < 0)

def variations(ch, xn, xd):
    s = [sign(peval(q, xn, xd)) for q in ch]
    s = [v for v in s if v]
    return sum(1 for a, b in zip(s, s[1:]) if a != b)

def squarefree(p):
    g = p[:]; h = deriv(p)
    while h:
        g, h = h, pdiv_rem(g, h)
    if len(g) <= 1:
        return p, []
    cont = 0
    for x in g: cont = math.gcd(cont, x)
    if g[0] < 0: cont = -cont
    g = [x // cont for x in g]
    # p / g by exact long division
    q, r = [], p[:]
    while len(r) >= len(g):
        f = r[0] // g[0]; assert r[0] % g[0] == 0
        q.append(f)
        for i in range(len(g)): r[i] -= f*g[i]
        r.pop(0)
    assert not any(r)
    return q, g

def roots(p, bits=60, bound=1 << 20):
    """isolate distinct real roots of integer p in (-bound, bound); dyadic brackets [lo, hi]/2^bits"""
    sq, _ = squarefree(p)
    ch = sturm_chain(sq)
    D = 1 << bits
    out = []
    def rec(lo, hi):  # numerators over D, half-open (lo, hi]
        n = variations(ch, lo, D) - variations(ch, hi, D)
        if n == 0: return
        if n == 1 and hi - lo <= 1:
            out.append((lo, hi)); return
        mid = (lo + hi) // 2
        rec(lo, mid); rec(mid, hi)
    rec(-bound*D, bound*D)
    return sorted(out), sq

def fmt(lo, hi, bits, scale=10**9):
    a = lo * scale >> bits
    return "%d.%09d" % (a // scale, a % scale) if a >= 0 else "-%d.%09d" % ((-a) // scale, (-a) % scale)

def main(path):
    n, w, mag, rows = load(path)
    flat = [i for i in range(n) if mag[i] == 0]
    pep = [[rows[a].get(b, 0) for b in flat] for a in flat]
    mids = sorted({mag[i] for a in flat for i in rows[a] if mag[i] > 0})
    L = 1
    for m in mids: L = L * m // math.gcd(L, m)
    K = [[sum(rows[a][i] * rows[i].get(b, 0) * (L // mag[i]) for i in rows[a] if mag[i] > 0) for b in flat] for a in flat]
    p = charpoly(K)
    rts, sq = roots(p)
    bits = 60
    # eigenvalues of K are roots / L
    kap = [((lo, hi)) for lo, hi in rts]
    kap.sort(key=lambda x: -x[0])
    name = path.rsplit("/", 1)[-1]
    pdiag = sorted({pep[i][i] for i in range(len(flat))})
    poff = sum(1 for i in range(len(flat)) for j in range(len(flat)) if i != j and pep[i][j])
    print("SECOND_ORDER|matrix=%s|flat=%d|PEP_diag=%s|PEP_offdiag_nonzero=%d|intermediate_M=%s|L=%d|charpoly_LK=%s|distinct_roots=%d|json=0"
          % (name, len(flat), ",".join(map(str, pdiag)), poff, ",".join(map(str, mids)), L, ",".join(map(str, p)), len(rts)))
    for lo, hi in kap:
        # kappa = root/L ; print with 9 decimals (floor of lower bracket)
        print("SECOND_ORDER|matrix=%s|kappa_in=[%s,%s]|json=0" % (name, fmt(lo, bits + 0, 1, 10**9) if False else fmt_div(lo, bits, L), fmt_div(hi, bits, L)))
    k0, k1 = kap[0], kap[1]
    # gap coefficient bracket: (k0_lo - k1_hi, k0_hi - k1_lo)/L
    print("SECOND_ORDER|matrix=%s|gap_coeff_in=[%s,%s]|law=gap~coeff*JE^2/JB|json=0"
          % (name, fmt_div(k0[0] - k1[1], bits, L), fmt_div(k0[1] - k1[0], bits, L)))

def fmt_div(x, bits, L, scale=10**9):
    a = (x * scale) // (L << bits)
    return ("-" if a < 0 else "") + "%d.%09d" % (abs(a) // scale, abs(a) % scale)

if __name__ == "__main__":
    for path in sys.argv[1:]:
        main(path)
