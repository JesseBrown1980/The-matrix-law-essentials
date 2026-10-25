#!/usr/bin/env python3
"""Exact magnetic-end (second-order) law for the Lx x Ly S3 torus WITHOUT
building the full gauge-orbit matrix: only the 8 flat vacua and their
one-link electric neighbours are needed, since M is diagonal.

Same construction as torusl (Stage 2C generalised): comb tree, free-loop
encoding, global-conjugation orbits; electric per link sum_t (1 - L_t);
magnetic per plaquette 2 - chi_std. Flat vacua are built from the 18
commuting holonomy pairs (A, B): seam x-links = A, seam y-links = B.
K_ab = sum_j E_aj E_jb / M_j with E_jb = w_b E_bj / w_j (Gauss-weighted
self-adjointness). Integer-only; eigenvalues via second_order.py's integer
charpoly + Sturm isolation. Output json=0."""
import math, sys
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import second_order as SO

PERM = [(0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0)]
MUL = [[PERM.index(tuple(PERM[a][PERM[b][i]] for i in range(3))) for b in range(6)] for a in range(6)]
INV = [next(b for b in range(6) if MUL[a][b] == 0) for a in range(6)]
TRANS = (1, 2, 5)
def chi(a): return 2 if a == 0 else (-1 if a in (3, 4) else 0)

class Lat:
    def __init__(s, lx, ly):
        s.lx, s.ly, s.n = lx, ly, lx*ly
        site = s.site
        s.edges = [None]*(2*s.n)
        for x in range(lx):
            for y in range(ly):
                a = site(x, y)
                s.edges[a] = (a, site(x+1, y)); s.edges[s.n+a] = (a, site(x, y+1))
        s.tree = [s.n + site(0, y) for y in range(ly-1)] + [site(x, y) for x in range(lx-1) for y in range(ly)]
        ts = set(s.tree)
        s.free = [e for e in range(2*s.n) if e not in ts]
    def site(s, x, y): return (x % s.lx)*s.ly + (y % s.ly)
    def gauge_fix(s, links):
        h = [0]*s.n
        for e in s.tree:
            t, hd = s.edges[e]; h[hd] = MUL[h[t]][links[e]]
        return tuple(MUL[MUL[h[s.edges[e][0]]][links[e]]][INV[h[s.edges[e][1]]]] for e in s.free)
    def conj_codes(s, q):
        return [tuple(MUL[MUL[c][x]][INV[c]] for x in q) for c in range(6)]
    def key(s, q):
        # canonical = min base-6 little-endian code, as in torusl
        return min(sum(v * 6**i for i, v in enumerate(cq)) for cq in s.conj_codes(q))
    def weight(s, q):
        return len(set(s.conj_codes(q)))
    def rep(s, q):
        links = [0]*(2*s.n)
        for i, e in enumerate(s.free): links[e] = q[i]
        return links
    def magnetic(s, L):
        act = 0
        for x in range(s.lx):
            for y in range(s.ly):
                a, r, u = s.site(x, y), s.site(x+1, y), s.site(x, y+1)
                p = MUL[MUL[MUL[L[a]][L[s.n+r]]][INV[L[u]]]][INV[L[s.n+a]]]
                act += 2 - chi(p)
        return act

def law(lx, ly):
    lat = Lat(lx, ly)
    flats = {}
    for A in range(6):
        for B in range(6):
            if MUL[A][B] != MUL[B][A]: continue
            L = [0]*(2*lat.n)
            for y in range(ly): L[lat.site(lx-1, y)] = A           # x-links crossing the x seam
            for x in range(lx): L[lat.n + lat.site(x, ly-1)] = B    # y-links crossing the y seam
            assert lat.magnetic(L) == 0
            q = lat.gauge_fix(L)
            flats[lat.key(q)] = q
    keys = sorted(flats)
    nl = 2*lat.n
    E = {}      # E[a][j] for a flat, j neighbour key (off-diagonal counts, negative)
    info = {}   # j -> (weight, M)
    pdiag = {}
    for a in keys:
        q = flats[a]; raw = lat.rep(q)
        assert lat.magnetic(raw) == 0
        row = {}
        diag = 3*nl
        for e in range(nl):
            for t in TRANS:
                nx = raw[:]; nx[e] = MUL[t][nx[e]]
                qj = lat.gauge_fix(nx); j = lat.key(qj)
                if j == a:
                    diag -= 1; continue
                row[j] = row.get(j, 0) - 1
                if j not in info:
                    info[j] = (lat.weight(qj), lat.magnetic(lat.rep(qj)))
        E[a] = row; pdiag[a] = diag
    wf = {a: lat.weight(flats[a]) for a in keys}
    flat_flat = sum(1 for a in keys for j in E[a] if j in flats)
    mids = sorted({info[j][1] for a in keys for j in E[a] if j not in flats})
    Lc = 1
    for m in mids: Lc = Lc*m // math.gcd(Lc, m)
    # K_ab * Lc = sum_j E_aj * E_jb * Lc / M_j,  E_jb = w_b E_bj / w_j
    K = []
    for a in keys:
        row = []
        for b in keys:
            tot = 0
            for j, eaj in E[a].items():
                if j in flats or j not in E[b]: continue
                wj, mj = info[j]
                num = wf[b] * E[b][j]
                assert num % wj == 0
                tot += eaj * (num // wj) * (Lc // mj)
            row.append(tot)
        K.append(row)
    p = SO.charpoly(K)
    rts, _ = SO.roots(p)
    kap = sorted(rts, key=lambda x: -x[0])
    bits = 60
    # multiplicities for rational roots, if any
    out = ("FLAT_LAW|L=%dx%d|links=%d|flat_vacua=%d|PEP_diag=%s|flat_to_flat=%d|neighbours=%d|intermediate_M=%s|Lc=%d|charpoly=%s|distinct=%d|json=0"
           % (lx, ly, nl, len(keys), ",".join(map(str, sorted(set(pdiag.values())))), flat_flat, len(info), ",".join(map(str, mids)), Lc,
              ",".join(map(str, p)), len(kap)))
    print(out, flush=True)
    for lo, hi in kap:
        print("FLAT_LAW|L=%dx%d|kappa_in=[%s,%s]|json=0" % (lx, ly, SO.fmt_div(lo, bits, Lc), SO.fmt_div(hi, bits, Lc)), flush=True)
    k0, k1 = kap[0], kap[1]
    print("FLAT_LAW|L=%dx%d|gap_coeff_in=[%s,%s]|law=E_k~%d*JE-kappa_k*JE^2/JB|json=0"
          % (lx, ly, SO.fmt_div(k0[0]-k1[1], bits, Lc), SO.fmt_div(k0[1]-k1[0], bits, Lc), sorted(set(pdiag.values()))[0]), flush=True)

if __name__ == "__main__":
    for spec in sys.argv[1:]:
        lx, ly = map(int, spec.split("x"))
        law(lx, ly)
