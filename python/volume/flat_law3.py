#!/usr/bin/env python3
"""Third-order effective Hamiltonian on the 8 flat vacua (integer-only).

H = JE E + JB M, H0 = JB M (diagonal, 0 on P, >0 on Q), V = JE E.
With PVP = c*JE*I (checked: no flat->flat hops, equal diagonal c):
  H_eff = c JE  - (JE^2/JB) K2  + (JE^3/JB^2) H3
  K2_ab = sum_j E_aj E_jb / M_j
  H3_ab = sum_{j,k} E_aj E_jk E_kb / (M_j M_k)  -  c * sum_j E_aj E_jb / M_j^2
(Kato / Bloch effective Hamiltonian to third order; R = -Q/(JB M).)
All E entries are matrix rows (moves from each orbit representative),
exactly as torusl builds them. Uses the one- and two-step neighbourhoods
of the vacua only. Eigenvalues of H_eff(r) at an integer r come from the
integer matrix D*(c r^2 I - r K2 + H3) via integer charpoly + Sturm.
Output json=0."""
import math, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import second_order as SO
from flat_law import Lat, MUL, TRANS

def build(lx, ly):
    lat = Lat(lx, ly)
    flats = {}
    for A in range(6):
        for B in range(6):
            if MUL[A][B] != MUL[B][A]: continue
            L = [0]*(2*lat.n)
            for y in range(ly): L[lat.site(lx-1, y)] = A
            for x in range(lx): L[lat.n + lat.site(x, ly-1)] = B
            q = lat.gauge_fix(L); flats[lat.key(q)] = q
    nl = 2*lat.n
    reps, rows, mag = {}, {}, {}
    def row_of(key):
        if key in rows: return rows[key]
        raw = lat.rep(reps[key]); r = {key: 3*nl}
        for e in range(nl):
            for t in TRANS:
                nx = raw[:]; nx[e] = MUL[t][nx[e]]
                qj = lat.gauge_fix(nx); j = lat.key(qj)
                reps.setdefault(j, qj)
                r[j] = r.get(j, 0) - 1
        rows[key] = r
        return r
    for k, q in flats.items():
        reps[k] = q; mag[k] = 0
    fk = sorted(flats)
    N1 = set()
    for a in fk:
        for j in row_of(a):
            if j not in flats: N1.add(j)
    for j in N1:
        mag[j] = lat.magnetic(lat.rep(reps[j]))
        row_of(j)
    cset = {row_of(a)[a] for a in fk}
    assert len(cset) == 1, "PEP not scalar"
    c = cset.pop()
    assert all(j not in flats or j == a for a in fk for j in row_of(a)), "flat->flat hop"
    Ms = sorted({mag[j] for j in N1})
    Lm = 1
    for m in Ms: Lm = Lm*m // math.gcd(Lm, m)
    D = Lm*Lm
    # K2 scaled by D (so later D*(...) is integral), H3 scaled by D
    K2 = [[sum(row_of(a)[j]*row_of(j).get(b, 0)*(D//mag[j]) for j in row_of(a) if j in N1) for b in fk] for a in fk]
    H3 = []
    for a in fk:
        ra = row_of(a); hrow = []
        for b in fk:
            tot = 0
            for j, eaj in ra.items():
                if j not in N1: continue
                for k, ejk in row_of(j).items():
                    if k not in N1: continue
                    ekb = row_of(k).get(b, 0)
                    if ekb: tot += eaj*ejk*ekb*(Lm//mag[j])*(Lm//mag[k])
                ejb = row_of(j).get(b, 0)
                if ejb: tot -= c*eaj*ejb*(D//(mag[j]*mag[j]))
            hrow.append(tot)
        H3.append(hrow)
    return lat, fk, c, Ms, D, K2, H3

def eig_brackets(Mint, bits=60, bound=1 << 70):
    p = SO.charpoly(Mint)
    rts, _ = SO.roots(p, bits=bits, bound=bound)
    return sorted(rts), p

def show(x, bits, den):
    return SO.fmt_div(x, bits, den)

def main(spec, rs):
    lx, ly = map(int, spec.split("x"))
    lat, fk, c, Ms, D, K2, H3 = build(lx, ly)
    # third-order spectrum of H3 alone (meaningful when K2 is scalar)
    k2s = sorted({K2[i][i] for i in range(8)}); k2off = sum(1 for i in range(8) for j in range(8) if i != j and K2[i][j])
    print("FLAT_LAW3|L=%s|c=%d|intermediate_M=%s|D=%d|K2_diag_xD=%s|K2_offdiag_nonzero=%d|json=0" % (spec, c, ",".join(map(str, Ms)), D, ",".join(map(str, k2s)), k2off), flush=True)
    rts, p = eig_brackets(H3)
    bits = 60
    for lo, hi in rts:
        print("FLAT_LAW3|L=%s|H3_eig_in=[%s,%s]|json=0" % (spec, show(lo, bits, D), show(hi, bits, D)), flush=True)
    if k2off == 0 and len(k2s) == 1 and len(rts) >= 2:
        g_lo, g_hi = rts[1][0] - rts[0][1], rts[1][1] - rts[0][0]
        print("FLAT_LAW3|L=%s|second_order_degenerate=1|gap_coeff3_in=[%s,%s]|law=gap~coeff*JE^3/JB^2|json=0" % (spec, show(g_lo, bits, D), show(g_hi, bits, D)), flush=True)
    for r in rs:
        Mr = [[D*c*r*r*(i == j) - r*K2[i][j] + H3[i][j] for j in range(8)] for i in range(8)]
        er, _ = eig_brackets(Mr)
        den = D*r*r
        e0, e1 = er[0], er[1]
        print("FLAT_LAW3|L=%s|r=%d|E0_eff=%s|E1_eff=%s|gap_eff=%s|json=0" % (spec, r, show(e0[0], bits, den), show(e1[0], bits, den), show(e1[0]-e0[0], bits, den)), flush=True)

if __name__ == "__main__":
    main(sys.argv[1], [int(x) for x in sys.argv[2:]])
