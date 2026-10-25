#!/usr/bin/env python3
"""Fourth-order splitting of the 8 flat vacua (integer-only), for tori where
second order is scalar (K2 = k I) and third order vanishes (H3 = 0).
Then every fourth-order Kato term except P V R V R V R V P is a multiple
of the identity on P (they are built from length-2 returns with all
intermediate M equal, and from length-3 paths, which vanish); checked:
K2 scalar, H3 == 0, all one-step intermediate M equal, PVR^3VP scalar.
  splitting term: -(JE^4/JB^3) Y,  Y_ab = sum E_aj E_jk E_kl E_lb / (M_j M_k M_l)
  (R = -Q/(JB M), three R's give the minus sign); ground = top eigenvalue of Y.
E_kl for two-step states k is taken from Gauss-weighted self-adjointness
E_kl = w_l E_lk / w_k (asserted integral), so only one-step rows are built.
Output json=0."""
import math, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import second_order as SO
import flat_law3 as F3

def main(spec):
    lx, ly = map(int, spec.split("x"))
    lat, fk, c, Ms, D, K2, H3 = F3.build(lx, ly)
    assert all(K2[i][j] == (K2[0][0] if i == j else 0) for i in range(8) for j in range(8)), "K2 not scalar"
    assert all(v == 0 for r in H3 for v in r), "H3 nonzero: third order already splits"
    assert len(Ms) == 1, "one-step intermediate M not all equal"
    # rebuild rows of flats and N1 (cheap) via the same code path
    from flat_law import MUL, TRANS
    flats = set(fk)
    reps = {}
    for A in range(6):
        for B in range(6):
            if MUL[A][B] != MUL[B][A]: continue
            L = [0]*(2*lat.n)
            for y in range(ly): L[lat.site(lx-1, y)] = A
            for x in range(lx): L[lat.n + lat.site(x, ly-1)] = B
            q = lat.gauge_fix(L); reps[lat.key(q)] = q
    nl = 2*lat.n
    rows = {}
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
    N1 = sorted({j for a in fk for j in row_of(a) if j not in flats})
    for j in N1: row_of(j)
    mag, w = {}, {}
    def mw(k):
        if k not in mag:
            mag[k] = lat.magnetic(lat.rep(reps[k])); w[k] = lat.weight(reps[k])
        return mag[k], w[k]
    for k in list(fk) + N1: mw(k)
    N2 = sorted({k for j in N1 for k in rows[j] if k not in flats})
    for k in N2: mw(k)
    Mj = Ms[0]
    M2 = sorted({mag[k] for k in N2})
    Lk = 1
    for m in M2: Lk = Lk*m // math.gcd(Lk, m)
    # u_a[k] = sum_j E_aj E_jk / M_j  (x Mj),   v_b[k] = sum_l E_kl E_lb / M_l (x Mj)
    def u(a):
        out = {}
        for j, eaj in rows[a].items():
            if j in flats: continue
            for k, ejk in rows[j].items():
                if k in flats: continue
                out[k] = out.get(k, 0) + eaj*ejk
        return out
    def v(b):
        out = {}
        for l in rows[b]:
            if l in flats: continue
            elb = rows[l].get(b, 0)
            if not elb: continue
            for k, elk in rows[l].items():
                if k in flats: continue
                num = w[l]*elk
                assert num % w[k] == 0
                out[k] = out.get(k, 0) + (num // w[k])*elb
        return out
    U = {a: u(a) for a in fk}; V = {b: v(b) for b in fk}
    # Y * (Mj^2 * Lk) = sum_k u_a[k] v_b[k] * (Lk / M_k)
    Y = [[sum(U[a][k]*V[b].get(k, 0)*(Lk // mag[k]) for k in U[a]) for b in fk] for a in fk]
    den = Mj*Mj*Lk
    # scalar check of the c^2 PVR^3VP-type term: sum_j E_aj E_ja / M_j^3 equal for all a
    r3 = {sum(rows[a][j]*rows[j].get(a, 0) for j in rows[a] if j not in flats) for a in fk}
    assert len(r3) == 1, "PVR^3VP not scalar"
    p = SO.charpoly(Y)
    rts, _ = SO.roots(p, bits=60, bound=1 << 80)
    bits = 60
    print("FLAT_LAW4|L=%s|c=%d|one_step=%d|two_step=%d|two_step_M=%s|den=%d|charpoly_Y=%s|distinct=%d|json=0"
          % (spec, c, len(N1), len(N2), ",".join(map(str, M2)), den, ",".join(map(str, p)), len(rts)), flush=True)
    for lo, hi in sorted(rts, key=lambda x: -x[0]):
        print("FLAT_LAW4|L=%s|Y_eig_in=[%s,%s]|json=0" % (spec, SO.fmt_div(lo, bits, den), SO.fmt_div(hi, bits, den)), flush=True)
    top = sorted(rts, key=lambda x: -x[0])
    if len(top) >= 2:
        print("FLAT_LAW4|L=%s|gap_coeff4_in=[%s,%s]|law=gap~coeff*JE^4/JB^3|json=0"
              % (spec, SO.fmt_div(top[0][0]-top[1][1], bits, den), SO.fmt_div(top[0][1]-top[1][0], bits, den)), flush=True)

if __name__ == "__main__":
    for s in sys.argv[1:]:
        main(s)
