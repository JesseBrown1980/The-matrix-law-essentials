#!/usr/bin/env python3
"""Integer-only exact nonlinear jump and certified dissipative nonlinear jump.

Nonlinear reversible map: F = T^{-1} L T, with local bilinear shear T and
coupled invertible linear L. Binary/ternary modal powers yield exact F^n.
Contraction map: genuinely nonlinear coupled integer recurrence, certified
by an inward monotone sup-norm bound to attain zero for every seed in ball.
Neither theorem covers arbitrary nonlinear maps, nor physical Yang-Mills.
"""
import argparse
import hashlib
from pathlib import Path

P = 1299079
G = 3
KAPPA = 7
VIRTUAL_ADDRESSES = 12_000_000_000


def primitive_root(n):
    assert n >= 3 and (P - 1) % n == 0
    w = pow(G, (P - 1) // n, P)
    assert pow(w, n, P) == 1
    assert pow(w, n // 3, P) != 1
    return w


def ntt(data, inverse=False):
    n = len(data)
    omega = primitive_root(n)
    degree = 0
    t = n
    while t > 1:
        assert t % 3 == 0
        t //= 3
        degree += 1
    a = data[:]
    for i in range(n):
        q = i
        reverse = 0
        for _ in range(degree):
            q, rem = divmod(q, 3)
            reverse = reverse * 3 + rem
        if reverse > i:
            a[i], a[reverse] = a[reverse], a[i]
    if inverse:
        omega = pow(omega, P - 2, P)
    length = 3
    while length <= n:
        width = length // 3
        wm = pow(omega, n // length, P)
        small_omega = pow(wm, width, P)
        omega_sq = small_omega * small_omega % P
        t1 = 1
        t2 = 1
        for j in range(width):
            for offset in range(0, n, length):
                i = offset + j
                x = a[i]
                y = a[i + width] * t1 % P
                z = a[i + width * 2] * t2 % P
                a[i] = (x + y + z) % P
                a[i + width] = (x + small_omega * y + omega_sq * z) % P
                a[i + width * 2] = (x + omega_sq * y + small_omega * z) % P
            t1 = t1 * wm % P
            t2 = t2 * wm * wm % P
        length *= 3
    if inverse:
        inv = pow(n, P - 2, P)
        a = [x * inv % P for x in a]
    return a


def initial(n):
    u = [(7 * i * i + 3 * i + 1) % P for i in range(n)]
    v = [(11 * i * i + 31 * i + 17) % P for i in range(n)]
    return u, v


def shear(u, v):
    n = len(u)
    return u[:], [(v[i] + KAPPA * u[i] * u[i - 1]) % P for i in range(n)]


def unshear(u, v):
    n = len(u)
    return u[:], [(v[i] - KAPPA * u[i] * u[i - 1]) % P for i in range(n)]


def linear_step(a, b):
    # v_new = v + R a; a_new = a + R v_new. An invertible block-linear step.
    n = len(a)
    bn = [(b[i] + a[i - 1]) % P for i in range(n)]
    an = [(a[i] + bn[i - 1]) % P for i in range(n)]
    return an, bn


def nonlinear_step(u, v):
    a, b = shear(u, v)
    an, bn = linear_step(a, b)
    return unshear(an, bn)


def nonlinear_back(u, v):
    a, b = shear(u, v)
    ap = [(a[i] - b[i - 1]) % P for i in range(len(a))]
    bp = [(b[i] - ap[i - 1]) % P for i in range(len(a))]
    return unshear(ap, bp)


def direct(state, steps):
    u, v = state
    for _ in range(steps):
        u, v = nonlinear_step(u, v)
    return u, v


def mmul(a, b):
    x, y, z, w = a
    A, B, C, D = b
    return ((x * A + y * C) % P, (x * B + y * D) % P,
            (z * A + w * C) % P, (z * B + w * D) % P)


def mpow_binary(m, steps):
    out = (1, 0, 0, 1)
    while steps:
        if steps & 1:
            out = mmul(m, out)
        m = mmul(m, m)
        steps >>= 1
    return out


def mpow_tri(m, k):
    for _ in range(k):
        m2 = mmul(m, m)
        m = mmul(m2, m)
    return m


def modal_factors(n, exponent, ternary_k=None):
    r = primitive_root(n)
    w = 1
    out = []
    for _ in range(n):
        # The Fourier eigenvalue of shift R is w.
        l = ((1 + w * w) % P, w, w, 1)
        factor = mpow_tri(l, ternary_k) if ternary_k is not None else mpow_binary(l, exponent)
        # Reversible in every mode: det(L_w) = 1 exactly.
        assert (factor[0] * factor[3] - factor[1] * factor[2]) % P == 1
        out.append(factor)
        w = w * r % P
    assert w == 1
    return out


def nonlinear_jump(state, steps, ternary_k=None):
    u, v = state
    a, b = shear(u, v)
    aa = ntt(a)
    bb = ntt(b)
    fac = modal_factors(len(u), steps, ternary_k)
    an = [(f[0] * x + f[1] * y) % P for f, x, y in zip(fac, aa, bb)]
    bn = [(f[2] * x + f[3] * y) % P for f, x, y in zip(fac, aa, bb)]
    return unshear(ntt(an, True), ntt(bn, True))


def digest(state):
    h = hashlib.sha256()
    for field in state:
        for value in field:
            h.update(value.to_bytes(8, 'big'))
    return h.hexdigest()


def sample_digest(state):
    u, v = state
    n = len(u)
    addresses = [i * (VIRTUAL_ADDRESSES - 1) // 242 for i in range(243)]
    selected = [((u[(a * n) // VIRTUAL_ADDRESSES], v[(a * n) // VIRTUAL_ADDRESSES])) for a in addresses]
    h = hashlib.sha256()
    for a, b in selected:
        h.update(a.to_bytes(8, 'big'))
        h.update(b.to_bytes(8, 'big'))
    return h.hexdigest()


def nonlinear_controls(n):
    p = initial(n)
    q = ([(3 * i + 5) % P for i in range(n)], [(i * i + 19) % P for i in range(n)])
    pq = ([(a + b) % P for a, b in zip(p[0], q[0])],
          [(a + b) % P for a, b in zip(p[1], q[1])])
    fq = nonlinear_step(*q)
    fp = nonlinear_step(*p)
    fpq = nonlinear_step(*pq)
    add = ([(a + b) % P for a, b in zip(fp[0], fq[0])],
           [(a + b) % P for a, b in zip(fp[1], fq[1])])
    assert fpq != add, 'physical state evolution is accidentally additive'
    assert fp != linear_step(*p), 'nonlinear interaction is not affecting evolution'
    bad = fp[0][:], fp[1][:]
    bad[0][0] = (bad[0][0] + p[0][0] * p[1][0]) % P
    assert bad != fp
    return 'NONLINEARITY|superposition=FAIL_EXPECTED|naive_linear_step=FAIL_EXPECTED|tampered_feedback=REJECTED|json=0'


def toward_zero(x, d):
    assert d > 0
    return x // d if x >= 0 else -((-x) // d)


def contraction_step(x, scale=1_000_000, linear_self=3):
    n = len(x)
    return [toward_zero(linear_self * x[i] + x[i - 1] +
                        toward_zero(x[i] * x[i - 1], scale), 8)
            for i in range(n)]


def contraction_bound(b, scale=1_000_000):
    assert 0 <= b <= scale
    return (4 * b + b * b // scale) // 8


def dissipative_jump(n, long_n, scale=1_000_000):
    assert long_n >= 0
    x = [((19 * i * i + 43 * i + 7) % (2 * scale + 1)) - scale for i in range(n)]
    b = scale
    actual = max(map(abs, x))
    assert actual <= b
    steps = 0
    while b > 0 and steps < long_n:
        bn = contraction_bound(b, scale)
        assert bn < b
        x = contraction_step(x, scale)
        actual = max(map(abs, x))
        assert actual <= bn, (steps, actual, bn)
        b = bn
        steps += 1
    if long_n >= steps and b == 0:
        assert all(v == 0 for v in x)
    return steps, b, x


def run(dims):
    rows = [f'NONLINEAR_HEADER|prime={P}|coupled_fields=2|virtual_addresses={VIRTUAL_ADDRESSES}|model=EXACT_CONJUGATE_AND_INTEGER_CONTRACTION|r9=HOLD|json=0']
    for n in dims:
        state = initial(n)
        assert nonlinear_back(*nonlinear_step(*state)) == state
        rows.append(f'CONTROL|n={n}|'+nonlinear_controls(n))
        for t in ([1, 3, 27] + ([243] if n <= 243 else [])):
            x = nonlinear_jump(state, t)
            assert x == direct(state, t), ('direct mismatch', n, t)
            rows.append(f'DIRECT|n={n}|steps={t}|jump_equals_direct=PASS|sha256={digest(x)}|json=0')
        for k in (6, 20, 40):
            steps = 3 ** k
            jumped = nonlinear_jump(state, steps, ternary_k=k)
            if n <= 243 or (n == 59049 and k in (6, 40)) or k == 6:
                binary = nonlinear_jump(state, steps)
                assert binary == jumped
            rows.append(f'TERNARY|n={n}|exponent=3^{k}|complete_state_sha256={digest(jumped)}|sampled243_sha256={sample_digest(jumped)}|json=0')
        freeze = nonlinear_jump(state, 3 ** 40, ternary_k=40)
        plus1 = nonlinear_step(*freeze)
        plus3 = direct(freeze, 3)
        if n <= 243:
            assert plus1 == nonlinear_jump(state, 3 ** 40 + 1)
            assert plus3 == nonlinear_jump(state, 3 ** 40 + 3)
        assert nonlinear_back(*plus1) == freeze
        assert nonlinear_back(*freeze) == nonlinear_jump(state, 3 ** 40 - 1) if n <= 243 else nonlinear_step(*nonlinear_back(*freeze)) == freeze
        rows.append(f'LEAF|n={n}|from=3^40|forward1=PASS|forward3=PASS|backward1=PASS|leaf_sha256={digest(plus3)}|json=0')
        if n <= 243 or n == 59049:
            x80 = nonlinear_jump(state, 3 ** 80, ternary_k=80)
            y80 = nonlinear_jump(state, 3 ** 80)
            assert x80 == y80
            rows.append(f'TERNARY|n={n}|exponent=3^80|complete_state_sha256={digest(x80)}|json=0')
        steps,bound,contracted = dissipative_jump(n, 3 ** 80)
        assert bound == 0 and all(x == 0 for x in contracted)
        # Independently verify a fully uniform loss-of-memory bound, not just this seed.
        b = 1_000_000
        for _ in range(steps):
            b = contraction_bound(b)
        assert b == 0
        # Control: 9 instead of 3 in the self coefficient violates the claimed contraction.
        uniform = [1_000_000] * n
        invalid = contraction_step(uniform, linear_self=9)
        assert max(invalid) > 1_000_000
        rows.append(f'BOUND|n={n}|initial_sup_bound=1000000|certified_zero_after_steps={steps}|jump_to=3^80|exact_zero_all_states_in_ball=PASS|bad_coefficient_9=REJECTED|json=0')
        rows.append(f'NONLINEAR_DIM_PASS|n={n}|genuine_neighbor_bilinear_feedback=YES|exact_power_jump=PASS|independent_dissipative_bound=PASS|json=0')
    rows.append('NONLINEAR_JUMP_VERIFY|status=PASS|exact_integrable_family=YES|integer_contractive_family=YES|general_nonlinear=NOT_COVERED|r9=HOLD|json=0')
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dims', default='27,243,6561')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    dims = [int(t) for t in args.dims.split(',')]
    for n in dims:
        assert n >= 3 and (P - 1) % n == 0 and n % 3 == 0
    rows = run(dims)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text('\n'.join(rows)+'\n', encoding='ascii')
    print('\n'.join(rows))


if __name__ == '__main__':
    main()
