#!/usr/bin/env python3
"""Power-of-three jump: the state at step 3^k of one fixed linear step, exact integers mod p.

The law (docs/LAW-POWER-OF-THREE-JUMP.md). A step that is the same at every
site of a ring of D sites,

    x'[i] = sum over the stencil of c_s * x[(i - s) mod D]   (mod p),

reaches step n = 3^k two ways, each in work that grows with k = log_3 n:
  cube : M <- M*M*M, k times. D^3 products per matrix product; this route
         serves any fixed linear map, translation-invariant or not.
  ntt  : when D divides p - 1, the integer NTT of size D diagonalizes M with
         eigenvalues lambda_m = sum_s c_s w^(s*m) (w a primitive D-th root of
         unity mod p). The jump raises each eigenvalue to n: about D powers of
         at most 39 products each, so the work stays flat as k grows.
Leaves: one step from the jump equals the jump to n + 1, and the inverse
eigenvalues give the state at n - 1, which steps onto the jump. A zero
eigenvalue refuses the backward leaf: two different starts then share every
later state (JUMP3_ZERO_MODE rows). The SHA-256 of each jumped state is its
fingerprint for later comparison.

p = 1000081, the kernel prime, with p - 1 = 2^4 * 3^3 * 5 * 463. mulmod
fields are the modular products the code performed, counted as it runs
(direct=NOT_RUN rows give n*D*S, the products those steps would take).
Wall-clock time appears only under --timing, outside the byte comparison.

usage: jump3_ntt.py --self-test [--timing]
       jump3_ntt.py --ring D --stencil s:c,s:c --ks k,k [--direct-max K] [--case NAME] [--step TEXT]
Output: HBP rows, json=0. crates/jump3 (binary jump3) prints the same bytes.
"""
from hashlib import sha256
from time import monotonic_ns
import argparse
import sys

P = 1000081
K_MAX = 80                     # 3^80 < 2^127, the Rust u128 range


class Ops:
    """Counter of modular products performed."""

    def __init__(self):
        self.n = 0


def factor(n):
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            out.append((d, e))
        d += 1
    if n > 1:
        out.append((n, 1))
    return out


def factor_text(n):
    return "*".join("%d^%d" % (q, e) if e > 1 else "%d" % q for q, e in factor(n))


def generator(p=P):
    """Smallest primitive root mod p."""
    qs = [q for q, _ in factor(p - 1)]
    g = 2
    while any(pow(g, (p - 1) // q, p) == 1 for q in qs):
        g += 1
    return g


def root_of_unity(d, p=P):
    """Primitive d-th root of unity mod p, or None when d does not divide p - 1."""
    if d < 2 or (p - 1) % d:
        return None
    w = pow(generator(p), (p - 1) // d, p)
    if pow(w, d, p) != 1 or any(pow(w, d // q, p) == 1 for q, _ in factor(d)):
        raise AssertionError("root of unity check failed")
    return w


def powmod(b, e, ops):
    """Right-to-left binary power mod P, every product counted."""
    r = 1
    b %= P
    while e > 0:
        if e & 1:
            r = r * b % P
            ops.n += 1
        e >>= 1
        if e > 0:
            b = b * b % P
            ops.n += 1
    return r


def lam_power(lam, n, ops):
    """lambda^n mod P for n >= 0; a zero eigenvalue gives 0 for every n >= 1."""
    if lam == 0:
        return 1 if n == 0 else 0
    return powmod(lam, n % (P - 1), ops)


class Ring:
    """A step that is the same at every site of a ring of d sites."""

    def __init__(self, d, stencil):
        if d < 2 or not stencil:
            raise ValueError("ring of at least 2 sites and a nonempty stencil required")
        self.d = d
        self.signed = [(s, c) for s, c in stencil]
        self.stencil = [(s % d, c % P) for s, c in stencil]
        self.w = root_of_unity(d)
        if self.w is not None:
            self.tw = [1] * d
            for e in range(1, d):
                self.tw[e] = self.tw[e - 1] * self.w % P
            self.tw_inv = [self.tw[(d - e) % d] for e in range(d)]
            self.inv_d = pow(d, P - 2, P)

    def stencil_text(self):
        return ",".join("%d:%d" % (s, c) for s, c in self.signed)


def start_state(d):
    return [(7 * i * i + 3 * i + 1) % P for i in range(d)]


def state_sha(x):
    """SHA-256 of the components as 8-byte big-endian integers (the seat's encoding)."""
    return sha256(b"".join(v.to_bytes(8, "big") for v in x)).hexdigest()


def step(ring, x, ops):
    d = ring.d
    acc = [0] * d
    for s, c in ring.stencil:
        xs = x[d - s:] + x[:d - s]           # xs[i] = x[(i - s) mod d]
        acc = [a + c * v for a, v in zip(acc, xs)]
    ops.n += d * len(ring.stencil)
    return [a % P for a in acc]


def circulant(ring):
    d = ring.d
    m = [[0] * d for _ in range(d)]
    for i in range(d):
        for s, c in ring.stencil:
            j = (i - s) % d
            m[i][j] = (m[i][j] + c) % P
    return m


def matmul(a, b, ops):
    d = len(a)
    bt = list(zip(*b))
    ops.n += d * d * d
    return [[sum(u * v for u, v in zip(row, col)) % P for col in bt] for row in a]


def matvec(a, x, ops):
    ops.n += len(a) * len(x)
    return [sum(u * v for u, v in zip(row, x)) % P for row in a]


def cube(a, ops):
    """a^3 mod P: two matrix products."""
    return matmul(matmul(a, a, ops), a, ops)


def cube_ladder(m, kmax, ops):
    """[(m^(3^r), products so far) for r = 0..kmax]: one cache of cubings serves every k."""
    out = [(m, ops.n)]
    for _ in range(kmax):
        m = cube(m, ops)
        out.append((m, ops.n))
    return out


def ntt(x, tw, ops):
    d = len(x)
    ops.n += d * d
    return [sum(x[j] * tw[j * m % d] for j in range(d)) % P for m in range(d)]


def eigenvalues(ring, ops):
    d = ring.d
    ops.n += d * len(ring.stencil)
    return [sum(c * ring.tw[s * m % d] for s, c in ring.stencil) % P for m in range(d)]


def ntt_apply(ring, x, factors, ops):
    x_hat = ntt(x, ring.tw, ops)
    y_hat = [a * f % P for a, f in zip(x_hat, factors)]
    ops.n += ring.d
    y = ntt(y_hat, ring.tw_inv, ops)
    ops.n += ring.d
    return [v * ring.inv_d % P for v in y]


def ntt_jump(ring, lam, x, n, ops):
    """M^n x through the eigenvalues: D powers, two transforms."""
    return ntt_apply(ring, x, [lam_power(v, n, ops) for v in lam], ops)


def ntt_unstep(ring, lam, x, ops):
    """M^-1 x through the inverse eigenvalues; None when an eigenvalue is zero."""
    if any(v == 0 for v in lam):
        return None
    return ntt_apply(ring, x, [powmod(v, P - 2, ops) for v in lam], ops)


def mark(v):
    return "1" if v else "0"


def run_case(name, d, stencil, label, ks, direct_max, out, flags):
    """Append the JUMP3_CASE, JUMP3 and JUMP3_ZERO_MODE rows of one case to out."""
    if any(k < 1 or k > K_MAX for k in ks):
        raise ValueError("k must lie in 1..%d" % K_MAX)
    ring = Ring(d, stencil)
    x0 = start_state(d)
    s_len = len(ring.stencil)
    ops = Ops()
    has_ntt = ring.w is not None
    if has_ntt:
        lam = eigenvalues(ring, ops)
        zeros = [m for m in range(d) if lam[m] == 0]
        root, ntt_state = str(ring.w), "DEFINED"
        zm, zl = str(len(zeros)), ",".join(str(m) for m in zeros) or "none"
        backward = "REFUSED_ZERO_MODE" if zeros else "DEFINED"
    else:
        lam, zeros = None, []
        root, ntt_state = "none", "REFUSED_RING_DOES_NOT_DIVIDE_P_MINUS_1"
        zm, zl, backward = "NOT_DEFINED", "none", "NOT_DEFINED"
    out.append("JUMP3_CASE|case=%s|ring=%d|step=%s|stencil=%s|root=%s|ntt=%s|eigen_mulmods=%d|"
               "zero_modes=%s|zero_mode_list=%s|backward=%s|start_sha256=%s|json=0"
               % (name, d, label, ring.stencil_text(), root, ntt_state, ops.n, zm, zl, backward,
                  state_sha(x0)))
    # direct stepping, one step at a time, to 3^k and 3^k + 1 for every k <= direct_max
    targets = set()
    for k in ks:
        if k <= direct_max:
            targets.add(3 ** k)
            targets.add(3 ** k + 1)
    direct = {}
    if targets:
        dops = Ops()
        x = x0
        for t in range(1, max(targets) + 1):
            x = step(ring, x, dops)
            if t in targets:
                direct[t] = (x, dops.n)
    ladder = cube_ladder(circulant(ring), max(ks), Ops())
    scratch = Ops()
    for k in ks:
        n = 3 ** k
        mk, lad = ladder[k]
        cops = Ops()
        cops.n = lad
        xc = matvec(mk, x0, cops)
        if has_ntt:
            nops = Ops()
            xn = ntt_jump(ring, lam, x0, n, nops)
            ntt_mul, c_eq = str(nops.n), mark(xc == xn)
            flags["cube_ntt"] &= xc == xn
        else:
            ntt_mul, c_eq = "NOT_DEFINED", "NOT_DEFINED"
        if n in direct:
            d_run, d_mul, d_eq = "RUN", str(direct[n][1]), mark(direct[n][0] == xc)
            flags["direct"] &= direct[n][0] == xc
        else:
            d_run, d_mul, d_eq = "NOT_RUN", str(n * d * s_len), "NOT_RUN"
        fwd = step(ring, xc, scratch)
        if has_ntt:
            ok = fwd == ntt_jump(ring, lam, x0, n + 1, scratch)
            leaf_f = mark(ok)
            flags["leaf_f"] &= ok
        elif n + 1 in direct:
            ok = fwd == direct[n + 1][0]
            leaf_f = mark(ok)
            flags["leaf_f"] &= ok
        else:
            leaf_f = "NOT_DEFINED"
        if has_ntt and not zeros:
            b = ntt_unstep(ring, lam, xc, scratch)
            ok = (step(ring, b, scratch) == xc and b == ntt_jump(ring, lam, x0, n - 1, scratch)
                  and ntt_unstep(ring, lam, step(ring, x0, scratch), scratch) == x0)
            leaf_b = mark(ok)
            flags["leaf_b"] &= ok
        elif has_ntt:
            leaf_b = "REFUSED_ZERO_MODE"
        else:
            leaf_b = "NOT_DEFINED"
        at_start = xc == x0
        if name == "shift":
            flags["shift"] &= at_start == (k >= 3)
        flags["rows"] += 1
        out.append("JUMP3|case=%s|k=%d|n=%d|cubings=%d|cube_mulmods=%d|ntt_mulmods=%s|cube_equals_ntt=%s|"
                   "direct=%s|direct_mulmods=%s|direct_equals_jump=%s|equals_start=%s|leaf_forward=%s|"
                   "leaf_back=%s|sha256=%s|json=0"
                   % (name, k, n, k, cops.n, ntt_mul, c_eq, d_run, d_mul, d_eq, mark(at_start), leaf_f,
                      leaf_b, state_sha(xc)))
    if has_ntt and zeros:
        k = max(ks)
        n = 3 ** k
        mk = ladder[k][0]
        jump_k = matvec(mk, x0, scratch)
        for m in zeros:
            v = [ring.tw_inv[j * m % d] for j in range(d)]          # eigenvector w^(-j*m), eigenvalue 0
            xb = [(a + b) % P for a, b in zip(x0, v)]
            image_b = matvec(mk, xb, scratch)
            ntt_a = ntt_jump(ring, lam, x0, n, scratch)
            checks = (xb != x0,
                      step(ring, x0, scratch) == step(ring, xb, scratch),
                      jump_k == image_b,
                      ntt_a == ntt_jump(ring, lam, xb, n, scratch),
                      image_b == ntt_a)
            for c in checks:
                flags["zero"] &= c
            flags["zero_rows"] += 1
            out.append("JUMP3_ZERO_MODE|case=%s|mode=%d|k=%d|preimages_distinct=%s|images_equal_after_1_step=%s|"
                       "images_equal_after_3^%d_cube=%s|images_equal_after_3^%d_ntt=%s|image_equals_case_jump=%s|"
                       "preimage_b_sha256=%s|image_sha256=%s|json=0"
                       % (name, m, k, mark(checks[0]), mark(checks[1]), k, mark(checks[2]), k, mark(checks[3]),
                          mark(checks[4]), state_sha(xb), state_sha(image_b)))


CASES = (
    ("signed_return", 27, ((0, 1), (1, 1), (2, -1)), "I+R-R^2", (1, 2, 3, 6, 9, 10, 20, 40, 80), 10),
    ("curved_sum", 27, ((0, 1), (1, 1), (2, 1)), "I+R+R^2", (1, 2, 3, 9, 40), 9),
    ("shift", 27, ((1, 1),), "R", (1, 2, 3, 40), 3),
    ("ring45", 45, ((0, 2), (1, -1), (44, 3)), "2I-R+3R^-1", (1, 2, 3, 6, 9, 20), 9),
    ("ring81", 81, ((0, 1), (1, 1), (80, -1)), "I+R-R^-1", (1, 2, 3), 3),
)


def new_flags():
    return {"cube_ntt": True, "direct": True, "leaf_f": True, "leaf_b": True, "zero": True, "shift": True,
            "rows": 0, "zero_rows": 0}


def head_row():
    return ("JUMP3_HEAD|law=power_of_three_jump|p=%d|p_minus_1_factors=%s|generator=%d|"
            "start=7i^2+3i+1_mod_p|state_hash=sha256_of_8_byte_big_endian_components|json=0"
            % (P, factor_text(P - 1), generator()))


def verify_row(flags, cases):
    ok = all(flags[f] for f in ("cube_ntt", "direct", "leaf_f", "leaf_b", "zero", "shift"))
    return ("JUMP3_VERIFY|status=%s|cases=%d|jump_rows=%d|cube_equals_ntt_all=%s|direct_equals_jump_all=%s|"
            "leaf_forward_all=%s|leaf_back_all=%s|zero_mode_rows=%d|zero_mode_all=%s|shift_returns_at_27=%s|json=0"
            % ("PASS" if ok else "FAIL", cases, flags["rows"], mark(flags["cube_ntt"]), mark(flags["direct"]),
               mark(flags["leaf_f"]), mark(flags["leaf_b"]), flags["zero_rows"], mark(flags["zero"]),
               mark(flags["shift"]))), ok


def timing_rows():
    """Wall-clock nanoseconds for the signed_return case; outside the byte comparison."""
    name, d, stencil, label, ks, direct_max = CASES[0]
    ring = Ring(d, stencil)
    lam = eigenvalues(ring, Ops())
    x0 = start_state(d)
    rows = []
    t0 = monotonic_ns()
    x = x0
    for _ in range(3 ** direct_max):
        x = step(ring, x, Ops())
    direct_ns = monotonic_ns() - t0
    for k in ks:
        n = 3 ** k
        t0 = monotonic_ns()
        mk = cube_ladder(circulant(ring), k, Ops())[k][0]
        matvec(mk, x0, Ops())
        cube_ns = monotonic_ns() - t0
        t0 = monotonic_ns()
        ntt_jump(ring, lam, x0, n, Ops())
        ntt_ns = monotonic_ns() - t0
        need = direct_ns * n // 3 ** direct_max
        rows.append("JUMP3_TIMING|lane=python|case=%s|k=%d|n=%d|cube_ns=%d|ntt_ns=%d|direct_ns_for_3^%d=%d|"
                    "direct_ns_needed=%d|direct_years_needed=%d|basis=direct_rate_at_3^%d_steps|json=0"
                    % (name, k, n, cube_ns, ntt_ns, direct_max, direct_ns, need,
                       need // (365 * 24 * 3600 * 10 ** 9), direct_max))
    return rows


def self_test(timing):
    out = [head_row()]
    flags = new_flags()
    for case in CASES:
        run_case(*case, out=out, flags=flags)
    row, ok = verify_row(flags, len(CASES))
    out.append(row)
    if timing:
        out.extend(timing_rows())
    sys.stdout.write("".join(r + "\n" for r in out))
    return 0 if ok else 1


def parse_stencil(text):
    out = []
    for part in text.split(","):
        s, c = part.split(":")
        out.append((int(s), int(c)))
    return out


def main():
    ap = argparse.ArgumentParser(description="power-of-three jump, exact integers mod 1000081")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--timing", action="store_true")
    ap.add_argument("--ring", type=int)
    ap.add_argument("--stencil")
    ap.add_argument("--ks")
    ap.add_argument("--direct-max", type=int, default=0)
    ap.add_argument("--case", default="custom")
    ap.add_argument("--step", default="stencil")
    a = ap.parse_args()
    if a.self_test:
        return self_test(a.timing)
    if a.ring is None or a.stencil is None or a.ks is None:
        ap.print_help()
        return 2
    out = [head_row()]
    flags = new_flags()
    run_case(a.case, a.ring, parse_stencil(a.stencil), a.step, [int(k) for k in a.ks.split(",")],
             a.direct_max, out, flags)
    row, ok = verify_row(flags, 1)
    out.append(row)
    sys.stdout.write("".join(r + "\n" for r in out))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
