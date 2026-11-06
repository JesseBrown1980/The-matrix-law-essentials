# The jump to step 3^k and its hash, three independent ways, exact integers mod p = 1000081 (the kernel's prime; 27 | p-1).
# State: a 27-site ring. One step: x <- M x with M = I + R - R^2 (the signed curved return j' of curve_ABC), R the cyclic shift.
# (a) step 3^k times; (b) cube M k times; (c) the integer NTT of size 27: eigenvalues raised to 3^k. Then the leaves at 3^k - 1 and 3^k + 1.
import hashlib
import time

P = 1000081
D = 27
C = {0: 1, 1: 1, 2: P - 1}          # M = I + R - R^2: x'[i] = x[i] + x[i-1] - x[i-2]


def step(x):
    return [(x[i] + x[i - 1] - x[i - 2]) % P for i in range(D)]


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(D)) % P for j in range(D)] for i in range(D)]


def matvec(A, x):
    return [sum(A[i][k] * x[k] for k in range(D)) % P for i in range(D)]


def h(x):
    return hashlib.sha256(b"".join(v.to_bytes(8, "big") for v in x)).hexdigest()


M = [[0] * D for _ in range(D)]
for i in range(D):
    for s, c in C.items():
        M[i][(i - s) % D] = (M[i][(i - s) % D] + c) % P

# a primitive 27th root of unity mod P
fac = [2, 3, 5, 463]
g = next(a for a in range(2, P) if all(pow(a, (P - 1) // q, P) != 1 for q in fac))
w = pow(g, (P - 1) // 27, P)
assert pow(w, 27, P) == 1 and pow(w, 9, P) != 1
inv27 = pow(27, P - 2, P)


def ntt(x, root):
    return [sum(x[j] * pow(root, j * m, P) for j in range(D)) % P for m in range(D)]


lam = [sum(c * pow(w, s * m, P) for s, c in C.items()) % P for m in range(D)]   # eigenvalue of mode m
invertible = all(v != 0 for v in lam)


def jump_ntt(x, n):
    X = ntt(x, w)
    Y = [X[m] * pow(lam[m], n % (P - 1), P) % P for m in range(D)] if invertible else None
    winv = pow(w, P - 2, P)
    return [v * inv27 % P for v in ntt(Y, winv)]


x0 = [(7 * i * i + 3 * i + 1) % P for i in range(D)]
print("JUMP3_HEAD|p=%d|sites=%d|M=I+R-R^2|eigenvalues_nonzero=%d|root27=%d|json=0" % (P, D, int(invertible), w))
for k in (3, 6, 9, 10, 20, 40):
    n = 3 ** k
    t0 = time.perf_counter_ns()
    Mk = M
    for _ in range(k):
        Mk = matmul(matmul(Mk, Mk), Mk)
    xb = matvec(Mk, x0)
    tb = time.perf_counter_ns() - t0
    t0 = time.perf_counter_ns()
    xc = jump_ntt(x0, n)
    tc = time.perf_counter_ns() - t0
    row = "JUMP3|k=%d|n=3^%d=%d|cubings=%d|cube_ms=%d|ntt_ms=%d|cube_equals_ntt=%d|hash=%s" % (k, k, n, k, tb // 10 ** 6, tc // 10 ** 6, int(xb == xc), h(xc)[:16])
    if k <= 10:
        t0 = time.perf_counter_ns()
        xa = x0
        for _ in range(n):
            xa = step(xa)
        ta = time.perf_counter_ns() - t0
        row += "|steps=%d|step_ms=%d|step_equals_jump=%d" % (n, ta // 10 ** 6, int(xa == xc))
    else:
        row += "|steps_needed=%d|step_ms=NOT_RUN" % n
    fwd = step(xc)
    back = jump_ntt(x0, n - 1)
    row += "|leaf_forward_equals_jump_n+1=%d|leaf_back_steps_to_jump=%d|json=0" % (int(fwd == jump_ntt(x0, n + 1)), int(step(back) == xc))
    print(row)
