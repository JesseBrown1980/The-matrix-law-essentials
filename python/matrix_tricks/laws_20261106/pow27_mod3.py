# The two 27-step tests of the mod-3 collapse on a 27-site ring: (I + R + R^2)^27 = 0 and (I + R - R^2)^27 = I (mod 3). json=0.
D = 27


def circ(coef):
    M = [[0] * D for _ in range(D)]
    for i in range(D):
        for s, c in coef.items():
            M[i][(i - s) % D] = (M[i][(i - s) % D] + c) % 3
    return M


def mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(D)) % 3 for j in range(D)] for i in range(D)]


ident = [[int(i == j) for j in range(D)] for i in range(D)]
zero = [[0] * D for _ in range(D)]
for name, coef in (("I+R+R^2", {0: 1, 1: 1, 2: 1}), ("I+R-R^2", {0: 1, 1: 1, 2: 2})):
    P = circ(coef)
    for _ in range(3):
        P = mul(mul(P, P), P)
    print("POW27_MOD3|step=%s|ring=27|equals_zero=%d|equals_identity=%d|json=0" % (name, int(P == zero), int(P == ident)))
