#!/usr/bin/env python3
"""Independent period exponent and monotone-bound audit for nonlinear jumps.

All computations are integers, no solver imports, no timing/float claims.
The bound is a universal exponent, not necessarily the minimal orbit period.
"""
from math import gcd, isqrt

P = 1_299_079
D = 59_049
Q = 1_000_000


def matmul(a, b):
    return (
        (a[0] * b[0] + a[1] * b[2]) % P,
        (a[0] * b[1] + a[1] * b[3]) % P,
        (a[2] * b[0] + a[3] * b[2]) % P,
        (a[2] * b[1] + a[3] * b[3]) % P,
    )


def matpow(x, n):
    identity = (1, 0, 0, 1)
    result = identity
    while n:
        if n & 1:
            result = matmul(result, x)
        x = matmul(x, x)
        n >>= 1
    return result


def period_check():
    # Independent proof that the modulus really is a field prime.
    for divisor in range(2, isqrt(P) + 1):
        assert P % divisor != 0, ("composite modulus", divisor)
    assert P % 4 == 3
    assert (P - 1) % D == 0
    period_exponent = (P * P - 1) // 2
    assert gcd(P - 1, P + 1) == 2
    assert period_exponent == 843_803_124_120
    assert pow(3, 80) % period_exponent == 444_653_673_201
    w = 1
    r = pow(3, (P - 1) // D, P)
    assert pow(r, D, P) == 1
    assert pow(r, D // 3, P) != 1
    split = quadratic = defective = 0
    for index in range(D):
        w2 = w * w % P
        trace = (2 + w2) % P
        discriminant = (trace * trace - 4) % P
        assert discriminant == w2 * ((w2 + 4) % P) % P
        assert w != 0
        if discriminant == 0:
            defective += 1
        else:
            symbol = pow(discriminant, (P - 1) // 2, P)
            if symbol == 1:
                split += 1
            else:
                assert symbol == P - 1
                quadratic += 1
        matrix = ((1 + w2) % P, w, w, 1)
        assert (matrix[0] * matrix[3] - matrix[1] * matrix[2]) % P == 1
        assert matpow(matrix, period_exponent) == (1, 0, 0, 1), index
        w = w * r % P
    assert w == 1
    assert (split, quadratic, defective) == (29605, 29444, 0)
    # A one-step exponent sabotage must be refused (already at the uniform cat map).
    cat = (2, 1, 1, 1)
    assert matpow(cat, period_exponent) == (1, 0, 0, 1)
    assert matpow(cat, period_exponent - 1) != (1, 0, 0, 1)
    print(
        "NONLINEAR_PERIOD|d=59049|field=1299079|"
        "split_modes=29605|quadratic_extension_modes=29444|"
        "defective_modes=0|all_59049_mode_matrices_at_exponent_are_identity=PASS|"
        "universal_exponent=843803124120|"
        "3pow80_remainder=444653673201|"
        "minimal_period=NOT_CLAIMED|json=0"
    )


def scalar(b):
    return (4 * b + (b * b) // Q) // 8


def monotonicity_check():
    assert scalar(0) == 0
    previous = scalar(0)
    for b in range(1, Q + 1):
        current = scalar(b)
        assert previous <= current, (b, previous, current)
        assert current < b, (b, current)
        if b == Q:
            assert (10 * b + (b * b) // Q) // 8 > b
        previous = current
    seq = [Q]
    while seq[-1] > 0:
        seq.append(scalar(seq[-1]))
    assert len(seq) == 22 and seq[-1] == 0
    expected = [1000000,625000,361328,196983,103341,53005,26853,
                13516,6780,3395,1698,849,424,212,106,53,26,13,6,3,1,0]
    assert seq == expected
    print(
        "NONLINEAR_MONOTONICITY|q=1000000|"
        "all_b_from_0_to_1000000_tested=YES|"
        "nondecreasing_comparison_map=PASS|strict_for_positive_b=PASS|"
        "top_trajectory_reaches_zero_at=21|"
        "all_smaller_initial_supnorms_covered=PASS|json=0"
    )


if __name__ == "__main__":
    period_check()
    monotonicity_check()
    print("NONLINEAR_SCOPE_AUDIT|status=PASS|r9_launch_not_invoked=YES|json=0")
