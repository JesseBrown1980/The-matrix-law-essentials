"""Exact min-plus ternary jump; no floats, no external dependencies."""
from hashlib import sha256
from random import Random


def compose(a, b):
    n = len(a)
    out = [[None]*n for _ in range(n)]
    for i in range(n):
        for k, u in enumerate(a[i]):
            if u is None:
                continue
            for j, v in enumerate(b[k]):
                if v is not None:
                    w = u+v
                    if out[i][j] is None or w < out[i][j]:
                        out[i][j] = w
    return out


def cube(a):
    return compose(compose(a, a), a)


def reference(a, steps, reverse=False):
    # Independent one-source Bellman recurrence, optionally backward in time.
    n = len(a)
    rows = []
    for s in range(n):
        d = [None]*n
        d[s] = 0
        for _ in range(steps):
            nxt = []
            for j in range(n):
                vals = [d[i]+(a[j][i] if reverse else a[i][j])
                        for i in range(n) if d[i] is not None and
                        (a[j][i] if reverse else a[i][j]) is not None]
                nxt.append(min(vals) if vals else None)
            d = nxt
        rows.append(d)
    return [list(r) for r in zip(*rows)] if reverse else rows


def main():
    rng = Random(27181)
    count = 0
    digest = sha256()
    for n in (3, 9, 27):
        for case in range(12):
            a = [[None if rng.randrange(4)==0 else rng.randrange(-9, 28)
                  for _ in range(n)] for _ in range(n)]
            if case == 0:
                a[-1] = [None]*n  # unreachable paths
            k = a
            for steps in (1, 3, 9, 27, 81):
                assert k == reference(a, steps) == reference(a, steps, True)
                count += n*n
                digest.update(repr(k).encode('ascii'))
                if steps != 81:
                    k = cube(k)
    # Greedy first-step trap, forbidden paths, noncommuting blocks.
    a = [[None, 0, 1], [None, 9, None], [None, None, 0]]
    assert compose(a, a)[0][2] == 1
    assert compose(a, a)[0][0] is None
    b = [[0, None, None], [None, None, 0], [None, None, None]]
    assert compose(a, b) != compose(b, a)
    # Uniform cost graph: exact huge horizon by induction; no enumeration.
    k = [[7]*3 for _ in range(3)]
    for _ in range(40):
        k = cube(k)
    assert k == [[7*3**40]*3 for _ in range(3)]
    print('PASS endpoint_horizon_cases=%d forward_backward_comparisons=%d' % (count, 2*count))
    print('PASS controls=4 horizon=3^40 sha256='+digest.hexdigest())
    # Dense 27-state full-kernel arithmetic counts, not wall-clock estimates.
    print('dense_81_frames candidates_direct=%d candidates_jump=%d' % (80*27**3, 8*27**3))


if __name__ == '__main__':
    main()
