#!/usr/bin/env python3
"""Hash-keyed position memo on the 27-state cube, with the path-count recurrence and the invariant.

Positions are (a, b, c) in {0,1,2}^3, start 000; one turn adds 1 mod 3 to one
digit. Many move orders reach one position (transpositions, as in a chess
engine's table). Computing each position once and keying it by the SHA-256 of
its three digit bytes reuses every repeated arrival.

Counts per position differ (63, 81 or 90 at L = 6), so anything that depends
on how many orders reach a position uses the recurrence

    C_{n+1}(s) = sum over turns p -> s of C_n(p),   C_0 = [start],

which costs 27 * 3 updates per turn instead of 3^L paths. The digit sum mod 3
equals the number of turns mod 3, so after L turns the candidates are the 9
positions of one class (27 without the invariant), and the recurrence needs
only the 9 positions of the current class per turn.

The recurrence is one fixed linear map A on the 27 counts, so the counts after
L = 3^k turns come from k cubings of A (jump3_ntt.cube_ladder), exact mod p.
For L divisible by 3 a position with three different digits is reached by
exactly 3^(L-2) move orders (the character sum's sqrt(3)-size terms cancel).

A hash identifies a position only inside a space that is enumerated: the
hidden position of FINITE_RECOVER is found by hashing candidates in order,
with 27 candidates as the bound (9 with the invariant). An input of 64
unstructured bits costs 2^64 candidates by the same search.

usage: transposition_memo.py --self-test [--timing]
Output: HBP rows, json=0. crates/jump3 (binary transmemo) prints the same bytes.
"""
from hashlib import sha256
from time import monotonic_ns
import argparse
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jump3_ntt as J3  # noqa: E402

POS = [tuple(t) for t in itertools.product(range(3), repeat=3)]      # index a*9 + b*3 + c
INDEX = {p: i for i, p in enumerate(POS)}
KEY = [sha256(bytes(p)).hexdigest() for p in POS]
SLOT = {k: i for i, k in enumerate(KEY)}


def name(i):
    return "%d%d%d" % POS[i]


def turn(i, d):
    p = list(POS[i])
    p[d] = (p[d] + 1) % 3
    return INDEX[tuple(p)]


def digit_sum(i):
    return sum(POS[i])


def enumerate_orders(L):
    """Every one of the 3^L move orders; the end position looked up by its hash."""
    table = {}
    turns = 0
    for seq in itertools.product(range(3), repeat=L):
        p = [0, 0, 0]
        for d in seq:
            p[d] = (p[d] + 1) % 3
        turns += L
        key = sha256(bytes(p)).hexdigest()
        if key in table:
            table[key] += 1
        else:
            table[key] = 1
    counts = [0] * 27
    for key, c in table.items():
        counts[SLOT[key]] = c
    table_text = "".join("%s=%d\n" % (key, table[key]) for key in sorted(table))
    return counts, table, turns, sha256(table_text.encode("ascii")).hexdigest()


def recurrence(L, modulus=0):
    """C_L by C_{n+1}(s) = sum over p -> s of C_n(p) over all 27 positions; exact, or mod modulus."""
    c = [0] * 27
    c[0] = 1
    updates = 0
    for _ in range(L):
        nxt = [0] * 27
        for q in range(27):
            for d in range(3):
                nxt[turn(q, d)] += c[q]
                updates += 1
        c = [v % modulus for v in nxt] if modulus else nxt
    return c, updates


def recurrence_invariant(L):
    """The same recurrence over the 9 positions of the current class only."""
    c = [0] * 27
    c[0] = 1
    updates = 0
    for n in range(L):
        nxt = [0] * 27
        for q in range(27):
            if digit_sum(q) % 3 != n % 3:
                continue
            for d in range(3):
                nxt[turn(q, d)] += c[q]
                updates += 1
        c = nxt
    return c, updates


def adjacency():
    """A[s][q] = number of turns taking q to s, mod p."""
    a = [[0] * 27 for _ in range(27)]
    for q in range(27):
        for d in range(3):
            a[turn(q, d)][q] += 1
    return a


def mark(v):
    return "1" if v else "0"


def counts_text(counts):
    return ",".join("%s:%d" % (name(i), counts[i]) for i in range(27) if counts[i])


def exact_sha(counts):
    return sha256(",".join(str(v) for v in counts).encode("ascii")).hexdigest()


PERM_012 = [INDEX[p] for p in itertools.permutations((0, 1, 2))]


def self_test(timing):
    out = ["TRANSMEMO_HEAD|cube=3x3x3|positions=27|turn=add_1_mod_3_to_one_digit|start=000|"
           "position_key=sha256_of_3_digit_bytes|invariant=digit_sum_mod_3_equals_turns_mod_3|p=%d|json=0" % J3.P]
    ok = True
    rows_t = 0
    for L in (6, 9, 12):
        counts, table, turns, table_sha = enumerate_orders(L)
        paths = 3 ** L
        rec, upd = recurrence(L)
        inv, upd_inv = recurrence_invariant(L)
        cls = [i for i in range(27) if digit_sum(i) % 3 == L % 3]
        reached = [i for i in range(27) if counts[i]]
        nz = [v for v in counts if v]
        checks = (rec == counts, reached == cls, inv == rec, sum(counts) == paths)
        ok &= all(checks)
        rows_t += 1
        out.append("TRANSPOSITION|L=%d|paths=%d|enumeration_turns=%d|distinct_positions=%d|memo_keys=%d|"
                   "memo_hits=%d|reuse_avg=%d|paths_min=%d|paths_max=%d|counts=%s|table_sha256=%s|"
                   "recurrence_updates=%d|recurrence_equals_enumeration=%s|invariant_class=%d|"
                   "invariant_candidates=%d|reachable_equals_invariant_class=%s|invariant_updates=%d|"
                   "invariant_equals_recurrence=%s|json=0"
                   % (L, paths, turns, len(reached), len(table), paths - len(table), paths // len(table),
                      min(nz), max(nz), counts_text(counts), table_sha, upd, mark(checks[0]), L % 3, len(cls),
                      mark(checks[1]), upd_inv, mark(checks[2])))
    rows_p = 0
    for L in (27, 54):
        c, _ = recurrence(L)
        perm = [c[i] for i in PERM_012]
        checks = (sum(c) == 3 ** L, perm[0] == 3 ** (L - 2), len(set(perm)) == 1)
        ok &= all(checks) and len([v for v in c if v]) == 9
        rows_p += 1
        out.append("PATHCOUNT|L=%d|route=recurrence_exact|sum=%d|sum_equals_3^L=%s|support=%d|count_000=%d|"
                   "count_111=%d|count_222=%d|count_012=%d|count_012_equals_3^(L-2)=%s|perm_012_all_equal=%s|"
                   "counts_sha256=%s|json=0"
                   % (L, sum(c), mark(checks[0]), len([v for v in c if v]), c[INDEX[(0, 0, 0)]],
                      c[INDEX[(1, 1, 1)]], c[INDEX[(2, 2, 2)]], perm[0], mark(checks[1]), mark(checks[2]),
                      exact_sha(c)))
    ladder = J3.cube_ladder(adjacency(), 40, J3.Ops())
    e0 = [1] + [0] * 26
    rows_j = 0
    for k in (3, 4, 40):
        L = 3 ** k
        mk, lad = ladder[k]
        ops = J3.Ops()
        ops.n = lad
        c = J3.matvec(mk, e0, ops)
        if k <= 4:
            rec, _ = recurrence(L, J3.P)
            eq = mark(rec == c)
            ok &= rec == c
        else:
            eq = "NOT_RUN"
        s = sum(c) % J3.P
        want_sum = pow(3, L % (J3.P - 1), J3.P)
        want_012 = pow(3, (L - 2) % (J3.P - 1), J3.P)
        perm_ok = all(c[i] == want_012 for i in PERM_012)
        ok &= s == want_sum and perm_ok
        rows_j += 1
        out.append("PATHCOUNT_JUMP|k=%d|L=%d|route=cube_mod_p|cubings=%d|cube_mulmods=%d|equals_recurrence_mod_p=%s|"
                   "support_mod_p=%d|sum_mod_p=%d|sum_equals_3^L_mod_p=%s|count_012_mod_p=%d|"
                   "count_012_equals_3^(L-2)_mod_p=%s|counts_sha256=%s|json=0"
                   % (k, L, k, ops.n, eq, len([v for v in c if v]), s, mark(s == want_sum),
                      c[INDEX[(0, 1, 2)]], mark(perm_ok), J3.state_sha(c)))
    hidden = (0, 2, 1, 1, 0, 2, 2, 1, 0, 1)
    p = 0
    for d in hidden:
        p = turn(p, d)
    target = KEY[p]
    tried_full = 0
    found = None
    for i in range(27):
        tried_full += 1
        if sha256(bytes(POS[i])).hexdigest() == target:
            found = i
            break
    tried_inv = 0
    cls = [i for i in range(27) if digit_sum(i) % 3 == len(hidden) % 3]
    found_inv = None
    for i in cls:
        tried_inv += 1
        if sha256(bytes(POS[i])).hexdigest() == target:
            found_inv = i
            break
    rec_ok = found == p and found_inv == p
    ok &= rec_ok
    out.append("FINITE_RECOVER|moves_hidden=%d|hidden_moves=%s|position=%s|bound_full=27|tried_full=%d|"
               "bound_invariant=%d|tried_invariant=%d|recovered_equals_hidden=%s|json=0"
               % (len(hidden), ",".join(str(d) for d in hidden), name(p), tried_full, len(cls), tried_inv,
                  mark(rec_ok)))
    out.append("UNSTRUCTURED|input_bits=64|candidates=%d|structure=none|cost=one_sha256_per_candidate|json=0"
               % (2 ** 64))
    out.append("TRANSMEMO_VERIFY|status=%s|transposition_rows=%d|pathcount_rows=%d|pathcount_jump_rows=%d|json=0"
               % ("PASS" if ok else "FAIL", rows_t, rows_p, rows_j))
    if timing:
        n = 1 << 20
        t0 = monotonic_ns()
        for i in range(n):
            sha256(i.to_bytes(8, "big")).digest()
        ns = monotonic_ns() - t0
        rate = n * 10 ** 9 // ns
        out.append("UNSTRUCTURED_TIMING|lane=python|hashes=%d|ns=%d|hashes_per_second=%d|years_for_2^64=%d|json=0"
                   % (n, ns, rate, 2 ** 64 // rate // (365 * 24 * 3600)))
    sys.stdout.write("".join(r + "\n" for r in out))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description="hash-keyed transposition memo on the 27-state cube")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--timing", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test(a.timing)
    ap.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
