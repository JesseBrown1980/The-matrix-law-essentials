# When does knowing a hash give the answer? Three measured cases, exact integers, HBP rows (json=0).
# 1. A finite state space: recover a hidden state from its SHA-256 alone by hashing every candidate.
# 2. Transpositions (the chess-engine table): many move orders reach one position; compute each position once, look it up by hash.
# 3. An unstructured 64-bit input: the hash alone, with no structure, costs a search over the whole space.
import hashlib
import itertools
import time


def h(b):
    return hashlib.sha256(b).hexdigest()


def state_bytes(s):
    return bytes(s)


# 1. The 27-state ternary cube of the YM repo's curve_ABC A test: (a,b,c) in {0,1,2}^3, a turn adds 1 mod 3 to one digit.
hidden_moves = [0, 2, 1, 1, 0, 2, 2, 1, 0, 1]
s = [0, 0, 0]
for d in hidden_moves:
    s[d] = (s[d] + 1) % 3
target = h(state_bytes(s))
tried = 0
found = None
for cand in itertools.product(range(3), repeat=3):
    tried += 1
    if h(state_bytes(cand)) == target:
        found = cand
        break
print("FINITE_RECOVER|space=27|moves_hidden=%d|hashes_tried=%d|recovered=%s|equals_hidden=%d|json=0" % (len(hidden_moves), tried, ",".join(map(str, found)), int(list(found) == s)))

# 2. Transpositions: every sequence of L turns from (0,0,0); paths 3^L, distinct end positions at most 27.
for L in (6, 9, 12):
    table = {}
    paths = 0
    for seq in itertools.product(range(3), repeat=L):
        paths += 1
        st = [0, 0, 0]
        for d in seq:
            st[d] = (st[d] + 1) % 3
        table.setdefault(h(state_bytes(st)), tuple(st))
    print("TRANSPOSITION|L=%d|paths=%d|distinct_positions=%d|reuse_ratio=%d|json=0" % (L, paths, len(table), paths // len(table)))

# 3. An unstructured input: measure this machine's SHA-256 rate, then the cost of a search over 2^64 candidates.
n = 1 << 20
t0 = time.perf_counter_ns()
for i in range(n):
    hashlib.sha256(i.to_bytes(8, "big")).digest()
ns = time.perf_counter_ns() - t0
rate = n * 10 ** 9 // ns
years = (1 << 64) // rate // (365 * 24 * 3600)
print("UNSTRUCTURED_SEARCH|space=2^64|hashes_per_second_measured=%d|full_search_years_at_that_rate=%d|space=2^256_for_a_sha256_preimage_in_general|json=0" % (rate, years))
