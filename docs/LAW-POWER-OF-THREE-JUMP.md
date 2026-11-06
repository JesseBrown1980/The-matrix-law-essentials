# The power-of-three jump law

**Operator's aim (2026-11-06):** test with this law instead of running whole
matrices for hours.

**The law in plain words.** When one step of a system is one fixed linear map
M, the state after 3^k steps is M cubed k times, applied to the start. When
the step is also the same at every site of a ring of D sites, and D divides
p − 1, the integer NTT turns M into D separate numbers (its eigenvalues). The
jump to step n then raises each of those D numbers to the power n. That work
stays about the same whether n is 27 or 3^80. Every jumped state gets a
SHA-256, so a later run can be checked against it byte for byte.

Everything below is exact integer arithmetic modulo the kernel prime
**p = 1,000,081**, with p − 1 = 1,000,080 = 2⁴ · 3³ · 5 · 463 (smallest
primitive root 7). Every claim carries one tag:
- **MEASURED:** a value printed by a named run, kept in `data/jump3/`.
- **DERIVED:** follows from the stated algebra.
- **CONJECTURE:** expected, with no derivation in hand.
- **CANNOT_SEE:** outside what this lane read or ran.

## Exact form

A ring step with stencil {(s, c_s)} is

    x'[i] = Σ_s c_s · x[(i − s) mod D]   (mod p),

the circulant matrix M with M[i][(i − s) mod D] = c_s.

| route | what it computes | modular products |
|---|---|---|
| **direct** | n steps, one at a time | n · D · S (S = stencil size) |
| **cube** (any fixed linear map) | M^(3^k) = (…((M³)³)…)³ with k cubings, then one product with the start | 2kD³ + D² |
| **NTT** (circulant, D \| p − 1) | X[m] = Σ_j x[j]·w^(jm), λ_m = Σ_s c_s·w^(sm), x_n = D⁻¹ Σ_m λ_m^n X[m] w^(−jm) | 2D² + 2D + Σ_m (products in λ_m^n) ≤ 2D² + 41D |

Here w is a primitive D-th root of unity mod p, w = 7^((p−1)/D). The power
λ^n uses n mod (p − 1) when λ ≠ 0, so each power takes at most 39 products
(p − 1 < 2²⁰), whatever n is. A zero eigenvalue stays zero for every n ≥ 1.
The NTT route also reaches any step n directly, powers of three or not.

**Leaves.** Forward leaf: M · x_n = x_{n+1}. Backward leaf: the inverse
eigenvalues give M⁻¹ x_n = x_{n−1}, which steps back onto x_n. The backward
leaf needs λ_m ≠ 0 for every m.

**Zero mode.** If λ_m = 0, the vector v_m[j] = w^(−jm) satisfies M v_m = 0.
The starts x₀ and x₀ + v_m then have the same state, and the same SHA-256,
after every n ≥ 1. The backward leaf is refused for such a step.

**State hash.** SHA-256 of the D components, each as an 8-byte big-endian
integer. This is the encoding of the seat's scratch run `jump3_test.py`, so
its 16-hex prefixes are reproduced exactly.

## Conditions (DERIVED)

1. **One fixed linear map.** The same M at every step. A nonlinear rule
   sampled as an ensemble lies outside the law. Its linear part can be jumped.
2. **Dense cubing costs about D³ per cube and stores D² entries.** At
   D = 6,561 one matrix holds 43,046,721 entries, and one product takes
   282,429,536,481 multiply-adds. The seat records 10.77 million states for
   the polariton's full space at N = 6,561. A dense matrix on that space holds
   about 1.16·10¹⁴ entries, more than this machine can store. (That figure is
   the seat's; this lane kept out of the polariton folders by rule, so it is
   CANNOT_SEE here.)
3. **A step that is the same at every site diagonalizes under the integer NTT
   of size D, when D divides p − 1.** For p = 1,000,081 this allows rings of
   27, 45 and 54 sites (all measured), and every other divisor of 1,000,080.
   A ring of 81 sites needs 3⁴ | p − 1. Its NTT is refused (measured), and the
   cube route still serves it. A ring of 81 or 6,561 sites jumps by NTT under
   a prime p′ ≡ 1 (mod 81) or (mod 6,561). Dirichlet's theorem guarantees
   infinitely many such primes.
4. **Playing backwards needs every eigenvalue nonzero.** Unitary evolution
   qualifies: the shift R, a permutation, has eigenvalues w^m, and its
   backward leaf holds (measured). The step j = I + R + R² on a 27-ring has
   1 + z + z² = 0 exactly when z is a primitive cube root of unity, so it has
   zero modes m = 9 and 18 (measured).
5. **The hash fingerprints and verifies a computed state.** Recovering a
   state from its hash means hashing candidates, so the size of the space sets
   the work: 27 positions on the 3×3×3 cube, 9 with the digit-sum invariant,
   and 2⁶⁴ candidates for an input of 64 unstructured bits.

## Measured cases (`python/matrix_tricks/jump3_ntt.py --self-test`, `crates/jump3`)

Start state x₀[i] = 7i² + 3i + 1 mod p. Python and Rust print the same
36 rows (`data/jump3/JUMP3-SELFTEST.hbp`, sha256 `8e3a7a60…153941`).

| case | ring | step | k | cube = NTT | direct = jump | leaves | zero modes |
|---|---|---|---|---|---|---|---|
| signed_return | 27 | I + R − R² | 1, 2, 3, 6, 9, 10, 20, 40, 80 | 9/9 | 6/6 (k ≤ 10) | forward 9/9, back 9/9 | none |
| curved_sum | 27 | I + R + R² | 1, 2, 3, 9, 40 | 5/5 | 4/4 (k ≤ 9) | forward 5/5, back refused | 9, 18 |
| shift | 27 | R | 1, 2, 3, 40 | 4/4 | 3/3 (k ≤ 3) | forward 4/4, back 4/4 | none |
| ring45 | 45 | 2I − R + 3R⁻¹ | 1, 2, 3, 6, 9, 20 | 6/6 | 5/5 (k ≤ 9) | forward 6/6, back 6/6 | none |
| ring81 | 81 | I + R − R⁻¹ | 1, 2, 3 | NTT refused (81 ∤ p − 1) | 3/3 | forward 3/3 (against direct n + 1) | — |

**The seat's hashes, reproduced** (signed_return, sha256 first 16 hex):
k = 3 `7d05096afc24b08c`, k = 6 `cc3581ab163de3a7`, k = 9 `ee4b40fb753713ab`,
k = 10 `81323c66a0e8a92f`, k = 20 `8ef1ccc86cd3ff0c`, k = 40
`177066f133a5c9ef`. At k = 80 the state hash is `fbd38706e0cfdc81…e557` (new).

**Work, measured counts (signed_return, D = 27, S = 3):**

| k | n = 3^k | cube products | NTT products | direct products |
|---|---|---|---|---|
| 3 | 27 | 118,827 | 1,728 | 2,187 (run) |
| 10 | 59,049 | 394,389 | 2,160 | 4,782,969 (run) |
| 20 | 3,486,784,401 | 788,049 | 2,322 | 282,429,536,481 |
| 40 | 1.2158·10¹⁹ | 1,575,369 | 2,079 | 984,770,902,183,611,232,881 |
| 80 | 1.4781·10³⁸ | 3,150,009 | 2,241 | 11,972,515,182,562,019,788,602,740,026,717,047,105,681 |

The NTT products stay between 1,593 and 2,322 for every k (bound 2,565). The
cube products grow as 2kD³, which is linear in k, i.e. logarithmic in n.

**Wall clock** (`data/jump3/TIMING-20261106.hbp`, this machine, with one
other heavy process running; outside the byte comparison):

| lane | direct, 3^10 steps | NTT jump to 3^40 | cube to 3^40 | direct to 3^40 at the measured rate |
|---|---|---|---|---|
| Rust 1.81 | 19.4 ms | 9.5 µs | 2.7 ms | 126,437 years (extrapolated) |
| Python 3.12 | 591.8 ms | 0.26 ms | 155 ms | 3,863,505 years (extrapolated) |

The extrapolations assume a constant cost per step (CONJECTURE).

**Zero mode, measured** (curved_sum, mode 9 and mode 18): x₀ and x₀ + v_m are
different starts. Their states agree after 1 step, and after 3⁴⁰ steps by
both routes. Both images hash to `4696b834afc58ab9…6bac`, the case's own
k = 40 hash.

**Any ring, any stencil:**

```bash
python3 python/matrix_tricks/jump3_ntt.py --ring 54 --stencil 0:1,1:-2,3:5 --ks 1,4,30 --direct-max 4
target/release/jump3                      --ring 54 --stencil 0:1,1:-2,3:5 --ks 1,4,30 --direct-max 4
```

This example (`data/jump3/JUMP3-CLI-RING54.hbp`) jumps 3³⁰ steps of
x[i] − 2x[i−1] + 5x[i−3] on 54 sites in 7,398 products, against
33,354,363,399,333,138 products stepping. The two lanes print the same bytes.

## Transpositions: the hash-keyed memo (`transposition_memo.py`, `transmemo`)

The cube has positions (a, b, c) ∈ {0,1,2}³, start 000, and one turn adds 1
mod 3 to one digit. Each end position is computed once and keyed by the
SHA-256 of its three digit bytes. Every later arrival is a memo hit.

| L | move orders | distinct positions | memo hits | orders per position (min, max) | recurrence updates (27 / 9 positions) |
|---|---|---|---|---|---|
| 6 | 729 | 9 | 720 | 63, 90 | 486 / 162 |
| 9 | 19,683 | 9 | 19,674 | 2,106, 2,268 | 729 / 243 |
| 12 | 531,441 | 9 | 531,432 | 58,806, 59,535 | 972 / 324 |

- **Correction beside the seat's record.** The average reuse (81, 2,187,
  59,049) hides unequal counts per position. At L = 6, 000 has 63 orders,
  111 and 222 have 90, and the six orderings of 012 have 81 each (GPT-6's
  table; measured here).
- **Counting rule.** Any result that depends on path counts uses the
  recurrence C_{n+1}(s) = Σ_{p→s} C_n(p). It equals the full enumeration at
  L = 6, 9 and 12 (measured).
- **The invariant.** The digit sum mod 3 equals L mod 3. The reachable
  positions are exactly the 9 of that class, so the recurrence needs 9
  positions per turn instead of 27 (measured).
- **Three different digits (DERIVED, measured).** For L divisible by 3, a
  position with three different digits is reached by exactly 3^(L−2) move
  orders: in the character sum, every term of size √3 cancels. Measured at
  L = 6, 9, 12, 27 and 54. At L = 54, 000 has 6,461,081,889,221,589,567,275,583
  orders and 111 has 6,461,081,889,229,215,164,760,570.
- **Jumping the counts.** The recurrence is one fixed linear map on 27 counts,
  so the counts after 3⁴⁰ turns come from 40 cubings (1,575,369 products,
  exact mod p). They check against the recurrence mod p at 3³ and 3⁴ turns,
  against Σ = 3^L mod p, and against 3^(L−2) mod p on the six 012 orderings.
  9 positions are nonzero mod p (measured).
- **Recovery from a hash.** 10 hidden turns end at 010. Hashing candidates in
  order finds it at the 4th of 27, or the 2nd of 9 with the invariant. The
  seat's "4 hashes" is this place in the enumeration order, and 27 (or 9) is
  the bound.

## How it answers the operator's aim

The seat's measurement is the model for it. Stepping 3⁴⁰ times is 10²¹
products, about 126,000 years in Rust at the measured rate. The jump takes
2,079 products and 9.5 µs, and gives the same state, checked against direct
stepping wherever direct stepping is affordable. To use it as a test:

1. **Write the step as one fixed linear map.** If it is the same at every
   site of a ring whose size divides p − 1, jump by NTT to any step n. The
   cost is flat in n. If it is linear but site-dependent and D is a few
   hundred or less, cube: 2kD³ products to 3^k.
2. **Verify by hash at small k, then jump.** Agreement of the jump with
   direct stepping at k ≤ 10 checks the implementation. M^(3^k) = (M^(3^(k−1)))³
   is exact algebra for every k, so the large-k states carry the same
   identity. Their SHA-256 pins them for every later run.
3. **Land between the powers.** Use the leaves (n ± 1), the NTT at any n, or
   the triadic composition of `time_jump_ladder.py` (n written in base 3).
4. **Count paths by recurrence, then jump the recurrence.** When many orders
   reach one state, keep one memo entry per state, keyed by its hash. Count
   orders with C_{n+1}(s) = Σ_{p→s} C_n(p), reduced by every invariant.
5. **Check the conditions before trusting a jump.** One linear map
   (condition 1); dense storage D² that fits (2); D | p − 1, or a prime chosen
   for D (3); nonzero eigenvalues before playing backwards (4); a hash
   compares computed states, and recovery needs an enumerated space (5).

**Open for the seat (CANNOT_SEE in this lane).** This lane left the polariton
bridge's frozen blind runs and their folders untouched. Whether the
polariton's step, or a part of it, is one translation-invariant linear map on
a ring is the seat's call. Rule 3 then picks the prime: 6,561 = 3⁸ needs
p′ ≡ 1 (mod 6,561).

## Predictions first

`data/jump3/PREDICTIONS-20261106.hbp` (commit `b4df487`, before any of the
code) fixes 366 PREDICT rows from closed formulas and the seat's scratch
runs. They break down as 44 known from the seat, 318 derived and 4
conjectures. `data/jump3/check_predictions.py` re-reads them against the
outputs: 366 met, 0 missed (`data/jump3/PREDICTIONS-CHECK-20261106.hbp`). The
4 conjectures were ring45 having no zero mode, and 9 positions nonzero mod p
at 3³, 3⁴ and 3⁴⁰ turns. Correction beside `b4df487`: its commit message
says 369 PREDICT rows. The file holds 366 PREDICT rows and 3 PREDICT_NOTE
rows under its header.

## Reproduce

```bash
python3 python/matrix_tricks/jump3_ntt.py --self-test          # 36 rows, JUMP3_VERIFY|status=PASS
python3 python/matrix_tricks/transposition_memo.py --self-test # 12 rows, TRANSMEMO_VERIFY|status=PASS
cargo build --release -p jump3
target/release/jump3 --self-test | cmp - data/jump3/JUMP3-SELFTEST.hbp
target/release/transmemo --self-test | cmp - data/jump3/TRANSMEMO-SELFTEST.hbp
./verify.sh                                                     # runs all of the above (section 5e-3)
```

Add `--timing` for wall-clock rows. Integer-only throughout: the timing uses
`monotonic_ns` and `Instant` nanoseconds.

Theory and architecture: Jesse Daniel Brown.
