# The matrix law essentials

Integer-only tools, tricks and results for the finite S₃ lattice gauge
Hamiltonian H = JE·E + JB·M on periodic tori. r = JB/JE. This repository has
everything needed to rebuild the operators, re-run the certificates, and
reproduce the measurements and exact laws below. Each result is labelled:

- **Certified**: a complete integer proof receipt (interval Sturm inertia), reproducible here.
- **Reported**: certificates in the published ledger that are pending independent review.
- **Measured**: Euclidean-evolution measurements. They are precise but are *not* certificates.
- **Exact (perturbative)**: closed-form weak-coupling laws, derived in exact integer arithmetic.

Source of record for the 2×2 programme, its PR stack and its reviews:
[Finite-universe Yang–Mills (HEB-C)](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C).
Nothing here is a continuum Yang–Mills result.

```
./verify.sh      # builds, rebuilds 2x2 byte-exact, builds 2x3, certificates, exact laws (2nd/3rd/4th order), integer audit
```

## The model

Kogut–Susskind Hamiltonian with finite group S₃. Links carry S₃ elements.
States are gauge orbits: a comb spanning tree is fixed to the identity, and
the k = LxLy + 1 free loop variables are taken modulo global conjugation.
Each orbit carries a weight (its orbit size: 1, 2, 3 or 6). The electric term
on each link is Σ over the 3 transpositions t of (1 − L_t); its single-link
spectrum is {0, 3, 6}. The magnetic term on each plaquette is 2 − χ_std(U_p),
with values in {0, 2, 3}. E is a Gauss-weighted self-adjoint graph Laplacian
and M is diagonal. Dimension (Burnside):

N = (6^k + 3·2^k + 2·3^k) / 6

| torus | states | nonzeros | flat vacua | operator sha256 |
|---|---|---|---|---|
| 2×2 | 1,393 | 32,385 | 8 | `137285e0…1d2c` (byte-identical to the pinned HEB-C operator) |
| 2×3 | 47,449 | 1,726,321 | 8 | `2bd842b0…66de` |
| 3×3 | 10,097,891 | — | 8 | (not built; see exact laws) |

Every torus has exactly 8 flat vacua: the commuting holonomy pairs of S₃ up to conjugation.

## The tricks

| trick | what it buys |
|---|---|
| **Scaling invariance** count(kJE, kJB, kt) = count(JE, JB, t) | One certificate covers a whole ray. Micro-unit thresholds come free by running at JE = 64,000 (`ladder2.py`). |
| **Monotone cell lemma** gap ≥ E₁(a) − E₀(b) on [a, b] | Eigenvalues of E + rM never decrease in r, since M ⪰ 0, so two certificates cover a continuous interval. |
| **Z₂×Z₂ translation sectors** (385 + 336 + 336 + 336) | The four sector counts add to the full count. About 150 s → 6.9 s per certificate; every 10th anchor is spot-checked on the full matrix. |
| **Predict-then-verify** | The second-order law picks the next jump, one interval-Sturm certificate verifies it, and a refuted jump is halved. Only verified results are recorded. |
| **SHA-keyed certificate ledger** | Key = sha256(mode, matrix shas, JE, JB, t, digits). Reruns resume from the ledger. Equivalent matrices hash differently, so the ledger reuses results but never predicts them. |
| **Fail-closed thresholds** | A threshold exactly on an eigenvalue returns "uncertified" (r = 0: use 5999, not 6000). |
| **Matrix jump** P = T^(2^J), T = c − H | Each squaring doubles the Euclidean time, so one multiply advances 2^J steps (J = 26 is about 67 million steps, ~3 s per 385-state sector). |
| **Gershgorin shift** | Puts c rigorously above the spectrum. Assuming "E₀ < 24" was a real bug: it fails for non-vacuum sectors at large r (kept as evidence in `data/sweep/sweep_bugged_shift.hbp`). |
| **Deflated Rayleigh energies** | E₁ comes from evolution with the ground vector projected out (Hotelling), so gap = E₁ − E₀ with no ratio fitting and quadratic error. |
| **Integer Chebyshev filter** | Three-term recurrence on ((a+b) − 2H)/(b − a) with doubling blocks: about √ as many steps as plain stepping. Matches the matrix jump to the micro-unit at r = 9 on 2×2. |
| **Vacuum-neighbourhood perturbation** | The weak-coupling law needs only the 8 vacua and their one-, two- and three-step neighbours, never the full matrix, so it reaches 3×3, 4×4 and beyond. |

## Results: 2×2 torus

**Certified** (PRs #9 and #14 in the source repository): E₁ > 13.240 at
(1, 1) and E₁ > 50.840 at (3, 5); the PR #14 anchors; the magnetic-end law
gap → (√14 − √6)·JE²/JB with a rigorous bound for r ≥ 10⁴ (PR #11);
coverage cells on 0 ≤ r ≤ 2 (PR #13) and a global positive bound (PR #12).

**Reported**: all 17 rungs, from `data/ledger/ledger.hbp` (846 certificates):

| r | E₀ < | E₁ > | gap/JE > |
|---|---|---|---|
| 0 | 0.001 | 5.999 | 5.998 |
| 1 | 7.750 | 13.187 | 5.437 |
| 2 | 14.375 | 18.125 | 3.750 |
| 3 | 19.313 | 20.625 | 1.312 |
| 5 | 21.885 | 22.210 | 0.325 |
| 6 | 22.283 | 22.537 | 0.254 |
| 8 | 22.746 | 22.921 | 0.175 |
| 9 | 22.895 | 23.045 | 0.150 |
| 15 | 23.347 | 23.432 | 0.085 |
| 21 | 23.536 | 23.593 | 0.057 |
| 27 | 23.639 | 23.683 | 0.044 |
| 45 | 23.787 | 23.812 | 0.025 |
| 63 | 23.846 | 23.865 | 0.019 |
| 81 | 23.882 | 23.895 | 0.013 |
| 135 | 23.930 | 23.937 | 0.007 |
| 189 | 23.951 | 23.955 | 0.004 |
| 243 | 23.961 | 23.965 | 0.004 |

The adaptive ladder certifies continuous cells from r = 3 upward; it was past
r ≈ 37 at this push. The micro-unit ladder (`ladder2.py`) continues toward 10⁴.

**Measured** (matrix jump on the four sectors, `data/sweep/sweep.hbp`): all 17
rungs land inside their brackets. Examples: r = 1 gap 5.608493, r = 3 gap
1.473689, r = 243 gap 0.005318. r·gap → 1.2923 against √14 − √6 = 1.29217.
For r > 0, E₀ and E₁ both lie in the vacuum sector (++), and the +− and −+
sectors are exactly degenerate.

## Results: 2×3 torus (first volume data)

**Measured** (Chebyshev evolution on all 47,449 states, `data/volume/`):

| r | 2×2 gap | 2×3 gap |
|---|---|---|
| 0 | 6.000000 | 6.000000 |
| 1 | 5.608493 | 5.617344 |
| 2 | 4.002279 | 4.165752 |
| 3 | 1.473689 | 1.231379 |
| 5 | 0.341237 | 0.640450 |
| 6 | 0.257248 | 0.525434 |
| 8 | 0.177326 | 0.386040 |
| 9 | 0.154383 | 0.341110 |
| 15 | 0.088346 | 0.201684 |
| 21 | 0.062323 | 0.143471 |
| 27 | 0.048228 | 0.111398 |

E₀ per site agrees between the two volumes (r = 1: 1.9085 against 1.9089). The
gaps cross between r = 3 and r = 5. At weak coupling the measured 2×3 gap
approaches the exact law 3/r with a residual that falls off as 1/r³:

| r | 2×3 gap | (gap − 3/r)·r³ |
|---|---|---|
| 9 | 0.341110 | 5.669 |
| 15 | 0.201684 | 5.684 |
| 21 | 0.143471 | 5.685 |
| 27 | 0.111398 | 5.647 |

This supports g(r) = 3/r + O(r⁻³) over the measured range. It is not a
rigorous bound on the remainder. The rung receipts (E₀, E₁, gap, steps,
settled status, tolerance, runtime, `json=0`) are the original measurement
lines in `data/volume/volume_t23.hbp` (sha256 `53b29aeb…`, checked by `verify.sh`).

## Results: exact weak-coupling laws

These come from integer characteristic polynomials and Sturm root isolation
(`python/volume/`). The full-matrix and neighbourhood-only computations agree
coefficient for coefficient on 2×2 and 2×3.

| torus | vacua split at | gap at large r = JB/JE |
|---|---|---|
| 2×2 | 2nd order | (√14 − √6)·JE²/JB |
| 2×3, 2×4, 2×5, 2×6 | 2nd order | exactly 3, 4, 5, 6 × JE²/JB |
| 3×3 | 3rd order | (9/8)(√14 − √6)·JE³/JB² |
| 3×4, 3×5 | 3rd order | exactly 3 and 15/4 × JE³/JB² |
| 4×4 | 4th order | (5/4)(√14 − √6)·JE⁴/JB³ |

On 2×3 the characteristic polynomial factors exactly:
(x−54)(x−42)(x−36)⁴(x−30)(x−18) for 4K. On 3×3, 16·H₃ has characteristic
polynomial x⁴(x⁴ − 6480x² + 8,817,984), with roots ±18√14, ±18√6, 0 (×4).

What this says: at weak coupling the finite-torus gap is a tunnelling
splitting among the 8 vacua. It enters at the order set by the shortest
non-contractible loop, JE^L/JB^(L−1) on an L×L torus, so at fixed r it
shrinks steeply as the lattice grows. Square tori share one vacuum-mixing
pattern (√14 and √6) with coefficients 1, 9/8, 5/4 for L = 2, 3, 4.

On 3×3 (and on 3×4, 3×5, 4×4, 4×5, 5×5) the 8×8 second-order matrix is
verified to be exactly scalar. Its gap coefficient is therefore 0, but that is
a second-order statement, not a vanishing gap: the splitting starts at third
order (fourth order on 4×4). `flat_law.py` reports
`second_order_degenerate=1|second_order_gap_coeff=0|next_order_required=3`,
and `flat_law.hbp` ends with `status=COMPLETE`. The first run, which stopped
with an `IndexError` at this degenerate case, is kept with an explanatory
header as `data/volume/flat_law_first_run_INDEXERROR_degenerate_3x3.hbp`; its
2×2–2×6 lines are valid.

**Pending checks**: the 3rd- and 4th-order formulas (3×3, 4×4) are exact
perturbative calculations pending independent mathematical review, and have
not yet been compared with direct numerics on those lattices. The 2×3 numerics confirm
the 2nd-order law. Every exact law needs independent mathematical review.

## Not established

- A gap bound uniform in the lattice volume. The data so far point the other way at weak coupling: the tunnelling gap shrinks with size.
- Any continuum Yang–Mills statement, the 4D mass gap, or the Clay problem.
- Authenticated provenance: receipts prove integrity relative to trusted inputs, not who produced them.

## Layout

```
crates/torusl     Lx×Ly gauge-orbit operator builder (Rust 1.81, no external crates)
                  + onfly: matrix-free H application for lattices too large to store (3×3)
crates/sturmrs    interval Sturm / LDLᵀ inertia verifier (num-bigint)
crates/evolvers   evolvers (sparse stepping), jump (matrix jump), cheb (Chebyshev)
python/ladder     ladder.py (milli), ladder2.py (micro, scaling trick), rungs.py
python/testground HTTP-0 page and routes (127.0.0.1): /evolve, /svg, /jump, /cheb
python/sweep      17-rung sector sweep
python/volume     volume.py, second_order.py, flat_law.py, flat_law3.py, flat_law4.py
python/audit      integer-only AST audit
data/             operators, ledgers, sweeps, volume runs, exact-law receipts (json=0)
```

Rules throughout: integers only (no float, Decimal, Fraction or float timing),
outward-rounded integer division in certificates, `json=0` on every receipt
line, Rust 1.81 with `clippy -D warnings`.
