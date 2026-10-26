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

## Exact matrix jumps: guaranteed integer equality on a fixed matrix

The reusable [exact matrix-power jump](python/matrix_tricks/exact_jump.py) is
now part of the checked toolset. For any explicitly supplied small integer
square matrix `S`, the calculation `S^N v` is **bit-for-bit identical**
to `N` consecutive applications, using exponentiation by squaring with
no floating-point arithmetic, normalization or lossy amplitude shifts.

The 3×3 S₃ eight-vacuum **third-order effective** matrix passes:
64 ordinary steps = 6 squarings; 4,096 steps = 12 squarings and also
two exact 2,048-step jumps composed. Pinned SHA-256 state hashes are checked
by `verify.sh`. The guarantee is equality, **not** bounded execution time,
computational feasibility for a dense 10-million-state matrix, or physical
time evolution. The pilot's 206-state and 4,454-state vacuum-ball Ritz
measurements are separate approximations, **not full-system gap certificates**.

Methods, step-count complexity, exact hashes, scope and vacuum-ball
measurements: [docs/EXACT_MATRIX_JUMP_AND_VACUUM_BALL.md](docs/EXACT_MATRIX_JUMP_AND_VACUUM_BALL.md).
Run: `python3 python/matrix_tricks/exact_jump.py --self-test`. Reproduce the actual 3×3 radius-1/radius-2 principal-ball Ritz brackets with `python3 python/matrix_tricks/vacuum_ball_ritz.py 1 27 --order 1` and `python3 python/matrix_tricks/vacuum_ball_ritz.py 2 27 --order 2`.

### Cosmological time-jumps, triadic ladders, and hash-altering in-between sampling

The [time-jump ladder](python/matrix_tricks/time_jump_ladder.py) tool provides:
- **Cosmological Time-Jump**: 60 squarings advance 2^60 steps (>36.5 billion years at 1 s/step) in ~5 ms, converging to the exact analytical ground state of the effective 8x8 vacuum matrix.
- **Triadic Ladder Skips**: Base-3 progression (1, 3, 9, 27, 81, 243 ... x 3^n) and intermediate harmonic rungs (2, 5, 6, 8) evaluate ANY intermediate checkpoint t in at most 2 * ceil(log3(t)) matrix-vector multiplies, without simulating intervening steps.
- **Skip with Hash Change**: Advances directly to target epoch T_skip, content-addresses the jumped state with SHA-256, alters the operator hash (shift, coupling, or perturbation), and measures local in-between dynamics along the ladder (dt in {1, 2, 3, 5, 6, 8, 9, 27}) directly from the jumped state.

Run:
```bash
python3 python/matrix_tricks/time_jump_ladder.py --self-test
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
| **Vacuum-ball compression** (`ball`) | H restricted to every gauge orbit within K moves of the 8 vacua, written in the standard matrix format so the evolvers run unchanged. By interlacing its energies are upper bounds that fall as K grows. 3×3: 206 / 4,454 / 59,246 / 518,177 states for K = 1–4, instead of 10,097,891. |
| **Block cluster eigensolver** (`cluster`) | Resolves all 8 near-degenerate vacuum levels at once: Chebyshev-filtered subspace iteration on a 12-vector block, fixed-point Gram–Schmidt with re-randomisation of collapsed vectors, Rayleigh–Ritz by integer Jacobi rotations. Each level gets a rigorous residual interval computed in exact big integers, and interval-Sturm counts between the intervals certify which levels are E₀ and E₁ (`cluster_certify.py`). |
| **Long jump + play-forward** (`jumpto`) | One cache of squarings T^(2^j) reaches any step N ≤ 2^(J+1) − 1, and the state can be played step by step from there. **Exact mode** (modulo a prime, exact algebraic replay only, no real amplitudes): the full 1,393-state 2×2 operator jumped to 10¹⁰ steps (33 squarings, ~4 min), with jump(1024) = 1,024 sequential steps, jump(5·10⁹)² = jump(10¹⁰), 27 replayed frames = jump(10¹⁰ + 27), and an independent Python sparse-step reference matching the state fingerprints at 1,024 steps. **Fixed-point mode** measures the energy at 10⁰ … 10¹⁰ steps (3×3 ball, r = 81: relaxes over 10⁴–10⁷ steps); at 10¹⁰ an exact residual radius places an eigenvalue in [53.832773, 53.832776], the certified E₀ interval of that ball. One step is v ← (c − H)v, a forward-Euler imaginary-time step 1 − H/c: not e^(−τH) and not years. A single relaxed vector suppresses the excited levels, so splittings are measured with `cluster`, and the ball remains a truncation of the full 3×3 operator. |
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

## Results: vacuum-ball compression (Rust `ball`)

**Validation on 2×3**, where the full answer is known (Chebyshev; full lattice: r = 9 gap 0.341110, r = 27 gap 0.111398):

| K | ball states | r = 9: E₀ | r = 9: gap | r = 27: gap |
|---|---|---|---|---|
| 1 | 122 | 34.557113 | 0.309247 | 0.110133 |
| 2 | 1,520 | 34.474701 | 0.339132 | 0.111391 |
| 3 | 10,332 | 34.470864 | 0.341032 | 0.111399 |
| 4 | 34,890 | 34.470733 | 0.341110 | 0.111399 |
| full | 47,449 | 34.470732 | 0.341110 | 0.111398 |

E₀ falls monotonically toward the full value (interlacing). At r = 9, K = 4 is
exact to 1 micro-unit, and at r = 27, K = 3 already is.

**3×3, K = 1 (206 states)**: converged ball eigenvalues by the dense matrix jump
(dE = 0), next to the 16-dimensional Krylov Ritz values of the same ball
(`python/matrix_tricks/vacuum_ball_ritz.py`) and the third-order law:

| r | ball K = 1 gap | Krylov-16 Ritz gap | (9/8)(√14 − √6)/r² |
|---|---|---|---|
| 9 | 0.016339 | – | 0.017947 |
| 27 | 0.001983 | 0.002005 | 0.001994 |
| 81 | 0.000221 | – | 0.000222 |

The ball's energies sit at or below the Krylov values (r = 27: E₀ 53.496482 ≤
53.496532, E₁ 53.498465 ≤ 53.498537), as interlacing requires for a subspace of
the same ball. The ball contains the length-3 tunnelling paths, and its gap
approaches the third-order law as r grows: 9%, 0.6% and 0.3% away at r = 9, 27
and 81. That is a check of the 3×3 law against the actual Hamiltonian
compressed to the ball, not against the 8×8 effective matrix.

**Block eigensolver results (3×3)** (`data/ball/cluster_runs.hbp`). Each
level's interval contains an eigenvalue of the compressed operator (exact
residual bound). Where marked certified, interval-Sturm counts (exactly 1
eigenvalue below t₁, exactly 2 below t₂) prove the intervals hold E₀ and E₁,
so the gap enclosure is rigorous **for that compressed operator**:

| r | K = 1 (206), certified | K = 2 (4,454), certified | K = 3 (59,246), ordering unverified | third-order law |
|---|---|---|---|---|
| 9 | gap ∈ [0.016336, 0.016342] | gap ∈ [0.019504, 0.019510] | gap ∈ [0.019802, 0.019808] | 0.017947 |
| 27 | gap ∈ [0.001980, 0.001986] | gap ∈ [0.002009, 0.002015] | gap ∈ [0.002010, 0.002016] | 0.001994 |
| 81 | gap ∈ [0.000218, 0.000224] | gap ∈ [0.000218, 0.000226] | gap ∈ [0.000216, 0.000228] | 0.000222 |

Read across: K = 1 undershoots, and from K = 2 on the gap settles slightly
above the third-order law (about +10% at r = 9, +1% at r = 27, within the
interval at r = 81), consistent with a positive next-order correction that
fades as r grows. At r = 27 the full 8-level cluster reproduces the
third-order spectrum level by level: offsets −0.00575, −0.00377, ≈ 0 (×4),
+0.00369, +0.00565 against (9/8)·{−√14, −√6, 0, √6, √14}/r².

Scope: these are certified for the compressed operators. The ball energies are
upper bounds on the full 3×3 E₀ and E₁ (interlacing), but the full-lattice gap
also needs a bound on the states outside the ball. That remains open.

**Superseded single-vector runs, kept as evidence**: Chebyshev runs on the 3×3 K = 2 and
K = 3 balls. Interlacing is violated (r = 9: E₀ at K = 3 is above K = 2), the
gap is negative at r = 81 with K = 1, and the final-block drift is large. The 8
vacuum levels lie within about 10⁻³, with a four-fold cluster next to E₁, so a
single-vector filter cannot resolve them; the block solver above does.
Receipts for every build and run are in `data/ball/ball_runs.hbp`.

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
                  + ball: vacuum-ball principal compression (upper bounds by interlacing)
crates/sturmrs    interval Sturm / LDLᵀ inertia verifier (num-bigint)
crates/evolvers   evolvers (sparse stepping), jump (matrix jump), cheb (Chebyshev), cluster (block eigensolver), jumpto (long jump + play)
python/ladder     ladder.py (milli), ladder2.py (micro, scaling trick), rungs.py
python/testground HTTP-0 page and routes (127.0.0.1): /evolve, /svg, /jump, /cheb
python/sweep      17-rung sector sweep
python/volume     volume.py, second_order.py, flat_law.py, flat_law3.py, flat_law4.py
python/matrix_tricks exact_jump.py (exact jumping) + vacuum_ball_ritz.py (3x3 reduced Ritz)
python/audit      integer-only AST audit
data/             operators, ledgers, sweeps, volume runs, exact-law receipts (json=0)
```

Rules throughout: integers only (no float, Decimal, Fraction or float timing),
outward-rounded integer division in certificates, `json=0` on every receipt
line, Rust 1.81 with `clippy -D warnings`.
