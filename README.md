# The matrix law essentials

Integer-only tools, matrix tricks and certified results for the finite S₃
lattice gauge Hamiltonian H = J_E·E + J_B·M on periodic tori, with r = J_B/J_E.
Everything here rebuilds the operators, re-runs the certificates and reproduces
the measurements and exact laws below from this repository:

```
./verify.sh      # builds, rebuilds every operator byte-exact, replays every certificate, exact laws, integer audit
```

Results carry one of four labels:
- **Certified:** integer interval-inertia proof receipts, replayed by `verify.sh`.
- **Reported:** ledger certificates under independent review.
- **Measured:** precise Euclidean-evolution or eigensolver values.
- **Exact:** closed-form weak-coupling laws in exact integer arithmetic.

## November 12: state-dependent nonlinear jumps, exact and certified-bounded

Two distinct tested results extend the [fixed-source affine jump](docs/LAW-AFFINE-FORCED-RADIX3-JUMP-20261112.md) to **genuinely state-dependent nearest-neighbor interactions**:

1. **Exact reversible nonlinear jump:** the two-field update `F=T⁻¹LT`, where `T(u,v)=(u,v+7u_i u_(i-1))` and `L` is a coupled invertible sparse step, is nonlinear in the original variables. A **proven invertible change of coordinates** gives `F^n=T⁻¹L^nT`, allowing the power-of-three Fourier jump, forward/backward frozen leaves and direct-state matching at **27, 243, 6,561 and 59,049 spatial sites** (largest **118,098 interacting field components**), including `3^80` on 59,049 sites.
2. **Certified long-time bound for a different irreversible nonlinear map:** a bilinear neighbor-dependent signed-integer recurrence satisfies a verified uniform sup-norm comparison. For **every** initial state with `max |x_i|<=1,000,000`, its entire state is **exactly zero after at most 21 steps** and remains zero, including at `3^80`. This applies to any ring length and does not enumerate the state space.

The two results are not a general nonlinear fast-forward algorithm: the first requires a **specific proven conjugacy**; the second relies on **strict global contraction**. Deliberately invalid linear and nonlinear controls must fail. A separate independent checker uses explicit nonlinear equations and naive DFT instead of the solver's FFT, plus audits every one of the million possible scalar comparison radii.

[Mathematical proof and reproducibility](docs/LAW-NONLINEAR-EXACT-AND-BOUNDED-JUMPS-20261112.md) · [Exact code](python/matrix_tricks/nonlinear_exact_and_bound.py) · [Independent verifier](python/matrix_tricks/nonlinear_independent_audit.py) · [Pre-hosted seal](data/jump3/NONLINEAR-JUMP-FREEZE-20261112.hbp).

All calculations are integer-only mathematical models, **not** physical Yang–Mills, gluon, proton or gravity dynamics. The protected coupling `r=9` remains on HOLD.

## November 12: exact affine-source jump for much larger Rainbon matrices

The original [pure linear jump](docs/LAW-POWER-OF-THREE-JUMP.md) does **not** account for a source injected at every microscopic step. We have added its **exact affine extension**, where `x_next=M*x+f` and the source `f` is constant *within each composed interval*. In Fourier space `X_n=lambda^n*X_0+S_n(lambda)*F`; piecewise constant changes at frozen checkpoints also compose exactly. A **state-dependent nonlinear** force is a distinct problem and is deliberately rejected by the negative controls.

The first independent **local** experiment built and audited complete **6,561-, 19,683- and 59,049-row** sparse matrices `I+R−R²` over `GF(1,299,079)`. The largest matrix has **177,147 nonzero entries** (not a dense 3.49-billion-entry allocation); its actual CSR binary is about **1.65 MB**. Radix-three NTT reaches `3^40`, `3^80`, binary `2^40`, forward freeze/play and a **changed-source composite**; short direct forced stepping and an independent affine semigroup agree exactly. One-, two- and three-source readouts sample 243 addresses of a 12-billion-coordinate **virtual** screen. That does not simulate 12 billion independently interacting particles.

- [Mathematical theorem, evidence and limitations](docs/LAW-AFFINE-FORCED-RADIX3-JUMP-20261112.md)
- [Runnable exact Python implementation](python/matrix_tricks/affine_mega_rainbon.py) and [separate complete-CSR/direct-step auditor](python/matrix_tricks/affine_mega_audit.py)
- [Post-local-exploration frozen hosted expectations](data/jump3/AFFINE-MEGA-FREEZE-20261112.hbp) and [automated large-matrix replays](.github/workflows/affine-mega-rainbon.yml)

**Scope:** exact modular algebra and synthetic Rainbon-source readouts. Not a nonlinear Yang–Mills, gravity, radiation or physically validated gluon/glueball/proton simulation. The protected finite-model coupling `r=9` remains HOLD.

## Latest cross-scale results: 43-million jumps, molecule slices and Shadow Cats (November 6–7, 2026)

The matrix-law toolkit now has **three distinct kinds of evidence**: the certified *finite* `S₃` Hamiltonian calculations below, exact or high-precision **algorithm tests**, and deliberately simplified **physical models**. They must not be treated as interchangeable claims.

| Track | Status and specific result | Reproduce / scope |
|---|---|---|
| **Power-of-three quantum-wave jump** | [**PR #7, MERGED**](https://github.com/JesseBrown1980/The-matrix-law-essentials/pull/7). Exactly **`3^16 = 43,046,721`** model time subdivisions reached by **16 cubings** of Fourier free-flight factors. **80/80** independent modal comparisons passed, with maximum reported high-precision amplitude difference `6.10e-68`. | [Hosted PASS](https://github.com/JesseBrown1980/The-matrix-law-essentials/actions/runs/42329111050) · [model note](docs/VIENNA-NANOPARTICLE-INTERFERENCE-20261106.md) · [59-row measurement receipt](data/vienna_nanoparticle/MEASURED-20261106.hbp). |
| **Spatial slices / synthetic observations** | The finite, truncated three-grating wave model gives the same first Fourier harmonic at virtual **43,046,721** sample sites as the analytical spectral readout **without allocating a 43-million-element state**; a 9-point slice fails by `0.00692102` at one chosen grating phase. Synthetic detector profiles are committed. | [Model camera frames](data/vienna_nanoparticle/VIRTUAL-43M-DETECTOR-SNAPSHOTS.svg). Inspired by [Pedalino et al., *Nature* (2026)](https://doi.org/10.1038/s41586-025-09917-9): real **sodium nanoparticles**, not viruses, and **not quantum cloning**. The raw measured interference curve was **not** reproduced. |
| **Molecular interactions / ensemble stress tests** | [PR #6, DRAFT](https://github.com/JesseBrown1980/The-matrix-law-essentials/pull/6). CO Morse bond stretch at three temporal/spatial slices and exact harmonic time jumps; tests of **1–729** phase-balanced and coupled-model coordinates, plus a directly tested two-molecule case. Two opposite-phase errors cancel in the signed mean even when the RMS per molecule remains nonzero. | [Original hosted PASS](https://github.com/JesseBrown1980/The-matrix-law-essentials/actions/runs/42327155410) · [ensemble hosted PASS](https://github.com/JesseBrown1980/The-matrix-law-essentials/actions/runs/42327754135) · [two-molecule hosted PASS](https://github.com/JesseBrown1980/The-matrix-law-essentials/actions/runs/42328026246). A **classical fitted CO force model** and an uncalibrated interaction toy, not ab initio molecular proof. |
| **Finite three-scale `S₃` gauge transfer** | [YM4 PR #33, DRAFT](https://github.com/JesseBrown1980/yang-mills-4-dimensional-continum-proof/pull/33). In a solvable **1+1D** `S₃` gauge model, exact 9→3→1 spatial blocking reproduces selected gauge-invariant correlation functions when the induced interaction is retained. | [Hosted PASS](https://github.com/JesseBrown1980/yang-mills-4-dimensional-continum-proof/actions/runs/42326608078). Not a 4D SU(2) continuum result. |
| **Retained `C0/C1/C2` Shadow Cat residue** | [YM4 PR #34, DRAFT](https://github.com/JesseBrown1980/yang-mills-4-dimensional-continum-proof/pull/34). `(C0+C1+C2) mod 3` is preserved under 27→9→3→1 blocking and ternary permutation jumps to `3^16`. An exhaustive **1,296-case** `S₃` test shows a *candidate* signed holonomy flips under gauge conjugation, whereas its Wilson character is invariant. | [Hosted PR #34 proof-boundary gate PASS](https://github.com/JesseBrown1980/yang-mills-4-dimensional-continum-proof/actions/runs/42330584710) · [finite bridge PR #125, DRAFT](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C/pull/125). A retained ternary residual is **not** itself a mass-gap eigenvalue. |

**Interpretation:** Repeated cubing certifies the *chosen* transfer rule, and finite Fourier bandlimits certify the stated virtual observation. Neither establishes that spacetime is discrete, that a nonlinear molecular system can always be fast-forwarded, or that the Clay four-dimensional Yang–Mills existence and physical mass gap are proved. The continuum work still requires a gauge-invariant excitation observable and a lower physical gap bound uniform in lattice spacing; see [the 4D companion proof ledger](https://github.com/JesseBrown1980/yang-mills-4-dimensional-continum-proof/blob/main/PROBLEM.md).

These exploratory experiments run in their **own GitHub Actions workflows** and are **not claimed to be covered by `./verify.sh`**. The original strict integer-only certificate pipeline remains the separately verified core.

The execution backend and the larger-volume runs live in
[Finite-universe Yang–Mills (HEB-C)](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C),
which uses these tools unchanged.

## Results: 3×3 torus, full zero-momentum operator (certified)

The four C₂ᵥ sectors of the complete 10,097,891-state operator hold 288,645,
279,386, 277,090 and 277,090 states. Their eigenvalue counts certify, at each
coupling:

| r | Gap above the 8-state cluster, E₉ − E₈ | Cluster width, E₈ − E₁ | First splitting, E₂ − E₁ |
|---|---|---|---|
| 9 | [23.938578, 23.939630] | [0.112613, 0.112636] | [0.019821, 0.019847] |
| 15 | [43.371514, 43.372787] | [0.038477, 0.038513] | [0.006655, 0.006695] |
| 27 | [80.161963, 80.163325] | [0.011609, 0.011692] | [0.001967, 0.002059] |
| 45 | [134.512046, 134.514099] | [0.004088, 0.004256] | [0.000635, 0.000805] |
| 81 | [242.733980, 242.735244] | [0.001140, 0.001432] | [0.000074, 0.000370] |
| 243 | [728.911221, 728.913186] | ≤ 0.000671 | below interval resolution |

All in units of J_E. The first splitting follows the exact third-order law
(9/8)(√14 − √6)/r²: 0.001994, 0.000718 and 0.000222 at r = 27, 45 and 81 lie
inside the certified intervals.

**Every coupling r ≥ 9: E₉ − E₈ ≥ 23.357514 J_E.**
- **Cells between the points.** Every eigenvalue of E + rM rises with r
  (M ⪰ 0), so on [r_a, r_b] the gap is at least t_B(r_a) − t₁(r_b).
- **Beyond the last point.** The eight flat states carry M = 0, electric
  diagonal 54 and no mutual coupling (`flatcomp`), so E₈ ≤ 54 for every r.
  That gives E₉ − E₈ ≥ 728.856 for r ≥ 243.

**How it is computed.**
- **Magnetic-shell Schur reduction (`shellschur`).** Rows that are strictly
  diagonally dominant form a positive-definite block. Haynsworth additivity
  leaves only a small core to factor exactly: 7–1,802 states instead of about
  280,000. Interval Gauss–Seidel encloses the Schur complement, and interval
  LDLᵀ counts its negative eigenvalues.
- **Where each piece ran.** r = 9 was split across 60 GitRAM cells, the other
  couplings ran through a 30-cell ladder, and r = 243 ran locally.
- **Receipts:** `data/volume/ym3x3/`.

## Results: 2×2 torus

**Certified:**
- E₁ > 13.240 at (1, 1) and E₁ > 50.840 at (3, 5).
- The translation-sector anchors.
- The magnetic-end law gap → (√14 − √6)·J_E²/J_B, with a rigorous bound for r ≥ 10⁴.
- Coverage cells on 0 ≤ r ≤ 2 and a global positive bound.

All come from the HEB-C PR stack (#9, #11, #12, #13, #14).

**Reported:** all 17 rungs, from `data/ledger/ledger.hbp` (846 certificates):

| r | E₀ < | E₁ > | gap/J_E > |
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

The micro-unit ladder (`ladder2.py`) carries the certified cells on to
r = 10,022 (ledger on HEB-C branch `integer-sturm-ladder`).

**Measured** (matrix jump on the four sectors, `data/sweep/sweep.hbp`): all 17
rungs land inside their brackets. Examples: r = 1 gap 5.608493; r = 3 gap
1.473689; r = 243 gap 0.005318, where r·gap = 1.2923 against √14 − √6 =
1.29217. For r > 0, E₀ and E₁ both lie in the vacuum sector (++), and the +−
and −+ sectors are exactly degenerate.

## Results: 2×3 torus

**Measured** (Chebyshev evolution on all 47,449 states, `data/volume/`):

| r | 2×2 gap | 2×3 gap |
|---|---|---|
| 0 | 6.000000 | 6.000000 |
| 1 | 5.608493 | 5.617344 |
| 2 | 4.002279 | 4.165752 |
| 3 | 1.473689 | 1.231379 |
| 5 | 0.341237 | 0.640450 |
| 9 | 0.154383 | 0.341110 |
| 27 | 0.048228 | 0.111398 |

E₀ per site agrees across the two volumes (r = 1: 1.9085 and 1.9089). At weak
coupling the 2×3 gap approaches the exact law 3/r with a residual falling as
1/r³: (gap − 3/r)·r³ = 5.669, 5.684, 5.685 and 5.647 at r = 9, 15, 21 and 27.

## Results: exact weak-coupling laws

From integer characteristic polynomials and Sturm root isolation
(`python/volume/`). The full-matrix and neighbourhood-only computations agree
coefficient for coefficient on 2×2 and 2×3.

| torus | vacua split at | gap at large r = J_B/J_E |
|---|---|---|
| 2×2 | 2nd order | (√14 − √6)·J_E²/J_B |
| 2×3, 2×4, 2×5, 2×6 | 2nd order | exactly 3, 4, 5, 6 × J_E²/J_B |
| 3×3 | 3rd order | (9/8)(√14 − √6)·J_E³/J_B² |
| 3×4, 3×5 | 3rd order | exactly 3 and 15/4 × J_E³/J_B² |
| 4×4 | 4th order | (5/4)(√14 − √6)·J_E⁴/J_B³ |

On 2×3, 4K has characteristic polynomial (x−54)(x−42)(x−36)⁴(x−30)(x−18). On
3×3, 16·H₃ has x⁴(x⁴ − 6480x² + 8,817,984), with roots ±18√14, ±18√6 and 0
(×4). The splitting among the 8 vacua enters at the order set by the shortest
non-contractible loop. Square tori share one mixing pattern (√14, √6), with
coefficients 1, 9/8 and 5/4 for L = 2, 3 and 4. On 3×3 and larger square tori
the second-order 8×8 matrix is exactly scalar and the splitting starts at third
order (`flat_law.py`: `second_order_degenerate=1|next_order_required=3`). The
3×3 certificates above confirm the third-order law directly.

## Results: vacuum-ball compression (measured and certified on the ball)

`ball` restricts H to every gauge orbit within K moves of the 8 vacua. By
interlacing, its energies are upper bounds that fall as K grows. On 2×3 the
ball reaches the full answer:

| K | states | r = 9: E₀ | r = 9: gap | r = 27: gap |
|---|---|---|---|---|
| 1 | 122 | 34.557113 | 0.309247 | 0.110133 |
| 2 | 1,520 | 34.474701 | 0.339132 | 0.111391 |
| 3 | 10,332 | 34.470864 | 0.341032 | 0.111399 |
| 4 | 34,890 | 34.470733 | 0.341110 | 0.111399 |
| full | 47,449 | 34.470732 | 0.341110 | 0.111398 |

On the 3×3 balls (206 / 4,454 / 59,246 states for K = 1, 2, 3), the block
eigensolver with interval-Sturm ordering gives certified gap enclosures for
each compressed operator. They converge onto the full-operator values above:
for example r = 27, K = 2 gives [0.002009, 0.002015]. Receipts:
`data/ball/cluster_runs.hbp`.

## The matrix tricks

| trick | what it buys |
|---|---|
| **Scaling invariance** count(kJ_E, kJ_B, kt) = count(J_E, J_B, t) | One certificate covers a whole ray; micro-unit thresholds come from running at J_E = 1000 or 64,000. |
| **Monotone cells** | M ⪰ 0, so eigenvalues of E + rM rise with r: two endpoint certificates cover a continuous interval. |
| **Flat-state ray** | E₈ ≤ 3 × (number of links) for every r (`flatcomp`), so one point certificate extends to r → ∞. |
| **Symmetry sectors** (`ymsec`) | Zero-momentum sectors of the complete operator under translations ⋊ C₂ᵥ, with ±1 reflection characters. Checks: group relations, commutation with H on every gauge orbit (2×2, 2×3), exact W-self-adjointness, the 7/0/0/1 flat split, and spectral moments Tr(Hⁿ) equal to the unreduced operator. [Argument](docs/YM_SYMMETRY_ARGUMENT.md). |
| **Magnetic-shell Schur** (`shellschur`) | Factor only the non-dominant core. Column-split and merge for GitRAM, with every core entry required exactly once. |
| **Fill measurement** (`fillest`, `ddcore`) | Symbolic fill-in and core sizes before any exact arithmetic. |
| **Block cluster eigensolver** (`cluster`) | All 8 near-degenerate vacuum levels at once: Chebyshev-filtered subspace iteration, integer Jacobi Rayleigh–Ritz, exact residual intervals. |
| **Matrix jump** P = T^(2^J), T = c − H | Each squaring doubles the Euclidean time; J = 26 is about 67 million steps in seconds. Gershgorin shift places c above the spectrum. |
| **Long jump + play-forward** (`jumpto`) | One cache of squarings reaches any step N ≤ 2^(J+1) − 1. Exact mode modulo a prime: the 1,393-state operator jumped to 10¹⁰ steps, matching sequential stepping and an independent reference. |
| **Exact matrix-power jump** ([`exact_jump.py`](python/matrix_tricks/exact_jump.py)) | S^N v bit-for-bit equal to N sequential steps by squaring. [Methods](docs/EXACT_MATRIX_JUMP_AND_VACUUM_BALL.md). |
| **Time-jump ladder** ([`time_jump_ladder.py`](python/matrix_tricks/time_jump_ladder.py)) | 2^60 steps in 60 squarings. Any checkpoint t in at most 2·⌈log₃ t⌉ multiplies along the triadic ladder, with SHA-256 content addressing of jumped states. |
| **Power-of-three jump** ([`jump3_ntt.py`](python/matrix_tricks/jump3_ntt.py), `jump3`) | A step that is the same at every site of a ring reaches step 3^k by k cubings, or by the integer NTT when the ring size divides p − 1 = 1,000,080: 3^40 steps in 2,079 modular products and 3^80 in 2,241, where stepping to 3^40 one step at a time takes 9.8·10²⁰ products; mod p the states repeat with a period dividing p − 1, and the 3^40 state equals the state at step 38,961. Leaves to 3^k ± 1, zero-mode refusal of the backward leaf, SHA-256 of every jumped state; Python and Rust byte for byte. [Law](docs/LAW-POWER-OF-THREE-JUMP.md). |
| **Transposition memo** ([`transposition_memo.py`](python/matrix_tricks/transposition_memo.py), `transmemo`) | Each end position stored once, keyed by its SHA-256 (3^L − 9 memo hits); path counts by `C_{n+1}(s) = Σ_{p→s} C_n(p)`, 27·3 updates per turn instead of 3^L orders; the digit-sum invariant leaves 9 of 27 candidates; counts after 3^40 turns by 40 cubings, exact mod p. |
| **Integer Chebyshev filter** (`cheb`) | Three-term recurrence with doubling blocks, about √ as many steps as plain stepping. |
| **Deflated Rayleigh energies** | E₁ from evolution with the ground vector projected out, so gap = E₁ − E₀ with quadratic error. |
| **Vacuum-ball compression** (`ball`) | Principal restriction near the vacua, upper bounds by interlacing (3×3: 206 → 518,177 states for K = 1 → 4). |
| **Matrix-free rows** (`onfly`) | H applied without storing it, for lattices too large to hold. |
| **Predict-then-verify ladder + SHA-keyed ledger** | The second-order law picks the next jump, one certificate verifies it; reruns resume from the ledger. |
| **Symmetry-sector compression** (`symsec`) | The same reduction for the three-preon cloud (`crates/rainbon`): 151–200× fewer states, exact against the full lattice. |

## Laws of November 6, 2026

The operator stated three laws on November 6, 2026. Each check below is an exact-integer Python script in [`python/matrix_tricks/laws_20261106/`](python/matrix_tricks/laws_20261106/), with its recorded output beside it. Status: measured by the seat's scripts. The power-of-three law's full code, Rust crate and independent attack are in [`docs/LAW-POWER-OF-THREE-JUMP.md`](docs/LAW-POWER-OF-THREE-JUMP.md).

### 1. The power-of-three jump

A fixed linear step M reaches step 3^k in k cubings, skipping the 3^k steps. Some steps are the same at every site of a ring whose size divides p − 1. The integer NTT splits such a step into independent modes, and the jump then costs one power per mode.

- **Measured** (`jump3_jump.py`): a 27-site ring with M = I + R − R² mod 1,000,081, where 27 divides 1,000,080. Cubing, the NTT route and direct stepping give the same state at step 3^k, with the same SHA-256 at k = 3, 6, 9, 10, 20 and 40.
- **Timings:**
  - 3^40 ≈ 1.2 × 10^19 steps take 40 cubings, 223 ms.
  - The NTT route takes under 1 ms.
  - Stepping through 3^10 = 59,049 steps takes 235 ms.
- **Leaves:** one step forward from the jump equals the jump to 3^k + 1. The state at 3^k − 1, from the inverse eigenvalues, steps onto the jump.
- **Scope:** a fixed linear step. Dense cubing costs about D³ per cube. Going backward needs every eigenvalue nonzero.

### 2. Pascal's pyramid and the mod-3 collapse

- **Paths.** On the 27-cube each turn adds 1 mod 3 to one digit. The 3^L move orders reach 9 positions, because the digit sum mod 3 equals L mod 3.
  - Each position's count is a sum of trinomial coefficients L!/(i! j! k!), grouped by (i, j, k) mod 3 (`pascal_pyramid.py`, `hash_answer.py`).
  - The counts are 63, 81 and 90 at L = 6; 2,106, 2,187 and 2,268 at L = 9; and 58,806, 59,049 and 59,535 at L = 12.
  - Computing each position once under its hash replaces 81, 2,187 and 59,049 paths per position on average.
  - Results that depend on the counts keep the recurrence C_{n+1}(s) = Σ_{p→s} C_n(p) in exact integers.
- **The collapse.** For commuting operators mod 3, (A + B + C)^(3^k) = A^(3^k) + B^(3^k) + C^(3^k), by Frobenius and Lucas's theorem.
  - **Measured** on the 27-ring for k = 1, 2, 3: (I + R ± R²)^(3^k) = I + R^(3^k) ± R^(2·3^k).
  - The 27-step tests (`pow27_mod3.py`): (I + R + R²)^27 = 0 and (I + R − R²)^27 = I (mod 3).
  - **Scope:** arithmetic mod 3 with commuting operators.
- **What the hash does.** In a finite space, a state's SHA-256 alone identifies the state by hashing the candidates: at most 27 here, or 9 with the invariant. An input with no structure, 64 bits, would take 2^64 candidates. At the measured 765,864 hashes per second that is about 764,000 years.

### 3. Freeze and play

A state frozen with its HBI → HBP → SHA → SH → HASH chain can be played forward from later, and backward where the step is reversible. A tree of frozen slices on a ternary schedule lets any later read start from the nearest leaf.

**Four guarantees, stated separately:**
- **Integrity:** every frozen state's bytes match their recorded digest.
- **Provenance:** the root seal matches a seal pinned independently, in a committed GitHub record or an earlier frozen, trusted record.
- **Replay:** resuming from a frozen state reproduces the reference calculation byte for byte.
- **Scope:** backward replay holds for reversible dynamics with every needed part of the state retained.

**Status: FILED.** Final acceptance stays OPEN until both conditions are met:
- **The two MEDIUM holds resolve.**
  - H1: play checks the outside anchor, which the negative control in [`freeze_play/ANCHOR-CONTROL.hbp`](python/matrix_tricks/laws_20261106/freeze_play/ANCHOR-CONTROL.hbp) reads. That control replaces a leaf and its whole internal chain: it plays against the internal chain alone (measured) and must be refused against the anchor.
  - H2: a back-read check independent of the engine's own printed bound.
- **The expanded chain passes.**

The record: all eight holds with their dispositions are in [`freeze_play/HOLDS-DISPOSITIONS.hbp`](python/matrix_tricks/laws_20261106/freeze_play/HOLDS-DISPOSITIONS.hbp). Both lanes' reports are in [`BUILD-REPORT.md`](python/matrix_tricks/laws_20261106/freeze_play/BUILD-REPORT.md) and [`CHECK-REPORT.md`](python/matrix_tricks/laws_20261106/freeze_play/CHECK-REPORT.md).

**Measured** on the polariton quantum engine: a copy at QNN `e4fc6ba` with `polariton_freeze_play.rs` (sha256 `4a903f22…`), freezing a leaf every 9 windows on the tree of 9, 27, 81 and 243. An independent check passed it with holds. The records are [`freeze_play/PREDICTION.hbp`](python/matrix_tricks/laws_20261106/freeze_play/PREDICTION.hbp), committed before the code, and [`freeze_play/RESULTS.hbp`](python/matrix_tricks/laws_20261106/freeze_play/RESULTS.hbp).

- **Play equals the uninterrupted run, byte for byte:** 14 of 14 plays, from leaves 1 to 263 windows before the read, including N = 6561 at K = 1200 and 1600. 31 of 31 leaves written again by a play are the same bytes.
- **Any changed byte is refused:** a byte in the state, the record, the index, the SHA, the SH or the seal, or a truncation, each exits with no row, 10 of 10.
- **A later read from the nearest leaf** (one engine process at a time):
  - at N = 6561, K = 1200: 6.75× faster (15.5 s against 104.9 s) and 39.2× (2.2 s against 87.2 s);
  - at K = 1600: 7.37× and 38.4×.
- **Backward play:** 11 of 11 back reads lie inside the forward read's certified interval.
- **The anchor:** a seal proves a leaf's bytes are intact; it does not prove who wrote them. A leaf with a changed accumulator and its chain resealed plays, and so does one paired with a resealed tree (measured). Play therefore trusts a seal kept outside the leaves: one pinned in a committed record on GitHub, or in the frozen prediction.
- **Projection (DERIVED) to the full space at N = 6561, 10,768,242 states:**
  - leaves of 172 MB each, 2.4 GB to t = 2;
  - t = 2 from t = 0 in about 1.4 h;
  - any later read in at most 8 windows, about 5 min.

**Reproduce and verify:** `bash python/matrix_tricks/laws_20261106/verify_laws.sh`, which `verify.sh` also runs.
- **The chain it checks:** the HBI index `LAWS.hbi`, one row per file with bytes, sha256 and its 16-hex SH; the HBP records `*.out.hbp`; the SHA list `SHA256SUMS`; and the HASH seal `LAWS.seal`, the sha256 of `SHA256SUMS` followed by `LAWS.hbi`.
- **The replay:** it reruns every check and compares the output byte for byte, with the timing and hash-rate fields masked. Any changed byte fails it.

## The model

Kogut–Susskind Hamiltonian with finite group S₃.
- **States.** Links carry S₃ elements. States are gauge orbits: a comb spanning
  tree is fixed to the identity, and the k = LxLy + 1 free loop variables are
  taken modulo global conjugation. Each orbit carries its orbit size (1, 2, 3
  or 6) as a weight.
- **Electric term.** Per link, Σ over the 3 transpositions t of (1 − L_t),
  with single-link spectrum {0, 3, 6}.
- **Magnetic term.** Per plaquette, 2 − χ_std(U_p), with values {0, 2, 3}.
- **Dimension (Burnside):** N = (6^k + 3·2^k + 2·3^k)/6.

| torus | states | flat vacua | operator sha256 |
|---|---|---|---|
| 2×2 | 1,393 | 8 | `137285e0…1d2c` (byte-identical to the pinned HEB-C operator) |
| 2×3 | 47,449 | 8 | `2bd842b0…66de` |
| 3×3 | 10,097,891 | 8 | sectors: `ad0168f9…`, `67992a04…`, `6c53ff4a…`, `d8e07b51…` |

Every torus has exactly 8 flat vacua, the commuting holonomy pairs of S₃ up to
conjugation.

## Layout

```
crates/torusl     operator builder; ymsec (symmetry sectors), onfly (matrix-free), ball (vacuum-ball)
crates/sturmrs    interval Sturm / LDLᵀ inertia (sturmrs), sturmfast, shellschur (magnetic-shell Schur)
crates/evolvers   jump, cheb, cluster, jumpto, fillest, ddcore, flatcomp
crates/jump3      power-of-three jump (jump3) and transposition memo (transmemo)
crates/rainbon    three-preon cloud model: cloudb, cloud, symsec
python/ladder     ladder.py (milli), ladder2.py (micro), rungs.py
python/volume     exact laws, sector validation, trace moments, shell_gap_cert.py, monotone_cells.py
python/matrix_tricks exact_jump.py, vacuum_ball_ritz.py, time_jump_ladder.py, jump3_ntt.py, transposition_memo.py
python/testground HTTP-0 page and routes (127.0.0.1)
python/audit      integer-only AST audit
data/             operators, ledgers, sweeps, volume runs, 3×3 certificates, jump3 predictions and receipts (json=0)
```

Rules for the certified finite-`S₃` core (not the separately labeled exploratory wave/molecule scripts):
- integers only in the proof/certificate path: no float, Decimal, Fraction or float timing;
- outward-rounded integer division in certificates;
- `json=0` on every receipt line;
- Rust 1.81 with `clippy -D warnings`.

Theory and architecture: Jesse Daniel Brown.
