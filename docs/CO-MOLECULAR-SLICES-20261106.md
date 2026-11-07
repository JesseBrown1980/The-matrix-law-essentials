# CO molecular interaction through three ternary observation slices

**Status:** Registered predictions before code at `2e4b7e84f4360e6233af54b9e5aaf2ea6ea7a36f`. Hosted execution: see [PR #6](https://github.com/JesseBrown1980/The-matrix-law-essentials/pull/6) and its CI artifact. This is a physical-parameterized **classical** diatomic force-model experiment, not measured molecule trajectories or quantum Yang-Mills.

## What is calibrated and from where

Morse bond model for `12C16O`, with parameters from a CO and N2 spectroscopy-based collision study [Vibrational Energy Transfer in CO+N2 Collisions, Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC8659027/): `r_e=1.128 Å`, `w_e=2169.81 cm^-1`, `beta=2.298 Å^-1`, `D_e=11.07 eV`. This is a fitted potential, *not* an independent ab-initio molecular computation. The initial relative bond displacement is `+0.050 Å` and the initial velocity is zero. Model motion is **classical and one-dimensional**, neglects rotation, translation, environmental interactions, and quantum vibrational levels.

Write `x=r-r_e` in Å and `tau=t/fs`. With `omega=2*pi*c*w_e*(1e-15 s/fs)`, the calibrated equation is

`V(x)=D_e (1-exp(-beta x))^2`
`x''(tau)=-omega^2/beta * exp(-beta x) * (1-exp(-beta x))`.

This eliminates the reduced mass via the harmonic curvature relation `mu=2 D_e beta^2 / omega^2` with consistent unit conversion. The initial energy, energy drift and molecular bond coordinate are all explicit. A fourth-order Runge–Kutta evolution at `h/9` independently approximates the continuum differential equation; the real continuum analytic Morse solution is not claimed. Velocity-Verlet runs at `h,3h,9h` are compared at the same time coordinates.

## Spatial and temporal observation slices

The user-defined **nominal** Planck coordinates use `t_P=5.391247e-44 s` and `l_P=1.616255e-35 m`; only the units, not a fundamental physical Planck clock, are assumed. Reference constants: [NIST 2022 CODATA](https://www.nist.gov/publications/codata-recommended-values-fundamental-physical-constants-2022).

| Level | Temporal sampling | Time per step (fs) | Spatial readout bin width (Å) |
| --- | --- | --- | --- |
| Fine | 3^57 t_P | 0.0846448906955 | 3^45 l_P = 0.000477492268353 |
| Middle | 3^58 t_P | 0.253934672086 | 3^46 l_P = 0.00143247680506 |
| Coarse | 3^59 t_P | 0.761804016259 | 3^47 l_P = 0.00429743041517 |

The spatial bins quantify the *same* molecular bond coordinate readout. They do not form a genuine 3D molecular simulation cell and do not automatically supply a renormalized effective force law.

The common physical observation times are `n*h` for `n=9,27,81,243,729,2187` fine ticks. The sampled bond distance in each bin is rounded to the closest spatial lattice coordinate, with exact bound `|r_rnd-r| ≤ Δx/2`. The time evolution **does not** round the bond coordinate before each force evaluation; this cleanly separates time-integration error from spatial *measurement* resolution.

## Exact law versus numerical dynamics

**Linear harmonic limit:** Taylor expanding the Morse force at equilibrium gives `x''=-omega^2 x`. In normalized state coordinates `(x,v/omega)`, the velocity-Verlet update has matrix
`M(s)=[[1-s^2/2,s],[-s(1-s^2/4),1-s^2/2]]` for `s=omega h`.

The identity `M(s)^9=(M(s)^3)^3` holds algebraically. The script checks its 3- and 9-step version in Python **exact rational Fractions** for the chosen decimal `s`. At the long `3^15=14348907` fine-tick jump (about 1.21 ns), it compares matrix cubing fifteen times with independent binary matrix exponentiation using **65-digit decimal precision**, not a claim of exact real arithmetic. The script also computes a separate high-precision cosine reference, so a numerically exact discrete jump cannot be mistaken for exact continuum dynamics.

**Naive re-timestep:** An independent velocity-Verlet step at `3h` or `9h` is *not* equal to `M(s)^3` or `M(s)^9` at the same physical elapsed time. The exact-block operator is a different effective evolution; the difference is the **relevant renormalization question** rather than a failure of the matrix identity.

**Nonlinear Morse:** State evolution is nonlinear; a fixed 2×2 matrix does not reproduce the full force law. The independent RK4 reference and three Verlet runs expose actual discretization error and the effect of resolution. The prediction `max_error(h) < max_error(3h) < max_error(9h)` is labelled **CONJECTURE** before testing; if not observed, it must be reported as a missed prediction. The 4 negative controls test force-sign, wrong blocking, hash tampering and off-grid observation time.

## Why this is not a continuum Yang-Mills result

The matrix-jump rule is a mathematical method. This CO model is **not** a gauge field, does not use the S3 Hamiltonian, and its spectroscopic parameters are themselves inferred from existing physical measurements. It is useful as a unit-aware, falsifiable calibration of the three-slice methodology. It cannot independently prove spacetime discreteness, molecule-level quantum dynamics, or a four-dimensional SU(2) physical mass gap. Before connecting to [YM4 PR #33](https://github.com/JesseBrown1980/yang-mills-4-dimensional-continum-proof/pull/33), the bridge must use the same gauge-group transfer operator and *gauge-invariant* observables through all three spatial blocks with controlled effective interactions.

## Reproduce

```bash
python3 python/matrix_tricks/co_molecular_slices_20261106.py > MEASURED-CO-MOLECULAR-SLICES.hbp
```

The hosted job `.github/workflows/co-molecular-slices.yml` runs the script, checks the preregistered gates and saves a SHA256SUMS-sealed HBP receipt.
