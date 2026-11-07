# Planck-coordinate to molecular-scale observation — 6 November 2026

**Scope:** An exact integer, *specified* 45-cell periodic transport model. This experiment uses the actual `jump3_ntt.py` ring/NTT/power-jump API in this repository; it is **not** a simulation of an actual molecule, a molecular vacuum, the 96³ D3Q27 engine, or a test of Planck discreteness in nature. The test keeps evidence and limitations together.

## Measurement contract

Use CODATA nominal numerical Planck values, approximated here as
`lP = 1.616255e-35 m` and `tP = 5.391247e-44 s`. Their nominal decimal digits are a **coordinate assignment**, not evidence that a Planck second is a measured fundamental simulator tick.

- One **coarse spatial cell** = `3^54 lP = 58149737003040059690390169 lP ≈ 0.939848031798485 nm`.
- One **coarse temporal tick** = `3^65 tP = 10301051460877537453973547267843 tP ≈ 0.555355127853016 ps`.
- Hypothetical microscopic clock: one underlying fine-cell motion every `3^11 = 177147` Planck time counts. Then `3^65` Planck time counts exactly correspond to one coarse cell of `3^54` Planck lengths.
- Effective modeled transport speed is approximately **1692.337 m/s**. This is a chosen calibration, not a measured molecular velocity.
- The 45-cell periodic *coarse* ring would span approximately **42.293 nm**, not 45 Planck lengths.

The clocked microscopic rule is formal: `clock` increments on every fine tick and a fine-cell position changes on clock rollover. It need not be instantiated as an enormous 45×3^54-site array. An independently stepped 3^3-clock / 3^4-cell miniature verifies the exact clock arithmetic.

## Model and comparison

For the ideal noiseless environment, `R` shifts a 45-entry vector by one cell at each **coarse** tick:

`x_(n+1)[j] = x_n[(j-1) mod 45]`.

The finite-field computations are in `F_1000081`, using the existing exact NTT jumper. The independent mathematical continuum comparator samples the transport PDE `∂_t u+v ∂_x u=0` on the 45 grid coordinates. With the **same prescribed constant-speed rigid translation** and data sampled after integer numbers of coarse ticks, `u(j dx,n dt)=u_0((j-n)dx)`; hence perfect agreement is an identity of the two specified transport laws.

**This exact match is not independent evidence of physical equivalence to molecular quantum dynamics.** It does show the jump computes an exact sampled observable in the particular model.

## Registered result and negative controls

Locally executed, using an independent implementation of the published jump3 API; **GitHub Actions will execute the exact upstream repository implementation before acceptance**:

| Arm / check | Observation | Interpretation |
|---|---|---|
| Unperturbed background | Zero remains zero | Ideal mathematical control |
| Isolated occupancy-2 pulse | One-cell translation at each coarse tick, preserved exactly | Sampled classical advection agrees |
| One-cell occupancy-1 disturbance added at cell 13 | Exactly one extra cell, translated to `(13+n) mod 45`; no amplification | Separates an injected disturbance from deterministic dynamics |
| Long exact jumps | Coarse steps `3^10`, `3^15` corresponding to `3^75`, `3^80` assumed Planck time counts; NTT state agrees with sampled advection | Exact *modular* jump, not molecular evolution |
| Clock off-grid controls | 3/3 rejected, no rounding | Mapping is not applied outside registered observations |
| Tampered-state hash | 1/1 rejected | SHA-256 commits the modeled state |
| Diffusion red control | Alternating 54-site mode becomes exactly zero in **one discrete diffusion step**; continuum heat semigroup retains nonzero amplitude `exp(-π²/4)>0` | Strict continuum-equivalence counterexample |

The negative diffusion comparator uses a *different* operator, deliberately: `[1/4,1/2,1/4]` on the ring; a continuum heat equation with `ν Δt/Δx² = 1/4` gives a positive Nyquist-mode amplitude. Discrete stability or exact integer conservation does not itself prove continuum equivalence.

The 45-site state hashes computed from 8-byte big-endian entries are unchanged under reprojection into the SVG camera stills. The Python snapshot renderer checks each state SHA-256 again before drawing. `shot_*.svg` figures and the HBP receipt are uploaded by the GitHub workflow as one named Actions artifact. **These are simulator screenshots, not physical camera recordings.**

## The vacuum issue

A run without injected random perturbations is a **noise-free computational reference**. That does not mean a physical quantum vacuum has no zero-point field fluctuations, interactions, boundaries, or measurement back-action; nor is environmental noise present in every real measurement by assumption. The experiment intentionally separates ideal-model dynamics from one deliberately applied perturbation. It cannot settle quantum-vacuum physics.

## Reproduction

```bash
python3 python/matrix_tricks/molecular_scale_jump.py > MOLECULAR-RESULTS.hbp
python3 python/matrix_tricks/render_molecular_shots.py MOLECULAR-RESULTS.hbp shots/
```

Check `.github/workflows/molecular-scale-jump.yml` for hosted direct-vs-NTT and SVG checks. No changes to matrix-jump law, Hamiltonian proofs, HTTP-0 physical engine, continuum Yang–Mills results or numerical acceptance status are implied.

**Next actual-physics test:** calibrate an observable of the real D3Q27 collision-and-streaming dynamics using a known `Δx`, `Δt`, ensemble, and transport coefficients; compare coarse-grained density/momentum correlations versus continuum hydrodynamic predictions across at least three lattice resolutions, including the shuffled-noise and finite-cutoff controls. The current idealized shift model supplies a baseline, not that evidence.
