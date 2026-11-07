# 2026 Vienna nanoparticle interference: three-slice wave model and 3^16 temporal jump

**Publication:** [Pedalino et al., *Nature* 649, 866–870 (2026)](https://doi.org/10.1038/s41586-025-09917-9), "Probing quantum mechanics with nanoparticle matter-wave interferometry", first published January 21, 2026. **Original open data and analysis package:** [Zenodo 17502163](https://doi.org/10.5281/zenodo.17502163), available from the authors. Our current retrieval did **not** obtain or execute that archive. The paper demonstrates quantum interference by **sodium nanoparticles** containing >7,000 atoms, not a biological virus, and **does not** implement quantum cloning.

## Registered scope, pre-execution

Predictions were frozen in two sequential commits **before** implementation and before extension: [original 12 rows](../data/vienna_nanoparticle/PREDICTIONS-20261106.hbp) at `9e688d621ce87a44a99c07a4bc973ab6ae24560c`; [43-million extension 8 rows](../data/vienna_nanoparticle/PREDICTIONS-43046721-20261106.hbp) at `1245a6b984c9c82b62b955b694626b1f1306bc08`. The follow-on executable is `python/matrix_tricks/vienna_nanoparticle_slices_20261106.py`. Record misses rather than altering those frozen predictions.

**Important difference:** this is a declared **toy coherent phase-only three-grating model**, not a first-principles reproduction of the Vienna experiment's Wigner–Weyl model. It does **not** model the actual UV laser photodepletion efficiencies, nanoparticle beam mass/speed distributions, gas collisions, dipole forces, detailed grating power response, or the authors' classical ray-optics comparison. It cannot be used to establish a fit to the paper's raw data. The toy model is a mathematical check of Fourier slice and operator-jump methods.

## Literature anchors versus chosen model parameters

| Quantity | Value | Status |
|---|---:|---|
| Cluster mass in representative interference scans | 172 kDa | **MEASURED (paper)** |
| Cluster diameter | ~8 nm | **MEASURED (paper)** |
| Grating period and separation | 133 nm, 0.983 m | **PAPER APPARATUS** |
| Reported visibility at one laser-power set | 0.10 ± 0.01 and 0.08 ± 0.01 | **MEASURED (paper)** |
| Middle laser power | 15.2 ± 0.3 mW | **MEASURED (paper)** |
| Representative beam speed | 160 m/s | **CHOSEN MODEL PARAMETER**; no claim that this is the scan's measured speed distribution |
| First and third amplitude masks | 0.85 + 0.15 cos(2πx/d) | **CHOSEN MODEL** |
| Second phase grating | exp(i φ cos(2πx/d)), φ=0,1,2,3 | **CHOSEN MODEL**; φ is not calibrated to P₂ |
| A single coherent pure matter wave | Yes | **CHOSEN MODEL**; real Talbot-Lau sources require mixed-state and decoherence treatment |

The de Broglie wavelength calculated from the chosen mass and speed is about 1.449968×10⁻¹⁴ m; its Talbot length is about 1.219958 m. The measured grating geometry gives L/L_T≈0.805766 for the chosen speed. One 0.983 m flight takes 0.00614375 s in this simplified geometry.

## Quantum phase evolution and spatial slices

The Fourier spectrum of the first amplitude grating has modes m=−1,0,+1. The second phase grating uses the identity
`exp(iφ cos θ) = Σ_k i^k J_k(φ) exp(ikθ)`, truncated at |k|≤18. Each of two free-flight segments multiplies Fourier mode m by `exp[−iπ(L/L_T)m²]`. This is coherent Schrödinger phase evolution for the chosen reduced model. The output density `|ψ(x)|²` is sampled at spatial grids N=9,27,81,243,729 over one grating period. Fringe visibility is defined from the first density harmonic and an assumed third-mask contrast.

The truncated Fourier wave has |m|≤19, and the intensity has harmonic frequencies |m|≤38. By discrete Fourier orthogonality, an N-point uniform grid with **N≥40** gives the **exact** first density Fourier harmonic (in mathematical exact arithmetic) of **this truncated model**. N=81,243,729 and N=43,046,721 therefore have identical spectral readouts to the finite Fourier sum; N=27 is checked numerically and may incur small nonzero aliasing, and N=9 is explicitly retained as a low-resolution red control. Any discrepancy of the actual physical apparatus is outside this finite-mode certificate.

### 43 million is a virtual spectral observation, not an allocation

The user-requested upper bound is 45 million; the greatest power of three below it is **3^16=43,046,721**, while `3^17=129,140,163` is too large. The test does **not** construct 43,046,721 complex amplitudes, simulate 43 million particles, or read 43 million physical detector samples. Instead it proves the no-alias statement from the **finite Fourier bandlimit** and computes the virtual first-harmonic readout directly. It records `stored_samples=0` on `VIENNA_VIRTUAL_SLICE`.

## 43-million-step matrix jump: physically calibrated test interval

For each Fourier momentum mode m, the free-flight evolution is a scalar diagonal transfer eigenvalue `z_m(δt) = exp[−iπ(L/L_T)m² δt/T]`, where `T=L/v`. Raising the **same** transfer factor to `3^16` via sixteen cubings reaches the full physical flight duration `T=0.00614375 s`, from an elementary substep of about 1.427×10⁻¹⁰ s. An independently computed one-shot phase serves as the reference; 75-digit Decimal complex arithmetic verifies agreement to the registered numerical precision. The experiment tests modes m=1,2,3,5,11 for every power `3^k`, k=1..16, totaling 80 comparisons.

This **is** a valid power-of-three matrix jump because the free-flight operator is diagonal in the Fourier basis and identical for every substep. It does **not** skip 43 million nonlinear interactions, model quantum cloning, or prove continuum Yang–Mills. The exact mathematical identity does not make finite-precision Decimal arithmetic exact real arithmetic, and the representation of physical time uses an assumed constant velocity.

## Negative controls and physical limits

The code rejects a wrong mass/Talbot ratio, a wrong spatial grating period, an omitted middle phase grating and aliased 9-point sampling. It also records four limits/red tests for incorrectly labelling 3^16 as 45 million, claiming the virtual grid is a massive allocated state, calling aliased 9-site data exact, and calling a toy model experimental reproduction.

**How this helps the continuum work:** the phase propagation provides a known test of exact operator jumps; the Fourier certificate provides a known test of coarse observables. A gauge-invariant spectral correlation test in **continuous SU(2) Yang–Mills** is still separate and unperformed. Even perfect agreement of this toy wave model with published fringe visibility would not prove fundamental spacetime discreteness or the 4D mass gap.

## Reproduce

```bash
python3 python/matrix_tricks/vienna_nanoparticle_slices_20261106.py > MEASURED-20261106.hbp
```

GitHub workflow `.github/workflows/vienna-nanoparticle-43m.yml` publishes the numerical receipt with SHA-256 verification as an artifact, and subsequently checks the complete committed receipt byte-for-byte when one is present.

## Actual hosted results — measured November 7, 2026 UTC

The original two prediction freezes preceded implementation and extension. [Hosted GitHub Actions run 42328971963](https://github.com/JesseBrown1980/The-matrix-law-essentials/actions/runs/42328971963) **PASS**, 59 HBP output records preserved in [MEASURED-20261106.hbp](../data/vienna_nanoparticle/MEASURED-20261106.hbp). The workflow also executes 8 red checks: 4 perturbations of the toy physics and 4 false interpretations/alias claims. A later replay re-compares every record byte-for-byte to the committed receipt.

| Middle-grating phase φ, chosen model radians | Toy visibility at 43,046,721 virtual samples | Physical comparison status |
| ---: | ---: | --- |
| 0 | 0.010525865726 | **NOT A FIT** to measured laser power |
| 1 | 0.095778254134 | **NOT A FIT** to measured laser power |
| 2 | 0.107343875060 | **NOT A FIT** to measured laser power |
| 3 | 0.039192841218 | **NOT A FIT** to measured laser power |

Every 43,046,721-point virtual readout is obtained from the truncated finite Fourier coefficients with `stored_samples=0`, not by sampling tens of millions of particles. The explicit 9-site test at φ=3 instead read `0.046113860849` versus the fine spectral value `0.039192841218`, a discrepancy of **0.006921019631** due to aliasing. This red check demonstrates the importance of finite spatial resolution.

The temporal jump covered **43,046,721 hypothetical substeps** over a single model flight of **0.00614375 s**, each about `1.427228336392916e-10 s`, with **16** power-three cubings. 5 Fourier momentum modes were verified independently at each k=1..16: **80 comparisons**. Maximum amplitude discrepancy versus the separately calculated whole-flight phase was `6.100616890818564e-68` in 75-digit Decimal arithmetic. That number is **not** an experimental accuracy and is not a calibrated physical prediction.

**Context for the real published data:** the *Nature* paper's reported two independent measurement visibilities, `0.10±0.01` and `0.08±0.01`, are at the specified UV grating laser powers; they **cannot** be compared to arbitrary φ values without a verified laser-power-to-phase map, photodepletion efficiency, mass and velocity distributions, and decoherence corrections. A numerical proximity between a toy φ prediction and a measured visibility is not an experimental fit.
