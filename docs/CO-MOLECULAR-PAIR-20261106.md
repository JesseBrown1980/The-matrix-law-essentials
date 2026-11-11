# Two CO molecules: whether the first pair erases discrete/continuous error

**Experiment:** An exact preregistered *classical harmonic CO stretch model*. This is not a measurement of two actual molecules. The same CO harmonic constants were used as in [the molecular slice study](CO-MOLECULAR-SLICES-20261106.md). The model is chosen to test the hypothesis "adding more than one molecule makes the errors disappear".

**Chronology:** predicted at commit `6734e0c0aba6cb79292b427ad6a6491e389e5534` before code. [Hosted run 42328026246](https://github.com/JesseBrown1980/The-matrix-law-essentials/actions/runs/42328026246) passed the preregistered checks. Its actual output is [the 16-row receipt](../data/molecular_co/MEASURED-CO-PAIR-20261106.hbp). Four false controls rejected. Tests at 2,187 fine ticks (~185.118 fs) and `3^15` fine ticks (~1.21456 ns).

## Both phase arrangements, and an actually coupled pair

The time-step matrix for one harmonic oscillator, `M_h`, is evaluated at `n` fine ticks. The independent continuum reference is a harmonic cosine. Let `e=x_h(n)-x_cont(n)` be their nonzero individual difference.

| Arm | At ~1.21456 ns: signed mean error | At ~1.21456 ns: individual RMS error |
|---|---:|---:|
| Two separate CO stretches, identical phase | +0.01696656816 Å | 0.01696656816 Å |
| Two separate CO stretches, opposite phase | 0 Å | 0.01696656816 Å |
| Two coupled stretches, opposite phase | 0 Å | 0.01874201562 Å |
| Two coupled stretches, coherent phase | +0.01696656816 Å | 0.01696656816 Å |

The **two actually coupled** trajectories use `x_1''=−ω²x_1−κ(x_1−x_2)` and `x_2''=−ω²x_2−κ(x_2−x_1)`, with `κ=ω²/9` a **hypothetical, not experimentally calibrated** CO–CO interaction. The coherent mode has `x_1=x_2`, so the relative coupling vanishes and both share the same time-discretization error. The anti-phase mode `x_1=−x_2` has `ω_anti²=ω²+2κ=(11/9)ω²`. Both coupled direct integration and normal-mode matrix powers matched at ~185 fs within `3e-64 Å` (65-digit Decimal precision). Their continuous reference errors remained nonzero.

### Why 2 is enough for mean cancellation but not for physical exactness

For two molecules with *identical phase*, the error vector is `(e,e)`, so `mean=e` and `RMS=|e|`.

For *opposite phase* the error vector is `(e,−e)`, so `mean=0` and `RMS=|e|`. The average is zero even when `e != 0`.

The same applies to the coupled anti-symmetric normal mode, with the caveat that its frequency and thus its discretization error change when the coupling is turned on. **A smooth or zero mean is not a demonstration of zero per-molecule error.**

The `3^15` matrix cubing and independent binary powering agree to `1.45e-60 Å` at 65-digit precision for the anti-symmetric mode. This verifies the *jump implementation*. The difference against the continuous harmonic evolution at that time is still `0.01874201562 Å`. Correctly jumping an approximate discrete operator does not make that operator an exact continuum evolution.

## Consequence for the real-world claim

This test validates the **mathematical possibility** of two-molecule cancellation. It **does not** establish that two physical CO molecules must be anti-phase, have this precise force law, or reproduce molecular quantum correlations. Real molecular averaging involves thermal fluctuations, electromagnetic coupling, nuclear/electronic quantum dynamics, and experimental uncertainty, none of which are included here. One molecule can already produce smooth data after averaging over time or detector bandwidth, and two molecules can produce coherent error.

**Therefore:** "some two-molecule averages can look continuous" is demonstrated for this model; "the real world becomes continuous-looking as soon as there are two molecules" is **not established**. Neither observation proves that spacetime is physically discrete or that a four-dimensional SU(2) Yang–Mills mass gap persists under continuum scaling.

## Reproduce

```bash
python3 python/matrix_tricks/co_molecular_pair_20261106.py > MEASURED-CO-PAIR.hbp
```

The workflow `.github/workflows/co-molecular-pair.yml` checks four red controls and byte-for-byte fidelity to the published hosted receipt.
