# Exact affine-source jumps on very large sparse ternary matrices

**Date:** November 12, 2026. **Verification class:** exploratory local calculation plus exact independent hosted reproduction. The source, reference hashes and dimensions are frozen *after* the local calculation and *before* the new hosted run, not represented as a blind discovery. The mathematical identity is standard linear algebra; the new repository instance applies it to previously unhandled **forced** Rainbon-style dynamical channels.

## The newly needed rule

Our existing [power-of-three jump](LAW-POWER-OF-THREE-JUMP.md) proves that a fixed linear step `x[n+1] = M*x[n]` can be jumped. In the RAINBOR [12-billion-virtual-address optical experiment](https://github.com/JesseBrown1980/RAINBOR/pull/1), inserting a `+7` local source after **every** step caused `M^27*x0` to disagree with the forced recurrence. This is not a defect of the previous jump: the evolution was changed.

The appropriate exact extension, for a source `f` **constant within the time interval**, is

```
x[n+1] = M*x[n] + f
x[n]   = M^n*x[0] + (I + M + ... + M^(n-1))*f.
```

The augmented matrix `[[M,f],[0,1]]` is linear on `(x,1)` and can be powered by binary/ternary squaring, even when `M` is not invertible. For a circulant stencil, the finite-field Fourier transform diagonalizes it. In each mode with eigenvalue `lam` and source coefficient `F`,

```
X[n] = lam^n*X[0] + S_n(lam)*F,
S_n(lam) = (lam^n - 1)/(lam - 1),   lam != 1,
S_n(1)   = n mod p.
```

This holds in **GF(p)**; division means modular inverse, and the `lam=1` branch matters. A second independent implementation uses exact affine-pair binary squaring:

```
(a,b) represents z -> a*z + b*F
(a,b) after (c,d) = (a*c, a*d+b).
```

The two methods must agree, and direct stepping is a third route at feasible step counts.

A changing *externally specified*, piecewise-constant source is still exactly composable. If source `f1` acts for `n1` steps and `f2` for `n2`,

```
X_final = lam^(n1+n2)*X0
        + lam^n2*S_n1(lam)*F1
        + S_n2(lam)*F2.
```

**A state-dependent nonlinear source `f(x)` is not covered.** Substituting it into a fixed-`f` jump produces a different recurrence. The negative control directly applies nonlinear local kicks for 27 steps and rejects that substitution.

## Bigger matrices built and tested

Use `python/matrix_tricks/affine_mega_rainbon.py`, with independent `affine_mega_audit.py` all-row CSR parsing and direct 27-step reconstruction.

The stencil is `M = I+R-R^2` on periodic rings. Prime `p=1,299,079` has `p-1=2*11*3^10` and primitive generator `3`. A radix-three NTT supports all three dimensions exactly. Every matrix is **materialized in CSR binary** (row pointers, column indices, coefficients); no dense `D^2` matrix is allocated.

| Dimension | Dense entries if materialized | Actual sparse nonzeros | CSR file bytes | Expected binary SHA-256 |
| ---: | ---: | ---: | ---: | --- |
| 6,561 = 3^8 | 43,046,721 | 19,683 | 183,738 | `ebe0ee2bc29bf93d…` |
| 19,683 = 3^9 | 387,420,489 | 59,049 | 551,154 | `95f12bb8addcbf88…` |
| 59,049 = 3^10 | 3,486,784,401 | 177,147 | 1,653,402 | `3297505eda8ac2…` |

**Measured local exact-algebra gates:**
- NTT forward-inverse identity and shift-eigenvalue identity at each scale;
- every CSR row and one full CSR-matrix-vector application checked;
- independent `x[n+1]=M*x[n]+f` direct-stepping at **27** steps for all three dimensions, and **729** steps for `D=6561`;
- two independent affine formulations (geometric inverse formula versus affine-pair exponentiation) agree for every short test; the long `3^40` test agrees at `D=6561`;
- `3^40` and `2^40` jumps for every dimension; `3^80` for `D=6561` and `D=59049`;
- freeze at `3^40`, play forward three forced steps, and verify against a fresh jump to `3^40+3`;
- one source applied for `3^40` steps, then a *different* source for `2^40`: independently computed spectral composition matches both staged jumps in every component;
- phase/sourcing null controls: source set to zero equals unforced evolution exactly, and state-dependent kicks falsify the constant-affine shortcut.

The complete **59,049-component state SHA-256** after `3^40` forced steps is `3fedf4d67518bee4742cbe49995e3c7248843647b7f5ca6a9f5c42cf8da10395`; at `3^80`, `11d5ea99a88f2001dfc3ea0c1cf743759ef461cf4e26fccf45ff0b9365741ac1`. The `3^40` source followed by `2^40` altered-source composite hashes to `94a7470283abd6a6fec516b7a7b4b97f5edd266757f6fcf22f69d8f3334fe7e0`.

The test also produces one-, two-, and three-source sampled field readouts at **243 distinct positions** mapped from the **12-billion-address virtual** detector. Those finite-field readouts are **not** the optical Fraunhofer intensities of the earlier RAINBOR HTML program; their source-count labels are diagnostic. Do not conflate them with a 12-billion-particle simulation.

## Complexity

The radix-three NTT performs `O(D log_3 D)` arithmetic per transform, with approximately `O(D)` state and sparse-matrix storage. The source jump costs additional modular powers/inverses for `D` spectral modes. Jumping to `3^80` does not require `3^80` time steps.

The finite-field evolution is periodic or ultimately periodic; a long modular jump is **not** a physical simulation of a universe evolving for the same number of seconds. Full independent nonlinear interactions, gravitational backreaction and real 4D Yang–Mills dynamics are not handled by this rule.

## Reproduction

```bash
python3 -B python/matrix_tricks/affine_mega_rainbon.py --out build/affine-mega
(cd build/affine-mega && sha256sum -c SHA256SUMS)
python3 -B python/matrix_tricks/affine_mega_audit.py build/affine-mega
```

Compare the result to `data/jump3/AFFINE-MEGA-FREEZE-20261112.hbp`. The canonical binary files are attached to the hosted workflow artifact rather than committed directly.

**Physics and proof boundary:** this is a rigorously reproducible **modular linear/affine signal-processing model**, not a Standard Model gluon, glueball or proton simulation, a density-dependent gravity theory, a physically temperature-switched black hole, or a continuum Yang–Mills mass-gap proof. SHA-256 fingerprints states and receipts; it does not create interaction forces. The protected finite Yang–Mills coupling `r=9` remains **HOLD**.
