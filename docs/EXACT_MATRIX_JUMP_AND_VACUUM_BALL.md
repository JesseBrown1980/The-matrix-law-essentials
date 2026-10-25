# Exact integer matrix jumps and the 3×3 vacuum-ball pilot

**Reproducible algebraic capability (guaranteed for the specified finite input):**
For an integer square matrix `S`, an integer vector `v` and nonnegative
integer `N`, binary exponentiation computes exactly `S^N v` through
integer-only matrix multiplication. A jump over `2^J` sequential applications
requires `J` matrix squarings and one matrix-vector action. No amplitude
normalization, lossy shift, floating point, division or rounding is used.

This guarantees **equal mathematical integer output** to sequential
applications, not a guaranteed wall-clock speedup: dense squaring needs
O(n^3 log N) scalar operations and intermediate integers can grow rapidly.
The reusable pilot in `python/matrix_tricks/exact_jump.py` deliberately
fails closed for dense matrices bigger than 16×16. It does not square the
47,449-state or 10,097,891-state full Hamiltonian.

## Run the accepted exact-jump control

```bash
python3 python/matrix_tricks/exact_jump.py --self-test
```

For the *independently reconstructed* eight-vacuum **3×3 third-order
effective** integer matrix `16H3`, define the exploratory map
`S=128 I - 16H3`. The test pins the effective source matrix SHA
`69147c597f981eb9af15dc7f0b31848d52ef7190e9cfa4007de75ce2e7a9c4d0`,
then checks:

- 64 sequential steps equal exactly 6 squarings, identical SHA-256
  `c11971469a50925b416ce6edf780a45ee2d216241fbd7514b5cb5ab4877db8b6`.
- 4096 steps via 12 squarings equal two independently requested 2048-step
  jumps composed; state components reach 31,168 bits, SHA-256
  `07bf219237a7100a67ac7d4d3b01472af5777cb30a90efc7d7269118ba479eaa`.
- Other matrices and irregular step counts (including 0, 17 and 65) agree
  exactly with sequential evaluation.
- Negative/noninteger steps and wrong dimensions are rejected.

Every machine result is `json=0`. The state hash is tied to the explicitly
pinned 8×8 effective matrix and its basis/encoding.

## The separate 3×3 vacuum-ball compression pilot

The full periodic 3×3 finite-S₃ gauge model has **10,097,891 orbit states**.
Principal compressions include all eight flat vacua and states reachable by
electric graph moves:

| Radius | States | Stored directed nonzeros | Trial Ritz dimension |
|---|---:|---:|---:|
| 1 | 206 (8+198) | 818 | 16 |
| 2 | 4,454 (8+198+4,248) | 33,218 | 24 |

Independent integer-only block-Krylov Ritz pilot enclosures (for `JE=1`):

| r = JB/JE | radius-1 Ritz gap interval | radius-2 Ritz gap interval |
|---|---|---|
| 9 | [0.016822,0.016824] | [0.019186,0.019188] |
| 27 | [0.002004,0.002006] | [0.002014,0.002016] |
| 81 | [0.000222,0.000224] | [0.000221,0.000223] |

These are **Ritz trial-space estimates**, not spectral-gap certificates of
the full 3×3 Hamiltonian. Individual Ritz eigenvalues upper-bound
corresponding full-Hamiltonian eigenvalues, but their difference does **not**
bound the full gap. An independent exact integer BFS-tree implementation
generated the pilot; the reusable exact-jump tool does not claim to build or
certify these trial matrices. The larger radius-3 compression and a rigorous
bound on omitted states remain future verification tasks.

**Scientific limits:** The jumping test is on `16H3` in the eight-vacuum
effective space, **not** 4096 unitary time steps of the complete 3×3
Hamiltonian. The `S=128I-16H3` filter is artificial; no real-world days,
real-time Yang–Mills field propagation, continuum theory, or Clay
Millennium result is thereby simulated or proven.

The perturbative coefficient
`(9/8)(sqrt(14)-sqrt(6))*JE^3/JB^2` is a distinct, exact effective
Hamiltonian computation still pending independent mathematical review.
The existing full-model 2×2 Sturm certificates and the running ladder are
untouched.
