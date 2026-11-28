# Least-action ternary jump (2026-11-28)

Status: exact mathematical identity; local integer-only tests PASS; dedicated hosted replay PASS at commit 83d5382 ([run 46150870015](https://github.com/JesseBrown1980/The-matrix-law-essentials/actions/runs/46150870015)); full repository acceptance separate.

For a finite state space with complete boundary state and additive transition costs, let K_n(i,j) be the minimum cost of an exactly n-step path, or +infinity if none exists. Splitting a path at step n proves K_(n+m)(i,j)=min_k(K_n(i,k)+K_m(k,j)). Conversely any finite minimizing pair concatenates to a legal path. Thus min-plus multiplication is associative, and K_(3n)=(K_n ⊗ K_n) ⊗ K_n. For time-dependent costs, compose chronological blocks; cubing only applies to identical blocks. Negative finite costs are allowed at fixed horizons. This is a classical min-plus semiring identity, not a mathematical priority claim.

Dense full-kernel cost O(d^3 log n) for a homogeneous d-state graph, versus O(d^3 n) sequential dense full-kernel composition; not a guarantee against sparse, single-source or specialized solvers. State dimension, densification, integer bit lengths, and witness output must be accounted for. A minimum is not a sum over paths, a quantum propagator, or a spectral-gap certificate. Forward/backward cost recurrence is not inverse dynamics. An argmin witness can be retained at each split; explicit full-path output costs at least n steps.

Reproduce: python3 python/matrix_tricks/least_action_jump.py (do not use -O).
Local receipt:
```
PASS endpoint_horizon_cases=49140 forward_backward_comparisons=98280
PASS controls=4 horizon=3^40 sha256=994ef829b8095f846f7565b97febc3add77a80b4838a3a911b248c6a9bd9257c
dense_81_frames candidates_direct=1574640 candidates_jump=157464
```
36 deterministic seeded matrices of dimensions 3,9,27, five horizons; separate forward and backward Bellman recurrences. Four controls: greedy trap, unreachable path, noncommuting blocks, uniform-cost 3^40 closed form. Receipt digest is over tested kernels, not a signed or externally timestamped commitment. No Rust or hosted execution claimed yet. SAT application requires complete separator assignments; Yang–Mills requires separate sum-product/transfer or valid operator-bound argument.

## Exact normalized shadow-boundary quotient (local probe)

If A[i,k] = c_i + R[q(i),k] for every k, including identical unreachable support, then (A tensor B)[i,j] = c_i + (R tensor B)[q(i),j]. Normalize each row by its finite minimum, treating all-unreachable rows separately with offset zero. Compare complete normalized rows, not their hash or scalar minimum alone. This is a lossless factorization, not an assertion that arbitrary CNF frontier kernels compress. Constructing the quotient still reads the original matrix; dense candidate slots decrease from d^3 to r*d^2, plus construction and output costs. Witness reconstruction must retain original boundary IDs and argmin splits.

Local integer-only result: 72 matrices, 132840 exact output entries compared with ordinary composition; 81-row constructed examples have 27 normalized classes, 531441 to 177147 dense candidate slots (excluding quotient construction/expansion). Equal scalar-minimum but different rows must not merge. SHA-256: 0925dd6a85765a1bfd4bbacefa51fb1611c806d8b103711b4a905ca35ba8822e. Reproduce: python3 python/matrix_tricks/boundary_quotient.py. This is new local evidence, not inherited hosted or native acceptance.
