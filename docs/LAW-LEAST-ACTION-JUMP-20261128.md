# Least-action ternary jump (2026-11-28)

Status: exact mathematical identity; local integer-only tests PASS; hosted verification pending.

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
