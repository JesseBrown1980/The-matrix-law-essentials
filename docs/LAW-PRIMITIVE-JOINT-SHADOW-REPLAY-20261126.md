# Primitive joint-shadow replay and the address-order obstruction

**Status:** exact finite algebraic identity; local integer audit in the P–NP research track. Native Rust verification of the new companion test is pending. This note concerns finite nonnegative integer constraint functions, not a Yang–Mills spectral certificate or a change to the existing Matrix Laws numbering.

## Law

Let a finite joint shadow function on binary boundary variables be represented by an ordered decision DAG. Each edge stores a nonnegative integer weight and a primitive child node; the terminal represents the constant one. For child edges (a,L) and (b,R), extract g=gcd(a,b) and store (a/g,L),(b/g,R), returning g times the new node. Equal child edges collapse; zero has one canonical representation. Intern identical normalized nodes at the same address level.

Induction proves that this normalization leaves every addressed function value unchanged, and that each nonzero primitive node has leaf-value gcd one. A query multiplies edge weights along the selected address. For count replay, cache the suffix sum of the primitive node once, and multiply by the incoming scale; skipped binary coordinates contribute a factor of two each. A local projection eliminates internal variables by summing the complete joint factor bucket over those variables. Reverse playback uses retained positive choices and original-factor verification. A single positive choice recovers one witness; it does not list all extensions.

The identity is multiplicative weighted decision-diagram normalization. It is a reusable exact matrix-toolkit technique, not a mathematical priority claim. It is compatible with freeze/play, but the new experiment did not invoke the original Fischer or Brown–Hilbert implementation.

## Why ordering matters

For the joint 2x2 count matrix W=[[2,1],[1,2]], the k-pair product has rank 2^k over Q across the cut separating all first coordinates from all second coordinates. The same holds for the equality matrix I_2. Under paired order, each pair adds just three primitive nodes to a common suffix: 3k nodes represent the complete function on 2k bits. Under separated order the distinct left assignments yield 2^k pairwise nonproportional residual row functions, forcing at least that many residual nodes in a scalar-weighted ordered decision DAG. Rank is invariant under invertible transformations applied separately to the two sides, not under arbitrary global coordinate reordering. The obstruction is specific to that cut and representation.

This explains an operational matrix rule: discover a coordinate order or decomposition that preserves coupled addresses, then cache primitive functions separately from incoming integer scales. Never infer joint compatibility from agreeing marginals; never infer polynomial cost from compression of one fixed family. Constructing and verifying the representation, intermediate storage, and integer bit complexity must all be charged.

## Evidence

The P–NP companion [PR #71](https://github.com/JesseBrown1980/P-NP/pull/71) records the [theorem, cost ledger, source hashes and native test](https://github.com/JesseBrown1980/P-NP/blob/113178b4b4fb08b59319c31caf3ab8c80f8eba71/native/PRIMITIVE_JOINT_SHADOW_GENERATOR.md). Its earlier non-affine weighted shadow gate passed [hosted run 45879385528](https://github.com/JesseBrown1980/P-NP/actions/runs/45879385528).

The new local exact audit: 44/56 pair configurations completed, 12 held under a 60,000 allocated-node cap in separated order. At k=32 (64 boundary bits), paired weighted diagram: 96 live nodes, 624 allocated intermediate nodes, versus 1,617 live / 2,769 allocated for the same in-house unnormalized multi-terminal ADD. Count = 6^32. A second source-incidence discovery audit with randomized variable names and a supplied boundary/internal partition completed 42/42, 14,208 address checks and 7,224 verified witnesses. Primitive DP visited 1,800 states versus 10,791 for scale-dependent memoization. All 64 general small-CNF diagram configurations matched reference DPLL. Counts are not runtime or byte savings; no optimized third-party comparison is claimed.

Local source hashes and independent receipts reside in the companion ledger. The Matrix Law Essentials core verify.sh, original law numbering, finite S3 certificates and protected r=9 workflow remain unchanged. P=NP remains open.

## Reusable acceptance conditions

1. Preserve every joint value, nonzero root scale, and boundary address.
2. Check the recurrence independently, including zero, skipped variables and corrupted weights.
3. Retain a positive recovery choice and validate original constraints.
4. Report node count, allocated intermediates, planning, evaluation, verification and integer bit sizes separately; report resource cap as HOLD.
5. Test paired and separated order and connected high-width controls before broader efficiency claims.
