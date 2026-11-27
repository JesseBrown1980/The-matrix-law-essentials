# Signed residual reuse across addresses

**Status:** exact finite bijection identity; companion local integer research audit passed. New native verification pending. This is an extension of conditional shadow component replay, not a new numbered law or a Yang–Mills certificate.

Let F and G be finite Boolean CNFs. A signed coordinate bijection T maps each source variable to one target variable, possibly complemented. If applying T to every literal gives exactly the target normalized clause set, then F(x)=G(T(x)) with the appropriate direction of the assignment map. Consequently their active-variable counts agree and positive target witnesses map back exactly. Free-variable cardinality factors are handled separately. Reuse is authorized by exact mapped-clause equality, not by a hash or incidence signature alone. The same principle applies to weighted factors if every corresponding factor value and its scope are certified.

The companion implementation uses occurrence profiles and bounded backtracking to find signed permutations. Matching failure or a search cap is a cache miss, not a proof of non-equivalence. Source identity, map, normalization convention, and witness replay must be retained. This is standard syntactic symmetry/isomorphism reuse applied to the shadow architecture, not a priority claim, nor a universal compact generator theorem.

## Evidence and limitations

[P–NP PR #71, pinned full ledger](https://github.com/JesseBrown1980/P-NP/blob/d5e9446a619b08cb54af0f07e3b230b49587109c/native/SIGNED_RESIDUAL_REUSE.md).

- 720 independently signed/readdressed template copies: 695 certified cache hits, 25 generators computed, 720 exact count agreements, 600 reconstructed witnesses. Counting calls fell from 13296 to 465, but mapping required 135735 probes and 155726 clause checks. One capped match fell back to direct computation. No runtime advantage established.
- 16 connected selector CNFs: exact count agreement with independent exhaustive template counts; signed reuse reduced calls 628 to 334 and cached states 306 to 159, with 152 map probes and 437 clause checks.
- 24 random CNFs: exhaustive bit-parallel counts agree, but only four fewer recursive calls (608 vs 604). Eight-round incidence refinement barely improved the repeated-template search and incurred extra work.
- 24 pairs of prism and complete bipartite K3,3 constraints have indistinguishable tested local signatures but different counts (13 and 15). All exact matches rejected; all 24 corrupted signed maps rejected. Marginal structural signatures are insufficient.

Prior companion native automatic replay PASS: [run 45886348707](https://github.com/JesseBrown1980/P-NP/actions/runs/45886348707). New native test covers mapping and corruption controls only; larger bank/selector experiments executed locally. Original Brown–Hilbert/Fischer engines not invoked. No alteration to matrix verify.sh, historical laws, or physical theorem status. P=NP remains open.
