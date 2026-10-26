# Why translations and reflections commute with the S₃ lattice Hamiltonian

This is the algebraic argument behind `ymsec`. The sampled and exhaustive
checks in `verify.sh` test it numerically: every gauge orbit on 2×2 and 2×3,
plus stride and targeted samples on 3×3.

## Setting

On a periodic Lx × Ly lattice, link variables U_μ(s) ∈ S₃ sit on oriented
links s → s + μ̂. A gauge transformation h acts as
U_μ(s) ↦ h(s) U_μ(s) h(s + μ̂)⁻¹, and physical states are functions of
gauge orbits. The Hamiltonian is

H = JE Σ_links E_ℓ + JB Σ_plaquettes (2 − χ(U_p)),

where:

- E_ℓ = Σ_{t ∈ T} (1 − L_t^ℓ) is the electric term on link ℓ;
- L_t^ℓ is left multiplication U_ℓ ↦ t U_ℓ;
- T = {(12), (13), (23)} is the conjugacy class of transpositions;
- χ is the character of the real 2-dimensional representation;
- U_p is the ordered product of the links around a plaquette.

## The symmetry operations

A point map π of the lattice (translation, or reflection x ↦ −x or y ↦ −y)
sends each oriented link either to an oriented link or to an oriented link
traversed backwards. The induced action on configurations is:

- **(Sπ U)_{π(ℓ)} = U_ℓ** if π preserves the orientation of ℓ;
- **(Sπ U)_{π(ℓ)} = U_ℓ⁻¹** if π reverses it. The reversed link of the reflected lattice is read in its own forward direction, which inverts the group element.

This is exactly what `apply` in `ymsec.rs` implements. Under σₓ, an x-link
from s to s + x̂ maps to the x-link from σₓ(s) − x̂ to σₓ(s), carrying U⁻¹, and
y-links keep their orientation. σᵧ is the same with the axes exchanged, and
C₂ = σₓσᵧ.

## 1. Gauge covariance

If U ~ U' (gauge-equivalent through h), then Sπ U ~ Sπ U' through h ∘ π⁻¹.
A reversed link carries U⁻¹, and (h(s) U h(s')⁻¹)⁻¹ = h(s') U⁻¹ h(s)⁻¹,
which is again a gauge transformation of the inverted link with the endpoints
exchanged. So Sπ maps gauge orbits to gauge orbits, and the sector
construction (gauge-fix, then canonicalise) is well defined.

## 2. The electric term is invariant

**Orientation preserved.** E_{π(ℓ)} acting on Sπ U is E_ℓ acting on U,
relabelled.

**Orientation reversed.** Left multiplication of the *reversed* variable
becomes right multiplication of the original one. In the reflected
configuration the link carries V = U⁻¹, and

- L_t acting on V gives V ↦ tV = tU⁻¹ = (U t⁻¹)⁻¹, which is a right multiplication R_{t⁻¹} on U.

So under Sπ the electric operator of a reversed link becomes
Σ_{t ∈ T} (1 − R_{t⁻¹}). Two facts make this equal to the original:

- **Transpositions are self-inverse,** so t⁻¹ = t and the sum is Σ_{t ∈ T} (1 − R_t).
- **Left and right sums over a whole conjugacy class agree.** For any function f of U, Σ_{t∈T} f(U t) = Σ_{t∈T} f((U t U⁻¹) U) = Σ_{t'∈T} f(t' U), because U T U⁻¹ = T.

Hence Σ_t R_t = Σ_t L_t on every function of U, and E_ℓ is invariant under
link reversal. Summing over all links, which π permutes, Σ_ℓ E_ℓ commutes
with Sπ.

## 3. The magnetic term is invariant

π maps each plaquette to a plaquette. The image's ordered link product is
either a cyclic relabelling of U_p, which is conjugate to it, or, when π
reverses the orientation of the plaquette (a reflection), the inverse of a
conjugate of U_p. Characters are class functions, and the real
2-dimensional character satisfies χ(g⁻¹) = χ(g). So 2 − χ(U_p) is
preserved, and π permutes the plaquette terms.

## 4. Consequences used in the construction

- **Commutation.** Sπ commutes with H on gauge-invariant states, for every element of G = (Z_{Lx} × Z_{Ly}) ⋊ C₂ᵥ.
- **Semidirect structure.** Reflections conjugate translations into inverse translations: σₓ t(a, b) σₓ = t(−a, b). On zero-momentum states (trivial translation character) the four C₂ᵥ characters therefore label a complete set of sectors.
- **Weights.** The conjugation weight of a gauge orbit, its number of distinct tree-gauge configurations, is preserved by Sπ. So the sector weight W = (gauge orbits in the G-orbit) × (conjugation weight) makes the sector operator W-self-adjoint. `ymsec` checks this exactly for every matrix element.

## Numerical cross-checks (in `verify.sh`)

- **Group relations,** including the semidirect ones, on raw configurations.
- **Commutation:** the move targets of g·x equal g applied to the move targets of x, with multiplicity, and the magnetic energies agree. This runs on **every gauge orbit** of 2×2 (1,393) and 2×3 (47,449), and on stride and targeted samples of 3×3 (flat vacua, nontrivial-stabiliser states, near-identity configurations).
- **Flat-vacuum decomposition** 7/0/0/1, as predicted from the reflection action on commuting holonomy pairs.
- **Spectral reconstruction.** On 2×2 the four sectors reproduce the independently built zero-momentum sector level for level, and **exactly reproduce its spectral moments** Tr(Hⁿ), n = 1 to 4, in integer arithmetic, with no eigensolver.
