# Law: Schur divided difference as an exact cross-state Gram jump

**Date:** 2026-11-20. **Type:** exact operator algebra, with independently executed integer-only toy in [finite Yang–Mills PR #348](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C/pull/348). This is a new application/toolkit combination in this repository, not a claim of original mathematical discovery. The actual original finite first-Wilson source remains **SOURCE-CONDITIONAL / FORMAL HOLD**; 4D continuum remains OPEN.

## One stored pencil, two energy leaves, and an exact Gram identity

Fix a **W-selfadjoint** physical finite gauge Hamiltonian decomposed as a flat P and a nonflat Q:
\[
H=\begin{pmatrix}A&B^\ast_W\\B&D\end{pmatrix},\quad
R(z)=(D-zI)^{-1},\quad S(z)=A-zI-B^\ast_W R(z)B.
\]
For two distinct energies \(z_1<z_2<\inf\sigma(D)\), the resolvents commute and obey
\[
R(z_2)-R(z_1)=(z_2-z_1)R(z_2)R(z_1).
\]
Consequently
\[
\boxed{I+B^\ast_W R(z_2)R(z_1)B=
       \frac{S(z_1)-S(z_2)}{z_2-z_1}.}
\]
The left matrix is the **entire physical P+Q cross Gram** of the dressed vectors \((p,-R(z)Bp)\), in P coordinates. It is recovered from just **two small Schur matrices**; no new Q-overlap source solve is needed when the two Schur matrices are already authenticated.

This is an exact *parameter jump* between two frozen energy checkpoints. It is compatible with [the frozen-trial replay law](LAW-FROZEN-TRIAL-PARAMETER-REPLAY-20261114.md): the same source/basis/W matrix must be pinned. It is different from the circulant integer NTT time-jump, nonlinear conjugacy, and unproven full-gauge evolution compression.

## Exact observable centering without a new Q Gram

For coordinate-diagonal O, the source-specific Feshbach observable is
\[
M_O(z_2,z_1)=O_P+B^\ast_W R(z_2)O_Q R(z_1)B.
\]
Linearity implies
\[
\boxed{M_{O-cI}=M_O-c\,\frac{S(z_1)-S(z_2)}{z_2-z_1}.}
\]
For true eigenstates \(S(z_i)p_i=0\), the correction has zero \(\langle p_2,\cdot p_1\rangle_W\) matrix element, exactly by orthogonality. For **fixed approximate trial vectors** its Gram correction is **not zero**; dropping it makes a wrong observable source. A separate rational three-state exact toy in #348 checks true and nonorthogonal cases.

In the original 3×4 finite S₃ source, the seven-loop Wilson character has global range \([-7,14]\). The old doubled integer shift \(O_9=2O-16I\) has \(\|O_9\|\le30\). The minimax center is \(7/2\), hence \(O_{\rm opt}=2O-7I=O_9+9I\) with \(\|O_{\rm opt}\|\le21\). This is a **30% reduction in generic norm** while preserving the *true* transition. Re-estimating first-order directional residuals under the new source center still requires their own mathematical checks.

## Exact polynomial error gates at the existing original Wilson energy windows

For \(\gamma=94.068200\), whole-Q upper \(468\), and accepted original energy upper endpoints \(E_1^{hi}=69.909586556\), \(E_2^{hi}=69.949433884\), set \(g_i=\gamma-E_i^{hi}\), \(h_i=468-E_i^{hi}\). A Chebyshev inverse polynomial degree \(n-1\) can be chosen at each actual energy, with error \(\epsilon_i\le 2((\sqrt h_i-\sqrt g_i)/(\sqrt h_i+\sqrt g_i))^n/g_i\); these endpoint choices dominate every energy in the certified windows.

With the accepted physical \(\|B\|^2=72\), an observable norm \(C\), and a two-factor resolvent telescope,
\[
\|\Delta M\|\le72C(\epsilon_2/g_1+\epsilon_1/g_2+\epsilon_1\epsilon_2).
\]
The independently run **integer-only outward-rational** verifier reports:

| Method | C | Error upper (parts per 10¹²) |
|---|---:|---:|
| Cheb n=14, old center |30|13,006,928,810|
| Cheb n=14, optimal center |21|9,104,850,167|
| Cheb n=15, old center |30|7,864,055,405|
| Cheb n=18, old center |30|1,739,136,977|
| Cheb n=27, old center |30|18,831,448|
| Geometric \(3^4=81\), old center |30|396,756,589|

The source-conditional original signed O9 acceptance margin is \(9,285,530,088/10^{12}\). Thus the centered *degree 14* error envelope is smaller than that reference margin. **This is an operator approximation/error statement, not the actual signed polynomial Wilson matrix element.** Any independent polynomial certificate must also calculate/bound the *source-authenticated signed action*, true projector angles and full-state normalization. Polynomial degrees 14–27 could touch a huge fraction of the original 2.18-billion-state Q space and are not automatically cheaper than reusing the existing two Q-dual residual source artifacts.

The exact radix-three geometric identity \(\sum_{j<3^k}T^j=\prod_{\ell<k}(I+T^{3^\ell}+T^{2\cdot3^\ell})\) gives the conservative 81-term comparison. A direct fast NTT jump on the nonlinear S₃ gauge-orbit graph has **no** established translation-invariant/circulant proof, and is not authorized.

## One pass over the *existing* certificate proofs

1. Freeze and authenticate original H, W, eight P flats, 264 P/Q reciprocal links, B*B=72I, Wilson O, positive Perron Q whole-space lift.
2. Recompute the 8×8 signed four-moment Schur matrices and exact two eigenbranch/projector enclosures with source-pinned frozen-trial parameter replay and first-principles matrix-derived 7D compressed separation.
3. Reproduce the two original Q-dual residual/action W-bilinears from the already authenticated K3/Q4/Q5 source and record both forward/adjoint directions, including exterior closure.
4. Close on one **strict positive integer subtraction** and complete P+Q norm denominator. That mathematical implication is already source-conditional; its primary source binding is the remaining acceptance work.

The engineering pattern is **three independently pinned leaves (physical source/floor, small Schur/projector, Wilson Q-action) followed by one exact fan-in**, with altered-premise negative controls and SHA receipts. The proof power comes from the matrix identities and authentic physical source, not from hash names or compression by itself. The main finite [research PR #348](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C/pull/348) has the runnable integer test and detailed acceptance notes.
