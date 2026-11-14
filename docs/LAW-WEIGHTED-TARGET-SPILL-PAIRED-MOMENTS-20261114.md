# Law: exact weighted target-spill reduction and paired whole-operator moments

**2026-11-14.** A reusable sparse-integer computation law for
large W-self-adjoint finite operators. This is a mathematical
bookkeeping identity and a source-completeness protocol, **not** a
claim of mathematical priority and **not** a mass-gap certificate.

## Exact identities

Let \(\mathcal Q\) be a finite indexed Hilbert space with positive
diagonal weight \(W\), and \(A=A^{\!*_W}\) be an exact real/integer
linear operator. Let \(R\) be a finite \(Q\times p\) source matrix.
For any nonnegative integer \(j\), define \(V_j=A^jR\). Then

\[
 \mu_{2j}=R^{\!*_W}A^{2j}R=V_j^{\!*_W}V_j,\qquad
 \mu_{2j+1}=R^{\!*_W}A^{2j+1}R
           =V_j^{\!*_W}V_{j+1}.
\]

The paired-vector formula alone uses \(n\) complete applications to
construct \(V_n\) and recover \(2n\) moments, but **the on-support
odd-moment quadratic form below saves the last complete application**.
Precisely, only \(V_0,\ldots,V_{n-1}\) are required:
\(\mu_{2j}=V_j^{\!*_W}V_j\) and
\(\mu_{2j+1}=V_j^{\!*_W} A V_j\). Therefore **\(n-1\)
complete vector expansions yield \(2n\) moments**, plus \(n\)
on-support quadratic/edge evaluations. For degree **11**, the
production target is **10 complete outward vector expansions** and
**22 exact 8×8 moments per frozen energy**, not 11 expansions.
This saves the last and potentially largest vector \(V_{11}\),
while still requiring exact source and edge verification for \(V_{10}\).
Independent direct-power checks may use additional applications
solely as verifier work; they do not count toward the production algorithm.

For a directed raw electric adjacency \(T\) with positive orbit
weights \(w_x\), suppose detailed balance holds for every pair:

\[
 w_x T_{xy}=w_y T_{yx}.
\]

For any partition \(S=\bigsqcup_a S_a\) of source indices, define
the target contribution \(Z_a(y)=\sum_{x\in S_a}T_{yx}R(x)\).
The **exact** target Gram is

\[
 \sum_y w_y\Big(\sum_a Z_a(y)\Big)
                 \Big(\sum_a Z_a(y)\Big)^\top.
\]

It is **generally not equal** to
\(\sum_a\sum_y w_y Z_a(y)Z_a(y)^\top\). The difference consists
of the cross-shard terms. Sorting and reducing records by
**destination key** before forming any weighted outer product
preserves the full operator. The source row's 72-move census and
the independent target's complete incoming sources are different
requirements; both must be verified.

For Q-sector powers of a full Hamiltonian, remove P destinations
**after** enumerating all full 72 raw electric moves and
**before** the next Q operator action. In general
\((QHQ)^2\ne QH^2Q\); unprojected flat-sector reentry is a source
error.

## Extra exact trick: odd moments by an on-support edge form

The odd identity also admits a no-next-vector evaluation:

\[
 \mu_{2j+1}=V_j^{\!*_W}AV_j.
\]

For the finite gauge Hamiltonian
\(X_z=b(72I-T+9M)-z_{\rm num}I\), and source values
\(V_j(x)\) with finite known support \(S_j\),

\[
\mu_{2j+1}
 =\sum_{x\in S_j} w_x[b(72+9m_x)-z_{\rm num}]
      V_j(x)V_j(x)^\top
 -b\sum_{x,y\in S_j}w_xT_{xy}
      V_j(x)V_j(y)^\top.
\]

All 72 outgoing raw moves per source still must be enumerated and
their destinations resolved to prove support completeness. But once
\(V_j(y)=0\) is certified outside \(S_j\), the **quadratic form**
needs no *stored* \(V_{j+1}\) values on the next outer shell.
The inner electric edge term may be evaluated exactly by taking
each *undirected* off-diagonal edge once:

\[
 -b\sum_{\{x,y\}\subseteq S_j}w_xT_{xy}
   \{V_j(x)V_j(y)^\top+V_j(y)V_j(x)^\top\},
\]
plus the self-loop term
\(-b\sum_x w_xT_{xx}V_j(x)V_j(x)^\top\).

**The hosted synthetic must use at least one source vector whose
complete \(A V\) extends to a Q state outside the support of \(V\).**
It then compares the odd moment computed from *support edges only*
to the direct whole-Q product and checks that no exterior vector
was formed in that shortcut. This directly exercises the claimed
outer-shell saving.

**Completeness is a separate proof obligation:** row-sum 72 and
detailed balance alone do not detect an edge removed in both
directions and replaced by self-loop counts. The production verifier
must reconstruct every self-loop and every off-diagonal target from
an independently specified raw 72-move generator. A deliberately
self-consistent edge-to-self-loop tampering case should be rejected
*against that immutable raw move census*. If the purported raw generator
is itself truncated, a self-consistent moment checker is insufficient.

The balanced edge weight is the same in either orientation, by
the proved \(W\)-detailed-balance law. This is a further rigorous
memory/cycle saving for **odd** moments at a known complete support.
It does **not** supply the outer-shell \(V_{j+1}\) needed by the
next even moment \(\mu_{2j+2}\).

In particular, after constructing complete \(V_1=X_zR_z\)
for real \(\mu_2\), the \(\mu_3=V_1^{\!*_W}X_zV_1\) route
may avoid retaining a shell6 vector. Still, every real shell5
row and internal edge must be verified, and the result remains
HOLD until that source computation occurs.

## Verified toolkit accelerations and limits

1. **Sparse frontier rather than the global graph.** Keep exact
   nonzero source vectors, enumerate all outgoing neighbors for each
   active source, and append every reached target. Any unaccounted
   boundary edge is a certificate HOLD, not an approximation.
2. **Integer lane packing.** An unsigned fixed-width packed
   source is valid only when all packed entries are nonnegative.
   The accepted \(R\)-source outer-boundary representation is
   nonnegative after a fixed sign, but \(V_j\) for \(j\ge1\)
   generally has **mixed signs**. Later passes must use an exact
   per-lane bias with proved source-entry and incoming-degree bounds,
   then subtract bias times incoming multiplicity on each target,
   or a separately verified signed-lane representation. Pure unsigned
   packing of mixed-sign \(V_j\) is **forbidden**. Width bounds,
   no cross-lane carry, correct signed unpack/re-pack and independent
   matrix-row checks are hard gates.
3. **Out-of-core radix bucketing.** Store target canonical address,
   source index, weighted reverse multiplicity and target weight in
   a compact fixed-width record. Bucket by deterministic target key.
   Cross-shard coalescence precedes all Gram calculations, permitting
   bounded RAM and streamed integer accumulators.
4. **Freeze and replay.** SHA-256 pins the exact input and compressed
   event bytes; hashes are addresses/integrity checks only. The
   source data and operator identity compute the answer.
5. **Energy parameter reuse.** For \(X_z=bQHQ-z_{\rm num}I\),
   the graph and magnetic data are shared across thresholds.
   When the frozen trial \(R_z\) **also varies with \(z\)**, its
   numerical moment arrays cannot simply be shared; each source
   must be carried separately.
6. **Proof-gated polynomial reconstruction.** After the \(2n\)
   exact moments exist, \(U_n\) is a fixed integer combination of
   them. In the *accepted K3 historical integer-source units*,
   \[
      U_n=U_0-\frac{1}{K_n^2}
                    \sum_{k=0}^{2n-1}c_k\mu_k,
   \]
   where \(C_{2n-1}(x)=(K_n^2-P_n(x)^2)/x\).
   **No extra scale \(b\)** multiplies the accepted-source moments.
   The source-unit conversion must be anchored to the original
   Schur bracket, with old/new Loewner nesting and inertia checks.

### Formal production gates

- Every source row emits the exact mandated raw move count; no
  dropped destination, duplicate edge, incorrect orbit weight or
  unverified Q projection. A separately specified raw generator
  determines the **actual off-diagonal targets and self-loop counts**.
  Row sums and W-balance alone cannot defeat the planted
  self-consistent edge-to-self-loop fault.
- All disjoint source partitions are present, their row ranges
  cover the exact census, and all target-address buckets are present
  before the final target Gram is accepted.
- The immediately prior exact moments reproduce byte-for-byte;
  new moments are checked for symmetry, block-Hankel PSD, and
  spectral localizers
  \(\mu_{k+1}-g\mu_k\succeq0\),
  \(h\mu_k-\mu_{k+1}\succeq0\) on an accepted global spectrum \([g,h]\).
  These consistency gates supplement independent raw move verification.
- Any source-derived numerical Schur enclosure must nest inside
  the previously accepted Schur enclosure and retain its certified
  energy-threshold inertia counts. The original \(\times10^9\)
  source-unit fault is a mandatory planted negative control.

### Inapplicable shortcut

Radix-three NTT time jumping and nonlinear conjugacy methods
in this repository depend on separately proved Fourier/circulant
or transformation identities. The 3×4 periodic S₃ **gauge-orbit
operator has no established Fourier-diagonal reduction** of that
kind. Applying a \(3^k\) jump instead of full Q operator actions
without such a theorem would invalidate this certificate.

## Concrete source application and current limits

[Finite Yang–Mills issue #137](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C/issues/137)
uses the method in [draft PR #194](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C/pull/194).
The first verified exhaustive source partition at 3×4 S₃,
r=9, has **248,152 shell4 orbits**, **17,866,944** raw electric
moves, **14,204,850** signed-source shell5 edge records
(**198,867,900 bytes**), and exact reproduction of both prior
accepted \(\mu_0,\mu_1\) partials.
[Hosted shard-0 PASS](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C/actions/runs/43604542743).

A parallel eight-source-shard production pipeline is
[run 43605546210](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C/actions/runs/43605546210).
**Its full μ2 was not complete at the moment this note was frozen.**
Even a completed μ2 would leave the degree-11 moments μ3..μ21,
normalized E1/E2 Wilson transition and independent full raw
source review open. The original first-excitation Wilson weight
remains **HOLD**. This method does not prove continuum
SU(2)/SU(3) Yang–Mills or a physical mass gap.
