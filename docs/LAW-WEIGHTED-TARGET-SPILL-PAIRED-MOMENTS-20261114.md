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

Thus \(n\) *complete* operator applications yield the \(2n\) exact
moments \(\mu_0,\ldots,\mu_{2n-1}\). With \(n=11\), eleven actions
yield 22 moment matrices (44 at two distinct frozen source/energy
pairs). Source columns may be batched; an operator action must not
silently drop newly reached states.

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

## Verified toolkit accelerations and limits

1. **Sparse frontier rather than the global graph.** Keep exact
   nonzero source vectors, enumerate all outgoing neighbors for each
   active source, and append every reached target. Any unaccounted
   boundary edge is a certificate HOLD, not an approximation.
2. **Integer lane packing.** When source entries and move
   multiplicities are nonnegative, pack \(p\) columns in separated
   integer bit fields large enough for the *maximum sum of incoming
   weights*. The width bound itself is part of the proof; unpack,
   re-pack, and compare against an independent row method.
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
  unverified Q projection.
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
