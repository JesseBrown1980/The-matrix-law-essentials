# Exact two-channel integer-lane packing for signed Krylov vectors

2026-11-14 follow-up to the [independent review of law #14](https://github.com/JesseBrown1980/The-matrix-law-essentials/pull/14#issuecomment-6879238689). The earlier packed-integer source trick applies directly to **nonnegative** source columns; later Krylov vectors \(V_j=X^jR\) generally have **mixed signs** because the operator includes negative electric off-diagonals.

## Exact identity with two unsigned lanes

For each source column vector \(v(x)\), use its positive and negative parts

\[
v^+(x)=\max(v(x),0), \quad v^-(x)=\max(-v(x),0), \quad
v(x)=v^+(x)-v^-(x).
\]

For a nonnegative integer adjacency \(T_{yx}\),

\[
 (Tv)(y)=\sum_x T_{yx}v^+(x)-\sum_x T_{yx}v^-(x).
\]

Pack each of the two **nonnegative** parts independently in eight
integer digit fields of width \(w\). For a complete 72-move row and
source absolute bound \(M=\max_{x,k}|v_k(x)|\), choose

\[
 2^w>72M.
\]

Every target's positive *and* negative channel has each digit at most
\(72M\), so ordinary integer addition of packed words is exact and
has **no cross-lane carries**. The two channels are reduced by
**canonical target key across all source shards** before unpacking.
The final signed result is the difference of independently unpacked
positive and negative lane values.

For \(X=b(72I-T+9M_\mathrm{mag})-z_\mathrm{num}I\),
apply the signed electric \(Tv\) as above and combine the diagonal
term \([b(72+9m_x)-z_\mathrm{num}]v(x)\) exactly in ordinary signed
integer arithmetic. The diagonal may be signed and cannot silently
enter the unsigned packing routine.

**After** this signed reconstruction, take W-weighted outer products.
Squaring positive and negative channels separately loses cross terms
and does not equal the signed-vector Gram. Likewise, squaring source
shards before merging by target omits cross-shard terms.

## Instrument checks

The independent synthetic verifier
\`python/matrix_tricks/signed_two_channel_packing.py\` uses signed
eight-column values, several source shards feeding common targets,
nontrivial weights \(W\), and a strict 72-incoming-move maximum.

It requires exact equality with a separately written scalar
accumulation for every target and for the complete weighted 8×8
Gram. It refuses inadequate width **inside the reducer** before any
word accumulation, refuses a target with more than 72 reverse moves,
detects a changed source sign, and rejects the missing cross-shard
Gram terms from premature outer products.

This technique is suitable for source iterations \(V_1\to V_2\to
\cdots\) **only after** complete raw-move enumeration, original
source W balance, valid state projection, and a certified bound
on all incoming contributions. If the actual source has signed
coefficients outside the assumed domain, or target degrees exceed
the bound, the implementation must produce a resource/proof HOLD,
not truncate or wrap.

## Degree-11 Yang–Mills context

The real 3×4 finite S₃ accepted-K3 source has
\(\mu_0,\mu_1,\mu_2\) certified conditionally on the archived source;
the [complete mu2 run 43606393196](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C/actions/runs/43606393196)
accepted the 18,598,040-state shell-5 frontier. The
[ten-outward-pass terminal odd-law correction #16](https://github.com/JesseBrown1980/The-matrix-law-essentials/pull/16)
is merged. The [V1 shell5 recovery draft #195](https://github.com/JesseBrown1980/Finite-universe-Yang-mills-existence-and-mass-gap-using-HEB-C/pull/195)
reuses old source-addressed event artifacts for the next on-support
moment. This signed-packing law is a reusable verified helper, not
a claim that real \(\mu_3\) or the first Wilson weight has been computed.
