# Law: hash-addressed exact frozen-trial parameter replay (2026-11-14)

**Classification: DERIVED, with an independently executable integer-only toy
model; hosted verification is a separate gate.** This is a **new technique in
Matrix Law Essentials**, not a claim that the elementary algebra is new to
mathematics. It is distinct from the time-evolution power-of-three jump,
the affine forcing jump, and freeze/play of successive system states.

## Exact identity

Let P and Q be fixed finite-dimensional inner-product spaces; their metrics
and operator representation must be **fixed and certified**. Let

    A(z) = QHQ - z I,               F(z) = PHP - z I,
    B = QHP,

with A(z) positive definite on the energy interval being certified and
B independent of z. Fix once and for all a trial map Y:P -> Q, and define

    R(z) = A(z)Y - B,
    U(z) = F(z) - B*Y - Y*B + Y*A(z)Y,
    G(z) = R(z)*R(z).

The star is the adjoint in the correct weighted inner products. Here
G is the Gram with exactly the weights defining the Q space.

At the certified anchor z0, store U0, G0, R0, Y and the moments

    H = Y*Y,
    J = R0*Y + Y*R0,
    C = I_P + H.

Writing delta = z-z0, the **exact** algebra gives

    U(z) = U0 - delta C,
    G(z) = G0 - delta J + delta^2 H,
    R(z) = R0 - delta Y.

No graph enumeration, iterative solve or new Jacobi trial is required for
this *fixed* Y at additional threshold values z. The corresponding
Schur complement is

    S(z) = F(z) - B* A(z)^(-1) B
         = U(z) - R(z)* A(z)^(-1) R(z).

If the whole-Q lower spectral floor certifies A(z) >= gamma(z) I>0, then

    U(z) - G(z)/gamma(z) <= S(z) <= U(z).

Thus every new threshold can be checked from **small exact matrices and a
certified floor**. Strictly positive finite-dimensional z-gaps and
symmetry-sector inertia remain separate, necessary proof steps.

## Exact exterior boundary correction

When the fixed Y is zero on an exterior shell E, then R_E(z)=R_E(z0).
The exterior Gram G_E is **constant as z changes**, even though G as a
whole changes quadratically. Suppose the true operator principal
compression satisfies A_EE(z)<=D_E(z), where D_E(z)>0 is diagonal.
For V=D_E^(-1)R_E, the variational Schur identity gives

    S(z) <= U(z)-R_E*D_E(z)^(-1)R_E.

This bound is for the **whole Q operator**: edges out of E to the next
shell enter the full residual resolvent identity, while the quadratic
form of the E-supported corrector uses only A_EE. Do not replace the
full-Q residual with its truncated-shell version when computing a
lower bound.

For periodic 3x4 finite S3 at r=9, draft #189 has the independently
reviewed operator estimate

    A_EE(z) <= diag(144 + 9m_e - z),

because the weighted electric adjacency has norm <=72. Its follow-on
Yang-Mills draft #190 computes m-conditioned complete shell4 Gram
matrices and the 8x8 anchor moments for reuse. That work still needs
hosted acceptance and separate source review; this law's toy verifier
is **not** a replay of that 2-million-orbit graph.

## Hash-addressed cache: what SHA does

Store the source matrix basis, P/Q selection, W weights, Y, coupling
constants, source floor proof identifier, all anchor moments, exterior
Grams and the allowable interval under a canonical byte encoding.
SHA-256 binds the exact bytes and protects against substituting an
unrelated certificate. A SHA value by itself cannot reconstruct a
missing matrix, certify an unseen orbit or predict an arbitrary
eigenvalue. The speedup comes from the **proved polynomial identity**,
not from changing or guessing a hash.

Changing any of Y, the operator H, the metric W, B, P/Q selection or the
energy-dependent model invalidates the cache. First regenerate and
re-certify the changed coefficient set. A new z value in the permitted
interval, with the same frozen Y and pencil, requires only small-matrix
work plus integrity checks.

The independent source verifier should reject:
- modified cached coefficients, metric, operator, trial Y or source seal;
- an incorrect sign in the linear residual moment J;
- a threshold with gamma(z)<=0;
- missing exterior edge contributions or a failed PSD/Loewner check;
- an extrapolation of fixed-trial envelopes as though Y had been
  reoptimized at the new z.

## Exact runnable small model

The path python/matrix_tricks/frozen_trial_parameter_replay.py uses a
3-state Q, 2-state P *synthetic* integer model, a fixed integer trial Y,
and six preselected thresholds z=0..5. It directly computes the small
Schur complement by a 3x3 adjugate and compares its lower and
exterior-corrected upper bounds using exact integer cross-products,
while separately comparing the polynomial cache to direct recomputation
at every threshold. Its input source SHA-256 is sealed in the hosted
freeze and in the emitted receipt.

Three negative controls refuse a changed source-byte hash, a flipped
sign for J, and a z threshold at the Q floor. This makes the method
runnable and falsifiable without floats or external Yang-Mills data.

    python3 -I -B python/matrix_tricks/frozen_trial_parameter_replay.py

The source model and full operator are tiny, so no meaningful wall-clock
speedup is claimed from this demonstration. What is **derived** is the
cost separation: an expensive operator/source build once, then exact
polynomial replay at other z values using only cached small matrices.

### Proof/evidence boundaries

- **DERIVED:** exact fixed-Y polynomial and full-Q Schur inequalities.
- **LOCAL TESTED:** toy model's 6 exact energy thresholds and 3 fault controls.
- **REQUIRES INDEPENDENT APPLICATION CERTIFICATE:** any concrete
  r=9 Yang-Mills shell4 gram, energy inertia counts, normalized
  E1->E2 Wilson transition or full spectral weight.
- **NOT APPLICABLE:** generic changes to the Hamiltonian or nonlinear
  state-dependent maps without a proved parameter reduction.

This is an exact *parameter* jump law, not an exact general
time-forward jump law.
