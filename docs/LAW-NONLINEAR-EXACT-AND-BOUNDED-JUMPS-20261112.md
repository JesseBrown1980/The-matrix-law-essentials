# Exact and certified long-time jumps for two genuinely nonlinear interacting maps

November 12, 2026. [Executable exact Python](../python/matrix_tricks/nonlinear_exact_and_bound.py), [separately written independent verifier](../python/matrix_tricks/nonlinear_independent_audit.py), [post-local pre-hosted frozen hashes](../data/jump3/NONLINEAR-JUMP-FREEZE-20261112.hbp).

**Classification:** two established, narrowly stated mathematical mechanisms. Neither is a general-purpose shortcut for arbitrary nonlinear evolution. Both are exact integer computations, not a quantum field or physical Yang–Mills simulation.

## 1. A fully coupled nonlinear recurrence with an exact arbitrary-time jump

Work over `GF(1,299,079)` on a periodic `D=3^k` ring, with two fields `u_i,v_i` and neighbor shift `(Ra)_i=a_(i-1)`. Define a *nonlinear, invertible, local shear*:

```
T(u,v) = (u, v + 7*u_i*u_(i-1))
T^(-1)(a,b) = (a, b - 7*a_i*a_(i-1))
```

Define the reversible **coupled** linear shears on transformed coordinates:

```
b_new = b + R*a
a_new = a + R*b_new
L(a,b) = (a_new,b_new).
```

The actual nonlinear recurrence on the physical `(u,v)` variables is

```
F(u,v) = T^(-1) L T(u,v).
```

Its first coordinate depends on the bilinear neighbor product `u_i*u_(i-1)`, and its second depends on new bilinear products of the changed first coordinate. Thus it is genuinely **state-dependent, nearest-neighbor, nonlinear and two-way coupled** in these variables. Its tested superposition control fails: `F(x+y) != F(x)+F(y)`. Its reverse one-step map is exactly computable.

**Exact theorem:**

```
F^n = T^(-1) L^n T    for every integer n >= 0.
```

Proof: multiply `n` copies of `T^(-1)LT` and cancel adjacent `TT^(-1)`; no approximation or convergence assumption is used. The Fourier transform decomposes `L` into exact 2×2 mode matrices

```
L(w) = [ 1+w^2   w ]
       [ w       1 ]
det L(w) = 1.
```

The determinant identity proves reversibility in each mode. Binary powering or repeated cubing computes `L(w)^n` exactly modulo the field prime and reconstructs the physical nonlinear state through `T^(-1)`.

### Period, mode census, and what is inherited from the linear system

The physical map `F=T^(-1)LT` is **dynamically conjugate** to `L`: orbit structure, periods, recurrence relations and invariants are transported by `T`. In the original `(u,v)` state coordinates, `F` is nonlinear, so it does **not** have an ordinary linear-matrix spectrum. The Fourier **mode matrices of `L`** (or appropriately defined conjugate Koopman operators) provide the spectral information. This exact class is *nonlinearity caused by a change of coordinates*, not evidence for a general nonlinear jump.

For the largest tested ring `D=59,049`, the field prime is `p=1,299,079`. Every mode is

```
L(w) = [[1+w², w], [w, 1]]
tr L(w) = 2+w²       det L(w) = 1
Delta(w) = (tr L)^2 - 4 = w²(w²+4).
```

Since `p=3 mod 4`, `-4` is not a square modulo `p`, while `w != 0`; therefore **no discriminant is zero**. Every mode has distinct eigenvalues over at most the quadratic extension. Exhaustively enumerating all `59,049` roots gives exactly **29,605 split modes** (both eigenvalues in `GF(p)`) and **29,444 nonsplit modes** (both in `GF(p²)`), with **zero defective modes**.

For a split determinant-one mode, each eigenvalue has order dividing `p-1`. For a nonsplit mode, Frobenius swaps the two roots, and their product is one; thus each eigenvalue obeys `lambda^(p+1)=1`. As all modes are diagonalizable, the full transformed operator, and hence the nonlinear map `F`, satisfy

```
F^M = identity,    M = lcm(p-1,p+1) = (p²-1)/2
M = 843,803,124,120.
3^80 mod M = 444,653,673,201.
```

Therefore the `3^80` snapshot is **exactly** the snapshot at step `444,653,673,201` modulo the certified universal exponent. **This does not prove that `M` is the minimal common period**, nor that a long modular jump represents a new physical time scale. The new independent [all-mode period audit](../python/matrix_tricks/nonlinear_period_monotonicity_audit.py) directly checks `L(w)^M=I` at all `59,049` modes, counts the discriminant classes, and checks the long-jump residue. The uniform mode `w=1` is the Arnold cat map `[[2,1],[1,1]] mod p`.

**Method 1 scope line:** An exact jump exists here **because** an explicit invertible conjugacy to a linear matrix was constructed. Its dynamical invariants and periodicity come from that matrix. The result says nothing about a nonlinear model for which no such conjugacy or another valid reduction theorem is known.

**Why this matters and what it does not mean:** real nonlinear feedback occurs at each microscopic step, yet it can be skipped **because a specific invertible conjugacy was proven**. This does not imply that a generic nonlinear force is integrable or that our full SU(3) Yang–Mills equations possess such a conjugacy.

## 2. A different, irreducibly dissipative nonlinear model with a uniform long-time certificate

Use signed integers, *not* modular arithmetic. Define `trunc0(z/d)` as integer division **toward zero**, and for every spatial site

```
x'_i = trunc0( [3*x_i + x_(i-1) + trunc0(x_i*x_(i-1)/Q)] / 8 )
Q = 1,000,000.
```

Here the local source depends on both the current state and its neighbor and is applied **every** time step. The map is irreversible because of integer truncation; it is not the nonlinear conjugate of the invertible Fourier evolution above.

Let `B_t=max_i |x_i(t)|` and assume `B_0<=Q`. Since

```
|trunc0(z/d)| = floor(|z|/d),
|trunc0(x_i*x_(i-1)/Q)| <= floor(B_t^2/Q),
```

the entire dynamical system obeys the **uniform integer comparison**

```
B_(t+1) <= floor( (4*B_t + floor(B_t^2/Q)) / 8 ).
```

For every integer `1<=B<=Q`, the right side is at most `floor(5B/8)<B`. The exact comparison trajectory from `B_0=1,000,000` is

```
1,000,000, 625,000, 361,328, 196,983, 103,341, 53,005, 26,853,
13,516, 6,780, 3,395, 1,698, 849, 424, 212, 106, 53, 26,
13, 6, 3, 1, 0.
```

The missing comparison step is explicit: for nonnegative integers `B`, `g(B)=floor((4B+floor(B²/Q))/8)` is **nondecreasing**, because squaring on nonnegative integers, multiplication by positive coefficients, addition and floor division by positive integers all preserve order. Thus if `B_0 <= Q`, induction yields `B_t <= g^t(B_0) <= g^t(Q)`. This justifies using the single top-radius trajectory to enclose **every smaller initial radius**. The [independent monotonicity audit](../python/matrix_tricks/nonlinear_period_monotonicity_audit.py) additionally checks this ordering and strict shrinkage at all one million positive integer radii.

Therefore **every** initial condition in the entire sup-norm ball `|x_i|<=1,000,000`, for **any finite ring length D**, satisfies `x(21)=0`, and because zero is fixed, `x(n)=0` for every `n>=21`, including `3^80`. This is an exact universal-in-that-ball *long-time certificate*, derived from a proved inequality; it does **not** require enumerating states.

A stronger self coefficient `9` instead of `3` breaks the comparison at `x_i=Q`, explicitly showing that the bound is not indiscriminately reusable.

**Method 2 scope line:** This is a **certificate of extinction**. The long-time state is known because it becomes exactly zero by step 21. It cannot certify arbitrary persistent nonlinear waves, oscillations, chaotic attractors, or non-dissipative gauge dynamics.

This dissipative example reaches a fixed point; it does not claim to skip complicated persistent nonlinear transients or oscillations. It demonstrates a different valid general mathematical route: a comparison theorem over the complete state set.

## Experiments, reproducibility and controls

- **Four dynamical ring dimensions:** 27, 243, 6,561 and 59,049 sites. The exact nonlinear system has **two fields** at each site, so the largest complete dynamical state has **118,098 coupled integer components**.
- **Reversible exact lane:** tested ordinary direct stepping for 1, 3 and 27 steps at every dimension (and 243 steps on the small rings); exact modal jumps to `3^6,3^20,3^40` at every size, plus `3^80` at dimensions 27, 243 and 59,049; verified forward/backward freeze/play leaves at `3^40`.
- **Independent code:** expands the nonlinear physical equations directly without calling the conjugacy implementation; checks all 27 direct steps at all four sizes. Its second check computes **naive O(D²) DFTs** rather than reusing the radix-three NTT for `D=27,243` at both `3^40` and `3^80`.
- **Long-time bound:** separately computes all **1,000,000** positive scalar comparison radii, proving strict contraction at every one, and checks the independent spatial update for 21 steps on each tested dimension.
- **Negative controls:** generic linear superposition fails in physical variables; applying the uncorrected linear step to physical variables differs; a modified feedback term differs; and the altered self coefficient `9` cannot inherit the contraction bound.
- **Sparse readout:** 243 positions are sampled from a **12,000,000,000-address virtual detector**, but there are only 118,098 independently evolved components at the largest actual test. These are not physical qubits, gluons, glueballs or protons.

Local tested canonical receipt SHA-256:

`a308a754634f7738e1f68e75b9de394e6a9e98201859baba517d067c32b3d7f6`.

At `D=59,049` (118,098 components), the `3^40` nonlinear jump complete-state SHA is

`0233978db42073ae87a8c22680bf3664e5ab2545200ac5031908d000cd7c864c`

and the `3^80` complete-state SHA is

`14c12e884522631511e4ab86f38cfd22e799dad53d068bfc6d2df83688b7f012`.

Run the hosted workflow (two Python runtimes) or reproduce locally:

```bash
python3 -B python/matrix_tricks/nonlinear_exact_and_bound.py \
  --dims 27,243,6561,59049 --output build/nonlinear/NONLINEAR-JUMP.hbp
python3 -B python/matrix_tricks/nonlinear_independent_audit.py \
  build/nonlinear/NONLINEAR-JUMP.hbp
python3 python/audit/integer_only_audit.py \
  python/matrix_tricks/nonlinear_exact_and_bound.py \
  python/matrix_tricks/nonlinear_independent_audit.py
sha256sum build/nonlinear/NONLINEAR-JUMP.hbp
```

**Yang–Mills scope line:** The protected `r=9` task is a **linear Hamiltonian spectral floor certification**, not a nonlinear classical time-evolution problem. Neither nonlinear jump method in this note accelerates that existing spectral certificate, the outstanding full-symmetry `r=12` lower-sector bounds, or the continuum SU(N) mass-gap problem.

**Operations, reported November 12, 2026 at 22:36 UTC:** The operator reported the protected `r=9` verification had **already launched**, with nine steps completed at `rc=0` and no exit/result yet. Previous wording here that the *launch remained on HOLD* described the earlier pre-launch gate and is superseded. This note does not launch, restart, cancel or certify that run. Its final verdict remains pending based on the reported checkpoint.

**Proof boundary:** The exact lane is integrable by construction; the bounded lane is strictly dissipative. Neither supplies a general nonlinear Yang–Mills, gravitational, or persistent many-body dynamics shortcut.
