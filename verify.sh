#!/usr/bin/env bash
# End-to-end reproduction, integer-only. Every check fails closed.
set -euo pipefail
cd "$(dirname "$0")"
cargo build --release -q
mkdir -p build
B=target/release
# 1. generalised builder reproduces the pinned 2x2 operator byte for byte
$B/torusl 2 2 build/torus2x2.matrix
test "$(sha256sum build/torus2x2.matrix | cut -d' ' -f1)" = 137285e0a1b55bdf0c8c56bec6bd69e0ab9138e6e7c05350d677f186e1981d2c
cmp build/torus2x2.matrix data/matrices/torus2x2.matrix
# 2. 2x3 operator: 47,449 states (Burnside), all integer checks, pinned sha
$B/torusl 2 3 build/torus2x3.matrix | tee build/torus2x3.line
grep -q "dim=47449|burnside=47449" build/torus2x3.line
grep -q "checks=PASS" build/torus2x3.line
test "$(sha256sum build/torus2x3.matrix | cut -d' ' -f1)" = 2bd842b0a2358c93dc6b4195bae4bddcbb322829a24e28d4857f9a3a2f3f66de
# 3. one interval-Sturm certificate: exactly one level below 13.240 at (JE,JB)=(1,1)
$B/sturmrs data/matrices/torus2x2.matrix 1 1 13240 30 | tee build/sturm.line
grep -q "|negative=1|" build/sturm.line
# 4. four sectors sum to the full count at an anchor (1,3,19.300)
# (the verifier exits 2 whenever its count is not exactly 1, by design; read the count, not the code)
neg() { { "$B/sturmrs" "$@" || true; } | tr '|' '\n' | grep '^negative=' | cut -d= -f2; }
n=0; for s in pp pm mp mm; do c=$(neg data/matrices/sector_$s.matrix 1 3 19300 30); n=$((n+c)); done
f=$(neg data/matrices/torus2x2.matrix 1 3 19300 30)
test "$n" = "$f" && test "$f" = 1
# 5. exact magnetic-end laws (second order): 2x2 -> sqrt14 - sqrt6, 2x3 -> 3
python3 python/volume/second_order.py build/torus2x2.matrix | grep -q "gap_coeff_in=\[1.292167643,1.292167643\]"
python3 python/volume/flat_law.py 2x3 | grep -q "gap_coeff_in=\[2.999999999,3.000000000\]"
# 5b. 3x3: second order verified scalar and degenerate (coefficient 0, NOT a vanishing gap),
#     third-order splitting (9/8)(sqrt14-sqrt6); 4x4: fourth-order (5/4)(sqrt14-sqrt6).
#     Exact perturbative laws, pending independent mathematical review.
python3 python/volume/flat_law.py 3x3 | grep -q "^FLAT_LAW|L=3x3|second_order_degenerate=1|second_order_gap_coeff=0|next_order_required=3|json=0$"
python3 python/volume/flat_law3.py 3x3 | grep -q "gap_coeff3_in=\[1.453688599,1.453688599\]"
python3 python/volume/flat_law4.py 4x4 | grep -q "gap_coeff4_in=\[1.615209554,1.615209554\]"
python3 - <<'PY'
import math
K = 10**30
d = math.isqrt(14*K*K) - math.isqrt(6*K*K)          # (sqrt14 - sqrt6) * 10^30 (floor-level)
assert (9*d//8) // 10**21 == 1453688599, "3x3 coefficient"
assert (5*d//4) // 10**21 == 1615209554, "4x4 coefficient"
print("CLOSED_FORM|3x3=(9/8)(sqrt14-sqrt6)|4x4=(5/4)(sqrt14-sqrt6)|status=PASS|review=pending|json=0")
PY
# 5c. archives: published 2x3 receipts unchanged; regenerated law archive complete
test "$(sha256sum data/volume/volume_t23.hbp | cut -d' ' -f1)" = 53b29aeb421eed669aae9c538d2ffa6696988ddb1f5b014d29f48ed0729403df
tail -1 data/volume/flat_law.hbp | grep -q "^FLAT_LAW|status=COMPLETE|"
# 5d. matrix-free 3x3 operator (10,097,891-state lattice, rows generated on demand):
#     support after 4 applications from the trivial vacuum = 35,847 (independently
#     reproduced count), exact W-self-adjointness at every step, pinned norm at (1,1)
$B/onfly 3 3 1 1 4 | tee build/onfly.out
tail -1 build/onfly.out | grep -q "|support=35847|norm_w=122922827987436|"
test "$(grep -c "sa_check=PASS" build/onfly.out)" = 3
# 5e. Exact integer power-jump: the 64-step state MUST match sequential
#     evaluation byte for byte, and the 4096-step state MUST match two
#     2048-step jumps. This is the effective 8-vacuum matrix, not full 3x3.
python3 python/matrix_tricks/exact_jump.py --self-test | tee build/exact_jump.hbp
grep -q '|steps=64|squarings=6|exact_equal_sequential=1|' build/exact_jump.hbp
grep -q '|steps=4096|squarings=12|state_bits=31168|' build/exact_jump.hbp
grep -q '|rounding=NONE|exact_composition=PASS|' build/exact_jump.hbp
grep -q '^EXACT_JUMP|result=PASS|' build/exact_jump.hbp
# 5f. Actual 3x3 vacuum-ball principal compression via the repository's
#     comb-tree lattice source: do NOT confuse Ritz-gap intervals with a
#     certified full 10,097,891-state gap.
python3 python/matrix_tricks/vacuum_ball_ritz.py 1 27 --order 1 | tee build/vacuum_ball_r1.hbp
python3 python/matrix_tricks/vacuum_ball_ritz.py 2 27 --order 2 | tee build/vacuum_ball_r2.hbp
grep -q '|states=206|layers=8,198|nnz=818|' build/vacuum_ball_r1.hbp
grep -q '|ritz_dim=16|E0_in=\[53496532,53496533\]/1000000|E1_in=\[53498537,53498538\]/1000000|gap_in=\[2004,2006\]/1000000|' build/vacuum_ball_r1.hbp
grep -q '|states=4454|layers=8,198,4248|nnz=33218|' build/vacuum_ball_r2.hbp
grep -q '|ritz_dim=24|E0_in=\[53493574,53493575\]/1000000|E1_in=\[53495589,53495590\]/1000000|gap_in=\[2014,2016\]/1000000|' build/vacuum_ball_r2.hbp
# 5g. Rust vacuum-ball compression (crates/torusl/src/bin/ball.rs): principal
#     compression of H onto all orbits within K moves of the 8 vacua (upper bounds by
#     interlacing), solved by the existing evolvers.
#     (i) validation on 2x3: K=4 (34,890 states) reproduces the full 47,449-state
#         Chebyshev measurement at r=9 (E0 34.470733 vs 34.470732, gap 0.341110).
mkdir -p build/ball
$B/ball 2 3 4 build/ball/b23_k4.matrix | grep -q "|states=34890|nnz=1023198|.*|w_self_adjoint=PASS|"
$B/cheb build/ball/b23_k4.matrix 1 9 4000 341110 | grep -q "|E0=34470733/1000000|dE0=0|E1=34811843/1000000|dE1=0|gap=341110/1000000|"
#     (ii) 3x3 K=1 (206 states): converged ball eigenvalues by dense matrix jump at r=27;
#         they must sit at or below the 16-dim Krylov Ritz values of section 5f.
$B/ball 3 3 1 build/ball/b33_k1.matrix | grep -q "|states=206|nnz=818|.*|w_self_adjoint=PASS|"
$B/jump build/ball/b33_k1.matrix 1 27 27 6 | tee build/ball/b33_k1_r27.line | grep -q "|E0=53496482/1000000|dE0=0|E1=53498465/1000000|dE1=0|gap=1983/1000000|"
test 53496482 -le 53496532 && test 53498465 -le 53498537
# 5h. Block cluster eigensolver (crates/evolvers/src/bin/cluster.rs) on 3x3 vacuum balls:
#     all 8 vacuum levels resolved at once, rigorous residual intervals (exact big
#     integers), ordering certified by interval-Sturm counts between the intervals.
#     Scope: the compressed operator (ball), not the full 3x3 lattice.
$B/ball 3 3 2 build/ball/b33_k2.matrix | grep -q "|states=4454|nnz=33218|.*|w_self_adjoint=PASS|"
python3 python/volume/cluster_certify.py build/ball/b33_k1.matrix 27 | tee build/ball/cert_k1_r27.hbp | grep -q "|E0_in=\[53496481,53496484\]/1000000|E1_in=\[53498464,53498467\]/1000000|gap_in=\[1980,1986\]/1000000|.*|count_t1=1|.*|count_t2=2|ordering=CERTIFIED|"
python3 python/volume/cluster_certify.py build/ball/b33_k2.matrix 27 | tee build/ball/cert_k2_r27.hbp | grep -q "|E0_in=\[53493381,53493384\]/1000000|E1_in=\[53495393,53495396\]/1000000|gap_in=\[2009,2015\]/1000000|.*|count_t1=1|.*|count_t2=2|ordering=CERTIFIED|"
# 6. integer-only audit of every Python tool
python3 python/audit/integer_only_audit.py python/ladder/*.py python/testground/*.py python/sweep/*.py python/volume/*.py python/matrix_tricks/*.py
echo "VERIFY|status=PASS|json=0"
