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
# 6. integer-only audit of every Python tool
python3 python/audit/integer_only_audit.py python/ladder/*.py python/testground/*.py python/sweep/*.py python/volume/*.py python/matrix_tricks/*.py
echo "VERIFY|status=PASS|json=0"
