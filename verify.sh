#!/usr/bin/env bash
# End-to-end reproduction, integer-only. Every check fails closed.
set -euo pipefail
cd "$(dirname "$0")"
# Match without closing the pipe early: grep -q exits at the first match, and a
# producer that is still writing then dies of SIGPIPE under pipefail. grep -c
# reads all of its input and still fails when nothing matches.
has() { grep -c "$@" > /dev/null; }
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
python3 python/volume/second_order.py build/torus2x2.matrix | has "gap_coeff_in=\[1.292167643,1.292167643\]"
python3 python/volume/flat_law.py 2x3 | has "gap_coeff_in=\[2.999999999,3.000000000\]"
# 5b. 3x3: second order verified scalar and degenerate (coefficient 0, NOT a vanishing gap),
#     third-order splitting (9/8)(sqrt14-sqrt6); 4x4: fourth-order (5/4)(sqrt14-sqrt6).
#     Exact perturbative laws, pending independent mathematical review.
python3 python/volume/flat_law.py 3x3 | has "^FLAT_LAW|L=3x3|second_order_degenerate=1|second_order_gap_coeff=0|next_order_required=3|json=0$"
python3 python/volume/flat_law3.py 3x3 | has "gap_coeff3_in=\[1.453688599,1.453688599\]"
python3 python/volume/flat_law4.py 4x4 | has "gap_coeff4_in=\[1.615209554,1.615209554\]"
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
tail -1 data/volume/flat_law.hbp | has "^FLAT_LAW|status=COMPLETE|"
# 5d. matrix-free 3x3 operator (10,097,891-state lattice, rows generated on demand):
#     support after 4 applications from the trivial vacuum = 35,847 (independently
#     reproduced count), exact W-self-adjointness at every step, pinned norm at (1,1)
$B/onfly 3 3 1 1 4 | tee build/onfly.out
tail -1 build/onfly.out | has "|support=35847|norm_w=122922827987436|"
test "$(grep -c "sa_check=PASS" build/onfly.out)" = 3
# 5e. Exact integer power-jump: the 64-step state MUST match sequential
#     evaluation byte for byte, and the 4096-step state MUST match two
#     2048-step jumps. This is the effective 8-vacuum matrix, not full 3x3.
python3 python/matrix_tricks/exact_jump.py --self-test | tee build/exact_jump.hbp
grep -q '|steps=64|squarings=6|exact_equal_sequential=1|' build/exact_jump.hbp
grep -q '|steps=4096|squarings=12|state_bits=31168|' build/exact_jump.hbp
grep -q '|rounding=NONE|exact_composition=PASS|' build/exact_jump.hbp
grep -q '^EXACT_JUMP|result=PASS|' build/exact_jump.hbp
# 5e-2. Cosmological time-jump (2^60 steps) and triadic ladder in-between sampling:
#       evaluates skips of 1, 2, 3, 5, 6, 8, 9, 27, 81, 243 without simulating
#       intervening steps, and verifies skip-then-measure with hash change.
python3 python/matrix_tricks/time_jump_ladder.py --self-test | tee build/time_jump_ladder.hbp
grep -q '^COSMOLOGICAL_JUMP|squarings=60|steps=2\^60|.*|ground_convergence=PASS|' build/time_jump_ladder.hbp
grep -q '^TRIADIC_LADDER|rungs=1,2,3,5,6,8,9,27,81,243|all_exact_match=PASS|' build/time_jump_ladder.hbp
grep -q '^INBETWEEN_SAMPLING|.*|exact_sequential_match=PASS|' build/time_jump_ladder.hbp
grep -q '^HASH_CHANGE_JUMP|.*|inbetween_measured_match=PASS|' build/time_jump_ladder.hbp
grep -q '^TIME_JUMP_VERIFY|status=PASS|' build/time_jump_ladder.hbp
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
$B/ball 2 3 4 build/ball/b23_k4.matrix | has "|states=34890|nnz=1023198|.*|w_self_adjoint=PASS|"
$B/cheb build/ball/b23_k4.matrix 1 9 4000 341110 | has "|E0=34470733/1000000|dE0=0|E1=34811843/1000000|dE1=0|gap=341110/1000000|"
#     (ii) 3x3 K=1 (206 states): converged ball eigenvalues by dense matrix jump at r=27;
#         they must sit at or below the 16-dim Krylov Ritz values of section 5f.
$B/ball 3 3 1 build/ball/b33_k1.matrix | has "|states=206|nnz=818|.*|w_self_adjoint=PASS|"
$B/jump build/ball/b33_k1.matrix 1 27 27 6 | tee build/ball/b33_k1_r27.line | has "|E0=53496482/1000000|dE0=0|E1=53498465/1000000|dE1=0|gap=1983/1000000|"
test 53496482 -le 53496532 && test 53498465 -le 53498537
# 5h. Block cluster eigensolver (crates/evolvers/src/bin/cluster.rs) on 3x3 vacuum balls:
#     all 8 vacuum levels resolved at once, rigorous residual intervals (exact big
#     integers), ordering certified by interval-Sturm counts between the intervals.
#     Scope: the compressed operator (ball), not the full 3x3 lattice.
$B/ball 3 3 2 build/ball/b33_k2.matrix | has "|states=4454|nnz=33218|.*|w_self_adjoint=PASS|"
python3 python/volume/cluster_certify.py build/ball/b33_k1.matrix 27 | tee build/ball/cert_k1_r27.hbp | has "|E0_in=\[53496481,53496484\]/1000000|E1_in=\[53498464,53498467\]/1000000|gap_in=\[1980,1986\]/1000000|.*|count_t1=1|.*|count_t2=2|ordering=CERTIFIED|"
python3 python/volume/cluster_certify.py build/ball/b33_k2.matrix 27 | tee build/ball/cert_k2_r27.hbp | has "|E0_in=\[53493381,53493384\]/1000000|E1_in=\[53495393,53495396\]/1000000|gap_in=\[2009,2015\]/1000000|.*|count_t1=1|.*|count_t2=2|ordering=CERTIFIED|"
# 5i. Long jumps of the FULL evolution operator T = c - H with play-forward
#     (crates/evolvers/src/bin/jumpto.rs). Exact mode keeps every component modulo
#     a prime (no rounding): jump(1024) == 1024 sequential steps, jump(5e9)
#     applied twice == jump(1e10), 27 frames played from 1e10 == jump(1e10+27).
#     (The 1,393-state 2x2 run, ~4 min, is stored in data/jump/jumpto_runs.hbp;
#     CI uses the 206-state 3x3 vacuum ball.) Fixed mode measures E at every
#     stage 1..1e10; at 1e10 an exact residual radius puts an eigenvalue in
#     [53.832773, 53.832776], the certified E0 interval of that ball.
#     Modular mode = exact algebraic replay only (no real amplitudes or energies).
#     One step = v <- (c - H) v, forward-Euler imaginary time, not exp(-tau H).
$B/jumpto build/ball/b33_k1.matrix 1 81 exact 998244353 10000000000 | has "|seq1024_equal=1|composition_5000000000+5000000000_equal=1|replay27_equal=1|.*|rounding=NONE|"
$B/jumpto build/ball/b33_k1.matrix 1 81 fixed 10 27 | has "1e7:53832774,1e8:53832774,1e9:53832774,1e10:53832774|play_from_1e10=+1:53832774,+2:53832774,+3:53832774,+27:53832774|E_at_1e10=53832774/1000000|eps_up=1/1000000|eigenvalue_in=\[53832773,53832776\]/1000000|"
#     independent sparse-step reference (separate Python code reading the matrix file):
#     its state fingerprints after 1024 steps must equal the jump's
ref=$(python3 python/matrix_tricks/sparse_mod_reference.py build/ball/b33_k1.matrix 1 81 998244353 1024 | tr '|' '\n' | grep -E '^(state_sum_mod|state_weighted_mod|v0_component0)=' | tr '\n' '|')
jmp=$($B/jumpto build/ball/b33_k1.matrix 1 81 exact 998244353 1024 | tr '|' '\n' | grep -E '^(state_sum_mod|state_weighted_mod|v0_component0)=' | tr '\n' '|')
test -n "$ref" && test "$ref" = "$jmp"
# 5j. Symmetry-sector compression (crates/rainbon/src/bin/symsec.rs): the
#     three-preon cloud model reduced by its 288 symmetries (48 cube operations x 6
#     preon permutations; +-1 characters for the P-wave sector) must reproduce the
#     full 46,656-state lattice at L=6 exactly: sector A (309 states) gives the
#     full E0 and <r^2> to every digit, and sector Pz (1,680 states) gives an E1
#     inside the block solver's certified interval (JE units x 9).
full=$($B/cloud 6 harmonic 9 1 600 10000000)
secA=$($B/symsec 6 harmonic A 9 1 800 10000000)
secP=$($B/symsec 6 harmonic Pz 9 1 800 10000000)
fld() { echo "$1" | tr '|' '\n' | grep "^$2=" | cut -d= -f2 | cut -d/ -f1; }
test "$(fld "$full" E0)" = "$(fld "$secA" E_lowest)" && test "$(fld "$full" E0)" = 23960538
test "$(fld "$full" r2_ground)" = "$(fld "$secA" r2)" && test "$(fld "$secA" r2)" = 2597261
echo "$secA" | has "|sector_states=309|w_self_adjoint=PASS|"
echo "$secP" | has "|sector_states=1680|w_self_adjoint=PASS|"
$B/cloudb 6 harmonic build/rb_h6.matrix > /dev/null
iv=$($B/cluster build/rb_h6.matrix 9 1 8 4 | grep "^CLUSTER_LEVEL|k=1|" | tr '|' '\n' | grep '^interval=' | cut -d= -f2 | tr -d '[]' | cut -d/ -f1)
lo=$(( ${iv%,*} * 9 )); hi=$(( ${iv#*,} * 9 + 9 )); e1=$(fld "$secP" E_lowest)
test "$e1" -ge "$lo" && test "$e1" -le "$hi"
# 5k. Zero-momentum C2v sectors of the FULL operator (crates/torusl/src/bin/ymsec.rs):
#     group relations (incl. semidirect), commutation with H, exact W-self-adjointness,
#     predicted flat-vacuum decomposition 7/0/0/1 per sector (from the reflection
#     action on commuting holonomy pairs), and spectral reconstruction:
#     2x2 union == PR #14 zero-momentum sector, 12 levels exactly; 2x3 union contains
#     the full operator's zero-momentum levels and its ground state exactly.
for L in "2 2" "2 3"; do
  i=0; for s in "1 1" "1 -1" "-1 1" "-1 -1"; do
    set -- $L $s; want=$(echo "7 0 0 1" | cut -d' ' -f$((i+1))); i=$((i+1))
    $B/ymsec $1 $2 $3 $4 build/ymsec_check.matrix | has "|flat_states=$want|.*|group_relations=PASS|commutes_with_H=PASS|.*|w_self_adjoint=PASS|"
  done
done
python3 python/volume/ymsec_validate.py 2x2 27 | has "|flat_total=8|levels_compared=12|max_abs_diff=0|.*|result=PASS|"
python3 python/volume/ymsec_validate.py 2x3 9 | has "|flat_total=8|ground_equal=1|union_levels_in_window=10|matched_with_multiplicity=10|.*|result=PASS|"
#     exhaustive commutation on every gauge orbit (2x2: 1,393; 2x3: 47,449)
YMSEC_CHECK_SAMPLES=0 $B/ymsec 2 2 1 1 build/ymsec_check.matrix | has "|gauge_orbits=1393|.*|group_relations=PASS|commutes_with_H=PASS|structural_samples=1393|"
YMSEC_CHECK_SAMPLES=0 $B/ymsec 2 3 1 1 build/ymsec_check.matrix | has "|gauge_orbits=47449|.*|group_relations=PASS|commutes_with_H=PASS|structural_samples=47449|"
#     sector completeness by exact spectral moments (no eigensolver): Tr(H^n), n=1..4
mkdir -p build/ymsec
for s in "1 1" "1 -1" "-1 1" "-1 -1"; do set -- $s; $B/ymsec 2 2 $1 $2 build/ymsec/m_$1_$2.matrix > /dev/null; done
for jb in 1 27; do
  python3 python/volume/trace_moments.py data/matrices/sector_pp.matrix 1 $jb 4 build/ymsec/m_1_1.matrix build/ymsec/m_1_-1.matrix build/ymsec/m_-1_1.matrix build/ymsec/m_-1_-1.matrix | has "^TRACE_MOMENTS|.*|result=PASS|"
done
#     sturmfast (O(log n) minimum-degree pivot search, same interval arithmetic)
#     must reproduce sturmrs counts exactly: four 2x2 sectors x three thresholds
negf() { { "$B/sturmfast" "$@" || true; } | tr '|' '\n' | grep '^negative=' | cut -d= -f2; }
for s in "1 1" "1 -1" "-1 1" "-1 -1"; do set -- $s
  for t in 13240 20000 60000; do
    a=$(neg build/ymsec/m_$1_$2.matrix 1 1 $t 30); b=$(negf build/ymsec/m_$1_$2.matrix 1 1 $t 30)
    [ -n "$a" ] && [ "$a" = "$b" ] || { echo "sturmfast mismatch $1,$2 t=$t: $a vs $b"; exit 1; }
  done
done
#     shellschur (magnetic-shell Schur + Haynsworth inertia) must reproduce sturmrs
#     counts exactly: full 2x2 at JB=27, and four 2x2 sectors at JB=9 and 27
negs() { { "$B/shellschur" "$@" || true; } | tr '|' '\n' | grep '^negative=' | cut -d= -f2; }
[ "$(negs data/matrices/torus2x2.matrix 1 27 64000)" = "$(neg data/matrices/torus2x2.matrix 1 27 64000 30)" ] || { echo "shellschur mismatch 2x2 JB=27"; exit 1; }
for s in "1 1" "1 -1" "-1 1" "-1 -1"; do set -- $s
  for jb in 9 27; do
    a=$(neg build/ymsec/m_$1_$2.matrix 1 $jb 64000 30); b=$(negs build/ymsec/m_$1_$2.matrix 1 $jb 64000 -1 15 10)
    [ -n "$a" ] && [ "$a" = "$b" ] || { echo "shellschur mismatch $1,$2 JB=$jb: $a vs $b"; exit 1; }
  done
done
# 6. integer-only audit of every Python tool
python3 python/audit/integer_only_audit.py python/ladder/*.py python/testground/*.py python/sweep/*.py python/volume/*.py python/matrix_tricks/*.py
echo "VERIFY|status=PASS|json=0"
