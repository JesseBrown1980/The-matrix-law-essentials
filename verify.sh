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
# 6. integer-only audit of every Python tool
python3 python/audit/integer_only_audit.py python/ladder/*.py python/testground/*.py python/sweep/*.py python/volume/*.py
echo "VERIFY|status=PASS|json=0"
