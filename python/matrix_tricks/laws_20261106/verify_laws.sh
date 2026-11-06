#!/usr/bin/env bash
# The laws of November 6, 2026: the HBI -> HBP -> SHA -> SH -> HASH chain, then an independent replay of every check.
# Fails closed on any changed byte. HBP = the *.out.hbp records; HBI = LAWS.hbi (one row per file: bytes, sha256, sh);
# SHA = SHA256SUMS; SH = the first 16 hex of each sha256; HASH = LAWS.seal = sha256(SHA256SUMS bytes || LAWS.hbi bytes).
set -euo pipefail
cd "$(dirname "$0")"
FILES="hash_answer.py jump3_jump.py pascal_pyramid.py pow27_mod3.py hash_answer.out.hbp jump3_jump.out.hbp pascal_pyramid.out.hbp pow27_mod3.out.hbp verify_laws.sh"
build_hbi() {
  for f in $FILES; do
    s=$(sha256sum "$f" | cut -d' ' -f1)
    printf 'HBI|name=%s|bytes=%s|sha256=%s|sh=%s|json=0\n' "$f" "$(wc -c < "$f" | tr -d ' ')" "$s" "$(printf '%s' "$s" | cut -c1-16)"
  done
}
build_seal() {
  printf 'HASH|seal=%s|over=SHA256SUMS+LAWS.hbi|json=0\n' "$(cat SHA256SUMS LAWS.hbi | sha256sum | cut -d' ' -f1)"
}
if [ "${1:-}" = "--build" ]; then
  sha256sum $FILES | sed 's/ \*/  /' > SHA256SUMS
  build_hbi > LAWS.hbi
  build_seal > LAWS.seal
  echo "LAWS_CHAIN|built=1|json=0"
  exit 0
fi
sha256sum --strict -c --quiet SHA256SUMS
cmp <(build_hbi) LAWS.hbi
cmp <(build_seal) LAWS.seal
T=$(mktemp -d)
mask() { sed -E 's/(cube_ms|ntt_ms|step_ms|hashes_per_second_measured|full_search_years_at_that_rate)=[0-9]+/\1=T/g'; }
for s in pow27_mod3 pascal_pyramid; do python3 -I -B -X utf8 "$s.py" | tr -d '\r' > "$T/$s.out"; cmp "$T/$s.out" "$s.out.hbp"; done
for s in jump3_jump hash_answer; do python3 -I -B -X utf8 "$s.py" | tr -d '\r' | mask > "$T/$s.out"; cmp "$T/$s.out" <(mask < "$s.out.hbp"); done
rm -r "$T"
echo "LAWS_VERIFY|chain=PASS|replay=PASS|seal=$(cut -d'|' -f2 LAWS.seal)|json=0"
