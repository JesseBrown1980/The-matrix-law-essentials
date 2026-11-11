# CO molecular toy tests: integer-only replay (November 11, 2026)

## Why this exists

The three CO tests of November 6 (three slices, the two-molecule case, and the 1-to-729 ensemble) were written with Python `Decimal` at 65 digits, `Fraction` and ordinary division. `verify.sh` audits every script in `python/matrix_tricks/` for integer-only arithmetic, so the tests could not join `main` as they stood: the audit reported 22, 10 and 15 violations. [The Vienna repair of November 8](VIENNA-INTEGER-REPLAY-20261108.md) met the same problem and recorded that the two open molecular branches had to be brought up to date. This is the same repair, applied to the CO tests. It changes how the numbers are computed and checked. It does not change the model, the predictions or the interpretation.

## Historical evidence retained

Each file is byte-identical to the original; the replay checks its Git blob before it computes anything.

| Retained file | Git blob |
|---|---|
| `data/molecular_co/evidence/ORIGINAL-DECIMAL-SLICES-SOURCE-20261106.txt` | `32c0972b0a712c0ad9f6cd1d3cf39132ff9d1c31` |
| `data/molecular_co/evidence/ORIGINAL-DECIMAL-PAIR-SOURCE-20261106.txt` | `cb349b486a718d30f8e84430e76a6a05816b48b0` |
| `data/molecular_co/evidence/ORIGINAL-DECIMAL-ENSEMBLE-SOURCE-20261106.txt` | `ae8a0a0564d7b7c03a8d43c3dee3778dd286b2f4` |
| `data/molecular_co/MEASURED-CO-MOLECULAR-SLICES-20261106.hbp` | `ede652ea6ce55b5587fe98345afae0de39325e94` |
| `data/molecular_co/MEASURED-CO-PAIR-20261106.hbp` | `25db2a47bbac5a6f88856c91331763bf559b98a0` |
| `data/molecular_co/MEASURED-CO-ENSEMBLES-20261106.hbp` | `68731f98004fad1e33ec64dbe154a781446353e1` |

- The three prediction files are untouched.
- The archived sources are text. No executable imports or runs them.
- Before archiving, each Decimal script was run once more on the merged tree and reproduced its receipt byte for byte.

## The replacement executable

`python/matrix_tricks/co_molecular_integer_replay_20261106.py --self-test` recomputes all 72 receipt rows in 100-digit scaled **integer** arithmetic:

- **Morse bond:** velocity Verlet at the three time slices, the independent RK4 reference at one ninth of the fine step, energy drift, and the three spatial bin readouts.
- **Harmonic jumps:** matrix powers by squaring and by 15 cubings to 3^15 = 14,348,907 ticks; the exact third and ninth powers as unscaled integer matrices.
- **Two molecules:** in-phase, opposite-phase and coupled normal modes, with the directly stepped pair.
- **Ensembles:** 1 to 729 phase-balanced coordinates, and the coupled array with every one of its 3 to 729 coordinates stepped.
- **Red controls:** all 12 of the originals.

## Measured agreement

Every field of every row is compared with the retained receipt. A field passes in exactly one of three ways:

| Class | Count | Rule |
|---|---|---|
| Identical text | 436 | The replayed field equals the receipt, character for character. This includes every 12-digit value and all nine SHA-256 digests of 40-digit states. |
| Rounding residue | 30 | Both values are below 10⁻⁴⁰. These are the differences that are zero in exact arithmetic (jump against sequential, triadic against binary, balanced means); the receipt holds Decimal's residue near 10⁻⁶⁰. |
| Engine label | 15 | `precision` and `exact_discrete_jump`, which name the arithmetic: `DECIMAL_65…` in the receipt, `INTEGER_FIXED_POINT_100…` in the replay. |
| Anything else | 0 | Refused. |

`--rows` prints the replayed rows themselves.

## Scope

- **Not byte-equivalent arithmetic.** 100-digit integers and 65-digit Decimal round differently below about 10⁻⁶⁰. The receipts are not overwritten.
- **Still a toy.** A classical fitted CO stretch and a coupling of ω²/9 chosen by hand. No measured CO–CO interaction, no quantum chemistry, no observed molecular motion.
- **Not part of the certified finite-S₃ core.** Nothing here bears on the Hamiltonian certificates, on spacetime discreteness or on the four-dimensional Yang–Mills problem.

## Acceptance

1. The strict integer audit reports zero violations on the replay.
2. All 481 fields fall in the three classes above, with zero mismatches.
3. All six retained Git blobs match.
4. `./verify.sh` is green with the replay inside it.
