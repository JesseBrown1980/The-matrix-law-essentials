#!/usr/bin/env python3
"""Predictions for the power-of-three jump tests, written before the code they check.

Every value below comes from a closed formula or from the seat's scratch runs of
2026-11-06 (jump3_test.py, hash_answer_test.py). Nothing here calls the code
under test. Output: HBP rows, json=0, one PREDICT row per checked field.

basis=KNOWN_FROM_SEAT  value printed by the seat's scratch run (rerun 2026-11-06, same bytes)
basis=DERIVED          value from the formula given in the row's note
basis=CONJECTURE       expected value, no derivation in hand

usage: make_predictions.py > PREDICTIONS-20261106.hbp
"""

D27 = 27
P = 1000081

ROWS = []


def pred(pid, row, where, field, op, value, basis, note):
    ROWS.append("PREDICT|id=%s|row=%s|where=%s|field=%s|op=%s|value=%s|basis=%s|note=%s|json=0"
                % (pid, row, where, field, op, value, basis, note))


def note(pid, claim, basis, checked_by):
    ROWS.append("PREDICT_NOTE|id=%s|claim=%s|basis=%s|checked_by=%s|json=0" % (pid, claim, basis, checked_by))


# ---- jump3_ntt.py --self-test ------------------------------------------------
# cases: name, ring D, stencil size S, k list, largest k stepped directly, NTT defined
CASES = (
    ("signed_return", 27, 3, (1, 2, 3, 6, 9, 10, 20, 40, 80), 10, 1),
    ("curved_sum", 27, 3, (1, 2, 3, 9, 40), 9, 1),
    ("shift", 27, 1, (1, 2, 3, 40), 3, 1),
    ("ring45", 45, 3, (1, 2, 3, 6, 9, 20), 9, 1),
    ("ring81", 81, 3, (1, 2, 3), 3, 0),
)
SEAT_SHA = {3: "7d05096afc24b08c", 6: "cc3581ab163de3a7", 9: "ee4b40fb753713ab",
            10: "81323c66a0e8a92f", 20: "8ef1ccc86cd3ff0c", 40: "177066f133a5c9ef"}

pred("J000", "JUMP3_HEAD", "law=power_of_three_jump", "p", "eq", P, "DERIVED", "kernel_prime")
pred("J001", "JUMP3_HEAD", "law=power_of_three_jump", "p_minus_1_factors", "eq", "2^4*3^3*5*463",
     "DERIVED", "1000080=16*27*5*463")
pred("J002", "JUMP3_CASE", "case=signed_return", "root", "eq", 951846, "KNOWN_FROM_SEAT",
     "jump3_test.py_root27")
pred("J003", "JUMP3_CASE", "case=signed_return", "zero_modes", "eq", 0, "KNOWN_FROM_SEAT",
     "jump3_test.py_eigenvalues_nonzero=1")
pred("J004", "JUMP3_CASE", "case=signed_return", "backward", "eq", "DEFINED", "KNOWN_FROM_SEAT",
     "every_eigenvalue_nonzero")
pred("J005", "JUMP3_CASE", "case=curved_sum", "zero_modes", "eq", 2, "DERIVED",
     "1+z+z^2=0_iff_z_primitive_cube_root;w^m_is_one_iff_m_in_{9,18}")
pred("J006", "JUMP3_CASE", "case=curved_sum", "zero_mode_list", "eq", "9,18", "DERIVED", "as_J005")
pred("J007", "JUMP3_CASE", "case=curved_sum", "backward", "eq", "REFUSED_ZERO_MODE", "DERIVED", "as_J005")
pred("J008", "JUMP3_CASE", "case=curved_sum", "root", "eq", 951846, "KNOWN_FROM_SEAT", "same_ring_27")
pred("J009", "JUMP3_CASE", "case=shift", "zero_modes", "eq", 0, "DERIVED", "eigenvalues_w^m_are_units")
pred("J010", "JUMP3_CASE", "case=ring45", "zero_modes", "eq", 0, "CONJECTURE",
     "lambda=2-z+3/z_vanishes_iff_z_in_{3,-1};-1_is_no_45th_root;3^45_mod_p_not_evaluated")
pred("J011", "JUMP3_CASE", "case=ring45", "ntt", "eq", "DEFINED", "DERIVED", "45_divides_1000080")
pred("J012", "JUMP3_CASE", "case=ring81", "ntt", "eq", "REFUSED_RING_DOES_NOT_DIVIDE_P_MINUS_1",
     "DERIVED", "3-adic_valuation_of_1000080_is_3")
for name, d, s, ks, dmax, has_ntt in CASES:
    pred("J0E" + name, "JUMP3_CASE", "case=" + name, "eigen_mulmods", "eq", d * s if has_ntt else 0,
         "DERIVED", "D*stencil_size")

pid = 100
for name, d, s, ks, dmax, has_ntt in CASES:
    bound = 2 * d * d + 2 * d + d * 39          # transforms + spectrum + scale + D powers of <= 39 products (e < 2^20)
    for k in ks:
        n = 3 ** k
        w = "case=%s&k=%d" % (name, k)
        pred("J%d" % pid, "JUMP3", w, "n", "eq", n, "DERIVED", "3^k"); pid += 1
        pred("J%d" % pid, "JUMP3", w, "cubings", "eq", k, "DERIVED", "one_cube_per_power"); pid += 1
        pred("J%d" % pid, "JUMP3", w, "cube_mulmods", "eq", 2 * k * d ** 3 + d * d, "DERIVED",
             "2kD^3+D^2"); pid += 1
        pred("J%d" % pid, "JUMP3", w, "direct_mulmods", "eq", n * d * s, "DERIVED", "3^k*D*S"); pid += 1
        if has_ntt:
            pred("J%d" % pid, "JUMP3", w, "ntt_mulmods", "le", bound, "DERIVED",
                 "2D^2+2D+39D_independent_of_k"); pid += 1
            known = name == "signed_return" and k in (3, 6, 9, 10, 20, 40)
            pred("J%d" % pid, "JUMP3", w, "cube_equals_ntt", "eq", 1,
                 "KNOWN_FROM_SEAT" if known else "DERIVED", "same_linear_map"); pid += 1
            pred("J%d" % pid, "JUMP3", w, "leaf_forward", "eq", 1,
                 "KNOWN_FROM_SEAT" if known else "DERIVED", "M*M^n=M^(n+1)"); pid += 1
            if name == "curved_sum":
                pred("J%d" % pid, "JUMP3", w, "leaf_back", "eq", "REFUSED_ZERO_MODE", "DERIVED", "as_J005"); pid += 1
            else:
                pred("J%d" % pid, "JUMP3", w, "leaf_back", "eq", 1,
                     "KNOWN_FROM_SEAT" if known else "DERIVED", "inverse_eigenvalues_exist"); pid += 1
        else:
            pred("J%d" % pid, "JUMP3", w, "ntt_mulmods", "eq", "NOT_DEFINED", "DERIVED", "as_J012"); pid += 1
            pred("J%d" % pid, "JUMP3", w, "leaf_forward", "eq", 1, "DERIVED", "checked_against_direct_n+1"); pid += 1
            pred("J%d" % pid, "JUMP3", w, "leaf_back", "eq", "NOT_DEFINED", "DERIVED", "as_J012"); pid += 1
        if k <= dmax:
            known = name == "signed_return" and k in (3, 6, 9, 10)
            pred("J%d" % pid, "JUMP3", w, "direct_equals_jump", "eq", 1,
                 "KNOWN_FROM_SEAT" if known else "DERIVED", "same_linear_map"); pid += 1
        else:
            pred("J%d" % pid, "JUMP3", w, "direct", "eq", "NOT_RUN", "DERIVED", "3^k_steps_skipped"); pid += 1
        if name == "shift":
            pred("J%d" % pid, "JUMP3", w, "equals_start", "eq", 1 if k >= 3 else 0, "DERIVED",
                 "R^27=I"); pid += 1
        if name == "signed_return" and k in SEAT_SHA:
            pred("J%d" % pid, "JUMP3", w, "sha256", "prefix", SEAT_SHA[k], "KNOWN_FROM_SEAT",
                 "jump3_test.py_16_hex"); pid += 1
for m in (9, 18):
    w = "case=curved_sum&mode=%d" % m
    for f in ("preimages_distinct", "images_equal_after_1_step", "images_equal_after_3^40_cube",
              "images_equal_after_3^40_ntt", "image_equals_case_jump"):
        pred("J%d" % pid, "JUMP3_ZERO_MODE", w, f, "eq", 1, "DERIVED", "lambda_m=0_kills_v_m"); pid += 1
pred("J%d" % pid, "JUMP3_VERIFY", "status=PASS", "jump_rows", "eq", 27, "DERIVED", "9+5+4+6+3"); pid += 1
pred("J%d" % pid, "JUMP3_VERIFY", "status=PASS", "zero_mode_rows", "eq", 2, "DERIVED", "modes_9_18"); pid += 1

# ---- transposition_memo.py --self-test ---------------------------------------
COUNTS = {
    6: "000:63,012:81,021:81,102:81,111:90,120:81,201:81,210:81,222:90",
    9: "000:2187,012:2187,021:2187,102:2187,111:2106,120:2187,201:2187,210:2187,222:2268",
    12: "000:59535,012:59049,021:59049,102:59049,111:58806,120:59049,201:59049,210:59049,222:58806",
}
MINMAX = {6: (63, 90), 9: (2106, 2268), 12: (58806, 59535)}
pid = 500
for L in (6, 9, 12):
    w = "L=%d" % L
    paths = 3 ** L
    for f, v, b, nt in (
            ("paths", paths, "KNOWN_FROM_SEAT", "3^L"),
            ("enumeration_turns", L * paths, "DERIVED", "L*3^L"),
            ("distinct_positions", 9, "KNOWN_FROM_SEAT", "digit_sum_mod_3_fixed"),
            ("memo_keys", 9, "DERIVED", "one_key_per_position"),
            ("memo_hits", paths - 9, "DERIVED", "3^L-9"),
            ("reuse_avg", paths // 9, "KNOWN_FROM_SEAT", "3^L/9"),
            ("paths_min", MINMAX[L][0], "DERIVED", "multinomial_sums"),
            ("paths_max", MINMAX[L][1], "DERIVED", "multinomial_sums"),
            ("counts", COUNTS[L], "DERIVED", "multinomial_sums;L=6_matches_GPT-6_table_63_81_90"),
            ("recurrence_updates", 81 * L, "DERIVED", "27_states*3_turns*L"),
            ("recurrence_equals_enumeration", 1, "DERIVED", "C_{n+1}(s)=sum_{p->s}C_n(p)"),
            ("invariant_class", L % 3, "DERIVED", "digit_sum=L_mod_3"),
            ("invariant_candidates", 9, "DERIVED", "27/3"),
            ("reachable_equals_invariant_class", 1, "DERIVED", "every_class_member_reached_for_L>=2"),
            ("invariant_updates", 27 * L, "DERIVED", "9_states*3_turns*L"),
            ("invariant_equals_recurrence", 1, "DERIVED", "off_class_counts_are_zero")):
        pred("T%d" % pid, "TRANSPOSITION", w, f, "eq", v, b, nt); pid += 1
for L in (27, 54):
    w = "L=%d" % L
    for f, v, b, nt in (
            ("sum", 3 ** L, "DERIVED", "3^L"),
            ("sum_equals_3^L", 1, "DERIVED", "every_move_order_ends_somewhere"),
            ("support", 9, "DERIVED", "invariant_class"),
            ("count_012", 3 ** (L - 2), "DERIVED", "character_sum:all_sqrt3_terms_cancel_for_distinct_digits_L=0_mod_3"),
            ("count_012_equals_3^(L-2)", 1, "DERIVED", "as_count_012"),
            ("perm_012_all_equal", 1, "DERIVED", "digit_permutation_symmetry")):
        pred("T%d" % pid, "PATHCOUNT", w, f, "eq", v, b, nt); pid += 1
for k in (3, 4, 40):
    w = "k=%d" % k
    pred("T%d" % pid, "PATHCOUNT_JUMP", w, "cubings", "eq", k, "DERIVED", "one_cube_per_power"); pid += 1
    pred("T%d" % pid, "PATHCOUNT_JUMP", w, "cube_mulmods", "eq", 2 * k * D27 ** 3 + D27 * D27, "DERIVED",
         "2kD^3+D^2"); pid += 1
    pred("T%d" % pid, "PATHCOUNT_JUMP", w, "sum_equals_3^L_mod_p", "eq", 1, "DERIVED", "3^L_mod_p"); pid += 1
    pred("T%d" % pid, "PATHCOUNT_JUMP", w, "count_012_equals_3^(L-2)_mod_p", "eq", 1, "DERIVED",
         "as_count_012"); pid += 1
    pred("T%d" % pid, "PATHCOUNT_JUMP", w, "support_mod_p", "eq", 9, "CONJECTURE",
         "class_counts_000_111_222_nonzero_mod_p_not_derived"); pid += 1
    if k <= 4:
        pred("T%d" % pid, "PATHCOUNT_JUMP", w, "equals_recurrence_mod_p", "eq", 1, "DERIVED",
             "same_linear_map"); pid += 1
    else:
        pred("T%d" % pid, "PATHCOUNT_JUMP", w, "equals_recurrence_mod_p", "eq", "NOT_RUN", "DERIVED",
             "3^40_layers_skipped"); pid += 1
w = "moves_hidden=10"
for f, v, b, nt in (("position", "010", "KNOWN_FROM_SEAT", "hash_answer_test.py_recovered=0,1,0"),
                    ("bound_full", 27, "DERIVED", "27_positions"),
                    ("tried_full", 4, "KNOWN_FROM_SEAT", "place_in_enumeration_order"),
                    ("bound_invariant", 9, "DERIVED", "digit_sum=10_mod_3=1"),
                    ("tried_invariant", 2, "DERIVED", "001_then_010"),
                    ("recovered_equals_hidden", 1, "KNOWN_FROM_SEAT", "hash_answer_test.py")):
    pred("T%d" % pid, "FINITE_RECOVER", w, f, "eq", v, b, nt); pid += 1
pred("T%d" % pid, "UNSTRUCTURED", "input_bits=64", "candidates", "eq", 2 ** 64, "DERIVED", "2^64"); pid += 1
pred("T%d" % pid, "TRANSMEMO_VERIFY", "status=PASS", "transposition_rows", "eq", 3, "DERIVED", "L=6,9,12"); pid += 1

note("B01", "python_and_rust_print_the_same_rows_byte_for_byte", "DERIVED", "verify.sh_cmp")
note("B02", "ntt_mulmods_does_not_grow_with_k", "DERIVED", "rows_J*_op=le")
note("B03", "timing_rows_are_outside_the_byte_comparison", "DERIVED", "--timing_flag_only")

print("PREDICTIONS|written_before=jump3_ntt.py,transposition_memo.py,crates/jump3|date=2026-11-06|rows=%d|"
      "seat_scratch=jump3_test.py,hash_answer_test.py|json=0" % len(ROWS))
for r in ROWS:
    print(r)
