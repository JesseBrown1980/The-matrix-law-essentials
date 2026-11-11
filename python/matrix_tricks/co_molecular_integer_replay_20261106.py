#!/usr/bin/env python3
"""Integer-only replay of the three CO molecular toy tests of November 6, 2026.

The original 65-digit Decimal sources are kept byte-identical, unexecuted, in
data/molecular_co/evidence/. This tool recomputes every row of the three
retained MEASURED receipts in 100-digit fixed-point INTEGER arithmetic and
compares field by field. A classical fitted CO stretch model and an
uncalibrated coupling toy: no observed molecular motion is replayed.
"""
from hashlib import sha1, sha256
from math import isqrt
from pathlib import Path
import sys

PREC = 100
S = 10**PREC
NOISE = 10**(PREC - 40)
ROOT = Path(__file__).resolve().parents[2]
ENGINE = "INTEGER_FIXED_POINT_100"
ENGINE_FIELDS = {
    "precision": ("DECIMAL_65_NOT_EXACT_RATIONAL", ENGINE + "_NOT_EXACT_RATIONAL"),
    "exact_discrete_jump": ("DECIMAL_MATRIX_POWER", ENGINE + "_MATRIX_POWER"),
}
EVIDENCE = (
    ("data/molecular_co/evidence/ORIGINAL-DECIMAL-SLICES-SOURCE-20261106.txt",
     "32c0972b0a712c0ad9f6cd1d3cf39132ff9d1c31"),
    ("data/molecular_co/evidence/ORIGINAL-DECIMAL-PAIR-SOURCE-20261106.txt",
     "cb349b486a718d30f8e84430e76a6a05816b48b0"),
    ("data/molecular_co/evidence/ORIGINAL-DECIMAL-ENSEMBLE-SOURCE-20261106.txt",
     "ae8a0a0564d7b7c03a8d43c3dee3778dd286b2f4"),
)
RECEIPTS = (
    ("slices", "data/molecular_co/MEASURED-CO-MOLECULAR-SLICES-20261106.hbp",
     "ede652ea6ce55b5587fe98345afae0de39325e94", 23),
    ("pair", "data/molecular_co/MEASURED-CO-PAIR-20261106.hbp",
     "25db2a47bbac5a6f88856c91331763bf559b98a0", 16),
    ("ensemble", "data/molecular_co/MEASURED-CO-ENSEMBLES-20261106.hbp",
     "68731f98004fad1e33ec64dbe154a781446353e1", 33),
)


def require(ok, which):
    if not ok:
        print("CO_INT_FAIL|check=" + str(which) + "|json=0", flush=True)
        raise SystemExit(1)


def d(text):
    """Parse a signed finite decimal or scientific text to a scaled integer."""
    s = str(text).strip().lower()
    sign = -1 if s.startswith("-") else 1
    if s[:1] in ("+", "-"):
        s = s[1:]
    exp = 0
    if "e" in s:
        s, x = s.split("e")
        exp = int(x)
    left, _, right = s.partition(".")
    digits = left + right
    require(digits.isdigit(), "decimal_digits")
    shift = PREC + exp - len(right)
    a = int(digits)
    return sign * (a * 10**shift if shift >= 0 else a // 10**(-shift))


def numeric(text):
    t = text.lstrip("+-")
    mantissa, _, exponent = t.upper().partition("E")
    return (mantissa.replace(".", "", 1).isdigit() and mantissa.count(".") <= 1
            and (exponent == "" or exponent.lstrip("+-").isdigit()))


def m(a, b):
    return (a * b) // S


def q(a, b):
    require(b != 0, "division_nonzero")
    return (a * S) // b


def root(a):
    require(a >= 0, "root_nonnegative")
    return isqrt(a * S)


def sci(value, places):
    """Scientific text with `places` decimals, rounded half to even."""
    if value == 0:
        return "0." + "0" * places + "E+0"
    a = abs(value)
    width = len(str(a))
    exponent = width - 1 - PREC
    drop = width - (places + 1)
    if drop > 0:
        unit = 10**drop
        lead, rest = divmod(a, unit)
        if 2 * rest > unit or (2 * rest == unit and lead % 2):
            lead += 1
        if lead == 10**(places + 1):
            lead //= 10
            exponent += 1
    else:
        lead = a * 10**(-drop)
    text = str(lead)
    return (("-" if value < 0 else "") + text[0] + "." + text[1:] + "E"
            + ("-" if exponent < 0 else "+") + str(abs(exponent)))


def fmt(value):
    return sci(value, 12)


ZERO = 0
ONE = S
PI = d("3.141592653589793238462643383279502884197169399375105820974944")
EQUILIBRIUM = d("1.128")
BETA = d("2.298")
DISSOCIATION_EV = d("11.07")
X0 = d("0.050")
# Planck time 5.391247e-44 s times 3^57, in femtoseconds.
H = (5391247 * 3**57 * S) // 10**35
# omega = 2 pi c * 2169.81 cm^-1, in radians per femtosecond.
OMEGA = (2 * PI * 29979245800 * 216981) // 10**17
OM2 = m(OMEGA, OMEGA)
STEP = m(OMEGA, H)
# Planck length 1.616255e-35 m times 3^45, 3^46, 3^47, in angstrom.
SPATIAL = [(1616255 * 3**k * S) // 10**31 for k in (45, 46, 47)]
TIMES = (9, 27, 81, 243, 729, 2187)
SCALES = (1, 3, 9)
N_LONG = 3**15
KAPPA = OM2 // 9
OMEGA_ANTI = root(OM2 + 2 * KAPPA)
NS = (1, 3, 9, 27, 81, 243, 729)


def row(rows, tag, **fields):
    rows.append(tag + "|" + "|".join(str(k) + "=" + str(v) for k, v in fields.items())
                + "|json=0")


def exponential(t):
    term = ONE
    acc = ONE
    for j in range(1, 400):
        term = (term * t) // (S * j)
        acc += term
        if abs(term) < 2:
            return acc
    require(False, "exponential_convergence")


def cosine(x):
    period = 2 * PI
    x = x % period
    if x > PI:
        x -= period
    x2 = x * x
    term = ONE
    acc = ONE
    for j in range(1, 400):
        term = -(term * x2) // (S * S * (2 * j - 1) * (2 * j))
        acc += term
        if abs(term) < 2:
            return acc
    require(False, "cosine_convergence")


def sine(x):
    return cosine(PI // 2 - x)


def acceleration(x):
    z = exponential(-m(BETA, x))
    return -m(m(q(OM2, BETA), z), ONE - z)


def energy(x, v):
    z = exponential(-m(BETA, x))
    w = q(m(BETA, v), OMEGA)
    return m(DISSOCIATION_EV, m(ONE - z, ONE - z) + m(w, w))


def verlet(x, v, h):
    half = v + m(h, acceleration(x)) // 2
    newx = x + m(h, half)
    return newx, half + m(h, acceleration(newx)) // 2


def rk4(x, v, h):
    k1x, k1v = v, acceleration(x)
    k2x, k2v = v + m(h, k1v) // 2, acceleration(x + m(h, k1x) // 2)
    k3x, k3v = v + m(h, k2v) // 2, acceleration(x + m(h, k2x) // 2)
    k4x, k4v = v + m(h, k3v), acceleration(x + m(h, k3x))
    return (x + m(h, k1x + 2 * k2x + 2 * k3x + k4x) // 6,
            v + m(h, k1v + 2 * k2v + 2 * k3v + k4v) // 6)


def matmul(a, b):
    return ((m(a[0][0], b[0][0]) + m(a[0][1], b[1][0]), m(a[0][0], b[0][1]) + m(a[0][1], b[1][1])),
            (m(a[1][0], b[0][0]) + m(a[1][1], b[1][0]), m(a[1][0], b[0][1]) + m(a[1][1], b[1][1])))


def exact_matmul(a, b):
    """Unscaled integer product: the denominators S multiply outside."""
    return ((a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]),
            (a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]))


def matvec(a, x):
    return (m(a[0][0], x[0]) + m(a[0][1], x[1]), m(a[1][0], x[0]) + m(a[1][1], x[1]))


def harmonic(s):
    s2 = m(s, s)
    c = ONE - s2 // 2
    return ((c, s), (-m(s, ONE - s2 // 4), c))


def matpow(a, n, product=matmul, identity=((S, 0), (0, S))):
    outm = identity
    while n:
        if n % 2:
            outm = product(outm, a)
        a = product(a, a)
        n //= 2
    return outm


def triadic_jump(a, k):
    for _ in range(k):
        a = matmul(matmul(a, a), a)
    return a


def scaled(a, factor):
    return tuple(tuple(x * factor for x in line) for line in a)


def sha_state(x, v):
    return sha256((sci(x, 40) + "|" + sci(v, 40)).encode()).hexdigest()


def direct_verlet_states():
    data = {}
    start = energy(X0, ZERO)
    for scale in SCALES:
        x, v = X0, ZERO
        maxenergy = ZERO
        collected = {}
        for tick in range(1, 2187 // scale + 1):
            x, v = verlet(x, v, H * scale)
            drift = abs(q(energy(x, v), start) - ONE)
            maxenergy = max(maxenergy, drift)
            if tick * scale in TIMES:
                collected[tick * scale] = (x, v)
        data[scale] = (collected, maxenergy)
    return data


def independent_reference():
    step = H // 9
    x, v = X0, ZERO
    samples = {}
    for i in range(1, 2187 * 9 + 1):
        x, v = rk4(x, v, step)
        if i % 9 == 0 and i // 9 in TIMES:
            samples[i // 9] = (x, v)
    return samples


def read_clock(n, scale):
    if n < 0 or n % scale:
        raise ValueError("CLOCK_OFF_SLICE")
    return n // scale


def nearest_site(length, spacing):
    sites, rest = divmod(length, spacing)
    if 2 * rest > spacing or (2 * rest == spacing and sites % 2):
        sites += 1
    return sites


def replay_slices(models, reference):
    rows = []
    row(rows, "CO_MODEL", name="12C16O_stretch", omega_e_cm1="2169.81",
        beta_inverse_angstrom="2.298", De_eV="11.07", r_eq_angstrom="1.128",
        initial_displacement_angstrom="0.050",
        scope="CLASSICAL_MORSE_NOT_QUANTUM_CHEMISTRY")
    row(rows, "CO_SCALES", space_fine_angstrom=fmt(SPATIAL[0]),
        space_middle_angstrom=fmt(SPATIAL[1]), space_coarse_angstrom=fmt(SPATIAL[2]),
        time_fine_fs=fmt(H), time_middle_fs=fmt(H * 3), time_coarse_fs=fmt(H * 9),
        assumed_planck_coordinates=1)
    require(acceleration(X0) < ZERO and energy(X0, ZERO) > ZERO, "morse_sign")
    errors = {}
    for scale in SCALES:
        states, drift = models[scale]
        errors[scale] = max(abs(states[n][0] - reference[n][0]) for n in TIMES)
        require(errors[scale] > ZERO, "morse_error_positive")
        row(rows, "CO_MORSE", scale=scale, dt_fs=fmt(H * scale), samples=len(TIMES),
            max_abs_bond_error_angstrom=fmt(errors[scale]),
            max_relative_energy_drift=fmt(drift),
            final_bond_angstrom=fmt(states[2187][0] + EQUILIBRIUM),
            final_ref_bond_angstrom=fmt(reference[2187][0] + EQUILIBRIUM),
            sha256=sha_state(*states[2187]))
    require(errors[1] < errors[3] < errors[9], "scale_error_order")
    row(rows, "CO_MORSE_COMPARE", error_order="FINE_LT_MIDDLE_LT_COARSE",
        max_ratio_middle_to_fine=fmt(q(errors[3], errors[1])),
        max_ratio_coarse_to_middle=fmt(q(errors[9], errors[3])), result="PASS")
    for n in TIMES:
        bond = reference[n][0] + EQUILIBRIUM
        labels = []
        for spacing in SPATIAL:
            sites = nearest_site(bond, spacing)
            require(2 * abs(sites * spacing - bond) <= spacing, "site_rounding")
            labels.append(str(sites))
        row(rows, "CO_SNAPSHOT", fine_ticks=n, time_fs=fmt(H * n),
            reference_bond_angstrom=fmt(bond), coarse_measured_bins=",".join(labels),
            reference_sha256=sha_state(*reference[n]))
    mf = harmonic(STEP)
    unit = ((1, 0), (0, 1))
    mf3 = matpow(mf, 3, exact_matmul, unit)
    mf9 = matpow(mf, 9, exact_matmul, unit)
    require(mf9 == exact_matmul(exact_matmul(mf3, mf3), mf3), "exact_ninth_power")
    require(mf3 != scaled(harmonic(3 * STEP), S**2), "naive_coarse_3")
    require(mf9 != scaled(harmonic(9 * STEP), S**8), "naive_coarse_9")
    row(rows, "CO_HARMONIC_EXACT", fractional_matrix_power="FINE_CUBED_AND_NINTH",
        exact_3_9="PASS", naive_coarse_equal="FALSE")
    for n in (9, 27, 81, 243):
        jumped = matvec(matpow(mf, n), (X0, ZERO))[0]
        state = (X0, ZERO)
        for _ in range(n):
            state = matvec(mf, state)
        error = abs(jumped - state[0])
        require(error < NOISE, "jump_equals_sequential")
        row(rows, "CO_HARMONIC", fine_ticks=n, block3="IDENTICAL_BY_POWER_LAW",
            block9="IDENTICAL_BY_POWER_LAW", jump_vs_sequential_angstrom=fmt(error),
            jump_vs_continuum_angstrom=fmt(abs(jumped - m(X0, cosine(STEP * n)))))
    jt = matvec(triadic_jump(mf, 15), (X0, ZERO))[0]
    jb = matvec(matpow(mf, N_LONG), (X0, ZERO))[0]
    require(abs(jt - jb) < NOISE, "triadic_equals_binary")
    row(rows, "CO_HARMONIC_LONG", k=15, fine_ticks=N_LONG, elapsed_fs=fmt(H * N_LONG),
        triadic_vs_binary_angstrom=fmt(abs(jt - jb)),
        difference_from_continuum_angstrom=fmt(abs(jt - m(X0, cosine(STEP * N_LONG)))),
        precision=ENGINE_FIELDS["precision"][1])
    negative = 0
    if acceleration(X0) < ZERO:
        negative += 1
        row(rows, "CO_RED", case="wrong_Morse_force_sign", verdict="REJECT")
    if scaled(mf3, S**6) != mf9:
        negative += 1
        row(rows, "CO_RED", case="wrong_3_vs_9_block", verdict="REJECT")
    if sha_state(X0, ZERO) != sha_state(X0 + d("0.000001"), ZERO):
        negative += 1
        row(rows, "CO_RED", case="tampered_state", verdict="REJECT")
    try:
        read_clock(10, 9)
    except ValueError:
        negative += 1
        row(rows, "CO_RED", case="off_slice_clock", verdict="REJECT")
    require(negative == 4, "slice_red_controls")
    row(rows, "CO_VERIFY", status="PASS", frozen_predictions=13, red_controls=negative,
        temporal_slices=3, spatial_slices=3,
        molecular_trajectory="CLASSICAL_MODEL_NOT_OBSERVED", quantum_molecule="NOT_TESTED",
        Yang_Mills_SU2="NOT_TESTED", continuum_equivalence="NOT_PROVED", json_zero=1)
    return rows


def integration_error(omega, n):
    discrete = matvec(matpow(harmonic(m(omega, H)), n), (X0, ZERO))[0]
    return discrete, discrete - m(X0, cosine(m(omega, H) * n))


def direct_coupled_pair(amp1, amp2, n):
    x1, x2 = amp1, amp2
    v1, v2 = ZERO, ZERO
    for _ in range(n):
        hv1 = v1 + m(H, -m(OM2, x1) - m(KAPPA, x1 - x2)) // 2
        hv2 = v2 + m(H, -m(OM2, x2) - m(KAPPA, x2 - x1)) // 2
        x1 = x1 + m(H, hv1)
        x2 = x2 + m(H, hv2)
        v1 = hv1 + m(H, -m(OM2, x1) - m(KAPPA, x1 - x2)) // 2
        v2 = hv2 + m(H, -m(OM2, x2) - m(KAPPA, x2 - x1)) // 2
    return x1, x2


def replay_pair():
    rows = []
    row(rows, "CO_PAIR_SETUP", molecules=2, model="CO_CLASSICAL_HARMONIC_STRETCH",
        coupling="OMEGA_SQUARED_DIVIDED_BY_9_UNCALIBRATED", coupling_measurement="NONE",
        prediction_commit="6734e0c0aba6cb79292b427ad6a6491e389e5534")
    require(abs(m(OMEGA_ANTI, OMEGA_ANTI) - OM2 - 2 * KAPPA) < NOISE, "anti_frequency")
    bound = d("1e-12")
    for n in (2187, N_LONG):
        elapsed = fmt(H * n)
        discrete, err = integration_error(OMEGA, n)
        require(abs(err) > bound, "single_error_nonzero")
        c_mean = (err + err) // 2
        c_rms = root((m(err, err) + m(err, err)) // 2)
        require(abs(c_mean - err) < NOISE and abs(c_rms - abs(err)) < NOISE, "in_phase")
        row(rows, "CO_PAIR_IN_PHASE", ticks=n, elapsed_fs=elapsed,
            mean_signed_error_angstrom=fmt(c_mean), individual_rmse_angstrom=fmt(c_rms),
            vanished="NO")
        a_mean = (err - err) // 2
        a_rms = root((m(err, err) + m(-err, -err)) // 2)
        require(a_mean == ZERO and abs(a_rms - abs(err)) < NOISE, "opposite_phase")
        row(rows, "CO_PAIR_OPPOSITE_PHASE", ticks=n, elapsed_fs=elapsed,
            mean_signed_error_angstrom=fmt(a_mean), individual_rmse_angstrom=fmt(a_rms),
            vanished_in_mean="YES", individual_error_vanished="NO")
        anti, aerr = integration_error(OMEGA_ANTI, n)
        require(abs(aerr) > bound, "anti_error_nonzero")
        row(rows, "CO_PAIR_COUPLED_ANTI", ticks=n, elapsed_fs=elapsed,
            collective_frequency_squared_ratio="11/9",
            mean_signed_error_angstrom=fmt(ZERO), individual_rmse_angstrom=fmt(abs(aerr)),
            result="NONZERO_MOLECULE_ERROR")
        row(rows, "CO_PAIR_COUPLED_COHERENT", ticks=n, mean_signed_error_angstrom=fmt(err),
            individual_rmse_angstrom=fmt(abs(err)),
            coupling_force="ZERO_BUT_TIME_ERROR_NONZERO")
        if n == 2187:
            anti1, anti2 = direct_coupled_pair(X0, -X0, n)
            coh1, coh2 = direct_coupled_pair(X0, X0, n)
            anti_gap = max(abs(anti1 - anti), abs(anti2 + anti))
            coherent_gap = max(abs(coh1 - discrete), abs(coh2 - discrete))
            require(anti_gap < NOISE and coherent_gap < NOISE, "direct_pair")
            row(rows, "CO_PAIR_DIRECT", ticks=n, antisymmetric_mismatch_angstrom=fmt(anti_gap),
                coherent_mismatch_angstrom=fmt(coherent_gap), validation="PASS")
    step = harmonic(m(OMEGA_ANTI, H))
    tri = matvec(triadic_jump(step, 15), (X0, ZERO))[0]
    binary = matvec(matpow(step, N_LONG), (X0, ZERO))[0]
    continuous = m(X0, cosine(m(OMEGA_ANTI, H) * N_LONG))
    require(abs(tri - binary) < NOISE and abs(binary - continuous) > bound, "pair_jump")
    row(rows, "CO_PAIR_JUMP", fine_ticks=N_LONG,
        triadic_binary_difference_angstrom=fmt(abs(tri - binary)),
        continuous_model_difference_angstrom=fmt(abs(binary - continuous)),
        discrete_step_execution="PASS", continuum_equivalence="NOT_EXACT")
    require(abs(c_mean) > bound, "red_coherent")
    row(rows, "CO_PAIR_RED", case="coherent_average_must_cancel", verdict="REJECT")
    require(a_rms > bound, "red_opposite")
    row(rows, "CO_PAIR_RED", case="opposite_phase_rmse_zero", verdict="REJECT")
    require(KAPPA != ZERO and KAPPA == OM2 // 9, "red_coupling")
    row(rows, "CO_PAIR_RED", case="coupling_is_measured_CO_CO_force", verdict="REJECT")
    require(abs(m(OMEGA_ANTI, OMEGA_ANTI) - OM2) > bound, "red_frequency")
    row(rows, "CO_PAIR_RED", case="wrong_coupling_frequency", verdict="REJECT")
    row(rows, "CO_PAIR_VERIFY", status="PASS", registered_predictions=8, molecules=2,
        red_controls=4, phase_average_cancellation=1, individual_error_survives=1,
        real_CO_CO_measurement="NOT_PROVIDED", quantum_molecule="NOT_SIMULATED",
        YM4_continuum="NOT_PROVED", json_zero=1)
    return rows


def analytic_error_matrix(n):
    jumped = matpow(harmonic(STEP), n)
    angle = STEP * n
    return jumped[0][0] - cosine(angle), jumped[0][1] - sine(angle)


def balanced_ensemble(count, delta11, delta12):
    vals = []
    for i in range(count):
        phase = (2 * PI * i) // count
        vals.append(m(delta11, m(X0, cosine(phase))) - m(delta12, m(X0, sine(phase))))
    return sum(vals) // count, root(sum(m(v, v) for v in vals) // count)


def accelerations(xs):
    nn = len(xs)
    return [-(m(OM2, xs[i]) + m(KAPPA, 2 * xs[i] - xs[i - 1] - xs[(i + 1) % nn]))
            for i in range(nn)]


def direct_coupled(count, ticks):
    """Every one of `count` coupled coordinates is stepped, none is copied."""
    require(count >= 3 and count % 3 == 0, "period_three_array")
    x = [X0 if i % 3 == 0 else -X0 // 2 for i in range(count)]
    v = [ZERO] * count
    for _ in range(ticks):
        a = accelerations(x)
        vhalf = [v[i] + m(H, a[i]) // 2 for i in range(count)]
        x = [x[i] + m(H, vhalf[i]) for i in range(count)]
        a = accelerations(x)
        v = [vhalf[i] + m(H, a[i]) // 2 for i in range(count)]
    return x


def replay_ensemble(models, reference):
    rows = []
    row(rows, "CO_ENSEMBLE_SETUP", molecule_count="1,3,9,27,81,243,729",
        short_fine_ticks=2187, long_fine_ticks=N_LONG,
        coupling="omega_squared_over_9_nearest_neighbor", fitted_CO_internal_bond_only=1,
        coupling_physically_calibrated=0, observed_physical_vacuum=0)
    one_error = models[1][0][2187][0] - reference[2187][0]
    require(abs(one_error) > d("1e-15"), "morse_error_nonzero")
    for count in NS:
        average = (one_error * count) // count
        rmse = root((m(one_error, one_error) * count) // count)
        require(abs(average - one_error) < NOISE and abs(rmse - abs(one_error)) < NOISE,
                "identical_copies")
        row(rows, "CO_ENSEMBLE_MORSE", molecules=count, common_time_fs=fmt(H * 2187),
            mean_signed_bond_error_angstrom=fmt(average),
            rms_individual_error_angstrom=fmt(rmse),
            independent_repeated_trajectories=1, does_not_vanish=1)
    results = {}
    for n in (2187, N_LONG):
        delta11, delta12 = analytic_error_matrix(n)
        expected_rms = m(X0, root((m(delta11, delta11) + m(delta12, delta12)) // 2))
        coherent = m(X0, delta11)
        require(expected_rms > d("1e-12"), "harmonic_error_nonzero")
        for count in NS:
            if count == 1:
                mean, rms, group = coherent, abs(coherent), "identical_phase"
            else:
                mean, rms = balanced_ensemble(count, delta11, delta12)
                require(abs(mean) < NOISE, "phase_mean_cancels")
                require(abs(rms - expected_rms) < NOISE, "rms_scale_unchanged")
                group = "balanced_phases"
            row(rows, "CO_ENSEMBLE_HARMONIC", fine_ticks=n, molecules=count, group=group,
                mean_signed_bond_error_angstrom=fmt(mean),
                rms_individual_error_angstrom=fmt(rms),
                continuum_comparison="ANALYTIC_COSINE",
                exact_discrete_jump=ENGINE_FIELDS["exact_discrete_jump"][1])
            results[n, count] = (mean, rms)
    omega_eff = m(OMEGA, root((4 * S) // 3))
    mode_amp = matvec(matpow(harmonic(m(omega_eff, H)), 2187), (X0, ZERO))[0]
    continuum_amp = m(X0, cosine(m(omega_eff, H) * 2187))
    coupled_error = abs(mode_amp - continuum_amp)
    require(coupled_error > d("1e-12"), "coupled_error_nonzero")
    for count in NS[1:]:
        xs = direct_coupled(count, 2187)
        pattern = [2 if i % 3 == 0 else -1 for i in range(count)]
        gaps = [abs(xs[i] - (pattern[i] * mode_amp) // 2) for i in range(count)]
        errs = [xs[i] - (pattern[i] * continuum_amp) // 2 for i in range(count)]
        max_continuum = max(abs(e) for e in errs)
        require(max(gaps) < NOISE, "coupled_mode_reproduced")
        require(abs(max_continuum - coupled_error) < NOISE, "coupled_continuum_error")
        row(rows, "CO_ENSEMBLE_COUPLED", molecules=count, step_count=2187,
            coupling="OMEGA_SQUARED_DIVIDED_BY_9", collective_mode="PERIOD_3",
            max_discrete_mode_mismatch_angstrom=fmt(max(gaps)),
            maximum_continuum_error_angstrom=fmt(max_continuum),
            rms_individual_error_angstrom=fmt(root(sum(m(e, e) for e in errs) // count)),
            result="COUPLED_ERROR_NONZERO")
    wrong_eff = m(OMEGA, root(2 * S))
    wrong_amp = matvec(matpow(harmonic(m(wrong_eff, H)), 2187), (X0, ZERO))[0]
    require(abs(wrong_amp - mode_amp) > d("1e-8"), "red_wrong_coupling")
    row(rows, "CO_ENSEMBLE_RED", case="wrong_coupling", verdict="REJECT")
    require(results[N_LONG, 729][1] > d("1e-4"), "red_false_zero")
    row(rows, "CO_ENSEMBLE_RED", case="false_zero_individual_RMSE", verdict="REJECT")
    require(KAPPA == OM2 // 9 and KAPPA != ZERO, "red_uncalibrated")
    row(rows, "CO_ENSEMBLE_RED", case="uncalibrated_physical_claim", verdict="REJECT")
    require(2188 % 9 != 0, "red_clock")
    row(rows, "CO_ENSEMBLE_RED", case="wrong_clock", verdict="REJECT")
    row(rows, "CO_ENSEMBLE_VERIFY", status="PASS", registered_predictions=11,
        molecules_max=729, red_controls=4, coherent_error_persists=1,
        balanced_mean_cancels=1, balanced_individual_rmse_persists=1,
        coupled_errors_persist=1, molecules_are_simulated_in_coupled_arm=1,
        not_quantum_molecular_chemistry=1, not_4D_Yang_Mills=1,
        physical_discreteness_not_shown=1, json_zero=1)
    return rows


def retained(relative, blob):
    data = ROOT.joinpath(relative).read_bytes()
    digest = sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
    require(digest == blob, "git_blob_" + relative)
    return data.decode()


def pairs(line):
    parts = line.split("|")
    return parts[0], [tuple(x.split("=", 1)) for x in parts[1:]]


def compare(name, old_rows, new_rows):
    """Identical text, or both below 1e-40, or a declared engine label."""
    require(len(old_rows) == len(new_rows), name + "_row_count")
    counts = {"identical": 0, "noise": 0, "engine": 0}
    for index, (old, new) in enumerate(zip(old_rows, new_rows)):
        where = name + "_row_" + str(index + 1)
        old_tag, old_fields = pairs(old)
        new_tag, new_fields = pairs(new)
        require(old_tag == new_tag and len(old_fields) == len(new_fields), where)
        for (old_key, old_value), (new_key, new_value) in zip(old_fields, new_fields):
            field = where + "_" + old_key
            require(old_key == new_key, field)
            if old_value == new_value:
                counts["identical"] += 1
            elif old_key in ENGINE_FIELDS:
                require((old_value, new_value) == ENGINE_FIELDS[old_key], field)
                counts["engine"] += 1
            else:
                require(numeric(old_value) and numeric(new_value), field)
                require(abs(d(old_value)) < NOISE and abs(d(new_value)) < NOISE, field)
                counts["noise"] += 1
    return counts


def main():
    require(len(sys.argv) == 2 and sys.argv[1] in ("--self-test", "--rows"),
            "self_test_invocation")
    for relative, blob in EVIDENCE:
        require("getcontext().prec=65" in retained(relative, blob), "historical_precision")
    models = direct_verlet_states()
    reference = independent_reference()
    replayed = {"slices": replay_slices(models, reference), "pair": replay_pair(),
                "ensemble": replay_ensemble(models, reference)}
    if sys.argv[1] == "--rows":
        for name in ("slices", "pair", "ensemble"):
            print("\n".join(replayed[name]))
        return
    total = {"identical": 0, "noise": 0, "engine": 0}
    for name, relative, blob, expected_rows in RECEIPTS:
        old_rows = retained(relative, blob).splitlines()[3:]
        require(len(old_rows) == expected_rows, name + "_receipt_rows")
        counts = compare(name, old_rows, replayed[name])
        for key in total:
            total[key] += counts[key]
        print("CO_INT_RECEIPT|name=" + name + "|rows=" + str(len(old_rows))
              + "|fields_identical=" + str(counts["identical"])
              + "|fields_below_1e-40_both=" + str(counts["noise"])
              + "|fields_engine_label=" + str(counts["engine"])
              + "|original_measured_git_blob=" + blob + "|original_measured_unchanged=1|json=0")
    print("CO_INT_SOURCE|" + "|".join(
        "old_decimal_git_blob_" + name + "=" + blob
        for name, (_, blob) in zip(("slices", "pair", "ensemble"), EVIDENCE))
        + "|executed=0|json=0")
    print("CO_INT_VERIFY|status=PASS|engine=INTEGER_FIXED_POINT"
          + "|precision_decimal_digits=" + str(PREC)
          + "|receipts=3|fields_identical=" + str(total["identical"])
          + "|fields_below_1e-40_both=" + str(total["noise"])
          + "|fields_engine_label=" + str(total["engine"])
          + "|fields_mismatched=0|coupled_coordinates_stepped_max=729"
          + "|model=CLASSICAL_CO_STRETCH_TOY|coupling_calibrated=NO"
          + "|experiment_replayed=NO|json=0")


if __name__ == "__main__":
    main()
