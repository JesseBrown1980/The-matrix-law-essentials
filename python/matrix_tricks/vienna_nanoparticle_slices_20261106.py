#!/usr/bin/env python3
"""Vienna 2026 nanoparticle interferometry: a transparent reduced-model test.

Sources: Nature 649, 866-870 (2026); DOI 10.1038/s41586-025-09917-9.
The authors' actual data and Wigner-Weyl solver were NOT obtained or run.
This code is an independently specified coherent, phase-only Talbot model,
NOT the experiment's quantum prediction or its force-inclusive classical
baseline. The two published visibilities are contextual observations only.

Uses Decimal complex amplitudes (no binary floats), direct spatial quadrature,
and three-power free-flight phase jumps at measured apparatus length/velocity.
Predictions were frozen first: data/vienna_nanoparticle/PREDICTIONS-20261106.hbp.
"""
from decimal import Decimal, getcontext
from hashlib import sha256
from pathlib import Path

getcontext().prec = 75
D = Decimal
PI = D("3.141592653589793238462643383279502884197169399375105820974944592307816406")
ZERO = D(0)
ONE = D(1)
TWO = D(2)
H = D("6.62607015e-34")
AMU = D("1.66053906660e-27")
MASS = D(172000) * AMU
SPEED = D(160)
PERIOD = D("133e-9")
FLIGHT_LENGTH = D("0.983")
LAMBDA = H / (MASS * SPEED)
TALBOT = PERIOD * PERIOD / LAMBDA
RATIO = FLIGHT_LENGTH / TALBOT
FLIGHT_TIME = FLIGHT_LENGTH / SPEED
CONTRAST = D("0.30")
BESSEL_CUTOFF = 18
SPATIAL_GRIDS = (9, 27, 81, 243, 729)
PHASES = (0, 1, 2, 3)


class Z:
    __slots__ = ("r", "i")

    def __init__(self, r=ZERO, i=ZERO):
        self.r = r if isinstance(r, Decimal) else D(r)
        self.i = i if isinstance(i, Decimal) else D(i)

    def __add__(self, other):
        b = other if isinstance(other, Z) else Z(other)
        return Z(self.r+b.r, self.i+b.i)

    def __sub__(self, other):
        b = other if isinstance(other, Z) else Z(other)
        return Z(self.r-b.r, self.i-b.i)

    def __mul__(self, other):
        b = other if isinstance(other, Z) else Z(other)
        return Z(self.r*b.r-self.i*b.i, self.r*b.i+self.i*b.r)

    def __truediv__(self, other):
        b = other if isinstance(other, Z) else Z(other)
        den = b.r*b.r+b.i*b.i
        return Z((self.r*b.r+self.i*b.i)/den,
                 (self.i*b.r-self.r*b.i)/den)

    def conj(self):
        return Z(self.r, -self.i)

    def norm2(self):
        return self.r*self.r+self.i*self.i


def cis(angle):
    """70+-digit cosine/sine with exact pi-based argument reduction."""
    two_pi = TWO*PI
    a = angle % two_pi
    if a > PI:
        a -= two_pi
    t_s = a
    t_c = ONE
    s = a
    c = ONE
    eps = D("1e-72")
    for j in range(1, 180):
        t_s = -t_s*a*a / D((2*j)*(2*j+1))
        t_c = -t_c*a*a / D((2*j-1)*(2*j))
        s += t_s
        c += t_c
        if abs(t_s) < eps and abs(t_c) < eps:
            return Z(c, s)
    raise AssertionError("TAYLOR_SERIES_DID_NOT_CONVERGE")


def bessel_j(k, phi):
    """Integer-order Bessel J_k for positive integer phi and k>=0."""
    if k < 0:
        raise ValueError("INVALID_BESSEL_ORDER")
    z = D(phi)/TWO
    term = ONE
    for j in range(1, k+1):
        term *= z/D(j)
    total = term
    for j in range(1, 155):
        term = -term*z*z/(D(j)*D(j+k))
        total += term
        if abs(term)<D("1e-72"):
            break
    else:
        raise AssertionError("BESSEL_SERIES_DID_NOT_CONVERGE")
    return total


def grating_fourier(k, phi):
    """Jacobi-Anger: exp(i phi cos theta) = sum i^k J_k(phi) exp(ik theta)."""
    j = bessel_j(abs(k), phi)
    rotation = abs(k)%4
    if rotation==0:
        return Z(j)
    if rotation==1:
        return Z(ZERO,j)
    if rotation==2:
        return Z(-j)
    return Z(ZERO,-j)


def g1_state():
    # Amplitude mask 0.85 + 0.15*cos(theta).
    return {-1: Z(CONTRAST/D(4)),
             0: Z(ONE-CONTRAST/TWO),
             1: Z(CONTRAST/D(4))}


def free_flight(modes, ratio=RATIO):
    return {k:a*cis(-PI*ratio*D(k*k)) for k,a in modes.items()}


def phase_grating(modes, phase):
    result = {}
    for n,a in modes.items():
        for k in range(-BESSEL_CUTOFF,BESSEL_CUTOFF+1):
            key = n+k
            result[key]=result.get(key,Z()) + a*grating_fourier(k,phase)
    return result


def output_modes(phase,ratio=RATIO):
    return free_flight(phase_grating(free_flight(g1_state(),ratio),phase),ratio)


def spectral_visibility(modes):
    p0 = ZERO
    p1 = Z()
    for k,a in modes.items():
        p0 += a.norm2()
        p1 = p1 + a*modes.get(k-1,Z()).conj()
    return (CONTRAST/TWO)*p1.norm2().sqrt()/((ONE-CONTRAST/TWO)*p0)


def sampled_visibility(modes,n):
    # Finite-slice reading of density's Fourier mode with a real-space
    # molecular-scale observation grid. Aliasing is intentional at n=9.
    root=cis(TWO*PI/D(n))
    phase=Z(ONE)
    i0=ZERO
    i1=Z()
    ordered=tuple(sorted(modes.items()))
    for _ in range(n):
        powers={0:Z(ONE)}
        for k in range(1,BESSEL_CUTOFF+3):
            powers[k]=powers[k-1]*phase
            powers[-k]=powers[k].conj()
        psi=Z()
        for k,a in ordered:
            psi=psi+a*powers[k]
        intensity=psi.norm2()
        i0 += intensity
        i1 = i1 + phase.conj()*intensity
        phase=phase*root
    assert i0>ZERO
    return (CONTRAST/TWO)*i1.norm2().sqrt()/((ONE-CONTRAST/TWO)*i0)


def classical_straight_ray_visibility():
    # Deliberately simplified reference, NOT the paper's classical
    # trajectories, which include the optical dipole force.
    a=ONE-CONTRAST/TWO
    b=CONTRAST/TWO
    rho0=a*a+(b*b)/TWO
    rho1=a*b
    return b*rho1/(a*rho0)


def cubed(z,k):
    for _ in range(k):
        z=(z*z)*z
    return z


def report(tag,**fields):
    print(tag+"|"+"|".join(str(k)+"="+str(v) for k,v in fields.items())+"|json=0")


def ff(v):
    return format(v,'.15E')


def reference_sha(modes):
    source=";".join(str(k)+":"+format(a.r,'.35E')+":"+format(a.i,'.35E')
                    for k,a in sorted(modes.items()))
    return sha256(source.encode("ascii")).hexdigest()


def read_registered_predictions():
    path=Path(__file__).resolve().parents[2] / "data/vienna_nanoparticle/PREDICTIONS-20261106.hbp"
    data=path.read_text(encoding="utf-8")
    rows=[x for x in data.splitlines() if x.startswith("PREDICT|id=")]
    assert len(rows)==12 and "execution_status=NOT_RUN" in data
    return rows


def main():
    prediction_rows=read_registered_predictions()
    report("VIENNA_INPUT",article_DOI="10.1038/s41586-025-09917-9",
           source_data_DOI="10.5281/zenodo.17502163",
           measured_visibility_1="0.10+-0.01",
           measured_visibility_2="0.08+-0.01",
           particles="SODIUM_NANOPARTICLES_NOT_VIRUSES",
           raw_experiment_replay="NOT_RUN")
    report("VIENNA_SCALES",mass_Da=172000,velocity_m_s=160,
           grating_period_nm=133,grating_separation_m="0.983",
           deBroglie_wavelength_m=ff(LAMBDA),
           Talbot_length_m=ff(TALBOT),
           L_over_Talbot=ff(RATIO),flight_time_s=ff(FLIGHT_TIME))
    assert abs(LAMBDA-D("1.4499683e-14"))<D("1e-20")
    assert abs(RATIO-D("0.805766"))<D("0.000001")
    classical=classical_straight_ray_visibility()
    assert abs(classical-D("0.03065"))<D("0.0001")
    predicted={}
    highest={}
    for phi in PHASES:
        modes=output_modes(phi)
        high=spectral_visibility(modes)
        highest[phi]=high
        row={}
        for size in SPATIAL_GRIDS:
            visibility=sampled_visibility(modes,size)
            delta=abs(visibility-high)
            row[size]=visibility
            if size>=27:
                assert delta<D("1e-10"), ("SLICE_MISMATCH",phi,size,delta)
            report("VIENNA_SLICE",G2_phase_radians=phi,sites=size,
                   quantum_model_visibility=ff(visibility),
                   spectral_reference_visibility=ff(high),
                   delta_visibility=ff(delta),
                   classical_no_force_visibility=ff(classical),
                   model="TOY_COHERENT_THREE_GRATING")
        predicted[phi]=row
        report("VIENNA_QUANTUM",phase_radians=phi,
               spectral_visibility=ff(high),
               first_grating_model_amplitude="0.30",
               third_grating_model_contrast="0.30",
               sha256=reference_sha(modes))
    assert abs(highest[1]-D("0.09574"))<D("0.0001")
    assert abs(highest[2]-D("0.10730"))<D("0.0001")
    assert abs(highest[3]-D("0.03918"))<D("0.0001")
    assert abs(highest[1]-highest[0])>D("0.02")
    assert abs(highest[3]-highest[2])>D("0.05")
    assert abs(predicted[3][9]-predicted[3][729])>D("0.001")
    report("VIENNA_LOW_RES_RED",G2_phase_radians=3,sites=9,
           high_reference_sites=729,
           aliased_visibility=ff(predicted[3][9]),
           fine_visibility=ff(predicted[3][729]),
           discrepancy=ff(abs(predicted[3][9]-predicted[3][729])),
           verdict="FAIL_AS_EXPECTED")
    # Powers of three jump across the time of flight L/v for several
    # Fourier momentum modes; free flight is a fixed linear operator.
    max_error=ZERO
    jumps=0
    for k in range(1,6):
        n=3**k
        for mode in (1,2,3,5,11):
            theta=-PI*RATIO*D(mode*mode)
            step=cis(theta/D(n))
            jumped=cubed(step,k)
            direct=cis(theta)
            e=(jumped-direct).norm2().sqrt()
            max_error=max(max_error,e)
            assert e<D("1e-40"), ("TRIADIC_PHASE_MISMATCH",k,mode,e)
            jumps+=1
        report("VIENNA_TIME_JUMP",coarse_ticks=n,
               one_step_s=ff(FLIGHT_TIME/D(n)),
               whole_flight_s=ff(FLIGHT_TIME),
               max_phase_amplitude_error=ff(max_error),
               exact_real_arithmetic="NO_DECIMAL_75")
    report("VIENNA_NORM",intensity_before_G2=ff(sum(a.norm2() for a in g1_state().values())),
           intensity_after_first_free_flight=ff(
               sum(a.norm2() for a in free_flight(g1_state()).values())),
           conservation="WITH_DECIMAL_PRECISION")
    red=0
    # Another mass shifts the Talbot phase and observable.
    wrong= spectral_visibility(output_modes(1,RATIO/TWO))
    assert abs(wrong-highest[1])>D("0.001")
    report("VIENNA_RED",case="wrong_mass_changes_Talbot_phase",verdict="REJECT")
    red+=1
    wrong=spectral_visibility(output_modes(1,RATIO*D(4)))
    assert abs(wrong-highest[1])>D("0.001")
    report("VIENNA_RED",case="wrong_grating_period",verdict="REJECT")
    red+=1
    assert abs(highest[1]-highest[0])>D("0.02")
    report("VIENNA_RED",case="remove_phase_grating",verdict="REJECT")
    red+=1
    assert abs(predicted[3][9]-predicted[3][729])>D("0.001")
    report("VIENNA_RED",case="nine_cell_aliasing",verdict="REJECT")
    red+=1
    assert red==4
    report("VIENNA_VERIFY",status="PASS",frozen_predictions=len(prediction_rows),
           spatial_slices="27,81,243",independent_reference=729,
           phase_cases=len(PHASES),temporal_ternary_checks=jumps,
           red_controls=red,
           source_paper_raw_data="NOT_RETRIEVED",
           actual_experiment_reproduced="NO",
           paper_classical_dipole_force="NOT_IMPLEMENTED",
           Yang_Mills_4D_mass_gap="NOT_TESTED",
           fundamental_spacetime_discreteness="NOT_PROVED",json_zero=1)


if __name__=="__main__":
    main()
