#!/usr/bin/env python3
"""Two CO molecular stretch degrees of freedom: phase cancellation versus real error.

FROZEN FIRST: data/molecular_co/PREDICT-CO-PAIR-20261106.hbp, commit 6734e0c.
No physical CO-CO coupling measurement: kappa = omega**2/9 chosen as a toy.
No binary floats, no JSON.
"""
from decimal import Decimal, getcontext
from co_molecular_slices_20261106 import (
    D, ZERO, ONE, TWO, H, OMEGA, X0, harmonic, matpow, matvec,
    triadic_jump, cosine, fmt
)
getcontext().prec=65

KAPPA=OMEGA*OMEGA/D(9)
OMEGA_ANTI=(OMEGA*OMEGA+TWO*KAPPA).sqrt()
N_SHORT=2187
N_LONG=3**15
THRESHOLD=D("1e-50")


def out(label,**fields):
    print(label+"|"+"|".join(str(k)+"="+str(v) for k,v in fields.items())+"|json=0")


def integration_error(omega,n):
    m=matpow(harmonic(omega*H),n)
    discrete=matvec(m,(X0,ZERO))[0]
    continuous=X0*cosine(omega*H*D(n))
    return discrete,continuous,discrete-continuous


def direct_coupled_pair(amp1,amp2,n):
    x1,x2=amp1,amp2
    v1,v2=ZERO,ZERO
    w2=OMEGA*OMEGA
    for _ in range(n):
        a1=-w2*x1-KAPPA*(x1-x2)
        a2=-w2*x2-KAPPA*(x2-x1)
        hv1=v1+H*a1/TWO
        hv2=v2+H*a2/TWO
        x1=x1+H*hv1
        x2=x2+H*hv2
        a1=-w2*x1-KAPPA*(x1-x2)
        a2=-w2*x2-KAPPA*(x2-x1)
        v1=hv1+H*a1/TWO
        v2=hv2+H*a2/TWO
    return x1,x2


def main():
    out("CO_PAIR_SETUP",molecules=2,model="CO_CLASSICAL_HARMONIC_STRETCH",
        coupling="OMEGA_SQUARED_DIVIDED_BY_9_UNCALIBRATED",
        coupling_measurement="NONE",
        prediction_commit="6734e0c0aba6cb79292b427ad6a6491e389e5534")
    assert OMEGA_ANTI*OMEGA_ANTI==OMEGA*OMEGA+TWO*KAPPA or abs(
        OMEGA_ANTI*OMEGA_ANTI-(OMEGA*OMEGA+TWO*KAPPA)) < D("1e-60")
    for n in (N_SHORT,N_LONG):
        q,ref,err=integration_error(OMEGA,n)
        assert abs(err)>D("1e-12")
        # Two independent molecules starting in phase: both have same error.
        coherent_errors=(err,err)
        c_mean=sum(coherent_errors,ZERO)/TWO
        c_rms=(sum((e*e for e in coherent_errors),ZERO)/TWO).sqrt()
        assert abs(c_mean-err)<THRESHOLD
        assert abs(c_rms-abs(err))<THRESHOLD
        out("CO_PAIR_IN_PHASE",ticks=n,elapsed_fs=fmt(H*D(n)),
            mean_signed_error_angstrom=fmt(c_mean),individual_rmse_angstrom=fmt(c_rms),
            vanished="NO")
        # Two independent harmonic vibrational phases 0 and pi: errors cancel,
        # but the individual absolute error is identical to one-molecule error.
        antiphase_errors=(err,-err)
        a_mean=sum(antiphase_errors,ZERO)/TWO
        a_rms=(sum((e*e for e in antiphase_errors),ZERO)/TWO).sqrt()
        assert a_mean==ZERO and abs(a_rms-abs(err))<THRESHOLD
        out("CO_PAIR_OPPOSITE_PHASE",ticks=n,elapsed_fs=fmt(H*D(n)),
            mean_signed_error_angstrom=fmt(a_mean),individual_rmse_angstrom=fmt(a_rms),
            vanished_in_mean="YES",individual_error_vanished="NO")
        # Genuine two-body coupled model: antisymmetric normal mode
        # frequency sqrt(omega**2 + 2*kappa).
        aq,aref,aerr=integration_error(OMEGA_ANTI,n)
        assert abs(aerr)>D("1e-12")
        out("CO_PAIR_COUPLED_ANTI",ticks=n,elapsed_fs=fmt(H*D(n)),
            collective_frequency_squared_ratio="11/9",
            mean_signed_error_angstrom=fmt(ZERO),
            individual_rmse_angstrom=fmt(abs(aerr)),
            result="NONZERO_MOLECULE_ERROR")
        # Symmetric normal mode is immune to pair coupling and retains its bias.
        out("CO_PAIR_COUPLED_COHERENT",ticks=n,mean_signed_error_angstrom=fmt(err),
            individual_rmse_angstrom=fmt(abs(err)),
            coupling_force="ZERO_BUT_TIME_ERROR_NONZERO")
        if n==N_SHORT:
            anti1,anti2=direct_coupled_pair(X0,-X0,n)
            coh1,coh2=direct_coupled_pair(X0,X0,n)
            assert abs(anti1-aq)<THRESHOLD and abs(anti2+aq)<THRESHOLD
            assert abs(coh1-q)<THRESHOLD and abs(coh2-q)<THRESHOLD
            out("CO_PAIR_DIRECT",ticks=n,
                antisymmetric_mismatch_angstrom=fmt(max(abs(anti1-aq),abs(anti2+aq))),
                coherent_mismatch_angstrom=fmt(max(abs(coh1-q),abs(coh2-q))),
                validation="PASS")
    # Long-time triadic matrix jump is a faithful execution of discrete
    # dynamics; that does not mean its discretization bias vanishes.
    n=N_LONG
    tri=matvec(triadic_jump(harmonic(OMEGA_ANTI*H),15),(X0,ZERO))[0]
    binary=matvec(matpow(harmonic(OMEGA_ANTI*H),n),(X0,ZERO))[0]
    continuous=X0*cosine(OMEGA_ANTI*H*D(n))
    assert abs(tri-binary)<D("1e-46")
    assert abs(binary-continuous)>D("1e-12")
    out("CO_PAIR_JUMP",fine_ticks=n,
        triadic_binary_difference_angstrom=fmt(abs(tri-binary)),
        continuous_model_difference_angstrom=fmt(abs(binary-continuous)),
        discrete_step_execution="PASS",continuum_equivalence="NOT_EXACT")
    count=0
    assert abs(c_mean)>D("1e-12")
    count+=1
    out("CO_PAIR_RED",case="coherent_average_must_cancel",verdict="REJECT")
    assert a_rms>D("1e-12")
    count+=1
    out("CO_PAIR_RED",case="opposite_phase_rmse_zero",verdict="REJECT")
    assert KAPPA!=ZERO and KAPPA==OMEGA*OMEGA/D(9)
    count+=1
    out("CO_PAIR_RED",case="coupling_is_measured_CO_CO_force",verdict="REJECT")
    assert abs(OMEGA_ANTI*OMEGA_ANTI-OMEGA*OMEGA)>D("1e-12")
    count+=1
    out("CO_PAIR_RED",case="wrong_coupling_frequency",verdict="REJECT")
    assert count==4
    out("CO_PAIR_VERIFY",status="PASS",registered_predictions=8,
        molecules=2,red_controls=count,
        phase_average_cancellation=1,individual_error_survives=1,
        real_CO_CO_measurement="NOT_PROVIDED",
        quantum_molecule="NOT_SIMULATED",YM4_continuum="NOT_PROVED",json_zero=1)


if __name__=="__main__":
    main()
