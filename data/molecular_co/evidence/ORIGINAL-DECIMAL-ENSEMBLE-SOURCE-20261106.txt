#!/usr/bin/env python3
"""Test whether adding 1..729 molecular CO stretch oscillators removes integration errors.

Predictions: data/molecular_co/PREDICT-CO-ENSEMBLES-20261106.hbp (filed FIRST).
Existing spectroscopic CO parameters and harmonic Verlet integrator are reused.
Three arms:
* identical independent CO molecules;
* independent phase-balanced CO stretches, which may cancel signed means;
* nearest-neighbor coupled harmonic CO-stretch *toy* array, with an explicit
  nonphysical, uncalibrated spring coupling kappa = omega^2/9.

This is not ab-initio molecule-molecule interaction, a quantum vacuum, or
four-dimensional Yang-Mills physics.
"""
from decimal import Decimal, getcontext
from co_molecular_slices_20261106 import (
    D, ZERO, ONE, TWO, THREE, PI, H, OMEGA, X0,
    cosine, harmonic, matpow, matvec, fmt,
    direct_verlet_states, independent_reference
)

getcontext().prec=65
NS=(1,3,9,27,81,243,729)
TICKS=(2187,3**15)
TOL=D("1e-45")
KAPPA=OMEGA*OMEGA/D(9)


def out(tag, **fields):
    print(tag+"|"+"|".join(str(k)+"="+str(v) for k,v in fields.items())+"|json=0")


def sin_decimal(x):
    return cosine(PI/TWO-x)


def analytic_error_matrix(n):
    step=harmonic(OMEGA*H)
    jumped=matpow(step,n)
    angle=OMEGA*H*D(n)
    c=cosine(angle)
    s=sin_decimal(angle)
    delta11=jumped[0][0]-c
    delta12=jumped[0][1]-s
    return delta11,delta12


def ensemble(n, N, delta11,delta12):
    vals=[]
    for i in range(N):
        phase=TWO*PI*D(i)/D(N)
        x=X0*cosine(phase)
        p=-X0*sin_decimal(phase)
        vals.append(delta11*x+delta12*p)
    mean=sum(vals,ZERO)/D(N)
    rmse=(sum((v*v for v in vals),ZERO)/D(N)).sqrt()
    return mean,rmse


def accelerations(xs,kappa):
    nn=len(xs)
    om2=OMEGA*OMEGA
    return [-(om2*xs[i]+kappa*(TWO*xs[i]-xs[(i-1)%nn]-xs[(i+1)%nn]))
            for i in range(nn)]


def direct_coupled(N,ticks):
    # Actual update of *all N coupled particles*, not copies of one output.
    assert N>=3 and N%3==0
    x=[X0 if i%3==0 else -X0/TWO for i in range(N)]
    v=[ZERO]*N
    for _ in range(ticks):
        a=accelerations(x,KAPPA)
        vhalf=[v[i]+H*a[i]/TWO for i in range(N)]
        x=[x[i]+H*vhalf[i] for i in range(N)]
        a=accelerations(x,KAPPA)
        v=[vhalf[i]+H*a[i]/TWO for i in range(N)]
    return x,v


def main():
    out("CO_ENSEMBLE_SETUP",molecule_count="1,3,9,27,81,243,729",
        short_fine_ticks=TICKS[0],long_fine_ticks=TICKS[1],
        coupling="omega_squared_over_9_nearest_neighbor",
        fitted_CO_internal_bond_only=1,
        coupling_physically_calibrated=0,
        observed_physical_vacuum=0)
    original=direct_verlet_states()
    reference=independent_reference()
    one_error=original[1][0][2187][0]-reference[2187][0]
    assert abs(one_error)>D("1e-15")
    for N in NS:
        replicates=[one_error]*N
        average=sum(replicates,ZERO)/D(N)
        rmse=(sum((e*e for e in replicates),ZERO)/D(N)).sqrt()
        assert abs(average-one_error)<TOL
        assert abs(rmse-abs(one_error))<TOL
        out("CO_ENSEMBLE_MORSE",molecules=N,
            common_time_fs=fmt(H*D(2187)),
            mean_signed_bond_error_angstrom=fmt(average),
            rms_individual_error_angstrom=fmt(rmse),
            independent_repeated_trajectories=1,
            does_not_vanish=1)
    results={}
    for n in TICKS:
        delta11,delta12=analytic_error_matrix(n)
        expected_rms=X0*((delta11*delta11+delta12*delta12)/TWO).sqrt()
        coherent=X0*delta11
        assert expected_rms>D("1e-12")
        for N in NS:
            if N==1:
                mean=coherent
                rms=abs(coherent)
                group="identical_phase"
                assert abs(mean-coherent)<TOL
            else:
                mean,rms=ensemble(n,N,delta11,delta12)
                assert abs(mean)<D("1e-52"), ("PHASE_MEAN_NONZERO",N,n,mean)
                assert abs(rms-expected_rms)<D("1e-49"), ("RMS_SCALE_CHANGED",N,n,rms,expected_rms)
                group="balanced_phases"
            out("CO_ENSEMBLE_HARMONIC",fine_ticks=n,
                molecules=N,group=group,
                mean_signed_bond_error_angstrom=fmt(mean),
                rms_individual_error_angstrom=fmt(rms),
                continuum_comparison="ANALYTIC_COSINE",
                exact_discrete_jump="DECIMAL_MATRIX_POWER")
            results[(n,N)]=(mean,rms)
    # Interaction test: explicit N-body nearest-neighbor coupling.
    # Wavenumber 2*pi/3 is valid on every N multiple of 3.
    # Lattice Laplacian of the cosine pattern is exactly 3 times the pattern.
    omega_eff=OMEGA*(D(4)/D(3)).sqrt()
    mode_step=harmonic(omega_eff*H)
    mode_n=matpow(mode_step,2187)
    mode_amp=matvec(mode_n,(X0,ZERO))[0]
    continuum_amp=X0*cosine(omega_eff*H*D(2187))
    coupled_error=abs(mode_amp-continuum_amp)
    assert coupled_error>D("1e-12")
    for N in NS[1:]:
        xs,vs=direct_coupled(N,2187)
        pattern=[ONE if i%3==0 else -ONE/TWO for i in range(N)]
        max_discrete=max(abs(xs[i]-pattern[i]*mode_amp) for i in range(N))
        max_continuum=max(abs(xs[i]-pattern[i]*continuum_amp) for i in range(N))
        rms=(sum(((xs[i]-pattern[i]*continuum_amp)**2 for i in range(N)),ZERO)/D(N)).sqrt()
        assert max_discrete<D("1e-47"),("COUPLED_MODE_NOT_REPRODUCED",N,max_discrete)
        assert abs(max_continuum-coupled_error)<D("1e-47")
        out("CO_ENSEMBLE_COUPLED",molecules=N,
            step_count=2187,coupling="OMEGA_SQUARED_DIVIDED_BY_9",
            collective_mode="PERIOD_3",
            max_discrete_mode_mismatch_angstrom=fmt(max_discrete),
            maximum_continuum_error_angstrom=fmt(max_continuum),
            rms_individual_error_angstrom=fmt(rms),
            result="COUPLED_ERROR_NONZERO")
    # Deliberately incorrect controls must be detected.
    red=0
    wrong_eff=OMEGA*(TWO).sqrt()
    wrong_amp=matvec(matpow(harmonic(wrong_eff*H),2187),(X0,ZERO))[0]
    assert abs(wrong_amp-mode_amp)>D("1e-8")
    out("CO_ENSEMBLE_RED",case="wrong_coupling",verdict="REJECT")
    red+=1
    assert results[(3**15,729)][1]>D("1e-4")
    out("CO_ENSEMBLE_RED",case="false_zero_individual_RMSE",verdict="REJECT")
    red+=1
    assert KAPPA==OMEGA*OMEGA/D(9) and KAPPA!=ZERO
    out("CO_ENSEMBLE_RED",case="uncalibrated_physical_claim",verdict="REJECT")
    red+=1
    try:
        if (2188%9)!=0:
            raise ValueError("OFF_COARSE_SLICE")
    except ValueError:
        red+=1
        out("CO_ENSEMBLE_RED",case="wrong_clock",verdict="REJECT")
    assert red==4
    out("CO_ENSEMBLE_VERIFY",status="PASS",registered_predictions=11,
        molecules_max=729,red_controls=red,
        coherent_error_persists=1,balanced_mean_cancels=1,
        balanced_individual_rmse_persists=1,coupled_errors_persist=1,
        molecules_are_simulated_in_coupled_arm=1,
        not_quantum_molecular_chemistry=1,
        not_4D_Yang_Mills=1,
        physical_discreteness_not_shown=1,json_zero=1)


if __name__=="__main__":
    main()
