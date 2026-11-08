#!/usr/bin/env python3
"""Planck-unit coordinate scaling of the repository's exact 45-ring jump kernel.

Mathematical model only: finite-field advection, a declared phase clock,
a noise-free background, one controlled disturbance, and a heat-equation
counterexample. This does not model molecular forces or a quantum vacuum.
The physical Planck-length/tick interpretation is HYPOTHETICAL.
All simulated states and checks use integers, no float or JSON.
"""
from __future__ import annotations

from hashlib import sha256
from jump3_ntt import P, Ops, Ring, eigenvalues, ntt_jump, start_state, state_sha, step

D = 45
DX_PLANCK = 3**54
DT_PLANCK = 3**65
CLOCK_PERIOD = 3**11
assert DX_PLANCK*CLOCK_PERIOD == DT_PLANCK
assert (P-1)%D == 0


def record(tag, **fields):
    print(tag+'|'+'|'.join(str(k)+'='+str(v) for k,v in fields.items())+'|json=0')


def continuum_sample(v,n):
    """Sample u_t+v u_x=0 with exact one-coarse-cell translation per tick."""
    return [v[(i-n)%len(v)] for i in range(len(v))]


def clock_read(ticks):
    if not isinstance(ticks,int) or ticks<0:
        raise ValueError('INVALID_TICK')
    cells,phase=divmod(ticks,CLOCK_PERIOD)
    if phase:
        raise ValueError('OFF_CLOCK_BOUNDARY')
    coarse,residue=divmod(cells,DX_PLANCK)
    if residue:
        raise ValueError('OFF_COARSE_GRID')
    return coarse


def test_microclock_small():
    """Independent direct fine-tick stepping for a 3^4/3^3 miniature."""
    length=45*3**4
    phase_period=3**3
    for ticks in (0,1,26,27,81,729,3**7,2*3**7,3**8):
        position=0
        phase=0
        for _ in range(ticks):
            phase+=1
            if phase==phase_period:
                phase=0
                position=(position+1)%length
        assert (position,phase)==((ticks//phase_period)%length,ticks%phase_period)
        record('CLOCK_MINI',ticks=ticks,position=position,phase=phase,pass_check=1)


def main():
    ring=Ring(D,[(1,1)])
    assert ring.w is not None
    lam=eigenvalues(ring,Ops())
    assert all(v!=0 for v in lam)
    record('MOLECULAR_SETUP',dimension=D,prime=P,dx_planck_lengths=DX_PLANCK,
           dt_planck_times=DT_PLANCK,clock_period_planck_ticks=CLOCK_PERIOD,
           physical_calibration='ASSUMED_NOT_MEASURED',
           vacuum_status='IDEAL_ZERO_BACKGROUND_NOT_QUANTUM_VACUUM')
    test_microclock_small()
    starting=start_state(D)
    vacuum=[0]*D
    pulse=[0]*D
    pulse[5]=2
    disturbed=pulse[:]
    disturbed[13]=1
    times=(0,1,2,3,5,9,27,81,243,3**10,3**15)
    snapshot_times=(0,1,9,27)
    for arm,v in (('reference_pattern',starting),('ideal_vacuum',vacuum),
                  ('isolated_pulse',pulse),('single_perturbation',disturbed)):
        for n in times:
            ticks=n*DT_PLANCK
            if ticks>3**80:
                raise AssertionError('UNAUTHORIZED_JUMP_RANGE')
            assert clock_read(ticks)==n
            ops=Ops()
            actual=ntt_jump(ring,lam,v,n,ops)
            expected=continuum_sample(v,n)
            assert actual==expected
            if n<=243:
                direct=v[:]
                for _ in range(n):
                    direct=step(ring,direct,Ops())
                assert direct==actual
                status='PASS'
            else:
                status='NOT_RUN'
            assert sum(actual)==sum(v)
            record('MOLECULAR_JUMP',arm=arm,coarse_steps=n,planck_ticks=ticks,
                   site_displacement=n%D,ntt_products=ops.n,
                   sampled_continuum_advection='PASS',
                   direct_vs_jump=status,conserved_sum=sum(actual),
                   sha256=state_sha(actual))
            if n in snapshot_times and arm in ('ideal_vacuum','isolated_pulse','single_perturbation'):
                record('MOLECULAR_FRAME',arm=arm,coarse_steps=n,
                       planck_ticks=ticks,values=','.join(str(x) for x in actual),
                       sha256=state_sha(actual))
    for n in (0,1,9,27,81,3**15):
        a=ntt_jump(ring,lam,pulse,n,Ops())
        b=ntt_jump(ring,lam,disturbed,n,Ops())
        dif=[(v-u)%P for u,v in zip(a,b)]
        assert dif.count(1)==1 and dif.count(0)==D-1 and dif[(13+n)%D]==1
        record('MOLECULAR_NOISE',coarse_steps=n,changed_cells=1,
               changed_index=(13+n)%D,amplitude=1,
               propagation='EXACT_WITH_NO_DAMPING')
    for tag,ticks in (('off_time',DT_PLANCK+1),
                       ('off_position',DT_PLANCK+CLOCK_PERIOD),
                       ('not_full_step',3**64)):
        try:
            clock_read(ticks)
        except ValueError as exc:
            record('MOLECULAR_RED',control=tag,verdict='REFUSED',reason=str(exc))
        else:
            raise AssertionError('bad clock accepted')
    test=ntt_jump(ring,lam,pulse,1,Ops())
    altered=test[:]
    altered[0]=(altered[0]+1)%P
    assert state_sha(test)!=state_sha(altered)
    record('MOLECULAR_RED',control='hash_tamper',verdict='REFUSED',reason='SHA_MISMATCH')
    # Discrete lazy diffusion: [1/4,1/2,1/4] annihilates the
    # 54-site Nyquist alternating mode in one tick. The matching
    # continuum heat semigroup at k=pi/dx preserves strictly positive
    # amplitude exp(-pi^2/4) for every finite positive timestep.
    inv4=(3*P+1)//4
    assert 4*inv4%P==1
    diffusion=Ring(54,[(0,2*inv4%P),(1,inv4),(53,inv4)])
    mode=[1 if i%2==0 else P-1 for i in range(54)]
    one=step(diffusion,mode,Ops())
    assert one==[0]*54
    dl=eigenvalues(diffusion,Ops())
    assert any(v==0 for v in dl)
    assert ntt_jump(diffusion,dl,mode,1,Ops())==one
    record('DIFFUSION_COUNTEREXAMPLE',discrete_one_tick_amplitude=0,
           continuum_one_tick_amplitude='exp(-pi_squared/4)_positive',
           result='NOT_EQUAL',inverse='REFUSED_ZERO_EIGENVALUE')
    record('MOLECULAR_VERIFY',status='PASS',reference_equation='SAMPLED_CONTINUUM_ADVECTION',
           physical_molecule='NOT_SIMULATED',physical_vacuum='NOT_SIMULATED',
           diffusion_counterexample='PASS',json_zero=1)


if __name__=='__main__':
    main()
