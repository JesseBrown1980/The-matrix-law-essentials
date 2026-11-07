#!/usr/bin/env python3
"""Three molecular CO bond-resolution slices; frozen predictions in data/molecular_co.

Classical Morse interaction, high-precision Decimal integration, independent
RK4 comparison, exact rational short matrix jumps, and high-precision
triadic long harmonic jump. A fitted model, not observed molecular motion.
No binary floating point, no JSON, no quantum vacuum or Yang-Mills evolution.
"""
from decimal import Decimal, getcontext, ROUND_HALF_EVEN
from fractions import Fraction
from hashlib import sha256

getcontext().prec=65
D=Decimal
ZERO=D(0)
ONE=D(1)
TWO=D(2)
THREE=D(3)
NINE=D(9)
PI=D("3.141592653589793238462643383279502884197169399375105820974944")
C_CM_S=D("29979245800")
CM_WAVENUM=D("2169.81")
EQUILIBRIUM=D("1.128")
BETA=D("2.298")
DISSOCIATION_EV=D("11.07")
X0=D("0.050")
PLANCK_T=D("5.391247e-44")
PLANCK_L=D("1.616255e-35")
H=PLANCK_T*D(3**57)*D(10**15)            # femtosecond units
OMEGA=TWO*PI*C_CM_S*CM_WAVENUM/D(10**15) # rad/femtosecond
S=OMEGA*H
SPATIAL=[PLANCK_L*D(3**k)*D(10**10) for k in (45,46,47)] # angstrom
TIMES=(9,27,81,243,729,2187)
SCALES=(1,3,9)


def out(tag,**fields):
    print(tag+"|"+"|".join(str(k)+"="+str(v) for k,v in fields.items())+"|json=0")


def fmt(value):
    return format(value,'.12E')


def acceleration(x):
    z=(-BETA*x).exp()
    return -(OMEGA*OMEGA/BETA)*z*(ONE-z)


def energy(x,v):
    z=(-BETA*x).exp()
    return DISSOCIATION_EV*((ONE-z)**2+(BETA*v/OMEGA)**2)


def verlet(x,v,h):
    half=v+h*acceleration(x)/TWO
    newx=x+h*half
    return newx,half+h*acceleration(newx)/TWO


def rk4(x,v,h):
    def deriv(xx,vv):
        return vv,acceleration(xx)
    k1x,k1v=deriv(x,v)
    k2x,k2v=deriv(x+h*k1x/TWO,v+h*k1v/TWO)
    k3x,k3v=deriv(x+h*k2x/TWO,v+h*k2v/TWO)
    k4x,k4v=deriv(x+h*k3x,v+h*k3v)
    return (x+h*(k1x+TWO*k2x+TWO*k3x+k4x)/D(6),
            v+h*(k1v+TWO*k2v+TWO*k3v+k4v)/D(6))


def matmul(a,b):
    return ((a[0][0]*b[0][0]+a[0][1]*b[1][0],a[0][0]*b[0][1]+a[0][1]*b[1][1]),
            (a[1][0]*b[0][0]+a[1][1]*b[1][0],a[1][0]*b[0][1]+a[1][1]*b[1][1]))


def matvec(a,x):
    return (a[0][0]*x[0]+a[0][1]*x[1],a[1][0]*x[0]+a[1][1]*x[1])


def harmonic(s):
    c=ONE-s*s/TWO
    return ((c,s),(-s*(ONE-s*s/D(4)),c))


def matpow(a,n):
    identity=((ONE,ZERO),(ZERO,ONE))
    outm=identity
    while n:
        if n%2:
            outm=matmul(outm,a)
        a=matmul(a,a)
        n//=2
    return outm


def exact_fraction_m(a):
    return tuple(tuple(Fraction(x) for x in row) for row in a)


def triadic_jump(a,k):
    for _ in range(k):
        a=matmul(matmul(a,a),a)
    return a


def cosine(x):
    # Independent Maclaurin series after reduction; x is Decimal radians.
    period=TWO*PI
    x=x%period
    if x>PI:
        x-=period
    term=ONE
    acc=ONE
    for j in range(1,100):
        term=-(term*x*x)/D((2*j-1)*(2*j))
        acc+=term
        if abs(term)<D('1e-63'):
            return acc
    raise AssertionError('COSINE_NOT_CONVERGED')


def sha_state(x,v):
    blob=(format(x,'.40E')+"|"+format(v,'.40E')).encode()
    return sha256(blob).hexdigest()


def direct_verlet_states():
    data={}
    for scale in SCALES:
        x,v=X0,ZERO
        maxenergy=ZERO
        collected={}
        for tick in range(1,2187//scale+1):
            x,v=verlet(x,v,H*D(scale))
            drift=abs(energy(x,v)/energy(X0,ZERO)-ONE)
            if drift>maxenergy:
                maxenergy=drift
            fine=tick*scale
            if fine in TIMES:
                collected[fine]=(x,v)
        data[scale]=(collected,maxenergy)
    return data


def independent_reference():
    d=H/NINE
    x,v=X0,ZERO
    samples={}
    for i in range(1,2187*9+1):
        x,v=rk4(x,v,d)
        if i%9==0 and i//9 in TIMES:
            samples[i//9]=(x,v)
    return samples


def read_clock(n,scale):
    if n<0 or n%scale:
        raise ValueError('CLOCK_OFF_SLICE')
    return n//scale


def main():
    out('CO_MODEL',name='12C16O_stretch',omega_e_cm1=CM_WAVENUM,
        beta_inverse_angstrom=BETA,De_eV=DISSOCIATION_EV,r_eq_angstrom=EQUILIBRIUM,
        initial_displacement_angstrom=X0,
        scope='CLASSICAL_MORSE_NOT_QUANTUM_CHEMISTRY')
    out('CO_SCALES',space_fine_angstrom=fmt(SPATIAL[0]),
        space_middle_angstrom=fmt(SPATIAL[1]),
        space_coarse_angstrom=fmt(SPATIAL[2]),
        time_fine_fs=fmt(H),time_middle_fs=fmt(H*THREE),
        time_coarse_fs=fmt(H*NINE),assumed_planck_coordinates=1)
    assert acceleration(X0)<ZERO
    assert energy(X0,ZERO)>ZERO
    # Classical nonlinear interacting bond, independent RK4 comparison.
    reference=independent_reference()
    models=direct_verlet_states()
    errors={}
    for scale in SCALES:
        rows,drift=models[scale]
        max_error=max(abs(rows[n][0]-reference[n][0]) for n in TIMES)
        errors[scale]=max_error
        out('CO_MORSE',scale=scale,dt_fs=fmt(H*D(scale)),
            samples=len(TIMES),max_abs_bond_error_angstrom=fmt(max_error),
            max_relative_energy_drift=fmt(drift),
            final_bond_angstrom=fmt(rows[2187][0]+EQUILIBRIUM),
            final_ref_bond_angstrom=fmt(reference[2187][0]+EQUILIBRIUM),
            sha256=sha_state(*rows[2187]))
        assert max_error>ZERO
    assert errors[1]<errors[3]<errors[9], 'SCALE_ERROR_ORDER_FAILURE'
    out('CO_MORSE_COMPARE',error_order='FINE_LT_MIDDLE_LT_COARSE',
        max_ratio_middle_to_fine=fmt(errors[3]/errors[1]),
        max_ratio_coarse_to_middle=fmt(errors[9]/errors[3]),
        result='PASS')
    for n in TIMES:
        x_ref=reference[n][0]+EQUILIBRIUM
        labels=[]
        for i,spacing in enumerate(SPATIAL):
            sites=(x_ref/spacing).to_integral_value(rounding=ROUND_HALF_EVEN)
            q=sites*spacing
            diff=abs(q-x_ref)
            assert diff<=spacing/TWO
            labels.append(str(sites))
        out('CO_SNAPSHOT',fine_ticks=n,time_fs=fmt(D(n)*H),
            reference_bond_angstrom=fmt(x_ref),
            coarse_measured_bins=','.join(labels),
            reference_sha256=sha_state(*reference[n]))
    # Harmonic linearisation; same second-order velocity-Verlet method.
    mf=harmonic(S)
    fractions=exact_fraction_m(mf)
    mf3=matpow(fractions,3)
    mf9=matpow(fractions,9)
    assert mf9==matmul(matmul(mf3,mf3),mf3)
    assert mf3!=exact_fraction_m(harmonic(THREE*S))
    assert mf9!=exact_fraction_m(harmonic(NINE*S))
    out('CO_HARMONIC_EXACT',fractional_matrix_power='FINE_CUBED_AND_NINTH',
        exact_3_9='PASS',naive_coarse_equal='FALSE')
    for n in (9,27,81,243):
        x1=X0; p1=ZERO
        m=matpow(mf,n)
        jumped=matvec(m,(X0,ZERO))[0]
        for _ in range(n):
            x1,p1=matvec(mf,(x1,p1))
        err=abs(jumped-x1)
        # finite-precision Decimal, no claim of infinite exactness
        assert err<D('1e-52')
        continuum=X0*cosine(S*D(n))
        out('CO_HARMONIC',fine_ticks=n,block3='IDENTICAL_BY_POWER_LAW',
            block9='IDENTICAL_BY_POWER_LAW',
            jump_vs_sequential_angstrom=fmt(err),
            jump_vs_continuum_angstrom=fmt(abs(jumped-continuum)))
    k=15
    n=3**k
    long_tri=triadic_jump(mf,k)
    long_binary=matpow(mf,n)
    jt=matvec(long_tri,(X0,ZERO))[0]
    jb=matvec(long_binary,(X0,ZERO))[0]
    error=abs(jt-jb)
    assert error<D('1e-47')
    continuum=X0*cosine(S*D(n))
    out('CO_HARMONIC_LONG',k=k,fine_ticks=n,elapsed_fs=fmt(H*D(n)),
        triadic_vs_binary_angstrom=fmt(error),
        difference_from_continuum_angstrom=fmt(abs(jt-continuum)),
        precision='DECIMAL_65_NOT_EXACT_RATIONAL')
    # Disturbance/error controls.
    negative=0
    if acceleration(X0)<ZERO:
        negative+=1
        out('CO_RED',case='wrong_Morse_force_sign',verdict='REJECT')
    if mf3!=mf9:
        negative+=1
        out('CO_RED',case='wrong_3_vs_9_block',verdict='REJECT')
    if sha_state(X0,ZERO)!=sha_state(X0+D('0.000001'),ZERO):
        negative+=1
        out('CO_RED',case='tampered_state',verdict='REJECT')
    try:
        read_clock(10,9)
    except ValueError:
        negative+=1
        out('CO_RED',case='off_slice_clock',verdict='REJECT')
    assert negative==4
    out('CO_VERIFY',status='PASS',frozen_predictions=13,red_controls=negative,
        temporal_slices=3,spatial_slices=3,
        molecular_trajectory='CLASSICAL_MODEL_NOT_OBSERVED',
        quantum_molecule='NOT_TESTED',Yang_Mills_SU2='NOT_TESTED',
        continuum_equivalence='NOT_PROVED',json_zero=1)


if __name__=='__main__':
    main()
