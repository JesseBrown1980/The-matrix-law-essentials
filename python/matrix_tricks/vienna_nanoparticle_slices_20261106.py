#!/usr/bin/env python3
"""Integer-only reproduction of the simplified Vienna nanoparticle toy model.

Historical 75-digit Decimal version remains byte-identical in
data/vienna_nanoparticle/evidence/ORIGINAL-DECIMAL-SOURCE-20261106.txt.
This test recomputes the same truncated phase-only grating amplitudes using
fixed-point INTEGER arithmetic, then checks every historic sampled visibility
to a declared 10^-12 absolute tolerance. The old hosted result is evidence,
not an integer-only receipt and not a replay of the measured experiment.
"""
from functools import lru_cache
from hashlib import sha1
from math import isqrt
from pathlib import Path
import sys

PREC=100
S=10**PREC
TOL=10**(PREC-12)
CIS_TOL=10**(PREC-40)
KMAX=18
MAX_POWER=16
VIRTUAL=3**MAX_POWER
GRIDS=(9,27,81,243,729)
PHASES=(0,1,2,3)
ROOT=Path(__file__).resolve().parents[2]

def require(ok,which):
    if not ok:
        print("VIENNA_INT_FAIL|check="+which+"|json=0",flush=True)
        raise SystemExit(1)

def d(s):
    """Parse a signed finite decimal/scientific text to 100-digit scaled int."""
    s=str(s).strip().lower()
    sign=-1 if s.startswith("-") else 1
    if s[:1] in ("+","-"):
        s=s[1:]
    exp=0
    if "e" in s:
        s,x=s.split("e")
        exp=int(x)
    if "." in s:
        left,right=s.split(".")
    else:
        left,right=s,""
    digits=(left+right)
    require(digits.isdigit(),"decimal_digits")
    shift=PREC+exp-len(right)
    a=int(digits)
    return sign*(a*10**shift if shift>=0 else a//(10**(-shift)))

def m(a,b):
    return (a*b)//S

def q(a,b):
    require(b!=0,"division_nonzero")
    return (a*S)//b

PI=d("3.141592653589793238462643383279502884197169399375105820974944592307816406")
H=d("6.62607015e-34")
AMU=d("1.66053906660e-27")
MASS=m(d("172000"),AMU)
SPEED=d("160")
PERIOD=d("133e-9")
LENGTH=d("0.983")
LAMBDA=q(H,m(MASS,SPEED))
TALBOT=q(m(PERIOD,PERIOD),LAMBDA)
RATIO=q(LENGTH,TALBOT)
FLIGHT_TIME=q(LENGTH,SPEED)
CONTRAST=d("0.30")
A=CONTRAST//2
B=S-A

def zadd(a,b):
    return a[0]+b[0],a[1]+b[1]

def zmul(a,b):
    return m(a[0],b[0])-m(a[1],b[1]),m(a[0],b[1])+m(a[1],b[0])

def zconj(a):
    return a[0],-a[1]

def znorm(a):
    return m(a[0],a[0])+m(a[1],a[1])

def zroot(a):
    return isqrt(znorm(a)*S)

def cis(angle):
    a=angle%(2*PI)
    if a>PI:
        a-=2*PI
    a2=a*a
    ts=a
    tc=S
    ss=ts
    cc=tc
    for j in range(1,190):
        ts=-(ts*a2)//(S*S*(2*j)*(2*j+1))
        tc=-(tc*a2)//(S*S*(2*j-1)*(2*j))
        ss+=ts
        cc+=tc
        if abs(ts)<2 and abs(tc)<2:
            return cc,ss
    require(False,"trig_convergence")

@lru_cache(None)
def bessel_j(k,phase):
    require(k>=0 and phase in PHASES,"Bessel_argument")
    half=d(phase)//2
    term=S
    for j in range(1,k+1):
        term=(term*half)//(S*j)
    out=term
    for j in range(1,160):
        term=-(term*half*half)//(S*S*j*(j+k))
        out+=term
        if abs(term)<2:
            return out
    require(False,"Bessel_convergence")

def grating(k,phase):
    v=bessel_j(abs(k),phase)
    rotate=abs(k)%4
    return ((v,0),(0,v),(-v,0),(0,-v))[rotate]

def initial_modes():
    return {-1:(CONTRAST//4,0),0:(B,0),1:(CONTRAST//4,0)}

def flight(modes,ratio=RATIO):
    return {k:zmul(v,cis(-m(PI,ratio)*k*k)) for k,v in modes.items()}

def phase_mask(modes,phase):
    out={}
    for n,v in modes.items():
        for k in range(-KMAX,KMAX+1):
            j=n+k
            out[j]=zadd(out.get(j,(0,0)),zmul(v,grating(k,phase)))
    return out

def output_modes(phase,ratio=RATIO):
    return flight(phase_mask(flight(initial_modes(),ratio),phase),ratio)

def spectral(modes):
    p0=0
    p1=(0,0)
    for k,v in modes.items():
        p0+=znorm(v)
        p1=zadd(p1,zmul(v,zconj(modes.get(k-1,(0,0)))))
    return q(m(A,zroot(p1)),m(B,p0))

def sampled(modes,n):
    rot=cis((2*PI)//n)
    phase=(S,0)
    i0=0
    i1=(0,0)
    ordered=sorted(modes.items())
    for _ in range(n):
        powers={0:(S,0)}
        for k in range(1,KMAX+3):
            powers[k]=zmul(powers[k-1],phase)
            powers[-k]=zconj(powers[k])
        psi=(0,0)
        for k,v in ordered:
            psi=zadd(psi,zmul(v,powers[k]))
        power=znorm(psi)
        i0+=power
        i1=zadd(i1,(m(phase[0],power),-m(phase[1],power)))
        phase=zmul(phase,rot)
    require(i0>0,"positive_grid_norm")
    return q(m(A,zroot(i1)),m(B,i0))

def g1_classical():
    rho0=m(B,B)+m(A,A)//2
    rho1=m(B,A)
    return q(m(A,rho1),m(B,rho0))

def old_file(path,blob_sha):
    raw=ROOT.joinpath(path).read_bytes()
    obj=b"blob "+str(len(raw)).encode()+bytes((0,))+raw
    require(sha1(obj).hexdigest()==blob_sha,"historical_git_blob_pin")
    return raw.decode("utf-8")

def fields(line):
    return {a:b for a,b in (x.split("=",1) for x in line.split("|")[1:]
             if "=" in x)}

def old_values():
    text=old_file(Path("data/vienna_nanoparticle/MEASURED-20261106.hbp"),
                  "ee5065794fe55a879b1fb89afdcae8b8636e501a")
    archived=old_file(Path("data/vienna_nanoparticle/evidence/"
                  "ORIGINAL-DECIMAL-SOURCE-20261106.txt"),
                  "cd153dfaae045ff85ee2026da82268ce2299d353")
    require("getcontext().prec = 75" in archived,"historical_precision_provenance")
    rows={}
    reference={}
    for line in text.splitlines():
        if line.startswith("VIENNA_SLICE|"):
            f=fields(line)
            rows[int(f["G2_phase_radians"]),int(f["sites"])]=d(f["quantum_model_visibility"])
        if line.startswith("VIENNA_QUANTUM|"):
            f=fields(line)
            reference[int(f["phase_radians"])]=d(f["spectral_visibility"])
    require(len(rows)==20 and len(reference)==4,"complete_archived_reference")
    return rows,reference

def cube_jump(z,k):
    for _ in range(k):
        z=zmul(zmul(z,z),z)
    return z

def main():
    require(len(sys.argv)==2 and sys.argv[1]=="--self-test","self_test_invocation")
    require(VIRTUAL==43046721,"virtual_grid")
    require(abs(LAMBDA-d("1.449968282802470e-14"))<TOL,"debroglie")
    require(abs(RATIO-d("0.8057656294843282"))<TOL,"Talbot_ratio")
    expected,ref=old_values()
    worst=0
    aliases=0
    classical=g1_classical()
    require(abs(classical-d("0.03066439522998296"))<TOL,
            "classical_straight_ray")
    for phase in PHASES:
        modes=output_modes(phase)
        require(max(abs(k) for k in modes)==19,"truncated_bandlimit")
        spec=spectral(modes)
        delta=abs(spec-ref[phase])
        worst=max(worst,delta)
        require(delta<TOL,"old_spectral_phase_"+str(phase))
        for n in GRIDS:
            got=sampled(modes,n)
            error=abs(got-expected[phase,n])
            worst=max(worst,error)
            require(error<TOL,"old_slice_phase_"+str(phase)+"_n"+str(n))
            if n>=27:
                require(abs(got-spec)<10**(PREC-10),"fine_grid_alias_bound")
        if phase==3:
            aliases=abs(sampled(modes,9)-sampled(modes,729))
        print("VIENNA_INT_PHASE|phase="+str(phase)
              +"|visibility_micro="+str((spec*10**6)//S)
              +"|reference_agree_1e12=1|json=0")
    require(aliases>d("0.001"),"nine_cell_alias_control")
    require(VIRTUAL>2*38+1,"virtual_first_harmonic_no_alias")
    max_jump=0
    checks=0
    for k in range(1,MAX_POWER+1):
        n=3**k
        for mode in (1,2,3,5,11):
            theta=-m(PI,RATIO)*mode*mode
            small=cis(theta//n)
            jumped=cube_jump(small,k)
            direct=cis(theta)
            er=isqrt(znorm((jumped[0]-direct[0],jumped[1]-direct[1]))*S)
            max_jump=max(max_jump,er)
            require(er<CIS_TOL,"ternary_phase_jump")
            checks+=1
    require(checks==80,"ternary_checks")
    print("VIENNA_INT_SOURCE|old_decimal_git_blob=cd153dfaae045ff85ee2026da82268ce2299d353"
          +"|original_measured_git_blob=ee5065794fe55a879b1fb89afdcae8b8636e501a"
          +"|original_measured_unchanged=1|json=0")
    print("VIENNA_INT_VERIFY|status=PASS|engine=INTEGER_FIXED_POINT"
          +"|precision_decimal_digits="+str(PREC)
          +"|referenced_slice_rows=20|referenced_spectral_rows=4"
          +"|absolute_tolerance=1e-12|time_jump_checks="+str(checks)
          +"|max_error_less_1e-12=1|virtual_grid=43046721"
          +"|grating_model=TOY_ONLY|experiment_replayed=NO|json=0")

if __name__=="__main__":
    main()
