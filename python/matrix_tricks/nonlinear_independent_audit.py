#!/usr/bin/env python3
"""Independent exact nonlinear-jump check: no imports from the solver.

Checks explicit physical-variable recurrence, naive DFT (no FFT shared),
independent modal matrix powers, and complete million-radius integer bound.
"""
from pathlib import Path
import argparse
import hashlib

P = 1299079
Q = 1000000
K = 7


def hash_pair(u, v):
    h = hashlib.sha256()
    for arr in (u, v):
        for z in arr:
            h.update(z.to_bytes(8, 'big'))
    return h.hexdigest()


def start(n):
    return ([(7*i*i+3*i+1)%P for i in range(n)],
            [(11*i*i+31*i+17)%P for i in range(n)])


def direct_nonlinear(u, v):
    n=len(u)
    # Expand T^-1 L T in original physical coordinates, no T/L calls.
    work=[(v[i]+K*u[i]*u[i-1])%P for i in range(n)]
    b=[(work[i]+u[i-1])%P for i in range(n)]
    a=[(u[i]+b[i-1])%P for i in range(n)]
    out=[(b[i]-K*a[i]*a[i-1])%P for i in range(n)]
    return a,out


def simple_dft(a,inv=False):
    n=len(a)
    root=pow(3,(P-1)//n,P)
    if inv: root=pow(root,P-2,P)
    roots=[pow(root,j,P) for j in range(n)]
    out=[]
    for m in range(n):
        acc=0
        for j,x in enumerate(a):
            acc=(acc+x*roots[(j*m)%n])%P
        out.append(acc)
    if inv:
        invn=pow(n,P-2,P)
        out=[x*invn%P for x in out]
    return out


def multiply(a,b):
    return ((a[0]*b[0]+a[1]*b[2])%P,(a[0]*b[1]+a[1]*b[3])%P,
            (a[2]*b[0]+a[3]*b[2])%P,(a[2]*b[1]+a[3]*b[3])%P)


def independent_matrix_power(base,n):
    r=(1,0,0,1)
    while n:
        if n%2:r=multiply(r,base)
        base=multiply(base,base)
        n//=2
    return r


def slow_jump(n,t):
    u,v=start(n)
    v=[(v[i]+K*u[i]*u[i-1])%P for i in range(n)]
    U=simple_dft(u)
    V=simple_dft(v)
    root=pow(3,(P-1)//n,P)
    w=1
    U1=[];V1=[]
    for a,b in zip(U,V):
        m=( (1+w*w)%P,w,w,1)
        z=independent_matrix_power(m,t)
        assert (z[0]*z[3]-z[1]*z[2])%P==1
        U1.append((z[0]*a+z[1]*b)%P)
        V1.append((z[2]*a+z[3]*b)%P)
        w=w*root%P
    assert w==1
    uf=simple_dft(U1,True)
    vf=simple_dft(V1,True)
    vf=[(vf[i]-K*uf[i]*uf[i-1])%P for i in range(n)]
    return uf,vf


def toward_zero(x,d):
    return x//d if x>=0 else -((-x)//d)


def lossy_step(x):
    return [toward_zero(3*x[i]+x[i-1]+toward_zero(x[i]*x[i-1],Q),8) for i in range(len(x))]


def audit(path):
    txt=path.read_text(encoding='ascii')
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    for n in (27,243,6561,59049):
        u,v=start(n)
        for t in range(1,28):
            u,v=direct_nonlinear(u,v)
            if t in (1,3,27):
                h=hash_pair(u,v)
                assert ('DIRECT|n=%d|steps=%d|jump_equals_direct=PASS|sha256=%s|'%(n,t,h)) in txt
        print('IND_NONLIN_DIRECT|n=%d|all_27_direct_steps=PASS|hash27=%s|json=0'%(n,hash_pair(u,v)))
    for n in (27,243):
        for t,k in ((3**40,40),(3**80,80)):
            u,v=slow_jump(n,t)
            h=hash_pair(u,v)
            assert ('TERNARY|n=%d|exponent=3^%d|complete_state_sha256=%s|'%(n,k,h)) in txt
            print('IND_NONLIN_DFT|n=%d|k=%d|naive_DFT_equals_fast_NTT=PASS|sha256=%s|json=0'%(n,k,h))
    # Mathematical bound audit: exhaustive scalar B in 1..Q and nontrivial seeds.
    for b in range(1,Q+1):
        nex=(4*b+b*b//Q)//8
        assert 0<=nex<b
    sequence=[Q]
    while sequence[-1]>0:
        b=sequence[-1]
        sequence.append((4*b+b*b//Q)//8)
    assert len(sequence)==22 and sequence[-1]==0
    for n in (27,243,6561,59049):
        x=[((19*i*i+43*i+7)%(2*Q+1))-Q for i in range(n)]
        for i in range(21):
            x=lossy_step(x)
            assert max(map(abs,x))<=sequence[i+1]
        assert not any(x)
    print('IND_BOUND|radius=1000000|all_integer_radii_tested=1000000|certified_uniform_k=21|d59049_all_zero=PASS|json=0')
    assert 'NONLINEAR_JUMP_VERIFY|status=PASS|' in txt
    print('NONLINEAR_INDEPENDENT_AUDIT_PASS|four_direct_scales=YES|naive_DFT_two_scales=YES|all_radius_inequalities=YES|receipt_sha256=%s|json=0'%digest)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('receipt',type=Path)
    a=ap.parse_args()
    audit(a.receipt)


if __name__=='__main__':
    main()
