#!/usr/bin/env python3
"""Exact affine force jumps on a radix-3 NTT ring, 3^8 to 3^10 states.

One step is x <- (I+R-R^2)x + b, all arithmetic modulo 1,299,079.
The external force b is constant within each composed segment. The
radix-three Cooley-Tukey NTT makes jumps at 3^40, 3^80 and 2^40
possible without materializing a dense matrix or executing all steps.

The 12B detector addresses are sparse virtual coordinates, not 12B
independently evolved quantum degrees of freedom. All slit, density and
temperature source settings are synthetic controls, not astrophysics.
"""
from hashlib import sha256
from struct import pack, unpack_from
from pathlib import Path
import argparse

P=1_299_079
GENERATOR=3
VIRTUAL=12_000_000_000

def root(n):
    assert n in (6561,19683,59049) and (P-1)%n==0
    assert all(pow(GENERATOR,(P-1)//q,P)!=1 for q in (2,3,11))
    w=pow(GENERATOR,(P-1)//n,P)
    assert pow(w,n,P)==1 and pow(w,n//3,P)!=1
    return w

def ntt(x, inverse=False):
    n=len(x); r=root(n); k=0; m=n
    while m>1:
        assert m%3==0;m//=3;k+=1
    assert m==1
    a=x[:]
    for i in range(n):
        z=i; rev=0
        for _ in range(k): z,d=divmod(z,3);rev=rev*3+d
        if i<rev: a[i],a[rev]=a[rev],a[i]
    if inverse: r=pow(r,P-2,P)
    length=3
    while length<=n:
        step=length//3
        wm=pow(r,n//length,P)
        omega=pow(wm,step,P); omega2=omega*omega%P
        tw=1;tw2=1
        for j in range(step):
            for base in range(0,n,length):
                a0=base+j;a1=a0+step;a2=a1+step
                u=a[a0];v=a[a1]*tw%P;w=a[a2]*tw2%P
                a[a0]=(u+v+w)%P
                a[a1]=(u+omega*v+omega2*w)%P
                a[a2]=(u+omega2*v+omega*w)%P
            tw=tw*wm%P;tw2=tw2*wm*wm%P
        length*=3
    if inverse:
        inv=pow(n,P-2,P)
        a=[v*inv%P for v in a]
    return a

def eigen(n):
    r=root(n);w=1;out=[]
    for i in range(n):
        out.append((1+w-w*w)%P)
        w=w*r%P
    assert w==1
    return out

def source(n,strength=27,shift=0,slits=3):
    v=[0]*n
    places=(n//2,) if slits==1 else (n//2-n//9,n//2+n//9) if slits==2 else (n//2-n//9,n//2,n//2+n//9)
    for i in places:v[(i+shift)%n]=(v[(i+shift)%n]+strength)%P
    return v

def start(n):return [(7*i*i+3*i+1)%P for i in range(n)]

def matvec(x):return [(x[i]+x[i-1]-x[i-2])%P for i in range(len(x))]

def direct(x,b,n):
    for _ in range(n): x=[(y+z)%P for y,z in zip(matvec(x),b)]
    return x

def pairpow(lam,t):
    # Independent 2x2 affine semigroup exponentiation in each Fourier mode.
    ra,rb=1,0; a,b=lam,1
    while t:
        if t&1:
            ra,rb=(a*ra)%P,(a*rb+b)%P
        a,b=(a*a)%P,(a*b+b)%P
        t>>=1
    return ra,rb

def powers(lam,t):
    val=pow(lam,t,P)
    if lam==1: geo=t%P
    else:geo=(val-1)*pow((lam-1)%P,P-2,P)%P
    return val,geo

def jump(x,b,t,method='geom'):
    assert len(x)==len(b)
    X=ntt(x);B=ntt(b);l=eigen(len(x))
    if method=='geom': fac=[powers(v,t) for v in l]
    else:fac=[pairpow(v,t) for v in l]
    return ntt([(u*A+v*C)%P for u,v,(A,C) in zip(X,B,fac)],True)

def hashvec(v): return sha256(b''.join(x.to_bytes(8,'big') for x in v)).hexdigest()

def write_csr(n,path):
    # Full CSR representation of the exact sparse matrix, not a dense N^2 block.
    rows=[0];col=[];data=[]
    for i in range(n):
        items=sorted(((i,1),((i-1)%n,1),((i-2)%n,P-1)))
        for j,c in items:col.append(j);data.append(c)
        rows.append(len(col))
    assert len(col)==3*n
    payload=bytearray(b'RAINBON-CSR-V1')
    payload.extend(pack('<III',n,P,3*n))
    for ar in (rows,col,data):
        for v in ar:payload.extend(pack('<I',v))
    path.write_bytes(payload)
    return sha256(payload).hexdigest(),len(payload)

def check_csr(n,path):
    b=path.read_bytes();off=len(b'RAINBON-CSR-V1')
    nn,p,nnz=unpack_from('<III',b,off);off+=12
    assert(nn,p,nnz)==(n,P,3*n)
    nrows=n+1
    offs=[unpack_from('<I',b,off+4*i)[0] for i in range(nrows)];off+=4*nrows
    cols=[unpack_from('<I',b,off+4*i)[0] for i in range(nnz)];off+=4*nnz
    vals=[unpack_from('<I',b,off+4*i)[0] for i in range(nnz)];off+=4*nnz
    assert off==len(b) and offs[-1]==nnz
    for i in range(n):
        assert offs[i+1]-offs[i]==3
        assert [(cols[j],vals[j]) for j in range(offs[i],offs[i+1])]==sorted(((i,1),((i-1)%n,1),((i-2)%n,P-1)))
    x=start(n);out=[sum(x[cols[j]]*vals[j] for j in range(offs[i],offs[i+1]))%P for i in range(n)]
    assert out==matvec(x)

def sampled(v):
    n=len(v);addr=[i*(VIRTUAL-1)//242 for i in range(243)]
    assert len(set(addr))==243
    values=[v[(a*n)//VIRTUAL] for a in addr]
    return hashvec(values),sum(values)%P,len(set((a*n)//VIRTUAL for a in addr))

def run(outdir):
    outdir.mkdir(parents=True,exist_ok=True)
    rows=['AFFINE_MEGA_HEADER|prime=1299079|primitive_root=3|matrix_stencil=I+R-R2|dtype=u64_or_mod_p|virtual_meta_pixel_addresses=12000000000|json=0']
    for n in (6561,19683,59049):
        path=outdir.joinpath(f'RAINBON-META-{n}-CSR.bin')
        dig,bs=write_csr(n,path);check_csr(n,path)
        rows.append(f'MATRIX|d={n}|dense_entries={n*n}|sparse_nnz={3*n}|stored_csr_bytes={bs}|sha256={dig}|csr_all_rows=PASS|json=0')
        x=start(n);f=source(n,27)
        assert ntt(ntt(x),True)==x
        assert ntt(matvec(x))==[(v*w)%P for v,w in zip(ntt(x),eigen(n))]
        for t in ((27,729) if n==6561 else (27,)):
            z=jump(x,f,t)
            assert z==direct(x,f,t)
            assert z==jump(x,f,t,'pair')
            rows.append(f'DIRECT|d={n}|t={t}|stepped_equals_NTT=PASS|independent_affine_pair=PASS|hash={hashvec(z)}|json=0')
        t3=3**40
        big=jump(x,f,t3)
        bsum=jump(x,f,2**40)
        if n==6561:assert big==jump(x,f,t3,'pair')
        rows.append(f'JUMP|d={n}|t=3^40|hash={hashvec(big)}|sample_sha256={sampled(big)[0]}|sample_unique={sampled(big)[2]}|json=0')
        rows.append(f'BINARY|d={n}|t=2^40|hash={hashvec(bsum)}|json=0')
        leaf=direct(big,f,3)
        assert leaf==jump(x,f,t3+3)
        rows.append(f'LEAF|d={n}|at=3^40|play_forward_three=PASS|hash={hashvec(leaf)}|json=0')
        if n in (6561,59049):
            giant=jump(x,f,3**80)
            rows.append(f'JUMP|d={n}|t=3^80|hash={hashvec(giant)}|sample_sha256={sampled(giant)[0]}|json=0')
        f2=source(n,243,shift=n//27)
        t2=2**40
        staged=jump(big,f2,t2)
        X=ntt(x);F1=ntt(f);F2=ntt(f2);ll=eigen(n)
        params=[(powers(a,t3),powers(a,t2)) for a in ll]
        spec=[(pp2*pp1*u + pp2*ss1*v + ss2*w)%P
              for u,v,w,((pp1,ss1),(pp2,ss2)) in zip(X,F1,F2,params)]
        independently=ntt(spec,True)
        assert independently==staged
        rows.append(f'COMPOSE|d={n}|t1=3^40|t2=2^40|force_switch=27_to_243|independent_spectral_compose=PASS|hash={hashvec(staged)}|json=0')
        gate=jump(x,[0]*n,t3)
        assert gate==jump(x,source(n,0),t3)
        assert big!=gate
        rows.append(f'SOURCE_GATE|d={n}|neutron_like_force=27|hot_cutoff_force=0|forced_and_cutoff_differ=PASS|cutoff_equals_vacuum=PASS|json=0')
        for slits in (1,2,3):
            q=jump([0]*n,source(n,27,slits=slits),t3)
            h,s,un=sampled(q)
            rows.append(f'SLIT_READOUT|d={n}|nslits={slits}|time=3^40|sample_count=243|sample_unique={un}|sum_mod={s}|sha256={h}|model=MODULAR_FIELD_NOT_QUANTUM|json=0')
        nonlinear=x[:]
        for _ in range(27):
            nonlinear=matvec(nonlinear)
            nonlinear[n//2]=(nonlinear[n//2]+(7 if nonlinear[n//2]&1 else 11))%P
        false_affine=jump(x,source(n,7,slits=1),27)
        assert nonlinear!=false_affine
        rows.append(f'NEGATIVE|d={n}|state_dependent_local_kicks=27|constant_source_jump_agrees=NO|correctly_refused=PASS|json=0')
        rows.append(f'DIM_PASS|d={n}|csr=PASS|ntt=PASS|affine=PASS|composition=PASS|json=0')
    rows.append('AFFINE_MEGA_PASS|ring_states=6561,19683,59049|dense_matrix_allocation=NO|csr_matrix_materialized=YES|affine_forcing_exact=YES|nonlinear_force_refused=YES|r9=HOLD|json=0')
    receipt=outdir.joinpath('AFFINE-MEGA.hbp');receipt.write_text('\n'.join(rows)+'\n',encoding='ascii')
    outdir.joinpath('SHA256SUMS').write_text(''.join(f'{sha256(p.read_bytes()).hexdigest()}  {p.name}\n' for p in sorted(outdir.iterdir()) if p.is_file() and (p.name=='AFFINE-MEGA.hbp' or p.name.startswith('RAINBON-META-'))),encoding='ascii')
    print('\n'.join(rows))

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,default=Path('build/affine-mega'))
    args=ap.parse_args()
    run(args.out)
