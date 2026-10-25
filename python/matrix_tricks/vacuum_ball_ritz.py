#!/usr/bin/env python3
"""Integer-only 3x3 S3 vacuum-ball Ritz measurements, not full-gap proofs.

Reuses the production S3 lattice model but builds only the principal
radius-k submatrix about ALL eight flat vacua. A block Krylov basis of
8*(order+1) vectors yields rigorous Ritz eigenvalue brackets via fraction-
free inertia. The difference of two upper bounds is NOT a full gap bound.
No floats, Decimal, Fraction, lossy rescaling, or full 10M-state allocation.
"""
import argparse
import sys
from collections import deque
from hashlib import sha256
from pathlib import Path
from time import monotonic_ns

sys.path.insert(0,str(Path(__file__).resolve().parents[1].joinpath("volume")))
from flat_law import Lat, MUL, TRANS

L=3
lat=Lat(L,L)
reps={}
for a in range(6):
    for b in range(6):
        if MUL[a][b]!=MUL[b][a]:
            continue
        u=[0]*(2*lat.n)
        for y in range(L):
            u[lat.site(L-1,y)]=a
        for x in range(L):
            u[lat.n+lat.site(x,L-1)]=b
        assert lat.magnetic(u)==0
        q=lat.gauge_fix(u)
        reps.setdefault(lat.key(q),q)
FLATS=tuple(sorted(reps))
assert len(FLATS)==8

row_cache={}
def electric_row(key):
    if key in row_cache:
        return row_cache[key]
    q=reps[key]
    u=lat.rep(q)
    values={key:3*2*lat.n}
    for e in range(2*lat.n):
        for t in TRANS:
            v=u[:]
            v[e]=MUL[t][v[e]]
            qj=lat.gauge_fix(v)
            j=lat.key(qj)
            reps.setdefault(j,qj)
            values[j]=values.get(j,0)-1
    values={k:v for k,v in values.items() if v!=0}
    assert sum(values.values())==0
    row_cache[key]=values
    return values

def build_ball(radius):
    seen=set(FLATS)
    layer=list(FLATS)
    counts=[len(layer)]
    for _ in range(radius):
        next_layer=[]
        for key in layer:
            for j in electric_row(key):
                if j not in seen:
                    seen.add(j)
                    next_layer.append(j)
        layer=next_layer
        counts.append(len(layer))
    keys=sorted(seen)
    weights={key:lat.weight(reps[key]) for key in keys}
    mag={key:lat.magnetic(lat.rep(reps[key])) for key in keys}
    rows={key:{j:v for j,v in electric_row(key).items() if j in seen} for key in keys}
    assert all(weights[i]*v==weights[j]*rows[j].get(i,0)
               for i,row in rows.items() for j,v in row.items())
    assert all(rows[i].get(i,0)==electric_row(i).get(i,0) for i in keys)
    assert all(mag[i]==0 for i in FLATS)
    return keys,weights,mag,rows,counts

def apply(v,mag,rows,r):
    return {i:y for i,row in rows.items()
            if (y:=sum(z*v.get(j,0) for j,z in row.items())+r*mag[i]*v.get(i,0))!=0}

def gram(v,u,weights):
    if len(v)>len(u):
        v,u=u,v
    return sum(weights[i]*x*u.get(i,0) for i,x in v.items())

def project(weights,mag,rows,r,order):
    basis=[{a:1} for a in FLATS]
    last=basis
    for _ in range(order):
        last=[apply(v,mag,rows,r) for v in last]
        basis.extend(last)
    hcols=[apply(v,mag,rows,r) for v in basis]
    a=[[gram(u,v,weights) for v in basis] for u in basis]
    b=[[gram(u,hv,weights) for hv in hcols] for u in basis]
    n=len(a)
    assert all(a[i][j]==a[j][i] and b[i][j]==b[j][i]
               for i in range(n) for j in range(n))
    return a,b

def inertia(a):
    """Exact Bareiss leading-minor sequence; fail closed on any zero pivot."""
    m=[row[:] for row in a]
    n=len(m)
    prev=1
    signs=[1]
    for k in range(n):
        pivot=m[k][k]
        if pivot==0:
            raise ValueError("zero pivot: uncertified Ritz threshold")
        signs.append(1 if pivot>0 else -1)
        for i in range(k+1,n):
            x=m[i][k]
            for j in range(i,n):
                top=pivot*m[i][j]-x*m[k][j]
                q,rem=divmod(top,prev)
                if rem:
                    raise ValueError("non-exact Bareiss division")
                m[i][j]=q
                m[j][i]=q
            m[i][k]=m[k][i]=0
        prev=pivot
    return sum(s!=t for s,t in zip(signs,signs[1:]))

def bracket(a,b,k,scale=1000000,upper=200):
    n=len(a)
    def count(t):
        c=[[scale*b[i][j]-t*a[i][j] for j in range(n)] for i in range(n)]
        return inertia(c)
    lo,hi=0,upper*scale
    assert count(lo)<=k and count(hi)>k
    while hi-lo>1:
        middle=(lo+hi)//2
        if count(middle)<=k:
            lo=middle
        else:
            hi=middle
    return lo,hi

def measurement(radius,r,order):
    start=monotonic_ns()
    keys,w,mag,rows,layers=build_ball(radius)
    t=f"VACUUM_BALL|L=3x3|radius={radius}|JE=1|JB={r}|order={order}"
    raw=bytearray(b"HEBC-3X3-COMB-VACUUM-BALL-V1\0")
    for i in keys:
        raw.extend(i.to_bytes(8,"little"))
        raw.extend(w[i].to_bytes(1,"little"))
        raw.extend(mag[i].to_bytes(2,"little"))
        for j,v in sorted(rows[i].items()):
            raw.extend(j.to_bytes(8,"little"))
            raw.extend(v.to_bytes(2,"little",signed=True))
    print(t+f"|states={len(keys)}|layers={','.join(map(str,layers))}|nnz={sum(map(len,rows.values()))}|matrix_sha={sha256(raw).hexdigest()}|basis=COMB_TREE|json=0",flush=True)
    a,b=project(w,mag,rows,r,order)
    if inertia(a)!=0:
        raise ValueError("nonpositive/degenerate Krylov Gram matrix")
    lo0,hi0=bracket(a,b,0)
    lo1,hi1=bracket(a,b,1)
    print(t+f"|ritz_dim={len(a)}|E0_in=[{lo0},{hi0}]/1000000|E1_in=[{lo1},{hi1}]/1000000|gap_in=[{lo1-hi0},{hi1-lo0}]/1000000|elapsed_ns={monotonic_ns()-start}|kind=RITZ_MEASUREMENT_NOT_FULL_GAP_CERTIFICATE|json=0",flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("radius",type=int)
    p.add_argument("r",type=int)
    p.add_argument("--order",type=int,default=1)
    a=p.parse_args()
    if not (1<=a.radius<=2 and a.r>0 and 1<=a.order<=a.radius):
        p.error("radius 1..2; r>0; order in 1..radius")
    measurement(a.radius,a.r,a.order)
