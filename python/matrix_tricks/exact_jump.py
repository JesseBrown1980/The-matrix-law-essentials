#!/usr/bin/env python3
"""Exact integer repeated-squaring for small fixed matrices.

For integer S, integer v and nonnegative N, jump(S,v,N) = S**N @ v.
No floating point, normalization, amplitude right shifts, or approximate
division. This guarantees algebraic equality, not speed at unlimited size.
Pilot is the 8x8 third-order effective matrix of 3x3 S3, not the full lattice.
"""
from hashlib import sha256
from time import monotonic_ns
import argparse

# Pinned independently reconstructed 3x3 finite S3 effective matrix 16H3,
# in the BFS tree basis. The SHA guard prevents accidental matrix replacement.
H316 = (
    (0,-54,0,-54,0,0,0,0),
    (-18,0,-36,0,-18,0,0,0),
    (0,-54,0,0,0,0,0,0),
    (-18,0,0,0,-18,-36,0,0),
    (0,-18,0,-18,0,0,0,0),
    (0,0,0,-54,0,0,0,0),
    (0,0,0,0,0,0,0,0),
    (0,0,0,0,0,0,0,0),
)
SOURCE_SHA="69147c597f981eb9af15dc7f0b31848d52ef7190e9cfa4007de75ce2e7a9c4d0"
SHA64="c11971469a50925b416ce6edf780a45ee2d216241fbd7514b5cb5ab4877db8b6"
SHA4096="07bf219237a7100a67ac7d4d3b01472af5777cb30a90efc7d7269118ba479eaa"

def dimensions(a,v):
    n=len(a)
    if not 1<=n<=16 or len(v)!=n or any(len(row)!=n for row in a):
        raise ValueError("square matrix size 1..16 and matching vector required")
    if any(type(x) is not int for row in a for x in row):
        raise TypeError("integer matrix coefficients required")
    if any(type(x) is not int for x in v):
        raise TypeError("integer amplitudes required")
    return n

def mv(a,v):
    dimensions(a,v)
    return [sum(x*y for x,y in zip(row,v)) for row in a]

def mm(a,b):
    n=dimensions(a,[0]*len(a))
    if len(b)!=n or any(len(row)!=n for row in b):
        raise ValueError("matrix size mismatch")
    return [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)]
            for i in range(n)]

def exact_jump(a,v,steps):
    """Binary exponentiation: exact integer S**steps v, including steps=0."""
    dimensions(a,v)
    if type(steps) is not int or steps<0:
        raise ValueError("steps must be a nonnegative integer")
    p=[list(row) for row in a]
    result=list(v)
    squarings=0
    while steps:
        if steps & 1:
            result=mv(p,result)
        steps >>= 1  # shifts STEP COUNTER, not integer state amplitudes
        if steps:
            p=mm(p,p)
            squarings+=1
    return result,squarings

def sequential(a,v,steps):
    if type(steps) is not int or steps<0:
        raise ValueError("invalid steps")
    result=list(v)
    for _ in range(steps):
        result=mv(a,result)
    return result

def state_sha(v):
    raw=bytearray(b"HEBC-S3-3X3-H3-EXACT-JUMP-V1\0")
    for x in v:
        signed=x.to_bytes(max(1,(x.bit_length()+8)//8),"little",signed=True)
        raw.extend(len(signed).to_bytes(4,"little"))
        raw.extend(signed)
    return sha256(raw).hexdigest()

def source_sha(a):
    raw=bytearray(b"HEBC-3X3-THIRD-ORDER-INDEPENDENT-V1\0")
    for row in a:
        for x in row:
            raw.extend(x.to_bytes(8,"little",signed=True))
    return sha256(raw).hexdigest()

def self_test():
    assert source_sha(H316)==SOURCE_SHA
    s=[[128*(i==j)-H316[i][j] for j in range(8)] for i in range(8)]
    v=[1,0,0,0,0,0,0,0]
    assert exact_jump(s,v,0)[0]==v
    assert exact_jump(s,v,1)[0]==mv(s,v)
    t0=monotonic_ns()
    plain=sequential(s,v,64)
    short,n6=exact_jump(s,v,64)
    assert n6==6 and short==plain and state_sha(short)==SHA64
    print("VACUUM_JUMP|L=3x3|matrix=16H3|shift=128|steps=64|squarings=6|"
          "exact_equal_sequential=1|state_sha="+state_sha(short)+
          "|kind=EFFECTIVE_MATRIX_ONLY|json=0",flush=True)
    long,n12=exact_jump(s,v,4096)
    half,_=exact_jump(s,v,2048)
    double,_=exact_jump(s,half,2048)
    assert long==double and n12==12 and state_sha(long)==SHA4096
    print("VACUUM_JUMP|L=3x3|matrix=16H3|shift=128|steps=4096|squarings=12|"
          "state_bits="+str(max(x.bit_length() for x in long))+
          "|state_sha="+state_sha(long)+
          "|division=NONE|rounding=NONE|exact_composition=PASS|elapsed_ns="+
          str(monotonic_ns()-t0)+
          "|kind=EFFECTIVE_MATRIX_ONLY_NOT_FULL_DYNAMICS|json=0",flush=True)
    a=[[2,-3],[5,7]]
    b=[11,-13]
    for n in (0,1,2,3,9,17,65):
        assert exact_jump(a,b,n)[0]==sequential(a,b,n)
    for bad in (-1,"1.5"):
        try:
            exact_jump(s,v,bad)
            raise AssertionError("invalid step accepted")
        except ValueError:
            pass
    print("EXACT_JUMP|result=PASS|matrix_limit=16|integer_identity=S^N_v|"
          "sequential=64|composition=2048+2048|long_jump=4096|"
          "physics_claim=EFFECTIVE_MATRIX_ONLY|json=0",flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:
        self_test()
    else:
        p.print_help()
