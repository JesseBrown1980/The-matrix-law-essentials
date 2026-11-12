#!/usr/bin/env python3
"""Independent CSR and short-time checker for affine megapixel matrix jumps.

Does not import the generator or any of its computational functions.
Checks all rows of each original binary and directly computes
27 forced steps using independently parsed CSR data.
"""
from pathlib import Path
from hashlib import sha256
from struct import unpack_from
import argparse

P=1299079
EXPECTED={
    6561:("ebe0ee2bc29bf93d8ea3c5968a58fe676c4d461bb001e20165089a28bbf7c59b",
          "5281e8e0a7c1d30c334b03d505b15ddb162e2816d260835c50ec353fdc6b039a"),
    19683:("95f12bb8addcbf88cedafb14237e3a6d3e3727b22ab0c9262ceb834c6838d8bc",
           "c5243acb2fd07fe5328f05ea712a8ddf5307815f630598e02f3c994ce758a835"),
    59049:("3297505eda8ac2efa603d22bc5fac0211b0676ddeec13fd140477ec0a2208474",
           "33c51951c33b968373661489be2d052a36b6f2062b2c13e061b4525ec90bd3c7"),
}
def checked_prime(p):
    assert p>3 and p%2
    d=3
    while d*d<=p:
        assert p%d!=0,("not prime",d)
        d+=2

def independent_csr(path,n,expected_hash):
    raw=path.read_bytes()
    assert sha256(raw).hexdigest()==expected_hash
    magic=b"RAINBON-CSR-V1"
    assert raw.startswith(magic)
    at=len(magic)
    got,p,nnz=unpack_from("<III",raw,at);at+=12
    assert (got,p,nnz)==(n,P,3*n)
    offsets=[unpack_from("<I",raw,at+4*i)[0] for i in range(n+1)]
    at+=4*(n+1)
    cols=[unpack_from("<I",raw,at+4*i)[0] for i in range(nnz)]
    at+=4*nnz
    vals=[unpack_from("<I",raw,at+4*i)[0] for i in range(nnz)]
    at+=4*nnz
    assert at==len(raw)
    assert offsets[0]==0 and offsets[-1]==3*n
    assert all(offsets[i+1]-offsets[i]==3 for i in range(n))
    assert all(0<=v<P for v in vals)
    assert all(0<=j<n for j in cols)
    for i in range(n):
        sl=slice(offsets[i],offsets[i+1])
        entries=list(zip(cols[sl],vals[sl]))
        assert entries==sorted(((i,1),((i-1)%n,1),((i-2)%n,P-1)))
    return offsets,cols,vals

def independent_direct(n,offs,cols,vals):
    state=[(7*i*i+3*i+1)%P for i in range(n)]
    f=[0]*n
    for j in (n//2-n//9,n//2,n//2+n//9):
        f[j]=(f[j]+27)%P
    for _ in range(27):
        nxt=[0]*n
        for i in range(n):
            acc=f[i]
            for j in range(offs[i],offs[i+1]):
                acc+=vals[j]*state[cols[j]]
            nxt[i]=acc%P
        state=nxt
    digest=sha256(b"".join(v.to_bytes(8,"big") for v in state)).hexdigest()
    return digest

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("resultdir",type=Path)
    a=ap.parse_args()
    checked_prime(P)
    assert pow(3,P-1,P)==1
    assert all(pow(3,(P-1)//q,P)!=1 for q in (2,3,11))
    receipt=a.resultdir.joinpath("AFFINE-MEGA.hbp").read_text()
    assert "AFFINE_MEGA_PASS|ring_states=6561,19683,59049|" in receipt
    for n,(matrix_hash,direct_hash) in EXPECTED.items():
        path=a.resultdir.joinpath("RAINBON-META-%d-CSR.bin"%n)
        offs,cols,vals=independent_csr(path,n,matrix_hash)
        got=independent_direct(n,offs,cols,vals)
        assert got==direct_hash,("direct 27 CSR replay mismatch",n,got)
        assert ("DIRECT|d=%d|t=27|stepped_equals_NTT=PASS|independent_affine_pair=PASS|hash=%s|"%(n,got)) in receipt
        print("INDEPENDENT_CSR|d=%d|rows=%d|nonzeros=%d|checksum=%s|all_rows=PASS|direct27_hash=%s|json=0"%
              (n,n,3*n,matrix_hash,got))
    needed={
        59049:(
          "3fedf4d67518bee4742cbe49995e3c7248843647b7f5ca6a9f5c42cf8da10395",
          "11d5ea99a88f2001dfc3ea0c1cf743759ef461cf4e26fccf45ff0b9365741ac1",
          "94a7470283abd6a6fec516b7a7b4b97f5edd266757f6fcf22f69d8f3334fe7e0"),
    }
    for n,(h40,h80,h_switch) in needed.items():
        assert ("JUMP|d=%d|t=3^40|hash=%s|"%(n,h40)) in receipt
        assert ("JUMP|d=%d|t=3^80|hash=%s|"%(n,h80)) in receipt
        assert ("COMPOSE|d=%d|t1=3^40|t2=2^40|force_switch=27_to_243|independent_spectral_compose=PASS|hash=%s|"%(n,h_switch)) in receipt
    assert sha256(a.resultdir.joinpath("AFFINE-MEGA.hbp").read_bytes()).hexdigest()=="43440bb7d611ed26c0598b0091d70bfbc75d0339b9d31d9630d30f233112fd9a"
    print("AFFINE_MEGA_INDEPENDENT_AUDIT_PASS|matrices=3|total_rows=85293|direct_short_run=PASS|jump_receipt_hash=PASS|json=0")

if __name__=="__main__":
    main()
