#!/usr/bin/env python3
"""Exact frozen-trial parameter replay for symmetric operator pencils.

No floats, rounded eigensolvers, dynamic imports, or dependence on the
large finite-S3 data set. A checked small integer model demonstrates the law.
"""
import hashlib

Q0 = ((9,-1,0),(-1,10,-2),(0,-2,12))
B = ((2,1),(1,2),(0,0))
F0 = ((5,0),(0,7))
Y = ((1,0),(0,1),(0,0))
Z0 = 1
GRID = (0,1,2,3,4,5)
GAMMA0 = 7
BOUNDARY_DIAG = 14


def require(test, why):
    if not test:
        raise ValueError(why)


def fixed_bytes():
    return ("Q=" + str(Q0) + "\nB=" + str(B) + "\nF=" + str(F0) +
            "\nY=" + str(Y) + "\nz0=" + str(Z0) + "\n").encode("ascii")


def source_sha():
    return hashlib.sha256(fixed_bytes()).hexdigest()


def transpose(a):
    return tuple(zip(*a))


def mul(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def add(a,b):
    return tuple(tuple(x+y for x,y in zip(aa,bb)) for aa,bb in zip(a,b))


def sub(a,b):
    return tuple(tuple(x-y for x,y in zip(aa,bb)) for aa,bb in zip(a,b))


def scale(a,n):
    return tuple(tuple(n*x for x in row) for row in a)


def eye(n):
    return tuple(tuple(int(i==j) for j in range(n)) for i in range(n))


def det2(m):
    return m[0][0]*m[1][1]-m[0][1]*m[1][0]


def psd2(m):
    return m[0][0]>=0 and m[1][1]>=0 and m[0][1]==m[1][0] and det2(m)>=0


def det3(a):
    return sum((1 if j%2==0 else -1)*a[0][j]*det2(
        tuple(tuple(a[i][k] for k in range(3) if k!=j)
              for i in (1,2))) for j in range(3))


def adjugate3(a):
    cof=tuple(tuple((1 if (i+j)%2==0 else -1)*det2(
        tuple(tuple(a[ii][jj] for jj in range(3) if jj!=j)
              for ii in range(3) if ii!=i)) for j in range(3)) for i in range(3))
    return transpose(cof)


def model_at(z):
    return sub(Q0,scale(eye(3),z)),sub(F0,scale(eye(2),z))


def source_form(z):
    A,F=model_at(z)
    AY=mul(A,Y)
    residual=sub(AY,B)
    U=sub(sub(add(F,mul(transpose(Y),AY)),mul(transpose(B),Y)),
          mul(transpose(Y),B))
    G=mul(transpose(residual),residual)
    return U,G,residual


def moments():
    U0,G0,R0=source_form(Z0)
    H=mul(transpose(Y),Y)
    J=add(mul(transpose(R0),Y),mul(transpose(Y),R0))
    C=add(eye(2),H)
    R4=R0[-1:]
    G4=mul(transpose(R4),R4)
    return U0,G0,C,H,J,G4


def replay(z,cache):
    U0,G0,C,H,J,G4=cache
    dz=z-Z0
    U=sub(U0,scale(C,dz))
    G=add(sub(G0,scale(J,dz)),scale(H,dz*dz))
    return U,G,G4


def exact_schur(z):
    A,F=model_at(z)
    d=det3(A)
    require(d>0,"operator_not_positive_definite")
    numerator=sub(scale(F,d),mul(mul(transpose(B),adjugate3(A)),B))
    return numerator,d


def bounds(z,U,G,G4):
    gamma=GAMMA0-z
    dh=BOUNDARY_DIAG-z
    require(gamma>0 and dh>0,"outside_fullQ_floor")
    # The true Schur form S is bounded between these scaled integer forms.
    Lnum=sub(scale(U,gamma),G)
    Unum=sub(scale(U,dh),G4)
    return Lnum,gamma,Unum,dh


def check_sandwich(z,U,G,G4):
    S,den=exact_schur(z)
    L,dl,Up,du=bounds(z,U,G,G4)
    # Prove every matrix residual PSD by exact integer cross-products.
    lowerdiff=sub(scale(S,dl),scale(L,den))
    upperdiff=sub(scale(Up,den),scale(S,du))
    require(psd2(lowerdiff),"wrong_fullQ_lower")
    require(psd2(upperdiff),"wrong_boundary_upper")
    require(psd2(sub(scale(U,du),Up)),"boundary_not_positive")
    return det2(S)


def verify_hash(expected,raw):
    require(hashlib.sha256(raw).hexdigest()==expected,"source_hash_mismatch")


def negative_controls(cache):
    pin=source_sha()
    verify_hash(pin,fixed_bytes())
    changed=fixed_bytes().replace(b"Y=((1, 0)",b"Y=((2, 0)")
    require(changed!=fixed_bytes(),"fault_not_planted")
    failed=0
    try:
        verify_hash(pin,changed)
    except ValueError:
        failed+=1
    wrong=list(cache)
    wrong[4]=scale(wrong[4],-1)
    if source_form(2)[1]!=replay(2,tuple(wrong))[1]:
        failed+=1
    try:
        U,G,G4=replay(7,cache)
        bounds(7,U,G,G4)
    except ValueError:
        failed+=1
    require(failed==3,"control_accepted_wrong_input")
    return failed


def emit(name,**fields):
    print(name+"|"+"|".join(k+"="+str(v) for k,v in fields.items())+"|json=0")


def main():
    cache=moments()
    require(cache[-1][1][1]>0,"missing_exterior_residual")
    faults=negative_controls(cache)
    emit("FROZEN_PARAM_FREEZE",source_sha256=source_sha(),
         model="symmetric_integer_3Q_2P",z0=Z0,queries=len(GRID))
    for z in GRID:
        U,G,G4=replay(z,cache)
        Ud,Gd,Rd=source_form(z)
        require(U==Ud and G==Gd,"polynomial_replay_disagrees_with_direct")
        require(mul(transpose(Rd[-1:]),Rd[-1:])==G4,
                "exterior_energy_independence_failed")
        det=check_sandwich(z,U,G,G4)
        emit("FROZEN_PARAM_POINT",z=z,fullQ_floor=GAMMA0-z,
             boundary_denom=BOUNDARY_DIAG-z,
             exact_S_numerator_determinant=det,
             polynomial_equals_direct=1,Loewner_sandwich=1)
    emit("FROZEN_PARAM_VERDICT",points=len(GRID),
         negative_controls_rejected=faults,
         polynomial_degree_U=1,polynomial_degree_G=2,
         exterior_residual_degree=0,whole_graph_reused=1,
         general_spectral_weight="NOT_PROVED")


if __name__=="__main__":
    main()
