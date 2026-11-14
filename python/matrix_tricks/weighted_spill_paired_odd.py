#!/usr/bin/env python3
"""Independent exact W-balanced undirected edge and paired moment law.

The synthetic 72-move weighted quotient includes a flat P state and
nontrivial W. It compares odd Rayleigh edge forms against complete-Q
matrix actions and even paired Gram at orders 0..7. It rejects
planted dropped cross terms, missing self-loops and flat reentry.
No actual Yang-Mills 3x4 moments are computed here.
"""
import hashlib

W=(1,2,3,6,6)
M=(0,1,1,2,2)
T=((62,4,6,0,0),
   (2,61,3,6,0),
   (2,2,64,2,2),
   (0,2,1,66,3),
   (0,0,1,3,68))
Q=(1,2,3,4)
N=8
R=((2,0,0,0,1,0,0,1),
   (0,3,1,0,0,1,0,0),
   (0,0,0,4,0,1,1,0),
   (1,0,2,0,0,0,3,1))


def require(c,msg):
    if not c:
        raise ValueError("WEIGHTED_MOMENT_REFUSE_"+msg)


def emit(name,**fields):
    print(name+"|"+"|".join(str(k)+"="+str(v) for k,v in fields.items())+
          "|json=0",flush=True)


def operator(z=5):
    return tuple(tuple((72+9*M[i]-z if i==j else 0)-T[i][j]
                       for j in Q) for i in Q)


def matrix_apply(A,V):
    return tuple(tuple(sum(A[i][j]*V[j][k] for j in range(len(V)))
                       for k in range(N)) for i in range(len(V)))


def gram(U,V):
    return tuple(tuple(sum(W[Q[x]]*U[x][i]*V[x][j]
                           for x in range(len(Q)))
                       for j in range(N)) for i in range(N))


def edge_odd(V,z=5,drop_self=False,drop_cross=False):
    """Independent undirected off-diagonal and diagonal self-loop form."""
    out=[[0]*N for _ in range(N)]
    for x,i in enumerate(Q):
        w=W[i]
        diag=72+9*M[i]-z
        for a in range(N):
            for b in range(N):
                out[a][b]+=w*diag*V[x][a]*V[x][b]
        if not drop_self:
            for a in range(N):
                for b in range(N):
                    out[a][b]-=w*T[i][i]*V[x][a]*V[x][b]
    if not drop_cross:
        for x,i in enumerate(Q):
            for y,j in enumerate(Q):
                if y<=x or not T[i][j]:
                    continue
                balance=W[i]*T[i][j]
                require(balance==W[j]*T[j][i],"wrong_W_balance")
                for a in range(N):
                    for b in range(N):
                        out[a][b]-=balance*(V[x][a]*V[y][b]+
                                            V[y][a]*V[x][b])
    return tuple(tuple(row) for row in out)


def hashm(Mm):
    return hashlib.sha256("|".join(str(v) for row in Mm for v in row)
                          .encode("ascii")).hexdigest()


def run():
    for i,row in enumerate(T):
        require(sum(row)==72,"full_source_row_not_72")
        for j,v in enumerate(row):
            require(W[i]*v==W[j]*T[j][i],"source_weight_balance")
    A=operator()
    assert len(R)==4
    V=R
    moments=[]
    for j in range(4):
        Vnext=matrix_apply(A,V)
        even=gram(V,V)
        odd=gram(V,Vnext)
        sourceedge=edge_odd(V)
        require(odd==sourceedge,"odd_edge_not_complete")
        require(even==gram(R, repeated(A,R,2*j)),
                "even_direct_moment_disagrees")
        require(odd==gram(R, repeated(A,R,2*j+1)),
                "odd_direct_moment_disagrees")
        moments.extend((even,odd))
        emit("WEIGHTED_MOMENT_POINT",degree=j,
             even_order=2*j,odd_order=2*j+1,
             even_sha256=hashm(even),odd_sha256=hashm(odd),
             independent_undirected_edge="PASS",
             fullQ_pairing="PASS")
        V=Vnext

    require(edge_odd(R,drop_self=True)!=edge_odd(R),
            "omitted_self_loops_not_detected")
    require(edge_odd(R,drop_cross=True)!=edge_odd(R),
            "omitted_cross_edges_not_detected")
    # A full (P + Q) propagation can leave Q and return through P.
    H=tuple(tuple((72+9*M[i]-5 if i==j else 0)-T[i][j]
                  for j in range(5)) for i in range(5))
    raw2=tuple(tuple(sum(H[i][k]*H[k][j] for k in range(5))
                     for j in Q) for i in Q)
    A2=tuple(tuple(sum(A[i][k]*A[k][j] for k in range(4))
                   for j in range(4)) for i in range(4))
    require(raw2!=A2,"flat_reentry_truncation_not_detected")
    # Changing one W weight while retaining the source T breaks balance.
    bad=(1,2,6,6,6)
    require(any(bad[i]*T[i][j]!=bad[j]*T[j][i]
                for i in range(5) for j in range(5)),
            "wrong_metric_not_detected")
    require(len(moments)==8,"four_pass_moment_count")
    emit("WEIGHTED_MOMENT_VERDICT",operator_rows=5,
         P_states=1,Q_states=4,source_columns=N,
         raw_moves_per_full_row=72,
         exact_passes=4,exact_moments=8,
         fullQ_weighted_pairing="PASS",
         on_support_undirected_odd_form="PASS",
         missing_self_loop="REJECTED",missing_cross_edge="REJECTED",
         unprojected_P_return="REJECTED",wrong_W="REJECTED",
         finite_3x4_actual_source="NOT_COMPUTED",
         Wilson_weight="HOLD")


def repeated(A,V,n):
    out=V
    for _ in range(n):
        out=matrix_apply(A,out)
    return out


if __name__=="__main__":
    run()
