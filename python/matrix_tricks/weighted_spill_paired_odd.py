#!/usr/bin/env python3
"""Independent exact W-balanced undirected edge and paired moment law.

The synthetic 72-move weighted quotient includes a flat P state and
nontrivial W. It compares odd Rayleigh edge forms against complete-Q
matrix actions and even paired Gram at orders 0..7. It rejects
planted dropped cross terms, missing self-loops and flat reentry.
The final odd moment uses no final outward operator action; the
toy includes a genuine exterior spill and an independent raw-move oracle.
No actual Yang-Mills 3x4 moments are computed here.
"""
import hashlib
from collections import Counter

W=(1,2,3,6,6)
M=(0,1,1,2,2)
T=((62,4,6,0,0),
   (2,61,3,6,0),
   (2,2,64,2,2),
   (0,2,1,66,3),
   (0,0,1,3,68))
Q=(1,2,3,4)
N=8
# Fixed raw source list: 72 labeled moves for each 5-state row.
# This toy oracle is frozen independently of the candidate matrix
# supplied to the runtime. The real production verifier rebuilds
# the raw moves from the actual periodic-S3 gauge/link rule.
RAW_MOVES=tuple(tuple(dst for dst,count in enumerate(row)
                      for _ in range(count)) for row in T)
R=((2,0,0,0,1,0,0,1),
   (0,3,1,0,0,1,0,0),
   (0,0,0,4,0,1,1,0),
   (0,0,0,0,0,0,0,0))


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


def audit_source(candidate):
    """Independent frozen 72-move oracle, not only row sum/W balance."""
    require(len(candidate)==len(T) and
            all(len(row)==len(T) for row in candidate),
            "candidate_matrix_shape")
    for i,row in enumerate(candidate):
        require(sum(row)==72,"candidate_row_sum_72")
        for j,weight in enumerate(row):
            require(W[i]*weight==W[j]*candidate[j][i],
                    "candidate_W_balance")
        raw=Counter(RAW_MOVES[i])
        proposed=Counter({j:n for j,n in enumerate(row) if n})
        require(raw==proposed,"raw_generator_disagrees")


def edge_odd(V,z=5,source_table=T):
    """True on-support, W-balanced edge form, with raw-source refusal."""
    audit_source(source_table)
    support=tuple(i for i,row in enumerate(V) if any(v!=0 for v in row))
    out=[[0]*N for _ in range(N)]
    for x in support:
        i=Q[x]
        w=W[i]
        diag=72+9*M[i]-z-source_table[i][i]
        for a in range(N):
            for b in range(N):
                out[a][b]+=w*diag*V[x][a]*V[x][b]
    for pos,x in enumerate(support):
        i=Q[x]
        for y in support[pos+1:]:
            j=Q[y]
            if not source_table[i][j]:
                continue
            balance=W[i]*source_table[i][j]
            require(balance==W[j]*source_table[j][i],
                    "wrong_W_balance")
            for a in range(N):
                for b in range(N):
                    out[a][b]-=balance*(V[x][a]*V[y][b]+
                                        V[y][a]*V[x][b])
    return tuple(tuple(row) for row in out)


def refuse_bad_table(candidate):
    """The gate must actively throw, not just yield a changed answer."""
    try:
        audit_source(candidate)
    except ValueError:
        return True
    return False


def clone_table():
    return [list(row) for row in T]

def hashm(Mm):
    return hashlib.sha256("|".join(str(v) for row in Mm for v in row)
                          .encode("ascii")).hexdigest()


def signed_bias_lane_trial():
    """Proof-sized offset lanes for mixed-sign higher Krylov vectors.

    In production the source bound and incoming degree must be proven
    from all 72 raw moves. Never directly unsigned-pack negative V_j.
    """
    v0=(-13,15,-2,0,7,-5,11,-1)
    v1=(9,-4,6,-12,0,3,-1,12)
    offset=16
    degree=72
    largest=15
    bound=degree*(offset+largest)
    width=bound.bit_length()
    require(bound < (1<<width),"safe_width_not_large_enough")
    mask=(1<<width)-1

    def encode(vec,bits,bias):
        require(all(-bias<=v<=largest for v in vec),
                "signed_source_bound")
        require(degree*(bias+largest) < (1<<bits),
                "signed_lane_carry_bound")
        return sum((v+bias)<<(bits*k) for k,v in enumerate(vec))

    packed0=encode(v0,width,offset)
    packed1=encode(v1,width,offset)
    total=2*packed0+3*packed1
    actual=tuple(((total>>(width*k))&mask)-5*offset
                 for k in range(N))
    expected=tuple(2*v0[k]+3*v1[k] for k in range(N))
    require(actual==expected,"independent_signed_lane_unpack")
    require(all(v!=0 for v in expected),"nontrivial_signed_lane_test")

    refused_narrow=False
    try:
        encode(v0,7,offset)
    except ValueError:
        refused_narrow=True
    require(refused_narrow,"narrow_signed_width_must_refuse")
    bad=tuple(((total>>(width*k))&mask)-5*(offset+1)
              for k in range(N))
    require(bad!=expected,"wrong_signed_offset_not_detected")
    emit("WEIGHTED_SIGNED_LANES",signed_sources=2,
         weighted_incoming_degree=5,
         bound_incoming_degree=degree,
         tested_eight_columns=N,
         exact_signed_bias_roundtrip="PASS",
         unsigned_negative_values="FORBIDDEN",
         insufficient_lane_width="REFUSED",
         incorrect_bias_subtraction="REJECTED")


def run():
    audit_source(T)
    signed_bias_lane_trial()
    A=operator()
    require(len(R)==4,"Q_source_census")
    V=R
    moments=[]
    production_outward_passes=0

    # j=0 support is exactly three of four Q states. V1 has
    # nonzero entries in that missing state: this is a real spill.
    support0=tuple(i for i,row in enumerate(V) if any(row))
    require(len(support0)==3 and 3 not in support0,
            "must_start_with_missing_exterior_Q_state")
    raw_next=matrix_apply(A,V)  # independent probe; NOT retained by odd form
    require(any(raw_next[3]) and not any(V[3]),
            "nontrivial_next_shell_spill")
    require(edge_odd(V)==gram(V,raw_next),
            "support_only_odd_not_equal_to_direct_spill")

    for j in range(4):
        even=gram(V,V)
        odd=edge_odd(V)  # no V_{j+1} constructed here
        require(even==gram(R,repeated(A,R,2*j)),
                "even_direct_moment_disagrees")
        require(odd==gram(R,repeated(A,R,2*j+1)),
                "odd_direct_moment_disagrees")
        moments.extend((even,odd))
        emit("WEIGHTED_MOMENT_POINT",degree=j,
             even_order=2*j,odd_order=2*j+1,
             even_sha256=hashm(even),odd_sha256=hashm(odd),
             independent_undirected_edge="PASS",
             fullQ_pairing="PASS",
             outward_vector_passes_so_far=production_outward_passes)
        if j<3:
            V=matrix_apply(A,V)
            production_outward_passes+=1
    require(production_outward_passes==3 and len(moments)==8,
            "eight_moments_need_only_three_outward_passes")

    # Planted fault A: omit a self-loop; this fails the RAW oracle.
    lost_self=clone_table()
    lost_self[1][1]-=1
    require(refuse_bad_table(lost_self),
            "self_loop_drop_not_refused")
    # Planted fault B: omit a cross edge.
    lost_cross=clone_table()
    lost_cross[1][3]-=1
    require(refuse_bad_table(lost_cross),
            "cross_edge_drop_not_refused")
    # Planted fault C: remove an edge in both directions, restore
    # their full 72 counts with fabricated self-loops. Row sums
    # and detailed balance still PASS: only the independent raw
    # generator rejects the self-consistent counterfeit operator.
    folded=clone_table()
    edge13=folded[1][3]
    edge31=folded[3][1]
    require(edge13>0 and edge31>0 and
            W[1]*edge13==W[3]*edge31,"test_weighted_edge_exists")
    folded[1][3]=0
    folded[3][1]=0
    folded[1][1]+=edge13
    folded[3][3]+=edge31
    require(all(sum(row)==72 for row in folded),
            "folded_row_sum_must_still_pass")
    require(all(W[i]*folded[i][j]==W[j]*folded[j][i]
                for i in range(len(T)) for j in range(len(T))),
            "folded_detailed_balance_must_still_pass")
    require(refuse_bad_table(folded),
            "self_consistent_edge_fold_not_refused")
    # Planted fault D: alter the weighted metric, keeping adjacency.
    wrong_w=(1,2,6,6,6)
    require(any(wrong_w[i]*T[i][j]!=wrong_w[j]*T[j][i]
                for i in range(len(T)) for j in range(len(T))),
            "wrong_metric_not_detected")

    # A full (P + Q) propagation can leave Q and return through P.
    H=tuple(tuple((72+9*M[i]-5 if i==j else 0)-T[i][j]
                  for j in range(5)) for i in range(5))
    raw2=tuple(tuple(sum(H[i][k]*H[k][j] for k in range(5))
                     for j in Q) for i in Q)
    A2=tuple(tuple(sum(A[i][k]*A[k][j] for k in range(4))
                   for j in range(4)) for i in range(4))
    require(raw2!=A2,"unprojected_P_reentry_not_detected")

    emit("WEIGHTED_MOMENT_VERDICT",operator_rows=5,
         P_states=1,Q_states=4,source_columns=N,
         raw_moves_per_full_row=72,
         production_outward_passes=3,exact_moments=8,
         saved_final_full_vector_action=1,
         signed_eight_lane_offset_control="PASS",
         actual_outer_Q_spill_witness="PASS",
         on_support_undirected_odd_form="PASS",
         raw_oracle_self_loop_drop="REFUSED",
         raw_oracle_cross_edge_drop="REFUSED",
         balanced_72_row_forged_self_loops="REFUSED",
         wrong_W="REJECTED",unprojected_P_return="REJECTED",
         finite_3x4_actual_source="NOT_COMPUTED",
         Wilson_weight="HOLD")

def repeated(A,V,n):
    out=V
    for _ in range(n):
        out=matrix_apply(A,out)
    return out


if __name__=="__main__":
    run()
