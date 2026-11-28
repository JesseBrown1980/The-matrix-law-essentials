from random import Random
from hashlib import sha256
from least_action_jump import compose


def quotient(a):
    reps, ids, offsets = [], [], []
    index = {}
    for row in a:
        base = min((v for v in row if v is not None), default=0)
        key = tuple(None if v is None else v-base for v in row)
        if key not in index:
            index[key] = len(reps)
            reps.append(list(key))
        ids.append(index[key]); offsets.append(base)
    return reps, ids, offsets


def compact_compose(a,b):
    reps, ids, offsets = quotient(a)
    n = len(b)
    out = []
    for row in reps:
        vals = []
        for j in range(n):
            vals.append(min((v+b[k][j] for k,v in enumerate(row) if v is not None and b[k][j] is not None), default=None))
        out.append(vals)
    return [[None if v is None else v+offsets[i] for v in out[r]] for i,r in enumerate(ids)]


def main():
    rng = Random(281126)
    digest=sha256(); checks=0; counts=[]
    for n in (3,9,27,81):
        for mode in ('structured','generic','unreachable'):
            for case in range(6):
                r = max(1,n//3)
                reps = [[None if rng.randrange(5)==0 else rng.randrange(-30,31) for _ in range(n)] for _ in range(r)]
                a = [[None if x is None else x+i-10 for x in reps[i%r]] for i in range(n)] if mode=='structured' else [[None if rng.randrange(5)==0 else rng.randrange(-30,31) for _ in range(n)] for _ in range(n)]
                if mode=='unreachable': a[0]=[None]*n; a[-1]=[None]*n
                b = [[None if rng.randrange(5)==0 else rng.randrange(-30,31) for _ in range(n)] for _ in range(n)]
                actual=compact_compose(a,b); expected=compose(a,b); assert actual==expected
                checks+=n*n; digest.update(repr(actual).encode('ascii'))
                if n==81 and mode=='structured': counts.append(len(quotient(a)[0]))
    a=[[0,1],[1,0]]; b=[[0,None],[None,0]]
    assert len(quotient(a)[0])==2 and compact_compose(a,b)==a
    print(f'PASS cases=72 entries={checks} structured_81_classes={counts} sha256={digest.hexdigest()}')
    print('dense_81 candidate slots baseline=531441 quotient=177147 excluding quotient construction and expansion')


if __name__=='__main__':main()
