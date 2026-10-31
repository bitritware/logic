#!/usr/bin/env python3
"""Classification of all 127 nonempty subsets of the seven-operation signature, subset by subset.
Method A (certificates): a subset is complete (no nullary constants) iff NOT is in it and (DIF, or SYN together with ANT) is in it;
  an incomplete subset lies in a proper subclone, certified by one of three invariants (zero-fixing, preserves {0,+1}, preserves rho).
Method B (Slupecki): a subset is complete iff it contains an onto, essentially binary operation and the unary closure of its operations
  (with the three constants added when they are nullary operations) is all 27 unary functions.
Writes classification_127.csv. Every item is an assertion. Pure Python, a few seconds."""
import itertools, csv
from itertools import product
D=(-1,0,1); sgn=lambda z:(z>0)-(z<0)
F={'SYN':lambda a,b:sgn(a+b),'DIF':lambda a,b:sgn(a-b),'OR':lambda a,b:sgn((a+b)**2),'XOR':lambda a,b:sgn((a-b)**2),
   'AND':lambda a,b:sgn(4*a*b),'NOT':lambda a:sgn(1-a),'ANT':lambda a:-a}
NAMES=['SYN','DIF','OR','XOR','AND','NOT','ANT']
ar=lambda n:2 if n in NAMES[:5] else 1
RHO=[(a,b) for a in D for b in D if a*b>=0]
def fixes0(n): return F[n](*([0]*ar(n)))==0
def pres_set(n,S=(0,1)): return all(F[n](*t) in S for t in product(S,repeat=ar(n)))
def pres_rho(n):
    for ts in product(RHO,repeat=ar(n)):
        cols=list(zip(*ts))
        if F[n](*cols[0])*F[n](*cols[1])<0: return False
    return True
def slupecki(n):
    if ar(n)!=2: return False
    f=F[n]; vals={f(a,b) for a in D for b in D}
    e1=any(f(a,b)!=f(a2,b) for a in D for a2 in D for b in D); e2=any(f(a,b)!=f(a,b2) for a in D for b in D for b2 in D)
    return len(vals)==3 and e1 and e2
def ucl(sub,consts):
    known={tuple(D)}
    if consts: known|={(c,c,c) for c in D}
    ch=True
    while ch:
        ch=False; cur=list(known)
        for n in sub:
            if ar(n)==1:
                for g in cur:
                    r=tuple(F[n](v) for v in g)
                    if r not in known: known.add(r); ch=True
            else:
                for g in cur:
                    for h in cur:
                        r=tuple(F[n](a,b) for a,b in zip(g,h))
                        if r not in known: known.add(r); ch=True
    return len(known)
rows=[]
for r in range(1,8):
    for sub in itertools.combinations(NAMES,r):
        S=set(sub)
        A=('NOT' in S) and ('DIF' in S or {'SYN','ANT'}<=S)
        inv=[]
        if all(fixes0(n) for n in sub): inv.append('zero-fixing')
        if all(pres_set(n) for n in sub): inv.append('preserves {0,+1}')
        if all(pres_rho(n) for n in sub): inv.append('preserves rho')
        sl=any(slupecki(n) for n in sub)
        B0=sl and ucl(sub,False)==27; B1=sl and ucl(sub,True)==27
        assert A==B0==(not inv),(sub,A,B0,inv)
        why1='' if B1 else ('no onto essentially binary operation' if not sl else 'unary closure has %d of 27 functions'%ucl(sub,True))
        rows.append((' '.join(sub),int(A),'; '.join(inv),int(B1),why1))
assert len(rows)==127 and sum(r[1] for r in rows)==40 and sum(r[3] for r in rows)==95
with open('classification_127.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['subset','complete_without_constants','certificate_if_incomplete','complete_with_constants','reason_if_incomplete_with_constants'])
    w.writerows(rows)
print('classification_127.csv written: 127 subsets, 40 complete without constants, 95 with; methods A and B agree on every subset')
print('ALL ASSERTIONS PASSED')
