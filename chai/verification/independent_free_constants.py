import itertools
D=(-1,0,1); sgn=lambda z:(z>0)-(z<0)
B={'SYN':lambda a,b:sgn(a+b),'DIF':lambda a,b:sgn(a-b),'OR':lambda a,b:sgn((a+b)**2),'XOR':lambda a,b:sgn((a-b)**2),'AND':lambda a,b:sgn(4*a*b)}
U={'NOT':lambda a:sgn(1-a),'ANT':lambda a:-a}
names=list(B)+list(U)
def unary_closure(S,consts):
    base={tuple(D)}
    if consts:
        for c in D: base.add((c,c,c))
    known=set(base); changed=True
    while changed:
        changed=False
        cur=list(known)
        for n in S:
            if n in U:
                for f in cur:
                    r=tuple(U[n](v) for v in f)
                    if r not in known: known.add(r); changed=True
            else:
                for f in cur:
                    for g in cur:
                        r=tuple(B[n](a,b) for a,b in zip(f,g))
                        if r not in known: known.add(r); changed=True
    return known
def complete(S,consts):
    if not any(n in B for n in S): return False
    return len(unary_closure(S,consts))==27
res={}
for consts in (False,True):
    cnt=0; mins=[]
    sets=[S for r in range(1,8) for S in itertools.combinations(names,r)]
    comp=[S for S in sets if complete(S,consts)]
    print("constants free" if consts else "no free constants","-> complete subsets:",len(comp))
    # minimal
    m=[S for S in comp if not any(set(T)<set(S) for T in comp)]
    print("  minimal:",m)

# per-subset comparison with the main method (classification_table.py writes classification_127.csv)
import csv
ref={tuple(r['subset'].split()):(int(r['complete_without_constants']),int(r['complete_with_constants'])) for r in csv.DictReader(open('classification_127.csv'))}
mism=0
for S in [S for r in range(1,8) for S in itertools.combinations(names,r)]:
    mine=(int(complete(S,False)),int(complete(S,True)))
    if mine!=ref[tuple(sorted(S,key=names.index))]: mism+=1; print('MISMATCH',S,mine)
assert mism==0
print('per-subset comparison with classification_127.csv: 127 subsets, mismatches 0')
