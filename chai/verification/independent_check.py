import itertools, json, re
D=(-1,0,1); sgn=lambda z:(z>0)-(z<0)
SYN=lambda a,b:sgn(a+b); DIF=lambda a,b:sgn(a-b); OR=lambda a,b:sgn((a+b)**2); XOR=lambda a,b:sgn((a-b)**2); AND=lambda a,b:sgn(4*a*b)
NOT=lambda a:sgn(1-a); ANT=lambda a:-a
# --- k_c explicit DIF terms of the paper
kc={-1:lambda u,v:DIF(u,DIF(v,DIF(v,u))),0:lambda u,v:DIF(v,DIF(v,DIF(u,v))),1:lambda u,v:DIF(DIF(u,v),DIF(DIF(u,v),v))}
for c,f in kc.items():
    for u in D:
        for v in D:
            assert f(u,v)==(u if (u!=0 and v==c*u) else 0),(c,u,v)
print("k_c terms OK")
# --- L3 strong connectives / implications / half adder terms in Table III of the paper
for a in D:
    for b in D:
        assert DIF(SYN(a,b),NOT(AND(a,b)))==max(-1,a+b-1)
        assert SYN(SYN(a,b),NOT(AND(a,b)))==min(1,a+b+1)
        assert SYN(NOT(a),SYN(b,AND(a,b)))==min(1,1+b-a)
        assert SYN(DIF(b,a),OR(a,b))==max(-a,b)
        assert SYN(a,DIF(b,SYN(a,b)))==(1 if a==b==1 else -1 if a==b==-1 else 0)
        s=a+b; s3=((s+1)%3)-1
        assert AND(SYN(a,b),DIF(XOR(a,b),AND(a,b)))==s3
print("Table II terms OK")
# --- independent unary min-size DP over the seven operators (tree size = #operator nodes)
ops_b={'SYN':SYN,'DIF':DIF,'OR':OR,'XOR':XOR,'AND':AND}; ops_u={'NOT':NOT,'ANT':ANT}
def minsizes(allowed,maxsize=12):
    best={tuple(D):0}; layers={0:{tuple(D)}}
    for n in range(1,maxsize+1):
        new=set()
        for u in allowed:
            if u in ops_u:
                for f in layers[n-1]:
                    r=tuple(ops_u[u](v) for v in f)
                    if r not in best: new.add(r)
        for b in allowed:
            if b in ops_b:
                for i in range(0,n):
                    j=n-1-i
                    if i in layers and j in layers:
                        for f in layers[i]:
                            for g in layers[j]:
                                r=tuple(ops_b[b](x,y) for x,y in zip(f,g))
                                if r not in best: new.add(r)
        for r in new: best[r]=n
        layers[n]=new
        if len(best)==27: break
    return best
b7=minsizes(list(ops_b)+list(ops_u))
from collections import Counter
print("all seven: unary size histogram",sorted(Counter(b7.values()).items()),"count",len(b7))
b2=minsizes(['DIF','NOT'],14)
print("{DIF,NOT}: max size",max(b2.values()),"count",len(b2))
# --- compare with author's catalogue terms: evaluate every listed term
cat=json.load(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'unary_catalogue.json')))
def ev(term,a):
    t=term.replace('DIF','DIF_').replace('SYN','SYN_').replace('NOT','NOT_').replace('ANT','ANT_').replace('OR','OR_').replace('XOR','XOR_').replace('AND','AND_')
    return eval(t,{'a':a,'DIF_':DIF,'SYN_':SYN,'NOT_':NOT,'ANT_':ANT,'OR_':OR,'XOR_':XOR,'AND_':AND})
bad=0
for vec,(s7,t7),(s2,t2) in cat:
    for term in (t7,t2):
        if tuple(ev(term,a) for a in D)!=tuple(vec): bad+=1
    if b7[tuple(vec)]!=s7: bad+=1
    if b2[tuple(vec)]!=s2: bad+=1
print("catalogue mismatches (value or minimal size):",bad)
