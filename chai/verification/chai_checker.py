#!/usr/bin/env python3
"""Chai checker. Default run: pure Python + numpy, about 15 seconds.
   python3 chai_checker.py          -> results_default.json, results_default.md
   python3 chai_checker.py --full   -> results.json, results.md; also runs the five exhaustive closures (about 4 minutes)
   python3 chai_checker.py --full-consts [--only=SYN+OR]  -> closures with free constants (slow, ~10-30 min each; one set per process is fine)
Every item is an assertion; the script stops at the first failure."""
import itertools, json, random, sys, time
from itertools import product
import numpy as np
D=(-1,0,1); sgn=lambda z:(z>0)-(z<0); iv=lambda c:1 if c else 0
SYN=lambda a,b:sgn(a+b); DIF=lambda a,b:sgn(a-b); OR=lambda a,b:sgn((a+b)**2)
XOR=lambda a,b:sgn((a-b)**2); AND=lambda a,b:sgn(4*a*b); NOT=lambda a:sgn(1-a); ANT=lambda a:-a
B2={'SYN':SYN,'DIF':DIF,'OR':OR,'XOR':XOR,'AND':AND}; U1={'NOT':NOT,'ANT':ANT}
R={}; T0=time.time()
def rec(k,v): R[k]=v; print(f'[{time.time()-T0:6.1f}s] {k}: {v}',flush=True)

# 1. identities -------------------------------------------------------------
assert all((abs(a+b),abs(a-b)) in {(0,0),(0,2),(2,0),(1,1)} for a in D for b in D)
assert all(AND(a,b)==DIF(OR(a,b),XOR(a,b)) for a in D for b in D)
assert all(OR(a,b)==sgn((a+b)**2) and OR(a,b)==abs(SYN(a,b)) and XOR(a,b)==abs(DIF(a,b)) for a in D for b in D)
Phi=lambda a,b:(SYN(a,b),DIF(a,b))
assert all(Phi(*Phi(a,b))==(a,b) for a in D for b in D)
fixed=[p for p in product(D,D) if Phi(*p)==p]; assert fixed==[(0,0)]
assert all(DIF(a,b)==SYN(a,ANT(b)) and XOR(a,b)==OR(a,ANT(b)) and AND(a,ANT(b))==ANT(AND(a,b)) for a in D for b in D)
assert all((SYN(a,ANT(a)),DIF(a,ANT(a)),AND(a,ANT(a)))==(0,a,-a*a) for a in D)
assert sgn(2**2-1**2)==1 and DIF(sgn(4),sgn(1))==0          # m>=2 counterexample (2,1)
assert all(SYN(SYN(a,b),XOR(a,b))==max(a,b) and DIF(SYN(a,b),XOR(a,b))==min(a,b) for a in D for b in D)
assert [NOT(NOT(a)) for a in D]==[0,0,1] and all(max(a,-a)==a*a for a in D)
na={n:sum(f(f(a,b),c)!=f(a,f(b,c)) for a in D for b in D for c in D) for n,f in [('SYN',SYN),('DIF',DIF),('OR',OR),('XOR',XOR)]}
assert na=={'SYN':4,'DIF':16,'OR':8,'XOR':8}
B01=(0,1)
assert all(f(a,b) in B01 for f in (SYN,OR,XOR,AND) for a in B01 for b in B01) and all(NOT(a)==1-a for a in B01)
assert SYN(-1,1)==0 and DIF(0,1)==-1
rec('identities',{'parity':True,'polarization':True,'involution_fixed_points':fixed,'nonassoc_failures':na})

# 2. the explicit terms of the paper -----------------------------------------
for a in D:
    assert SYN(a,OR(a,a))==iv(a==1) and XOR(a,OR(a,a))==iv(a==-1)
    assert NOT(DIF(a,DIF(a,NOT(a))))==a*a
    for b in D:
        for s in (1,-1): assert AND(iv(a==1),iv(b==s))==iv(a==1)*iv(b==s)
        assert DIF(iv(a==1),AND(iv(a==1),OR(b,b)))==iv(a==1)*iv(b==0)
for c in D:
    for u in D:
        for v in D:
            k=u if (u!=0 and v==c*u) else 0
            t={-1:DIF(u,DIF(v,DIF(v,u))),0:DIF(v,DIF(v,DIF(u,v))),1:DIF(DIF(u,v),DIF(DIF(u,v),v))}[c]
            assert t==k
rho=[(a,b) for a in D for b in D if a*b>=0]; assert len(rho)==7
for f in (OR,XOR,AND):
    assert all(f(a,c)*f(b,d)>=0 for (a,b),(c,d) in product(rho,rho))
assert all(ANT(a)*ANT(b)>=0 and NOT(a)*NOT(b)>=0 for a,b in rho)
assert (SYN(1,0),SYN(0,-1))==(1,-1)
TT={'EQ':(lambda a,b:NOT(XOR(a,b)),lambda a,b:iv(a==b)),'LE':(lambda a,b:NOT(DIF(a,b)),lambda a,b:iv(a<=b)),
 'LT':(lambda a,b:NOT(NOT(DIF(b,a))),lambda a,b:iv(a<b)),'MAX':(lambda a,b:SYN(SYN(a,b),XOR(a,b)),max),
 'MIN':(lambda a,b:DIF(SYN(a,b),XOR(a,b)),min),'imp_K3':(lambda a,b:SYN(DIF(b,a),OR(a,b)),lambda a,b:max(-a,b)),
 'imp_L3':(lambda a,b:SYN(NOT(a),SYN(b,AND(a,b))),lambda a,b:min(1,1+b-a)),
 'and_L3':(lambda a,b:DIF(SYN(a,b),NOT(AND(a,b))),lambda a,b:max(-1,a+b-1)),
 'or_L3':(lambda a,b:SYN(SYN(a,b),NOT(AND(a,b))),lambda a,b:min(1,a+b+1)),
 'CARRY':(lambda a,b:SYN(a,DIF(b,SYN(a,b))),lambda a,b:1 if a+b>=2 else (-1 if a+b<=-2 else 0)),
 'SUM':(lambda a,b:AND(SYN(a,b),DIF(XOR(a,b),AND(a,b))),lambda a,b:((a+b+1)%3)-1)}
for n,(t,f) in TT.items(): assert all(t(a,b)==f(a,b) for a in D for b in D),n
rec('paper_terms',{'Lemma2':True,'k_c':True,'Theorem5_invariants':True,'TableII_terms':list(TT)})

# 3. five binary operations generate Z_n; terms only over DIF,NOT -------------
# term evaluation over {DIF,NOT} alone
d_zero=lambda x:DIF(x,x); d_ant=lambda x:DIF(0,x); d_syn=lambda a,b:DIF(a,d_ant(b))
d_mag=lambda x:NOT(DIF(x,DIF(x,NOT(x)))); d_or=lambda a,b:d_mag(d_syn(a,b)); d_xor=lambda a,b:d_mag(DIF(a,b))
d_and=lambda a,b:DIF(d_or(a,b),d_xor(a,b))
assert NOT(0)==1 and DIF(0,1)==-1
assert all(d_syn(a,b)==SYN(a,b) and d_or(a,b)==OR(a,b) and d_xor(a,b)==XOR(a,b) and d_and(a,b)==AND(a,b) for a in D for b in D)
dp1=lambda x:d_syn(x,d_or(x,x)); dm1=lambda x:d_xor(x,d_or(x,x))
hh={1:lambda u,v:d_and(dp1(u),dp1(v)),-1:lambda u,v:d_and(dp1(u),dm1(v)),0:lambda u,v:DIF(dp1(u),d_and(dp1(u),d_or(v,v)))}
def delta(p):
    j=next(i for i,c in enumerate(p) if c); d0=dp1 if p[j]==1 else dm1
    def f(x):
        u=d0(x[j])
        for i,c in enumerate(p):
            if i!=j: u=hh[c](u,x[i])
        return u
    return f
def build_any(tab,n):                      # arbitrary F over D^n from {DIF,NOT} terms (Theorem 3b,c)
    pts=list(product(D,repeat=n)); ds={p:delta(p) for p in pts if any(p)}
    def e0(x):
        r=NOT(d_or(x[0],x[0]))
        for xi in x[1:]: r=d_and(r,NOT(d_or(xi,xi)))
        return r
    c0=tab[(0,)*n]
    const=lambda: c0
    def F(x):
        acc=d_zero(x[0])
        for p in pts:
            if any(p) and tab[p]!=0:
                t=ds[p](x); acc=d_syn(acc,t if tab[p]==1 else d_ant(t))
        z=d_and(e0(x),c0)                 # constant c0 is a term: 0, +1=NOT(0), -1=DIF(0,+1)
        return d_syn(acc,z)
    return F
# all 3^9=19683 binary functions, all 27 unary, random n=3,4
cnt=0
for vals in product(D,repeat=9):
    tab=dict(zip(product(D,D),vals)); F=build_any(tab,2)
    assert all(F(x)==tab[x] for x in tab); cnt+=1
for vals in product(D,repeat=3):
    tab={(a,):v for (a,),v in zip(product(D),vals)}; tab=dict(zip([(a,) for a in D],vals)); F=build_any(tab,1)
    assert all(F(x)==tab[x] for x in tab); cnt+=1
random.seed(7)
for n in (3,4):
    pts=list(product(D,repeat=n))
    for _ in range(25):
        tab={p:random.choice(D) for p in pts}; F=build_any(tab,n); assert all(F(x)==tab[x] for x in pts)
rec('DIF_NOT_complete_constructive',{'functions_rebuilt_n<=2':cnt,'random_n3_n4':50})
# Z_n from the five operations (no NOT): all Z_2, random Z_3, Z_4
def build_Z(tab,n):
    ds={}
    pts=list(product(D,repeat=n))
    def dl(p):
        j=next(i for i,c in enumerate(p) if c); d0=(lambda x:SYN(x,OR(x,x))) if p[j]==1 else (lambda x:XOR(x,OR(x,x)))
        h={1:lambda u,v:AND(SYN(u,OR(u,u)),SYN(v,OR(v,v))),-1:lambda u,v:AND(SYN(u,OR(u,u)),XOR(v,OR(v,v))),
           0:lambda u,v:DIF(SYN(u,OR(u,u)),AND(SYN(u,OR(u,u)),OR(v,v)))}
        def f(x):
            u=d0(x[j])
            for i,c in enumerate(p):
                if i!=j: u=h[c](u,x[i]) if False else {1:lambda:AND(u,SYN(x[i],OR(x[i],x[i]))),-1:lambda:AND(u,XOR(x[i],OR(x[i],x[i]))),0:lambda:DIF(u,AND(u,OR(x[i],x[i])))}[c]()
            return u
        return f
    ds={p:dl(p) for p in pts if any(p)}
    def F(x):
        acc=DIF(x[0],x[0])
        for p in pts:
            if any(p) and tab[p]!=0:
                t=ds[p](x); acc=SYN(acc,t if tab[p]==1 else DIF(0,t))
        return acc
    return F
cz=0
nz=[p for p in product(D,D) if any(p)]
for vals in product(D,repeat=8):
    tab=dict(zip(nz,vals)); tab[(0,0)]=0; F=build_Z(tab,2)
    assert all(F(x)==tab[x] for x in tab); cz+=1
for n in (1,3,4):
    pts=list(product(D,repeat=n))
    reps=27 if n==1 else 25
    for _ in range(reps):
        tab={p:random.choice(D) for p in pts}; tab[(0,)*n]=0; F=build_Z(tab,n); assert all(F(x)==tab[x] for x in pts)
rec('Z_n_from_five_operations',{'all_Z2_rebuilt':cz,'unary_and_random_Z3_Z4':'ok'})

# 4. DIF alone: all odd functions ------------------------------------------------
kk={-1:lambda u,v:DIF(u,DIF(v,DIF(v,u))),0:lambda u,v:DIF(v,DIF(v,DIF(u,v))),1:lambda u,v:DIF(DIF(u,v),DIF(DIF(u,v),v))}
def s_p(p):
    j=next(i for i,c in enumerate(p) if c)
    def f(x):
        u=x[j] if p[j]==1 else DIF(0,x[j])
        for i,c in enumerate(p):
            if i!=j: u=kk[c](u,x[i])
        return u
    return f
for n in range(1,5):
    for p in product(D,repeat=n):
        if any(p):
            f=s_p(p); assert all(f(x)==(1 if x==p else (-1 if x==tuple(-c for c in p) else 0)) for x in product(D,repeat=n))
def build_odd(tab,n):
    pts=list(product(D,repeat=n)); reps=[p for p in pts if any(p) and p>tuple(-c for c in p)]
    ss={p:s_p(p) for p in reps}
    def F(x):
        acc=DIF(x[0],x[0])
        for p in reps:
            if tab[p]!=0:
                t=ss[p](x); acc=DIF(acc,DIF(0,t if tab[p]==1 else DIF(0,t)))   # SYN(acc,t) written with DIF only
        return acc
    return F
no=0
cells=list(product(D,D)); reps=[p for p in cells if any(p) and p>tuple(-c for c in p)]; assert len(reps)==4
for vals in product(D,repeat=4):
    tab={}
    for p,v in zip(reps,vals): tab[p]=v; tab[tuple(-c for c in p)]=-v
    tab[(0,0)]=0; F=build_odd(tab,2); assert all(F(x)==tab[x] for x in cells); no+=1
assert no==81
for n in (1,3,4):
    pts=list(product(D,repeat=n)); rr=[p for p in pts if any(p) and p>tuple(-c for c in p)]
    for _ in range(25):
        tab={(0,)*n:0}
        for p in rr: v=random.choice(D); tab[p]=v; tab[tuple(-c for c in p)]=-v
        F=build_odd(tab,n); assert all(F(x)==tab[x] for x in pts)
# DIF-closure on binary slice: exactly 81, all odd (breadth-first on 9-tuples)
A=tuple(a for a,b in cells); Bv=tuple(b for a,b in cells); known={A,Bv}; fr=[A,Bv]
while fr:
    allf=list(known); new=[]
    for f in fr:
        for g in allf:
            for x,y in ((f,g),(g,f)):
                r=tuple(DIF(x[i],y[i]) for i in range(9))
                if r not in known: known.add(r); new.append(r)
    fr=new
odd=[f for f in known if all(f[cells.index((-a,-b))]==-f[cells.index((a,b))] for a,b in cells)]
assert len(known)==81==len(odd)
# <DIF, g> complete iff g(0)!=0 (via NOT term: DIF(+1,.)) : count unary g
good=sum(1 for g in product(D,repeat=3) if g[1]!=0); assert good==18
rec('DIF_odd',{'binary_closure_size':len(known),'all_odd':True,'odd_construction_81':no,'unary_g_with_g0_nonzero':good})

# 5. classification of the 127 subsets -----------------------------------------
names=['SYN','DIF','OR','XOR','AND','NOT','ANT']
F7={**B2,'NOT':NOT,'ANT':ANT}
def ar(n): return 2 if n in B2 else 1
def preserves_set(n,S):
    k=ar(n); return all(F7[n](*t) in S for t in product(S,repeat=k))
def preserves_rho(n):
    k=ar(n)
    for ts in product(rho,repeat=k):
        cols=list(zip(*ts)); r=(F7[n](*cols[0]),F7[n](*cols[1]))
        if r[0]*r[1]<0: return False
    return True
fixes0=lambda n:F7[n](*([0]*ar(n)))==0
mins=[{'DIF','NOT'},{'SYN','ANT','NOT'}]
ncomplete=0; table=[]
for r in range(1,8):
    for sub in itertools.combinations(names,r):
        S=set(sub); crit=('NOT' in S) and ('DIF' in S or {'SYN','ANT'}<=S)
        by_min=any(m<=S for m in mins)
        inv=[]
        if all(fixes0(n) for n in S): inv.append('zero-fixing')
        if all(preserves_set(n,(0,1)) for n in S): inv.append('preserves {0,+1}')
        if all(preserves_rho(n) for n in S): inv.append('preserves rho')
        assert crit==by_min
        assert crit==(len(inv)==0), (sub,inv)   # complete exactly when no proper-subclone certificate exists
        ncomplete+=crit; table.append((sub,crit,inv))
assert ncomplete==40 and len(table)==127
assert 2**5+2**4-2**3==40
c2=[s for s,c,i in table if c and len(s)==2]; c3=[s for s,c,i in table if c and len(s)==3]
assert c2==[('DIF','NOT')] and len(c3)==6
rec('classification',{'subsets':127,'complete':ncomplete,'minimal':[sorted(m) for m in mins],'complete_pairs':len(c2),'complete_triples':len(c3)})

# 5b. free constants: completeness when 0,+1,-1 are given as nullary operations -------------
def unary_closure(Fb,Fu,consts):
    known={tuple(D)}
    if consts: known|={(c,c,c) for c in D}
    ch=True
    while ch:
        ch=False; cur=list(known)
        for n,f in Fu.items():
            for g in cur:
                r=tuple(f(v) for v in g)
                if r not in known: known.add(r); ch=True
        for n,f in Fb.items():
            for g in cur:
                for h in cur:
                    r=tuple(f(a,b) for a,b in zip(g,h))
                    if r not in known: known.add(r); ch=True
    return known
def slupecki_function(f):                   # onto, and essentially binary
    vals={f(a,b) for a in D for b in D}
    ess1=any(f(a,b)!=f(a2,b) for a in D for a2 in D for b in D); ess2=any(f(a,b)!=f(a,b2) for a in D for b in D for b2 in D)
    return len(vals)==3 and ess1 and ess2
def count_complete(Fb,Fu,consts):
    allops={**Fb,**Fu}; nm=list(allops); res=[]
    for r in range(1,len(nm)+1):
        for sub in itertools.combinations(nm,r):
            fb={n:Fb[n] for n in sub if n in Fb}; fu={n:Fu[n] for n in sub if n in Fu}
            if any(slupecki_function(f) for f in fb.values()) and len(unary_closure(fb,fu,consts))==27: res.append(sub)
    return res
c0=count_complete(B2,U1,False); c1=count_complete(B2,U1,True)
assert len(c0)==40 and {frozenset(x) for x in c0}=={frozenset(s_) for s_,c,i in table if c}          # same 40 as the certificate method
assert len(c1)==95
m1=[S for S in c1 if not any(set(T)<set(S) for T in c1)]
assert sorted(map(sorted,m1))==sorted(map(sorted,[('DIF',),('SYN','OR'),('SYN','XOR'),('SYN','AND'),('SYN','NOT'),('SYN','ANT')]))
assert sum(1 for S in c1 if 'DIF' in S)==64 and sum(1 for S in c1 if 'DIF' not in S)==31   # 2^6 and 2^5-1
rec('free_constants',{'complete_without_constants':len(c0),'complete_with_constants':len(c1),'minimal_with_constants':sorted(map(sorted,m1))})

# 5c. regression test: the counting procedure must react to changes of the signature ----------
def variant_count(Fb,Fu): return len(count_complete(Fb,Fu,False))
V={}
V['original']=variant_count(B2,U1)
V['AND->min']=variant_count({**B2,'AND':lambda a,b:min(a,b)},U1)
V['SYN->F3 addition']=variant_count({**B2,'SYN':lambda a,b:((a+b+1)%3)-1},U1)
V['NOT->cyclic shift']=variant_count(B2,{**U1,'NOT':lambda a:((a+2)%3)-1})
V['DIF->L3 strong disjunction']=variant_count({**B2,'DIF':lambda a,b:min(1,a+b+1)},U1)
assert V['original']==40 and any(v!=40 for k,v in V.items() if k!='original')
rec('regression_variants',V)

# 6. difference lift ------------------------------------------------------------
for n in (1,2):
    pts=list(product(D,repeat=n)); tot=0
    for vals in product(D,repeat=len(pts)):
        z=sum(1 for v in vals if v==0); tot+=2**z
    assert tot==4**(3**n)==(2**(3**n))**2
# exact count of (U,V) for n=1
cnt={}
for U in product((0,1),repeat=3):
    for V in product((0,1),repeat=3):
        F=tuple(u-v for u,v in zip(U,V)); cnt[F]=cnt.get(F,0)+1
assert len(cnt)==27 and all(cnt[F]==2**F.count(0) for F in cnt) and sum(cnt.values())==64
rec('difference_lift',{'n1_masks':8,'n1_pairs':64,'n2_masks':512,'n2_pairs':262144})

# 7. symmetry orbits ------------------------------------------------------------
pts=list(product(D,D)); idx={p:i for i,p in enumerate(pts)}
gens=[tuple(idx[g(p)] for p in pts) for g in (lambda p:(p[1],p[0]),lambda p:(-p[0],p[1]),lambda p:(p[0],-p[1]))]
G={tuple(range(9))}; fr=list(G)
while fr:
    nw=[]
    for g in fr:
        for h in gens:
            c=tuple(g[h[i]] for i in range(9))
            if c not in G: G.add(c); nw.append(c)
    fr=nw
assert len(G)==8
def cyc(g):
    seen=set(); c=0
    for i in range(9):
        if i not in seen:
            c+=1; j=i
            while j not in seen: seen.add(j); j=g[j]
    return c
cc=sorted(cyc(g) for g in G); assert cc==[3,3,5,6,6,6,6,9]
assert sum(3**c for c in cc)//8==2862 and sum(3**c for c in cc)%8==0
seen=set(); hist={}
for F in product(D,repeat=9):
    if F in seen: continue
    orb={tuple(F[g[i]] for i in range(9)) for g in G}; seen|=orb; hist[len(orb)]=hist.get(len(orb),0)+1
assert hist=={1:27,2:54,4:675,8:2106}
# unary: input reflection x->-x and output ANT
orb1=set(); h1={}
for g in product(D,repeat=3):
    if g in orb1: continue
    o={g,tuple(g[2-i] for i in range(3)),tuple(-v for v in g),tuple(-g[2-i] for i in range(3))}; orb1|=o; h1[len(o)]=h1.get(len(o),0)+1
assert sum(h1.values())==10 and h1=={1:1,2:5,4:4}
rec('orbits',{'group_order':8,'cycle_counts':cc,'burnside':2862,'orbit_sizes':hist,'unary_orbits':sum(h1.values()),'unary_hist':h1})
cw=lambda f:'0t'+''.join('012'[{0:0,1:1,-1:2}[f(a,b)]] if False else str({0:0,1:1,-1:2}[f(a,b)]) for a in D for b in D)
w=lambda f:(lambda s:s[:5]+'_'+s[5:8]+'_'+s[8:])(cw(f))
words={'SYN':w(SYN),'DIF':w(DIF),'AND':w(AND),'MAX':w(max),'SUM3':w(lambda a,b:((a+b+1)%3)-1),'CARRY':w(TT['CARRY'][1])}
assert words=={'SYN':'0t220_201_011','DIF':'0t022_102_110','AND':'0t102_000_201','MAX':'0t201_001_111','SUM3':'0t120_201_012','CARRY':'0t200_000_001'},words
assert 2**14<19683<=2**15
rec('T9ID_words',words)

# 8. shortest terms (tree size = operator nodes) ---------------------------------
def dp(pts,terms0,ops2,ops1,maxk):
    best={}; levels=[]
    L0={}
    for name,vec in terms0: L0[vec]=name; best[vec]=(0,name)
    levels.append(list(L0))
    for k in range(1,maxk+1):
        new={}
        def offer(r,t):                      # deterministic tie-break: lexicographically smallest printed term
            if r in best: return
            if r not in new or t<new[r]: new[r]=t
        for u in ops1:
            f=U1[u]
            for x in levels[k-1]:
                offer(tuple(f(p) for p in x),'%s(%s)'%(u,best[x][1]))
        for o in ops2:
            f=B2[o]
            for i in range(k):
                for x in levels[i]:
                    for y in levels[k-1-i]:
                        offer(tuple(f(p,q) for p,q in zip(x,y)),'%s(%s,%s)'%(o,best[x][1],best[y][1]))
        for r,t in new.items(): best[r]=(k,t)
        levels.append(list(new))
    return best,levels
b7,_=dp(list(D),[('a',tuple(D))],list(B2),list(U1),5)
assert len(b7)==27; h={}
for v,(k,t) in b7.items(): h[k]=h.get(k,0)+1
assert [h[k] for k in range(5)]==[1,4,10,7,5]
bdn,_=dp(list(D),[('a',tuple(D))],['DIF'],['NOT'],11)
assert len(bdn)==27 and max(k for k,t in bdn.values())==8
b2,_=dp(cells,[('a',tuple(p[0] for p in cells)),('b',tuple(p[1] for p in cells))],list(B2),list(U1),5)
size={}
for n,(t,f) in TT.items():
    r=b2[tuple(f(a,b) for a,b in cells)]; size[n]=r[0]
assert size=={'EQ':2,'LE':2,'LT':3,'MAX':3,'MIN':3,'imp_K3':3,'imp_L3':4,'and_L3':4,'or_L3':4,'CARRY':3,'SUM':5},size
rec('shortest_terms',{'unary_histogram_sizes_0_4':[h[k] for k in range(5)],'unary_max_over_DIF_NOT':8,'binary_connective_sizes':size})
json.dump([[list(v),list(b7[v]),list(bdn[v])] for v in sorted(b7)],open('unary_catalogue.json','w'))

# 9. optional exhaustive closures --------------------------------------------------
FULL='--full' in sys.argv; CONSTS='--full-consts' in sys.argv
ONLY=[a.split('=',1)[1] for a in sys.argv if a.startswith('--only=')]
if FULL or CONSTS:
    pw=3**np.arange(9)
    T2={n:np.array([[f(a,b)+1 for b in D] for a in D],dtype=np.int8) for n,f in B2.items()}
    T1={n:np.array([f(a)+1 for a in D],dtype=np.int8) for n,f in U1.items()}
    dec=lambda c:np.stack([(np.asarray(c,dtype=np.int64)//3**k)%3 for k in range(9)],axis=1).astype(np.int8)
    enc=lambda arr:(arr.astype(np.int64)*pw).sum(axis=-1)
    pa=int(enc(np.array([a+1 for a,b in cells]))); pb=int(enc(np.array([b+1 for a,b in cells])))
    def closure(S):
        bo=[n for n in S if n in T2]; uo=[n for n in S if n in T1]
        seen=np.zeros(3**9,bool); seen[[pa,pb]]=True; fr=np.array([pa,pb]); al=fr.copy(); layers=[2]
        while True:
            new=set(); Fa=dec(al); Ff=dec(fr)
            for n in bo:
                T=T2[n]
                for X,Y in ((Ff,Fa),(Fa,Ff)):
                    for i0 in range(0,len(X),64):
                        r=T[X[i0:i0+64][:,None,:],Y[None,:,:]]; cs=np.unique(enc(r).ravel()); new.update(cs[~seen[cs]].tolist())
            for n in uo:
                cs=np.unique(enc(T1[n][Fa])); new.update(cs[~seen[cs]].tolist())
            if not new: break
            nc=np.array(sorted(new)); seen[nc]=True; al=np.concatenate([al,nc]); fr=nc; layers.append(len(nc))
        return int(seen.sum()),layers
    full={}
    for S in ((['SYN','DIF','OR','XOR','AND'],['DIF','NOT'],['SYN','ANT','NOT'],['ANT','OR','XOR','AND','NOT'],['SYN','OR','XOR','AND','NOT']) if FULL else ()):
        full['+'.join(S)]=closure(S); rec('closure '+'+'.join(S),full['+'.join(S)])
    def closure_c(S):                      # same as closure(), constants added to the start set
        bo=[n for n in S if n in T2]; uo=[n for n in S if n in T1]
        cst=[int(enc(np.array([c+1]*9))) for c in D]
        seen=np.zeros(3**9,bool); st=[pa,pb]+cst; seen[st]=True; fr=np.array(st); al=fr.copy(); layers=[len(st)]
        while True:
            new=set(); Fa=dec(al); Ff=dec(fr)
            for n in bo:
                T=T2[n]
                for X,Y in ((Ff,Fa),(Fa,Ff)):
                    for i0 in range(0,len(X),64):
                        r=T[X[i0:i0+64][:,None,:],Y[None,:,:]]; cs=np.unique(enc(r).ravel()); new.update(cs[~seen[cs]].tolist())
            for n in uo:
                cs=np.unique(enc(T1[n][Fa])); new.update(cs[~seen[cs]].tolist())
            if not new: break
            nc=np.array(sorted(new)); seen[nc]=True; al=np.concatenate([al,nc]); fr=nc; layers.append(len(nc))
        return int(seen.sum()),layers
    for S in [S for S in (['DIF'],['SYN','OR'],['SYN','XOR'],['SYN','AND'],['SYN','NOT'],['SYN','ANT']) if CONSTS and (not ONLY or '+'.join(S) in ONLY)]:
        full['consts+'+'+'.join(S)]=closure_c(S); rec('closure with free constants '+'+'.join(S),full['consts+'+'+'.join(S)])
        assert full['consts+'+'+'.join(S)][0]==19683
    if FULL:
        assert full['SYN+DIF+OR+XOR+AND'][0]==6561 and full['DIF+NOT'][0]==19683 and full['SYN+ANT+NOT'][0]==19683
        assert full['ANT+OR+XOR+AND+NOT'][0]<19683 and full['SYN+OR+XOR+AND+NOT'][0]<19683
import platform
rec('environment',{'python':platform.python_version(),'numpy':np.__version__})
rec('status','ALL ASSERTIONS PASSED')
tag='_default' if (not ONLY and not FULL) else ('' if not ONLY else '_consts_'+'_'.join(ONLY))   # results.json is written only by --full runs
json.dump(R,open('results%s.json'%tag,'w'),indent=1,default=str)
with open('results%s.md'%tag,'w') as f:
    f.write('# Chai checker results\n\n')
    for k,v in R.items(): f.write(f'- **{k}**: `{v}`\n')
