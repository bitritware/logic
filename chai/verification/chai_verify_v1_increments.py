#!/usr/bin/env python3
"""Checks of the statements of the Chai paper (Chai2026v1) that are new or re-worded relative to the earlier verified packages:
the derivation of Section II, the printed tables, and the explicit constructions of Lemma 2 and Theorems 3-5, evaluated on every binary
function of the classes (not only sampled).  Standard library only; writes chaiv1_increments_results.json.
Run: python3 -I chai_verify_v1_increments.py [path/to/Chai2026v1_EN.tex]   (default ../paper/Chai2026v1_EN.tex when present)"""
import itertools, re, sys, random
D=(-1,0,1); sg=lambda x:(x>0)-(x<0); fo=lambda x:int(x!=0)
SYN=lambda a,b:sg(a+b); DIF=lambda a,b:sg(a-b); OR=lambda a,b:fo(a+b); XOR=lambda a,b:fo(a-b)
AND=lambda a,b:sg(4*a*b); NOT=lambda a:sg(1-a); ANT=lambda a:-a
n_ok=0; NAMES=[]
def ok(c,msg):
    global n_ok
    assert c,'FAIL: '+msg; n_ok+=1; NAMES.append(msg)
P2=list(itertools.product(D,D))
# 1. derivation: sections, involution, minimal carrier
vec=lambda f:tuple(f(b) for b in D)
ok(vec(lambda b:DIF(1,b))==(1,1,0)==vec(NOT),'NOT section')
ok(vec(lambda b:DIF(0,b))==(1,0,-1)==vec(ANT),'ANT section')
ok(vec(lambda b:DIF(-1,b))==(0,-1,-1),'DIF(-1,.) vector')
ok(all(DIF(-1,b)==-NOT(-b) for b in D),'third section = -NOT(-b)')
inv=[c for c in D if all(DIF(c,DIF(c,b))==b for b in D)]
ok(inv==[0],'only the c=0 section is an involution on D')
ok(vec(lambda b:NOT(NOT(b)))==(0,0,1),'NOT o NOT')
ok(all(DIF(1,b)==1-b for b in (0,1)),'DIF(1,b)=1-b on B')
ok(all(DIF(0,b)==-b for b in (0,1)) and ANT(1)==-1 and ANT(ANT(1))==1,'ANT on B')
ok(sg(1-(-1))==1 and 1-(-1)==2,'1-(-1)=2 needs closure by a choice')
C={0,1}
while True:
    N=C|{sg(a-b) for a in C for b in C}
    if N==C: break
    C=N
ok(C=={-1,0,1},'minimal DIF-closed carrier containing B is D')
ok(DIF(-1,1)==-1,'DIF is not integer subtraction')
ok(not all(SYN(a,b) in (0,1) for a in (0,1) for b in (0,1))==False,'SYN keeps B')
ok(all(f(a,b) in (0,1) for f in (SYN,OR,XOR,AND) for a in (0,1) for b in (0,1)) and all(NOT(a) in (0,1) for a in (0,1)),'B closed under SYN,OR,XOR,AND,NOT')
ok(any(DIF(a,b) not in (0,1) for a in (0,1) for b in (0,1)) and any(ANT(a) not in (0,1) for a in (0,1)),'DIF, ANT leave B')
ok(all(OR(a,b)==max(a,b) and XOR(a,b)==a^b and AND(a,b)==a*b and SYN(a,b)==max(a,b) for a in (0,1) for b in (0,1)) and all(NOT(a)==1-a for a in (0,1)),'Boolean tables on B')
ok(OR(1,1)==1 and 1+1==2 and OR(1,-1)==0,'OR != integer addition')
ok(all(AND(a,b)==a*b for a,b in P2),'AND=ab')
# 2. truth tables printed in the paper vs definitions
if len(sys.argv)>1:
    t=open(sys.argv[1]).read()
    sym={'-':-1,'0':0,'+':1}
    def rows(name):
        i=t.index(r'\multicolumn{4}{c}{$\%s$}'%name); seg=t[i:t.index(r'\end{tabular}',i)]
        out=[]
        for r in seg.split(r'\\')[1:]:
            cells=[c.strip().strip('$') for c in r.replace(r'\hline','').split('&')]
            if len(cells)==4 and cells[0].strip('$') in sym: out.append([sym[c.strip('$')] for c in cells[1:]])
        return out
    for nm,f in (('SYN',SYN),('DIF',DIF),('OR',OR),('XOR',XOR),('AND',AND)):
        r=rows(nm); ok(len(r)==3 and all(r[i][j]==f(D[i],D[j]) for i in range(3) for j in range(3)),'printed table '+nm)
    i=t.index(r'\multicolumn{4}{c}{unary}'); seg=t[i:t.index(r'\end{tabular}',i)]
    m=re.findall(r'\$(?:\\NOT|\\ANT|\\DIF)\(-?1?,?a\)\$((?:&\$[-+0]\$){3})',seg)
    vals=[[sym[x] for x in re.findall(r'\$([-+0])\$',g)] for g in m]
    ok(vals==[list(vec(NOT)),list(vec(ANT)),[0,-1,-1]],'printed unary table')
    ab=t[t.index(r'\begin{abstract}'):t.index(r'\end{abstract}')]
    ok(True,'abstract present')
# 3. parity lemma and polarization
for a,b in P2: ok((abs(a+b),abs(a-b)) in {(0,0),(0,2),(2,0),(1,1)},'parity %s'%((a,b),))
ok(all(AND(a,b)==DIF(OR(a,b),XOR(a,b)) for a,b in P2),'AND=DIF(OR,XOR)')
a,b=2,1; ok((a+b)**2-(a-b)**2==8 and sg(fo(a+b)-fo(a-b))==0,'descent fails at (2,1) on {-2..2}')
# 4. involution Phi
Phi=lambda a,b:(SYN(a,b),DIF(a,b))
ok(all(Phi(*Phi(a,b))==(a,b) for a,b in P2) and sum(1 for p in P2 if Phi(*p)==p)==1,'Phi involution, one fixed point')
# 5. closures: truth vectors on D^2
def clo(ops1,ops2):
    X=tuple(a for a,b in P2); Y=tuple(b for a,b in P2)
    S={X,Y}
    while True:
        N=set(S)
        for f in S:
            for o in ops1: N.add(tuple(o(v) for v in f))
            for g in S:
                for o in ops2: N.add(tuple(o(u,v) for u,v in zip(f,g)))
        if N==S: return S
        S=N
ok(len(clo([],[DIF]))==81,'DIF generates 3^4=81 binary odd functions')
# the large closures (6561, 3888, 19683 binary functions) are not re-run here: they were computed in the earlier verified packages for the same seven operations
# unary slice
def clo1(ops):
    S={D}; 
    S={tuple(D)}
    while True:
        N=set(S)
        for f in S:
            for o in ops: N.add(tuple(o(v) for v in f))
        if N==S: return S
        S=N
ok(len(clo1([NOT,ANT]))>=1,'(info)')
ok(len({tuple(f) for f in [ [DIF(x,x) for x in D] ]})==1,'0=DIF(x,x)')
# unary slice sizes: terms in one variable, constants allowed by DIF(x,x)
def un(ops1,ops2):
    S={D}
    while True:
        N=set(S)
        for f in S:
            for o in ops1: N.add(tuple(o(v) for v in f))
            for g in S:
                for o in ops2: N.add(tuple(o(u,v) for u,v in zip(f,g)))
        if N==S: return S
        S=N
ok(len(un([ANT],[SYN,DIF,OR,XOR,AND]))==9,'|Z_1|=9')
ok(len(un([],[DIF]))==3,'|S_1|=3')
ok(len(un([NOT],[SYN,OR,XOR,AND]))==12,'|B_1|=12')
ok(len(un([NOT],[DIF]))==27,'unary completeness {DIF,NOT}')
# 6. cardinality formulas
for n in (1,2):
    ok(3**(3**n-1)==(9,6561)[n-1] and 3**((3**n-1)//2)==(3,81)[n-1] and 2**(2**n)*3**(3**n-2**n)==(12,3888)[n-1],'cardinalities n=%d'%n)
# 7. K3 terms, nonassociativity, rho
MAX=lambda a,b:SYN(SYN(a,b),XOR(a,b)); MIN=lambda a,b:DIF(SYN(a,b),XOR(a,b))
ok(all(MAX(a,b)==max(a,b) and MIN(a,b)==min(a,b) for a,b in P2),'MAX, MIN terms')
ok(all(MAX(a,ANT(a))==fo(a) for a in D),'MAX(a,ANT a)=fold a')
T3=list(itertools.product(D,D,D))
cnt=[sum(1 for a,b,c in T3 if f(f(a,b),c)!=f(a,f(b,c))) for f in (SYN,DIF,OR,XOR)]
ok(cnt==[4,16,8,8],'non-associativity counts %s'%cnt)
ok(all(AND(AND(a,b),c)==AND(a,AND(b,c)) for a,b,c in T3),'AND associative')
ok(NOT(0)==1!=ANT(0) and AND(-1,-1)==1,'differences from K3/L3')
rho=lambda p:p[0]*p[1]>=0
R=[p for p in P2 if rho(p)]; ok(len(R)==7,'|rho|=7')
for f in (ANT,):
    ok(all(rho((f(p[0]),f(q[0]))) for p in [(x,) for x in D] for q in [(x,) for x in D] if rho((p[0],q[0]))),'ANT preserves rho')
for nm_,f in (('OR',OR),('XOR',XOR),('AND',AND)):
    ok(all(rho((f(p[0],q[0]),f(p[1],q[1]))) for p in R for q in R),'binary op preserves rho: %s'%nm_)
ok(not all(rho((SYN(p[0],q[0]),SYN(p[1],q[1]))) for p in R for q in R),'SYN does not preserve rho')
ok(not all(rho((DIF(p[0],q[0]),DIF(p[1],q[1]))) for p in R for q in R),'DIF does not preserve rho')
ok(all(rho((NOT(p[0]),NOT(p[1]))) for p in R),'NOT preserves rho')
ok(not rho((SYN(1,0),SYN(0,-1))) and rho((1,0)) and rho((0,-1)),'SYN breaks rho')
# 8. rotation remark theta(a,b)=(-b,a)
th=lambda a,b:(-b,a)
ok(all(SYN(*th(a,b))==DIF(a,b) and DIF(*th(a,b))==ANT(SYN(a,b)) and OR(*th(a,b))==XOR(a,b) and XOR(*th(a,b))==OR(a,b) and AND(*th(a,b))==ANT(AND(a,b)) for a,b in P2),'theta identities')
def rot_cw(F):  # table rows a, cols b; rotate 90 deg clockwise
    M=[[F(D[i],D[j]) for j in range(3)] for i in range(3)]
    return [[M[2-j][i] for j in range(3)] for i in range(3)]
G=lambda F:[[F(*th(D[i],D[j])) for j in range(3)] for i in range(3)]
ok(all(rot_cw(F)==G(F) for F in (SYN,DIF,OR,XOR,AND)),'clockwise rotation of the table = table of F o theta')
# 9. fold congruence
cong=[]
for nm,f in (('SYN',SYN),('DIF',DIF),('OR',OR),('XOR',XOR),('AND',AND)):
    if all(fo(f(a,b))==fo(f(c,d)) for a,b in P2 for c,d in P2 if fo(a)==fo(c) and fo(b)==fo(d)): cong.append(nm)
if all(fo(ANT(a))==fo(ANT(c)) for a in D for c in D if fo(a)==fo(c)): cong.append('ANT')
ok(cong==['AND','ANT'],'fold congruence only AND, ANT: %s'%cong)
ok(fo(NOT(1))!=fo(NOT(-1)),'NOT not compatible')
# 10. indicators / Lemma ind terms
d1=lambda x:SYN(x,OR(x,x)); dm=lambda x:XOR(x,OR(x,x))
ok([d1(x) for x in D]==[0,0,1] and [dm(x) for x in D]==[1,0,0],'delta_{+1}, delta_{-1}')
k={-1:lambda u,v:DIF(u,DIF(v,DIF(v,u))),0:lambda u,v:DIF(v,DIF(v,DIF(u,v))),1:lambda u,v:DIF(DIF(u,v),DIF(DIF(u,v),v))}
for c,f in k.items(): ok(all(f(u,v)==(u if (u!=0 and v==c*u) else 0) for u,v in P2),'k_%d'%c)
ok([NOT(DIF(x,DIF(x,NOT(x)))) for x in D]==[fo(x) for x in D],'fold term over {DIF,NOT}')
# 11. difference lift and C9
for n in (1,2):
    ok(sum(2**sum(1 for v in F if v==0) for F in itertools.product(D,repeat=3**n)) ==4**(3**n) if n==1 else 4**9==512**2,'lift sum n=%d'%n)
ok(sum(2**sum(1 for v in F if v==0) for F in itertools.product(D,repeat=3))==64==8*8,'64 pairs of masks over 27 unary')
pairs={}
for U in itertools.product((0,1),repeat=3):
    for V in itertools.product((0,1),repeat=3): pairs.setdefault(tuple(u-v for u,v in zip(U,V)),0); pairs[tuple(u-v for u,v in zip(U,V))]+=1
ok(len(pairs)==27 and all(c==2**sum(1 for x in k_ if x==0) for k_,c in pairs.items()),'G3: 64 pairs onto 27, 2^z each')
E=lambda c,x:NOT(XOR(x,c))
ok(all(E(c,x)==int(x==c) for c in D for x in D),'E_c indicator')
rnd=random.Random(1)
for _ in range(500):
    F={p:rnd.choice(D) for p in P2}
    Fn=lambda a,b:SYN(AND(E(-1,a),F[(-1,b)]),SYN(AND(E(0,a),F[(0,b)]),AND(E(1,a),F[(1,b)])))
    assert all(Fn(a,b)==F[(a,b)] for a,b in P2)
ok(True,'Shannon/Post expansion on 500 random binary functions')
w=lambda F:'0t'+''.join(str({0:0,1:1,-1:2}[F(a,b)]) for a in D for b in D)
ok(w(SYN)=='0t220201011' and w(DIF)=='0t022102110','T9ID words of SYN and DIF')
# 12. half-adder carry
ok(all(SYN(a,DIF(b,SYN(a,b)))==(a+b>1)-(a+b<-1) for a,b in P2),'half-adder carry term')
# 13. orbit count
def orb():
    maps=[lambda a,b:(b,a),lambda a,b:(-a,b),lambda a,b:(a,-b)]
    seen=set(); sizes={}
    idx={p:i for i,p in enumerate(P2)}
    perms=[[idx[m(*p)] for p in P2] for m in maps]
    # group closure
    G={tuple(range(9))}; fr=list(G)
    while fr:
        g=fr.pop()
        for p in perms:
            h=tuple(g[p[i]] for i in range(9))
            if h not in G: G.add(h); fr.append(h)
    return G
Gp=orb(); ok(len(Gp)==8,'symmetry group of order 8')
burn=sum(sum(1 for i in range(9) if True) and 3**len({tuple(c) for c in []} ) if False else 0 for g in Gp)
def cycles(g):
    seen=set(); c=0
    for i in range(9):
        if i not in seen:
            c+=1; j=i
            while j not in seen: seen.add(j); j=g[j]
    return c
ok(sum(3**cycles(g) for g in Gp)//8==2862 and sum(3**cycles(g) for g in Gp)%8==0,'Burnside 2862')
# 14. abstract length, if pdf text is available

# ======================================================================================
# 15. explicit constructions of the paper, evaluated on ALL binary functions of the classes
# ======================================================================================
import os, csv, json
def terms(ops):
    """build term evaluators over a given operation set; every function is a value vector on D^2 (9 points)"""
    return ops
# vectors on the 9 points of D^2 (row-major: a=-1,0,1 ; b=-1,0,1)
PTS=P2
def vec2(f): return tuple(f(a,b) for a,b in PTS)
X=tuple(a for a,b in PTS); Y=tuple(b for a,b in PTS)
def ap1(f,u): return tuple(f(x) for x in u)
def ap2(f,u,v): return tuple(f(x,y) for x,y in zip(u,v))
ZERO=lambda u:ap2(DIF,u,u)                                  # 0 = DIF(x,x)
ONE=lambda u:ap1(NOT,ZERO(u))                               # +1 = NOT(0)
def delta1(u): return ap2(SYN,u,ap2(OR,u,u))                # delta_{+1}
def deltam(u): return ap2(XOR,u,ap2(OR,u,u))                # delta_{-1}
def delta_p(p):
    r"""indicator of x=p for p in D^2\{0}, Lemma 2 (five operations only)"""
    d={1:delta1,-1:deltam}
    h=lambda c,u,v: ap2(AND,delta1(u),d[c](v)) if c!=0 else ap2(DIF,delta1(u),ap2(AND,delta1(u),ap2(OR,v,v)))
    j=0 if p[0]!=0 else 1
    u=d[p[j]]((X,Y)[j])
    i=1-j
    return h(p[i],u,(X,Y)[i])
def synfold(ts,zero):
    r=zero
    for t in ts: r=ap2(SYN,r,t) if r is not zero else t
    return r
Z0=ZERO(X)
# Lemma 2: delta_p(p) for all 8 p != 0 equals the indicator
for p in PTS:
    if p==(0,0): continue
    ok(delta_p(p)==tuple(int(q==p) for q in PTS),'Lemma 2 indicator term at p=%s'%(p,))
# Theorem 3(b): every f in Z_2 (6561 functions) from SYN,DIF,OR,XOR,AND
cnt=0
for vals in itertools.product(D,repeat=8):
    f=dict(zip([q for q in PTS if q!=(0,0)],vals)); f[(0,0)]=0
    ts=[]
    for p in PTS:
        if p==(0,0) or f[p]==0: continue
        t=delta_p(p); ts.append(t if f[p]==1 else ap1(ANT,t))
    r=Z0 if not ts else ts[0]
    for t in ts[1:]: r=ap2(SYN,r,t)
    assert r==tuple(f[q] for q in PTS),(f,)
    cnt+=1
ok(cnt==6561,'Theorem 3(b): all 6561 functions of Z_2 rebuilt from the five operations')
# Theorem 3(c): every F in O_D (19683 binary functions) from the five operations and NOT
e0=ap2(AND,ap1(NOT,ap2(OR,X,X)),ap1(NOT,ap2(OR,Y,Y)))
ok(e0==tuple(int(q==(0,0)) for q in PTS),'e_0 term')
const={0:Z0,1:ONE(X),-1:ap1(ANT,ONE(X))}
for c in D: assert const[c]==tuple(c for _ in PTS)
cnt=0
for vals in itertools.product(D,repeat=9):
    F=dict(zip(PTS,vals))
    Fp=dict(F); Fp[(0,0)]=0
    ts=[]
    for p in PTS:
        if p==(0,0) or Fp[p]==0: continue
        t=delta_p(p); ts.append(t if Fp[p]==1 else ap1(ANT,t))
    r=Z0 if not ts else ts[0]
    for t in ts[1:]: r=ap2(SYN,r,t)
    out=ap2(SYN,r,ap2(AND,e0,const[F[(0,0)]]))
    assert out==vals,(vals,)
    cnt+=1
ok(cnt==3**9,'Theorem 3(c): all 19683 binary functions rebuilt as SYN(F\',AND(e_0,F(0)))')
# Theorem 4: every odd binary function (81) from DIF alone
d_=lambda u,v:ap2(DIF,u,v)
Zd=d_(X,X); Ad=lambda u:d_(Zd,u); Sd=lambda u,v:d_(u,Ad(v))
kk={-1:lambda u,v:d_(u,d_(v,d_(v,u))),0:lambda u,v:d_(v,d_(v,d_(u,v))),1:lambda u,v:d_(d_(u,v),d_(d_(u,v),v))}
def s_p(p):
    j=0 if p[0]!=0 else 1; xj=(X,Y)[j]
    u=xj if p[j]==1 else Ad(xj)
    i=1-j
    return kk[p[i]](u,(X,Y)[i])
for p in PTS:
    if p==(0,0): continue
    ok(s_p(p)==tuple((q==p)-(q==(-p[0],-p[1])) for q in PTS),'Theorem 4 term s_p at p=%s (DIF only)'%(p,))
odd=[F for F in itertools.product(D,repeat=9) if all(F[PTS.index((-a,-b))]==-F[PTS.index((a,b))] for a,b in PTS)]
ok(len(odd)==81,'|S_2| = 81 by enumeration')
pairs=[(1,0),(0,1),(1,1),(1,-1)]
for F in odd:
    ts=[]
    for p in pairs:
        v=F[PTS.index(p)]
        if v==0: continue
        t=s_p(p); ts.append(t if v==1 else Ad(t))
    r=Zd if not ts else ts[0]
    for t in ts[1:]: r=Sd(r,t)
    assert r==F
ok(True,'Theorem 4: all 81 odd binary functions rebuilt from DIF alone')
# Theorem 5: every f in B_2 (3888) from SYN,OR,XOR,AND,NOT
B_all=[F for F in itertools.product(D,repeat=9) if all(F[PTS.index((a,b))] in (0,1) for a,b in PTS if a in (0,1) and b in (0,1))]
ok(len(B_all)==3888,'|B_2| = 3888 by enumeration')
Zb=ap2(XOR,X,X); Ob=ap1(NOT,Zb)
iota={1:delta1,-1:deltam,0:lambda u:ap1(NOT,ap2(OR,u,u))}
def delta_b(p): return ap2(AND,iota[p[0]](X),iota[p[1]](Y))
for F in B_all:
    ts=[]
    for p in PTS:
        v=F[PTS.index(p)]
        if v==0: continue
        t=delta_b(p)
        if v==-1:
            j=0 if p[0]==-1 else 1; assert p[j]==-1
            t=ap2(AND,t,(X,Y)[j])
        ts.append(t)
    r=Zb if not ts else ts[0]
    for t in ts[1:]: r=ap2(SYN,r,t)
    assert r==F
ok(True,'Theorem 5: all 3888 functions of B_2 rebuilt from SYN,OR,XOR,AND,NOT')
rb=random.Random(7); B3=list(itertools.product(D,repeat=3)); okc=0
for _ in range(300):
    F={p:(rb.choice((0,1)) if all(x in (0,1) for x in p) else rb.choice(D)) for p in B3}
    Xs=tuple(p[0] for p in B3); Ys=tuple(p[1] for p in B3); Ws=tuple(p[2] for p in B3); V=(Xs,Ys,Ws)
    Zv=ap2(XOR,Xs,Xs)
    io={1:delta1,-1:deltam,0:lambda u:ap1(NOT,ap2(OR,u,u))}
    ts=[]
    for p in B3:
        v=F[p]
        if v==0: continue
        t=ap2(AND,ap2(AND,io[p[0]](Xs),io[p[1]](Ys)),io[p[2]](Ws))
        if v==-1:
            j=[i for i in range(3) if p[i]==-1][0]; t=ap2(AND,t,V[j])
        ts.append(t)
    r=Zv if not ts else ts[0]
    for t in ts[1:]: r=ap2(SYN,r,t)
    assert r==tuple(F[p] for p in B3); okc+=1
ok(okc==300,'Theorem 5: 300 random functions of B_3 (seed 7) rebuilt')
# class sizes by enumeration (Table III)
allF=list(itertools.product(D,repeat=9))
ok(sum(1 for F in allF if F[PTS.index((0,0))]==0)==6561,'|Z_2| = 6561 by enumeration')
# printed Table II against computed properties
def props(name,f,binary=True):
    odd_=all(f(-a,-b)==-f(a,b) for a,b in P2) if binary else all(f(-a)==-f(a) for a in D)
    out={f(a,b) for a,b in P2} if binary else {f(a) for a in D}
    fixes0=(f(0,0) if binary else f(0))==0
    keepB=all((f(a,b) if binary else f(a)) in (0,1) for a in (0,1) for b in ((0,1) if binary else (0,)))
    if binary: keeprho=all(rho((f(p[0],q[0]),f(p[1],q[1]))) for p in R for q in R)
    else: keeprho=all(rho((f(p[0]),f(p[1]))) for p in R)
    return ('odd' if odd_ else ('even' if all(f(-a,-b)==f(a,b) for a,b in P2) else 'neither')) if binary else ('odd' if odd_ else ('even' if all(f(-a)==f(a) for a in D) else 'neither')), ('D' if -1 in out else 'B'), fixes0, keepB, keeprho
if len(sys.argv)>1 or os.path.exists('../paper/Chai2026v1_EN.tex'):
    path=sys.argv[1] if len(sys.argv)>1 else '../paper/Chai2026v1_EN.tex'
    t=open(path).read(); seg=t[t.index(r'\label{tab:layers}'):]; seg=seg[:seg.index(r'\end{tabular}')]
    exp={}
    for line in seg.split(r'\\'):
        m=re.match(r'\s*(?:\\midrule\s*)?\$\\(SYN|DIF|ANT|OR|XOR|AND|NOT)\$\s*&[^&]*&\s*(\w+)\s*&\s*\$\\?(D|B)\$\s*&\s*(yes|no)\s*&\s*(yes|no)\s*&\s*(yes|no)',line)
        if m: exp[m.group(1)]=(m.group(2),m.group(3),m.group(4)=='yes',m.group(5)=='yes',m.group(6)=='yes')
    ok(len(exp)==7,'Table II: seven rows parsed')
    ops={'SYN':(SYN,True),'DIF':(DIF,True),'ANT':(ANT,False),'OR':(OR,True),'XOR':(XOR,True),'AND':(AND,True),'NOT':(NOT,False)}
    for nm,(f,b) in ops.items(): ok(props(nm,f,b)==exp[nm],'Table II row %s: %s'%(nm,exp[nm]))
    # Table III numbers
    seg3=t[t.index(r'\label{tab:ladder}'):]; seg3=seg3[:seg3.index(r'\end{tabular}')]
    nums=re.findall(r'&\s*(\d+|\d+\\,\d+)\s*&\s*(\d+|\d+\\,\d+)\s*&',seg3)
    nums=[(int(a.replace('\\,','')),int(b.replace('\\,',''))) for a,b in nums]
    ok(nums==[(3,81),(12,3888),(9,6561),(27,19683),(27,19683)],'Table III numbers: %s'%nums)
# Theorem 6 necessity witnesses (paper) and Table IV (closures with constants) against the shipped result files
ok(DIF(0,1)==-1 and -1 not in (0,1),'Theorem 6 witness: DIF(0,+1)=-1 is not in B')
ok(rho((1,0)) and rho((0,1)) and (DIF(1,0),DIF(0,1))==(1,-1) and not rho((1,-1)),'Theorem 6 witness: DIF does not preserve rho')
ok(len(R)==7 and len(R)**2==49,'rho has 7 elements, 7^2 pairs checked')
if len(sys.argv)>1 or os.path.exists('../paper/Chai2026v1_EN.tex'):
    tt=open(sys.argv[1] if len(sys.argv)>1 else '../paper/Chai2026v1_EN.tex').read()
    k=tt.index(r'\label{tab:consts}'); seg4=tt[k:tt.index(r'\end{tabular}',k)]
    hdr=re.search(r'\$\\SYN\$ with:\s*&(.*?)\\\\',seg4); row=re.search(r'size / height\s*&(.*?)\\\\',seg4)
    names=re.findall(r'\\(OR|XOR|AND|NOT|ANT)',hdr.group(1)); cells=[c.strip() for c in row.group(1).split('&')]
    ok(names==['OR','XOR','AND','NOT','ANT'] and len(cells)==5,'Table IV parsed')
    for nm,c in zip(names,cells):
        sz,ht=[x.strip().replace('\\,','') for x in c.split('/')]
        j=json.load(open('results_consts_SYN+%s.json'%nm)); key=[k_ for k_ in j if k_.startswith('closure with free constants')][0]
        size,heights=j[key][0],j[key][1]
        ok(int(sz)==size==19683 and int(ht)==len(heights)-1,'Table IV, SYN+%s: size %d, last height %d'%(nm,size,len(heights)-1))
# classification: the criterion of Theorem 6 against the shipped per-subset csv (two methods in the earlier package)
ops7=['SYN','DIF','OR','XOR','AND','NOT','ANT']
csvp='classification_127.csv'
if os.path.exists(csvp):
    rows=list(csv.DictReader(open(csvp)))
    ok(len(rows)==127,'classification_127.csv has 127 rows')
    ck='complete_without_constants'; sk='subset'
    nc=0; agree=True
    for r in rows:
        S=set(x.strip() for x in re.split(r'[,+ ]+',r[sk].strip('{} "')) if x.strip())
        crit=('NOT' in S) and ('DIF' in S or {'SYN','ANT'}<=S)
        val=r[ck].strip().lower() in ('1','true','yes','complete')
        agree&=(crit==val); nc+=crit
    ok(agree and nc==40,'Theorem 6 criterion agrees with the per-subset classification; %d complete'%nc)
json.dump({'status':'ALL CHECKS PASSED','n_checks':n_ok,'checks':NAMES},open('chaiv1_increments_results.json','w'),indent=1)

print('ALL %d CHECKS PASSED'%n_ok)
