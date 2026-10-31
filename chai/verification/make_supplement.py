#!/usr/bin/env python3
"""Builds chai2026v1_supplement.tex from the result files. Needs the results of the --full run (results.json),
results_consts_*.json, unary_catalogue.json and classification_127.csv; it stops with an error if any of them lacks data,
so an incomplete supplement cannot be produced. Needs also direction_rotation_results.json (verify_direction_rotation.py). Usage: python3 make_supplement.py && pdflatex chai2026v1_supplement.tex (twice)."""
import json, csv, glob, os, re, itertools
from itertools import product
HERE=os.path.dirname(os.path.abspath(__file__)); P=lambda f:os.path.join(HERE,f)
D=(-1,0,1); sgn=lambda z:(z>0)-(z<0)
SYN=lambda a,b:sgn(a+b); DIF=lambda a,b:sgn(a-b); OR=lambda a,b:sgn((a+b)**2); XOR=lambda a,b:sgn((a-b)**2)
AND=lambda a,b:sgn(4*a*b); NOT=lambda a:sgn(1-a); ANT=lambda a:-a
B2={'SYN':SYN,'DIF':DIF,'OR':OR,'XOR':XOR,'AND':AND}; U1={'NOT':NOT,'ANT':ANT}
sym={-1:'-',0:'0',1:'+'}; M=lambda v:'$%s$'%sym[v]
res=json.load(open(P('results.json')))
need=['closure SYN+DIF+OR+XOR+AND','closure DIF+NOT','closure SYN+ANT+NOT','closure ANT+OR+XOR+AND+NOT','closure SYN+OR+XOR+AND+NOT']
assert all(k in res for k in need) and res['status']=='ALL ASSERTIONS PASSED', 'results.json lacks the closure entries: run chai_checker.py --full first'
for g in sorted(glob.glob(P('results_consts_*.json'))):
    res.update({k:v for k,v in json.load(open(g)).items() if k.startswith('closure with free constants')})
CK=['DIF','SYN+OR','SYN+XOR','SYN+AND','SYN+NOT','SYN+ANT']
assert all('closure with free constants '+k in res for k in CK), 'missing results_consts_*.json'
cat=json.load(open(P('unary_catalogue.json')))
cls=list(csv.DictReader(open(P('classification_127.csv'))))
assert len(cls)==127
log=open(P('full_run.log')).read(); runtime=float(re.findall(r'\[\s*([0-9.]+)s\] status',log)[-1])

def esc(t): return t.replace('_','\\_')
def tex(t): return r'\texttt{%s}'%esc(t).replace(',',',\\allowbreak ').replace('(','(\\allowbreak ')
def ml(t): return re.sub(r'\b(SYN|DIF|OR|XOR|AND|NOT|ANT|MAX|MIN|MAG)\b',r'\\\1',t)
def tt3(f,title):
    rows=''.join(r'%s&%s&%s&%s\\'%(M(a),M(f(a,-1)),M(f(a,0)),M(f(a,1))) for a in D)
    return r'\begin{tabular}{c|ccc}\multicolumn{4}{c}{%s}\\\hline &%s&%s&%s\\\hline %s\end{tabular}'%(title,M(-1),M(0),M(1),rows)
def code(v): return '0t'+''.join({0:'0',1:'1',-1:'2'}[x] for x in v)
def word(f):
    s=code([f(a,b) for a in D for b in D])[2:]; return '0t'+s[:3]+'_'+s[3:6]+'_'+s[6:]
iv=lambda c:1 if c else 0

# ---- S1 laws
ops={**B2,'MAX':max,'MIN':min}
rows=[]
for n,f in ops.items():
    comm=all(f(a,b)==f(b,a) for a in D for b in D); idem=all(f(a,a)==a for a in D)
    idn=[sym[e] for e in D if all(f(e,a)==a and f(a,e)==a for a in D)]; ab=[sym[z] for z in D if all(f(z,a)==z and f(a,z)==z for a in D)]
    af=sum(f(f(a,b),c)!=f(a,f(b,c)) for a in D for b in D for c in D)
    rows.append(r'\textsf{%s}&%s&%d/27&%s&%s&%s\\'%(n,'yes' if comm else 'no',af,'yes' if idem else 'no',','.join(idn) or '--',','.join(ab) or '--'))
LAWS='\n'.join(rows)
ld=[(f,g) for f in ops for g in ops if all(ops[f](a,ops[g](b,c))==ops[g](ops[f](a,b),ops[f](a,c)) for a in D for b in D for c in D)]
LDT=', '.join(r'$(\mathsf{%s},\mathsf{%s})$'%p for p in ld)

# ---- S3 explicit terms and their full tables
d1=lambda x:SYN(x,OR(x,x)); dm=lambda x:XOR(x,OR(x,x)); mag=lambda x:NOT(DIF(x,DIF(x,NOT(x))))
for a in D: assert d1(a)==iv(a==1) and dm(a)==iv(a==-1) and mag(a)==a*a
h={1:lambda u,v:AND(d1(u),d1(v)),-1:lambda u,v:AND(d1(u),dm(v)),0:lambda u,v:DIF(d1(u),AND(d1(u),OR(v,v)))}
for u,v in product(D,D):
    assert h[1](u,v)==iv(u==1 and v==1) and h[-1](u,v)==iv(u==1 and v==-1) and h[0](u,v)==iv(u==1 and v==0)
kt={-1:lambda u,v:DIF(u,DIF(v,DIF(v,u))),0:lambda u,v:DIF(v,DIF(v,DIF(u,v))),1:lambda u,v:DIF(DIF(u,v),DIF(DIF(u,v),v))}
for c in D:
    for u,v in product(D,D): assert kt[c](u,v)==(u if (u!=0 and v==c*u) else 0)
tabs_h=r'\quad'.join(tt3(h[c],r'$h_{%s}(u,v)$'%{-1:'-1',0:'0',1:'+1'}[c]) for c in (-1,0,1))
tabs_k=r'\quad'.join(tt3(kt[c],r'$k_{%s}(u,v)$'%{-1:'-1',0:'0',1:'+1'}[c]) for c in (-1,0,1))
unary_small=r'\begin{tabular}{c|ccc}$x$&%s&%s&%s\\\hline $\delta_{+1}$&%s&%s&%s\\ $\delta_{-1}$&%s&%s&%s\\ $\mathsf{MAG}$&%s&%s&%s\end{tabular}'%tuple(
    [M(-1),M(0),M(1)]+[M(d1(a)) for a in D]+[M(dm(a)) for a in D]+[M(mag(a)) for a in D])

# ---- S4 heights
names={'SYN+DIF+OR+XOR+AND':r'$\SYN,\DIF,\OR,\XOR,\AND$ (zero-fixing clone $\Z_2$)','DIF+NOT':r'$\DIF,\NOT$','SYN+ANT+NOT':r'$\SYN,\ANT,\NOT$',
 'SYN+OR+XOR+AND+NOT':r'$\SYN,\OR,\XOR,\AND,\NOT$ ($M_1$, incomplete)','ANT+OR+XOR+AND+NOT':r'$\ANT,\OR,\XOR,\AND,\NOT$ ($M_2$, incomplete)'}
hr=[r'%s&$%d$&%s\\'%(lab,res['closure '+k][0],', '.join(map(str,res['closure '+k][1]))) for k,lab in names.items()]
cn={'DIF':r'$\DIF$','SYN+OR':r'$\SYN,\OR$','SYN+XOR':r'$\SYN,\XOR$','SYN+AND':r'$\SYN,\AND$','SYN+NOT':r'$\SYN,\NOT$','SYN+ANT':r'$\SYN,\ANT$'}
hr+=[r'%s + nullary constants&$%d$&%s\\'%(lab,res['closure with free constants '+k][0],', '.join(map(str,res['closure with free constants '+k][1]))) for k,lab in cn.items()]
HEIGHTS='\n'.join(hr)

# ---- S5 unary tables
def dp(maxk,ops2,ops1):
    best={tuple(D):(0,'a')}; levels=[[tuple(D)]]
    for k in range(1,maxk+1):
        new={}
        def offer(r,t):
            if r in best: return
            if r not in new or t<new[r]: new[r]=t
        for u in ops1:
            for x in levels[k-1]: offer(tuple(U1[u](p) for p in x),'%s(%s)'%(u,best[x][1]))
        for o in ops2:
            for i in range(k):
                for x in levels[i]:
                    for y in levels[k-1-i]: offer(tuple(B2[o](p,q) for p,q in zip(x,y)),'%s(%s,%s)'%(o,best[x][1],best[y][1]))
        for r,t in new.items(): best[r]=(k,t)
        levels.append(list(new))
        if len(best)==27: break
    return best
san=dp(14,['SYN'],['ANT','NOT']); assert len(san)==27
catd={tuple(v):(a,b) for v,a,b in cat}
rows=[]
for v in sorted(catd,key=lambda v:code(v)):
    a,b=catd[v]; c=san[v]
    rows.append(r'\texttt{%s}&$[%s]$&%d&%s&%d&%s&%d&%s\\'%(code(v),','.join(sym[x] for x in v),a[0],tex(a[1]),b[0],tex(b[1]),c[0],tex(c[1])))
UNARY='\n'.join(rows)
mx=lambda i:max(catd[v][i][0] for v in catd)
mx_san=max(c[0] for c in san.values())

# ---- S6 classification table
nm=['SYN','DIF','OR','XOR','AND','NOT','ANT']
def fm(s): return ', '.join('\\'+x for x in s.split())
rows=[]
for i,r in enumerate(cls,1):
    rows.append(r'%d&$%s$&%s&%s&%s&%s\\'%(i,fm(r['subset']),'yes' if r['complete_without_constants']=='1' else 'no',esc(r['certificate_if_incomplete']).replace('rho',r'$\varrho$').replace('{0,+1}',r'$\{0,+1\}$'),'yes' if r['complete_with_constants']=='1' else 'no',esc(r['reason_if_incomplete_with_constants'])))
CLS='\n'.join(rows)

# ---- S7 invariants
F7={**B2,**U1}; ar=lambda n:2 if n in B2 else 1
RHO=[(a,b) for a in D for b in D if a*b>=0]
def pset(n): return all(F7[n](*t) in (0,1) for t in product((0,1),repeat=ar(n)))
def prho(n): return all(F7[n](*[c[0] for c in ts])*F7[n](*[c[1] for c in ts])>=0 for ts in product(RHO,repeat=ar(n)))
rows=[r'$\%s$&%s&%s&%s\\'%(n,'yes' if F7[n](*([0]*ar(n)))==0 else 'no','yes' if pset(n) else 'no','yes' if prho(n) else 'no') for n in nm]
INV='\n'.join(rows)

# ---- S8 involution, lift, orbits
phi=lambda a,b:(SYN(a,b),DIF(a,b))
for a,b in product(D,D): assert phi(*phi(a,b))==(a,b)
PHI=r'\quad'.join([r'\begin{tabular}{cc}$(a,b)$&$\Phi(a,b)$\\\hline '+''.join(r'$(%s,%s)$&$(%s,%s)$\\'%(sym[a],sym[b],sym[phi(a,b)[0]],sym[phi(a,b)[1]]) for a,b in product(D,D))+r'\end{tabular}'])
zc={}
for v in catd: z=sum(1 for x in v if x==0); zc[z]=zc.get(z,0)+1
assert sum(n*2**z for z,n in zc.items())==64
LIFT=', '.join(r'$z=%d$: %d functions, $2^z=%d$ pairs each'%(z,zc[z],2**z) for z in sorted(zc))
orb=res['orbits']; ORB=', '.join('%s of size %s'%(orb['orbit_sizes'][k],k) for k in ('1','2','4','8'))

# ---- S9 C9 example
E=lambda c:'NOT(XOR(x,%s))'%{0:'DIF(x,x)',1:'NOT(DIF(x,x))',-1:'ANT(NOT(DIF(x,x)))'}[c]
env={'SYN':SYN,'DIF':DIF,'OR':OR,'XOR':XOR,'AND':AND,'NOT':NOT,'ANT':ANT}
ev=lambda t,**kw:eval(t,{'__builtins__':{}},{**env,**kw})
sub=lambda t,var:re.sub(r'\ba\b',var,t)
SUMF=lambda a,b:((a+b+1)%3)-1
Pc={c:'AND(%s,%s)'%(E(c).replace('x','a'),sub(catd[tuple(SUMF(c,b) for b in D)][0][1],'b')) for c in D}
term='SYN(%s,SYN(%s,%s))'%(Pc[-1],Pc[0],Pc[1])
for a,b in product(D,D): assert ev(term,a=a,b=b)==SUMF(a,b)
c9rows=[r'%s&$[%s]$&%s&%s\\'%(M(c),','.join(sym[SUMF(c,b)] for b in D),tex(catd[tuple(SUMF(c,b) for b in D)][0][1].replace('a','b')),tex(Pc[c])) for c in D]
C9='\n'.join(c9rows)

# ---- S10 T9ID
T9=[('SYN',SYN),('DIF',DIF),('OR',OR),('XOR',XOR),('AND',AND),('MAX',max),('MIN',min),('SUM (mod 3)',SUMF),('CARRY',lambda a,b:1 if a+b>=2 else (-1 if a+b<=-2 else 0))]
w=res['T9ID_words']
for n,f in T9:
    k={'SUM (mod 3)':'SUM3'}.get(n,n)
    if k in w: assert word(f)==w[k],(n,word(f),w[k])
T9R='\n'.join(r'\textsf{%s}&\texttt{%s}&%s\\'%(n.replace(' (mod 3)',r'\,(mod 3)'),esc(word(f)),' '.join('$%s$'%''.join(sym[f(a,b)] for b in D) for a in D)) for n,f in T9)


# ---- v10 additions: connectives (S3), BID states (S8), direction/rotation/folding (S14-S16)
dr=json.load(open(P('direction_rotation_results.json')))
assert dr['_summary']['status']=='ALL ASSERTIONS PASSED' and dr['_summary']['failed']==[], 'run verify_direction_rotation.py first'
DR_TESTS=dr['_summary']['tests']; DR_PY=dr['_summary']['python']
MAXf=max; MINf=min; MAG=lambda a:a*a
CONN=[('EQ','$[a=b]$','NOT(XOR(a,b))',lambda a,b:iv(a==b)),('LE','$[a\\le b]$','NOT(DIF(a,b))',lambda a,b:iv(a<=b)),('LT','$[a<b]$','NOT(NOT(DIF(b,a)))',lambda a,b:iv(a<b)),
 ('MAX','$\\max(a,b)$','SYN(SYN(a,b),XOR(a,b))',max),('MIN','$\\min(a,b)$','DIF(SYN(a,b),XOR(a,b))',min),
 ('imp$_{K_3}$','$\\max(-a,b)$','SYN(DIF(b,a),OR(a,b))',lambda a,b:max(-a,b)),
 ('imp$_{L_3}$','$\\min(1,1-a+b)$','SYN(NOT(a),SYN(b,AND(a,b)))',lambda a,b:min(1,1-a+b)),
 ('and$_{L_3}$','$\\max(-1,a+b-1)$','DIF(SYN(a,b),NOT(AND(a,b)))',lambda a,b:max(-1,a+b-1)),
 ('or$_{L_3}$','$\\min(1,a+b+1)$','SYN(SYN(a,b),NOT(AND(a,b)))',lambda a,b:min(1,a+b+1)),
 ('CARRY','half-adder carry','SYN(a,DIF(b,SYN(a,b)))',lambda a,b:1 if a+b>=2 else (-1 if a+b<=-2 else 0)),
 ('SUM','$a+b$ mod 3','AND(SYN(a,b),DIF(XOR(a,b),AND(a,b)))',lambda a,b:((a+b+1)%3)-1)]
sizes=res['shortest_terms']['binary_connective_sizes']; kmap={'imp$_{K_3}$':'imp_K3','imp$_{L_3}$':'imp_L3','and$_{L_3}$':'and_L3','or$_{L_3}$':'or_L3'}
crow=[]
for nme,mean,t,f in CONN:
    for a,b in product(D,D): assert ev(t,a=a,b=b)==f(a,b),(nme,a,b)
    k=kmap.get(nme,nme); assert k in sizes, k
    nodes=len(re.findall(r'\b(SYN|DIF|OR|XOR|AND|NOT|ANT)\(',t)); assert nodes==sizes[k],(nme,nodes,sizes[k])
    crow.append(r'%s&%s&%d&%s\\'%(nme if '$' in nme else r'\textsf{%s}'%nme,mean,nodes,tex(t)))
CONNT='\n'.join(crow)
bidrows=r'\\ '.join(r'$(%d,%d)$&%s&%d&%s'%(X,Y,M(X-Y),X+Y,{(0,0):'$S_0$',(0,1):'$S_{-1}$',(1,0):'$S_{+1}$',(1,1):r'$S_0$ (convention)'}[(X,Y)]) for X,Y in ((0,0),(0,1),(1,0),(1,1)))+r'\\'

def Tm(F): return tuple(tuple(F(a,b) for b in D) for a in D)
def negM(Mx): return tuple(tuple(-x for x in r_) for r_ in Mx)
def Rcw(Mx):
    (a,b,c),(d,e,f),(g,h,i)=Mx; return ((g,d,a),(h,e,b),(i,f,c))
def Rccw(Mx):
    (a,b,c),(d,e,f),(g,h,i)=Mx; return ((c,f,i),(b,e,h),(a,d,g))
def mat(Mx,title,w=''):
    rows=''.join(r'%s&%s&%s&%s\\'%(M(D[r]),M(Mx[r][0]),M(Mx[r][1]),M(Mx[r][2])) for r in range(3))
    return r'\begin{tabular}{c|ccc}\multicolumn{4}{c}{%s}\\\hline &%s&%s&%s\\\hline %s\end{tabular}'%(title,M(-1),M(0),M(1),rows)
TS,TD,TO,TX,TA=Tm(SYN),Tm(DIF),Tm(OR),Tm(XOR),Tm(AND)
assert Rcw(TS)==TD and Rcw(TD)==negM(TS) and Rcw(negM(TS))==negM(TD) and Rcw(negM(TD))==TS
assert Rcw(TO)==TX and Rcw(TX)==TO and Rcw(TA)==negM(TA) and Rccw(Rcw(TS))==TS
for f in (SYN,DIF,OR,XOR,AND):
    assert Rcw(Tm(f))==Tm(lambda a,b,f=f:f(-b,a))
CYC=r'$\to$'.join([r'\begin{minipage}[t]{3.1cm}\centering '+mat(m,t)+r'\end{minipage}' for m,t in ((TS,r'$T_{\mathsf{SYN}}$'),(TD,r'$T_{\mathsf{DIF}}=\mathcal R T_{\mathsf{SYN}}$'),(negM(TS),r'$-T_{\mathsf{SYN}}=\mathcal R T_{\mathsf{DIF}}$'),(negM(TD),r'$-T_{\mathsf{DIF}}$'))])
CYC2=r'\quad '.join([r'\begin{minipage}[t]{3.1cm}\centering '+mat(m,t)+r'\end{minipage}' for m,t in ((TO,r'$T_{\mathsf{OR}}$'),(TX,r'$T_{\mathsf{XOR}}=\mathcal R T_{\mathsf{OR}}$'),(TA,r'$T_{\mathsf{AND}}$'),(negM(TA),r'$\mathcal R T_{\mathsf{AND}}=-T_{\mathsf{AND}}$'))])
REFL=r'\quad '.join([r'\begin{minipage}[t]{3.1cm}\centering '+mat(m,t)+r'\end{minipage}' for m,t in ((Tm(lambda a,b:DIF(-b,a)),r'$\mathsf{DIF}(-b,a)$: rotation'),(Tm(lambda a,b:DIF(a,-b)),r'$\mathsf{DIF}(a,-b)$: reflection'),(Tm(lambda a,b:SYN(-b,a)),r'$\mathsf{SYN}(-b,a)$: rotation'),(Tm(lambda a,b:SYN(a,-b)),r'$\mathsf{SYN}(a,-b)$: reflection'))])
d_c4=dr['3.2 orbits of the 19683 binary functions under the rotation group C4: 4995 (Burnside (3^9+2*3^3+3^5)/4)']['detail']
d_d4=dr['3.2 orbits under the dihedral group D4 generated by rotation and swap: 2862 (= v9 Proposition 5)']['detail']
assert sum(d_c4.values())==4995 and sum(d_d4.values())==2862
C4S=', '.join('%s of size %s'%(d_c4[k],k) for k in sorted(d_c4,key=int)); D4S=', '.join('%s of size %s'%(d_d4[k],k) for k in sorted(d_d4,key=int))
# sign-linear functions
lin={}
for al,be in product(range(-12,13),repeat=2):
    v=tuple(sgn(al*a+be*b) for a in D for b in D); lin.setdefault(v,(al,be))
assert len(lin)==17
def wordv(v): s_=''.join({0:'0',1:'1',-1:'2'}[x] for x in v); return '0t'+s_[:3]+'_'+s_[3:6]+'_'+s_[6:]
LINROWS='\n'.join(r'$%s$&\texttt{%s}&%s\\'%('(%d,%d)'%lin[v],esc(wordv(v)),' '.join('$%s$'%''.join(sym[x] for x in v[3*r_:3*r_+3]) for r_ in range(3))) for v in sorted(lin,key=lambda v:(lin[v][1]**2+lin[v][0]**2,lin[v])))
# complex tables
def gm(x,y): return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def gtex(z): 
    r_,i_=z
    if i_==0: return '$%d$'%r_
    if r_==0: return '$%si$'%({1:'',-1:'-'}[i_])
    return '$%d%s i$'%(r_,'+' if i_>0 else '-')
def zm(Zm,title):
    rows=''.join(r'%s&%s&%s&%s\\'%(M(D[r]),gtex(Zm[r][0]),gtex(Zm[r][1]),gtex(Zm[r][2])) for r in range(3))
    return r'\begin{tabular}{c|ccc}\multicolumn{4}{c}{%s}\\\hline &%s&%s&%s\\\hline %s\end{tabular}'%(title,M(-1),M(0),M(1),rows)
Zc=tuple(tuple((SYN(a,b),DIF(a,b)) for b in D) for a in D)
iZc=tuple(tuple(gm((0,1),z) for z in r_) for r_ in Zc)
Zrot=tuple(tuple((Rccw(Tm(SYN))[r][c],Rccw(Tm(DIF))[r][c]) for c in range(3)) for r in range(3))   # R_ccw acts on real and imaginary parts separately
assert Zrot==iZc
for a,b in product(D,D): assert (SYN(a,b),DIF(a,b))==(sgn(a+b),sgn(a-b))==(lambda u:(sgn(u[0]),sgn(u[1])))(gm((1,1),(a,-b)))
ZT=r'\quad '.join([r'\begin{minipage}[t]{3.6cm}\centering '+zm(Zc,r'$Z=T_{\mathsf{SYN}}+iT_{\mathsf{DIF}}$')+r'\end{minipage}',r'\begin{minipage}[t]{3.6cm}\centering '+zm(iZc,r'$iZ=\mathcal R^{-1}Z$')+r'\end{minipage}'])
Wc=tuple(tuple((OR(a,b),XOR(a,b)) for b in D) for a in D)
iWb=tuple(tuple(gm((0,1),(z[0],-z[1])) for z in r_) for r_ in Wc)
Wrot=tuple(tuple((Rcw(Tm(OR))[r][c],Rcw(Tm(XOR))[r][c]) for c in range(3)) for r in range(3))
assert Wrot==iWb
WT=r'\quad '.join([r'\begin{minipage}[t]{3.6cm}\centering '+zm(Wc,r'$W=T_{\mathsf{OR}}+iT_{\mathsf{XOR}}$')+r'\end{minipage}',r'\begin{minipage}[t]{3.6cm}\centering '+zm(iWb,r'$i\overline W=\mathcal R W$')+r'\end{minipage}'])
BIDCD=r'\\ '.join(r'$(%d,%d)$&$%d$&$%d$&%s'%(X,Y,X+Y,X-Y,gtex((X+Y,X-Y))) for X,Y in ((0,0),(0,1),(1,0),(1,1)))+r'\\'
for X,Y in ((0,0),(0,1),(1,0),(1,1)): assert (X+Y,X-Y)==gm((1,1),(X,-Y))
# congruence compatibility of ker(MAG)
OPSd={'SYN':(2,SYN),'DIF':(2,DIF),'OR':(2,OR),'XOR':(2,XOR),'AND':(2,AND),'NOT':(1,NOT),'ANT':(1,ANT)}
kerM={-1:'n',0:'z',1:'n'}
def firstwit(ar,f):
    for xs in product(D,repeat=ar):
        for ys in product(D,repeat=ar):
            if all(kerM[x]==kerM[y] for x,y in zip(xs,ys)) and kerM[f(*xs)]!=kerM[f(*ys)]: return xs,ys
    return None
KR=[]
for n_,(ar_,f_) in OPSd.items():
    w_=firstwit(ar_,f_)
    if w_ is None: KR.append(r'$\%s$&yes&\\'%n_)
    else:
        xs,ys=w_; arg=lambda t:'('+','.join(sym[x] for x in t)+')'
        KR.append(r'$\%s$&no&$\MAG(\%s%s)=%d\ne%d=\MAG(\%s%s)$ although $\MAG$ agrees on the arguments\\'%(n_,n_,arg(xs),f_(*xs)**2,f_(*ys)**2,n_,arg(ys)))
assert [r_.endswith('yes&\\\\') for r_ in KR]==[False,False,False,False,True,False,True]
KRT='\n'.join(KR)
BOOLT=r'\\ '.join(r'$(%d,%d)$&%d&%d&%d&%d&%d&%s&%s'%(x,y,SYN(x,y),OR(x,y),XOR(x,y),AND(x,y),NOT(x),M(DIF(x,y)),M(MAG(DIF(x,y)))) for x,y in product((0,1),(0,1)))+r'\\'
Uu=(1,2,3); Aset={1,2}; Bset_={2,3}
ind_=lambda S_:[1 if u in S_ else 0 for u in Uu]
IA,IB=ind_(Aset),ind_(Bset_)
SETROWS=[('$A$',IA),('$B$',IB),(r'$\mathbf 1_{A\cup B}$ = OR',[OR(x,y) for x,y in zip(IA,IB)]),(r'$\mathbf 1_{A\triangle B}$ = XOR',[XOR(x,y) for x,y in zip(IA,IB)]),(r'$\mathbf 1_{A\cap B}$ = AND',[AND(x,y) for x,y in zip(IA,IB)]),(r'DIF(OR,XOR)',[DIF(OR(x,y),XOR(x,y)) for x,y in zip(IA,IB)])]
assert ind_(Aset|Bset_)==SETROWS[2][1] and ind_(Aset^Bset_)==SETROWS[3][1] and ind_(Aset&Bset_)==SETROWS[4][1]==SETROWS[5][1]
SETT=r'\\ '.join('%s&%s'%(a_,'&'.join(str(x) for x in v_)) for a_,v_ in SETROWS)+r'\\'

# ---- S11 traceability
TR=[('Lemma 1, Theorem 1 (parity, polarization)','chai_checker.py; check_c9_chain.py; chai_verify_v1_increments.py','identities.parity, identities.polarization; v1: parity and polarization at all nine pairs'),
('Theorem 2 (involution)','chai_checker.py; chai_verify_v1_increments.py','identities.involution_fixed_points'),
('Proposition 1 (minimal carrier), the three sections','chai_verify_v1_increments.py','chaiv1_increments_results.json (closure iteration, section vectors, involution test)'),
('Remark 2 (table symmetry), Proposition 5 (fold), S14--S16','verify_direction_rotation.py; chai_verify_v1_increments.py','direction_rotation_results.json (%d tests, negative controls included); v1: clockwise rotation of the five tables, fold congruence'%DR_TESTS),
('Lemma 2, Theorems 3 and 4 (explicit terms)','chai_checker.py; chai_verify_v1_increments.py','paper_terms; DIF_NOT_complete_constructive; Z_n_from_five_operations; DIF_odd; v1: all 6561 functions of Z_2, all 19683 binary functions (Thm 3c), all 81 odd functions from DIF alone'),
('Theorem 5 (Boolean layer), Table III','chai_verify_v1_increments.py','v1: all 3888 functions of B_2 and 300 random functions of B_3 rebuilt; class sizes and the printed numbers of Table III'),
('Table II (layers)','chai_verify_v1_increments.py','v1: parity, range, zero-fixing, B-keeping and rho-keeping computed for the seven primitives and compared with the printed rows'),
('Theorem 3(b): $\\Z_n$ from five operations','chai_checker.py --full','closure SYN+DIF+OR+XOR+AND (6561 functions)'),
('Theorem 6: 40 of 127','chai_checker.py; classification_table.py; independent_free_constants.py; chai_verify_v1_increments.py','classification; classification_127.csv (per-subset, two methods plus an independent script); v1: the criterion agrees with the csv'),
('Remark 3: 95 with constants','chai_checker.py; results_consts_*.json; verify_direction_rotation.py','free_constants; closure with free constants *; direction_rotation_results.json (terms, invariants, 127-subset criterion)'),
('Remark 1 (sign rule of AND)','verify_direction_rotation.py; chai_verify_v1_increments.py','direction_rotation_results.json; failure at (2,1) on {-2..2}'),
('Connectives as terms (S3), Proposition 3 ($K_3$)','chai_checker.py; chai_verify_v1_increments.py; this supplement','paper_terms; shortest_terms; every term evaluated by this generator; v1: MAX, MIN terms'),
('27 unary terms (S5)','chai_checker.py; independent_check.py','unary_catalogue.json; catalogue mismatches 0'),
('Difference lift (Prop. 6)','chai_checker.py; check_c9_chain.py; chai_verify_v1_increments.py','difference_lift; n=2 pairs enumerated'),
('2862 orbits (Prop. 7)','chai_checker.py; verify_direction_rotation.py; chai_verify_v1_increments.py','orbits; direction_rotation_results.json; v1: group of order 8, Burnside sum'),
('T9ID words','chai_checker.py; chai_verify_v1_increments.py; this supplement','T9ID_words'),
('Table IV (closures with the constants), Remark 3','chai_verify_v1_increments.py; chai_checker.py --full-consts','results_consts_*.json (size and heights of the five closures); v1: the printed table against these files'),
('XY3 and G3 tables, hexaogram and bi-trigram symbols (S17)','verify_direction_rotation.py; this supplement','direction_rotation_results.json (S17 tests); every row asserted by the generator'),
('C9 three-slice selector (19683 functions)','check_c9_chain.py; chai_verify_v1_increments.py','all binary functions rebuilt as $U-V$ and as the term of Section VI; v1: Shannon/Post expansion on 500 random functions')]
for k in ['identities','paper_terms','DIF_NOT_complete_constructive','Z_n_from_five_operations','DIF_odd','classification','free_constants','shortest_terms','difference_lift','orbits','T9ID_words']: assert k in res,k
TRR='\n'.join(r'%s&\texttt{%s}&%s\\'%(a,esc(b).replace('; ','; \\allowbreak '),esc(c).replace('$\\Z_n$','')) for a,b,c in TR)


SEC1416=r'''
\section*{S14. The quarter turn of a truth table}
\textbf{Convention.} The matrix $T_F$ of a binary $F$ has rows $a=-,0,+$ from top to bottom and columns $b=-,0,+$ from left to right. The clockwise quarter turn is
\[\mathcal R\begin{pmatrix}p&q&r\\s&t&u\\v&w&x\end{pmatrix}=\begin{pmatrix}v&s&p\\w&t&q\\x&u&r\end{pmatrix},\]
$\mathcal R^{-1}$ is the counterclockwise turn, and $-T$ negates every entry. Put $\rho(a,b)=(-b,a)$ (the quarter turn of the argument grid; $a+ib\mapsto i(a+ib)$).

\textbf{Proposition S14 (quarter turn).} For every $F:\D^2\to\D$, $\mathcal R(T_F)=T_{F\circ\rho}$. Hence $\mathcal RT_{\mathsf{SYN}}=T_{\mathsf{DIF}}$, $\mathcal RT_{\mathsf{DIF}}=-T_{\mathsf{SYN}}$, $\mathcal RT_{\mathsf{OR}}=T_{\mathsf{XOR}}$, $\mathcal RT_{\mathsf{XOR}}=T_{\mathsf{OR}}$, $\mathcal RT_{\mathsf{AND}}=-T_{\mathsf{AND}}$, $\mathcal R^4=\mathrm{id}$, and $\mathcal R^2T_F=-T_F$ iff $F$ is odd ($\mathcal R^2T_F=T_F$ iff $F$ is even).

\emph{Proof.} The entry of $\mathcal R(M)$ in row $a$ and column $b$ is the entry of $M$ in the row labelled $-b$ and the column labelled $a$; for $M=T_F$ this is $F(-b,a)$. Then $\mathsf{SYN}(-b,a)=\mathrm{sgn}(a-b)$, $\mathsf{DIF}(-b,a)=\mathrm{sgn}(-b-a)$, $\mathsf{OR}(-b,a)=\mathrm{sgn}((a-b)^2)$, $\mathsf{XOR}(-b,a)=\mathrm{sgn}((a+b)^2)$, $\mathsf{AND}(-b,a)=\mathrm{sgn}(-4ab)$. Since $\rho^2(a,b)=(-a,-b)$, $\mathcal R^2T_F=T_{F(-a,-b)}$. \hfill$\square$

The generator rotates the matrices below and asserts every equality; the program \texttt{verify\_direction\_rotation.py} checks the lemma $\mathcal R(T_F)=T_{F\circ\rho}$ for all $19\,683$ binary $F$.
\begin{center}\small '''+CYC+r'''\end{center}
\begin{center}\small '''+CYC2+r'''\end{center}
The identities of Proposition S14 for the sgn-defined operations hold on every symmetric grid $\{-m,\dots,m\}$ (checked for $m=2,3,5$; they follow from the same one-line computation), whereas the polarization identity of Theorem 1 does not (Remark 1): only the descent to signs is special to $\D$. The even operations $\mathsf{OR},\mathsf{XOR}$ have period $2$ under $\mathcal R$ and the odd operations $\mathsf{SYN},\mathsf{DIF}$ period $4$, because $\mathcal R^2$ is the substitution $(a,b)\mapsto(-a,-b)$. Among the $19\,683$ binary functions exactly $81=3^4$ are odd and $243=3^5$ are even.

\textbf{Three notions to keep apart.} (1)~turning the picture of the table, (2)~substituting $\rho$ in the arguments (equal to (1) by Proposition S14), (3)~identities such as $\mathsf{DIF}(a,b)=\mathsf{SYN}(a,\mathsf{ANT}(b))$, which negate one argument and \emph{reflect} the table (reverse the columns). For commutative $F$ the substitution $(-b,a)$ and $(a,-b)$ coincide, which is why rotation and reflection agree on $\mathsf{SYN}$ and $\mathsf{OR}$; $\mathsf{DIF}$ is antisymmetric, so rotation gives $-\mathsf{SYN}$ and reflection gives $+\mathsf{SYN}$:
\begin{center}\small '''+REFL+r'''\end{center}
\textbf{Orbits.} Under the rotation group $C_4=\langle\rho\rangle$ the $19\,683$ binary functions fall into $4995=\tfrac14(3^9+2\cdot3^3+3^5)$ orbits ('''+C4S+r'''); under the dihedral group $\langle\rho,\mathsf{SWAP}\rangle$ of order $8$ (the group of S8) into $2862$ orbits ('''+D4S+r'''), as in Section VI of the paper. The orbit of $\mathsf{DIF}$ is $\{\mathsf{DIF},\mathsf{SYN},-\mathsf{DIF},-\mathsf{SYN}\}$, that of $\mathsf{OR}$ is $\{\mathsf{OR},\mathsf{XOR}\}$ and that of $\mathsf{AND}$ is $\{\mathsf{AND},-\mathsf{AND}\}$ (all enumerated).

\textbf{Sign-linear binary functions (degree one).} The functions $\mathrm{sgn}(\alpha a+\beta b)$, $\alpha,\beta\in\mathbb Z$, are exactly $17$: the value at the eight non-zero points of the grid is $\pm\mathrm{sgn}\,\alpha$, $\pm\mathrm{sgn}\,\beta$, $\pm\mathrm{sgn}(\alpha+\beta)$, $\pm\mathrm{sgn}(\alpha-\beta)$, so a function is determined by the signs of these four numbers, and the directions of $(\alpha,\beta)$ form $8$ open sectors, $8$ rays and the origin. They form dihedral orbits of sizes $8$ (sectors), $4$ ($\mathsf{SYN},\mathsf{DIF},-\mathsf{SYN},-\mathsf{DIF}$), $4$ ($\pm a,\pm b$) and $1$ (zero). All are odd, hence lie in the clone generated by $\mathsf{DIF}$ (Theorem 4). The program enumerates $|\alpha|,|\beta|\le12$ and finds $17$ functions; the count for all integers follows from the sign argument. This is an elementary observation, not a result used in the paper.
\begin{center}\scriptsize\begin{tabular}{lcc}\toprule $(\alpha,\beta)$&T9ID word&rows $a=-,0,+$\\\midrule
'''+LINROWS+r'''
\bottomrule\end{tabular}\end{center}

\section*{S15. The complex table $Z=T_{\mathsf{SYN}}+i\,T_{\mathsf{DIF}}$}
Let $\mathrm{sgn}_c(x+iy)=\mathrm{sgn}\,x+i\,\mathrm{sgn}\,y$ on Gaussian integers, and let $\mathcal R$ act on a complex matrix by acting on its real and imaginary parts (the linear extension).

\textbf{Statement (paper, Section III-C).} $\mathcal R^{-1}Z=iZ$ and $\mathcal RZ=-iZ$; hence $\mathcal R^4Z=Z$ and $\mathcal R^2Z=-Z$.

\emph{Proof.} (a)~$\mathrm{sgn}_c(iu)=i\,\mathrm{sgn}_c(u)$: multiplication by $i$ sends $(x,y)$ to $(-y,x)$ and $\mathrm{sgn}$ is odd; applying it to $-iu$ gives $\mathrm{sgn}_c(-iu)=-i\,\mathrm{sgn}_c(u)$. (b)~$Z(a,b)=\mathrm{sgn}_c\big((a+b)+i(a-b)\big)=\mathrm{sgn}_c\big((1+i)\overline{(a+ib)}\big)$, because $(1+i)(a-ib)=(a+b)+i(a-b)$. (c)~With $w=a+ib$, $\rho$ is $w\mapsto iw$ and $\overline{iw}=-i\bar w$, so $Z(\rho(a,b))=\mathrm{sgn}_c\big(-i(1+i)\bar w\big)=-iZ(a,b)$. By Proposition S14 of the paper applied to the real and imaginary parts, $\mathcal RZ=Z\circ\rho=-iZ$, and $\mathcal R^{-1}Z=iZ$ follows. \hfill$\square$

The program checks (a) for all Gaussian integers with $|x|,|y|\le6$, (b) at the nine points, (c) for all $w$ with $|a|,|b|\le6$ (the formula of (b) extended to all integers), and the proposition entry by entry in exact integer arithmetic; wrong-sign versions ($\mathcal RZ=iZ$, $\mathcal R^{-1}Z=-iZ$) are checked to be \emph{false}.
\begin{center}\small '''+ZT+r'''\end{center}
\textbf{BID view.} At $(a,b)=(X,Y)\in\{0,1\}^2$ the unsaturated values are $c+id=(1+i)\overline{(X+iY)}$ with $c=X+Y$, $d=X-Y$ (the sum and difference views of a BID), and $Z(X,Y)=\mathrm{sgn}\,c+i\,\mathrm{sgn}\,d$ is their saturation. The identity $c^2-d^2=4XY$ is the BID form of the polarization $\mathsf{AND}=\mathsf{DIF}(\mathsf{OR},\mathsf{XOR})$.
\begin{center}\small\begin{tabular}{cccc}\toprule $(X,Y)$&$c$&$d$&$c+id$\\\midrule
'''+BIDCD+r'''
\bottomrule\end{tabular}\end{center}
\textbf{After folding.} The folded table $W=T_{\mathsf{OR}}+i\,T_{\mathsf{XOR}}$ satisfies $\mathcal RW=\mathcal R^{-1}W=i\overline W$: the turn acts as a reflection (order $2$) and not as multiplication by $i$.
\begin{center}\small '''+WT+r'''\end{center}
\textbf{What this does and does not say.} It says that the two real tables are the two components of one complex-valued table and that the quarter turn of the table acts on it by multiplication by $-i$ (clockwise) or $i$ (counterclockwise), because $\mathrm{sgn}$ commutes with the dihedral symmetry of $\mathbb Z^2$. It does not make $i$ the cause of the algebra, and it does not give the algebra an interpretation by complex numbers.

\section*{S16. Folding, congruences, Boolean restriction and sets}
\textbf{Proposition S16(ii): the kernel of $\mathsf{MAG}$.} Let $\theta=\{\{-1,+1\},\{0\}\}$ (an equivalence relation, not a congruence of the algebra) be the kernel of $\mathsf{MAG}(a)=a^2$. An operation respects $\theta$ if equivalent arguments give equivalent values. Checked over all pairs of $\theta$-equivalent argument tuples:
\begin{center}\small\begin{tabular}{lcp{10cm}}\toprule operation&respects $\theta$&first violation\\\midrule
'''+KRT+r'''
\bottomrule\end{tabular}\end{center}
So $\mathsf{MAG}$ is a homomorphism from the reduct $(\D;\mathsf{AND},\mathsf{ANT})$ onto $(\{0,1\};\min,\mathrm{id})$, and from no reduct containing another of the seven operations ($\mathsf{MAG}(\mathsf{AND}(a,b))=\min(\mathsf{MAG}\,a,\mathsf{MAG}\,b)$, $\mathsf{MAG}(\mathsf{ANT}\,a)=\mathsf{MAG}\,a$).

\textbf{Proposition S16(iii): no nontrivial congruence, proper nonempty subalgebra or nontrivial automorphism.} Constants $0,+1,-1$ and all $27$ unary functions are terms (Theorem 6 and its proof; S5). \emph{Congruences:} let $\theta'$ identify $x\neq y$. A unary term $u$ with $u(x)=0$, $u(y)=+1$ gives $0\ \theta'\ {+1}$; for any $p,q$ a unary term $v$ with $v(0)=p$, $v(+1)=q$ gives $p\ \theta'\ q$. \emph{Subalgebras:} a subalgebra containing $s$ contains $c(s)$ for every constant unary term $c$, hence all of $\D$. \emph{Automorphisms:} an automorphism $\varphi$ satisfies $\varphi(c(s))=c(\varphi(s))$ for the constant unary terms, hence fixes every element. The program confirms by exhaustion: of the five partitions of $\D$ only the two trivial ones are congruences, no nonempty proper subset is closed, and only the identity of the six permutations is an automorphism. (The same facts hold for any primal algebra; we prove them directly.)

The reducts behave differently: $\{0,+1\}$ is a subalgebra of $(\D;\mathsf{SYN},\mathsf{OR},\mathsf{XOR},\mathsf{AND},\mathsf{NOT})$, and not of the full algebra since $\mathsf{DIF}(0,+1)=-1$ and $\mathsf{ANT}(+1)=-1$.

\textbf{Boolean restriction.} On $\{0,1\}$, $\mathsf{SYN}=\mathsf{OR}=$ Boolean or, $\mathsf{XOR}=$ Boolean exclusive or, $\mathsf{AND}$ is Boolean and, $\mathsf{NOT}(x)=1-x$; $\mathsf{DIF}$ leaves the carrier. Two-valued inputs thus produce a signed result outside the carrier, and $\mathsf{MAG}\circ\mathsf{DIF}=\mathsf{XOR}$ folds the two directions into one value. This does not mean that Boolean logic contains $-1$.
\begin{center}\small\begin{tabular}{cccccccc}\toprule $(x,y)$&SYN&OR&XOR&AND&NOT$(x)$&DIF&MAG(DIF)\\\midrule
'''+BOOLT+r'''
\bottomrule\end{tabular}\end{center}
\textbf{Sets.} For subsets $A,B$ of a universe, the characteristic functions satisfy $\mathbf 1_{A\cup B}=\mathbf 1_A\vee\mathbf 1_B$, $\mathbf 1_{A\triangle B}=\mathbf 1_A\oplus\mathbf 1_B$, $\mathbf 1_{A\cap B}=\mathbf 1_A\wedge\mathbf 1_B$ and $\mathbf 1_{A\cap B}=\mathbf 1_{A\cup B}-\mathbf 1_{A\triangle B}$ (inclusion--exclusion); with $\mathsf{OR}=\mathsf{SYN}$, $\mathsf{XOR}$, $\mathsf{AND}$ on $\{0,1\}$ this is the polarization identity of Theorem 1, which extends it to signed values. Example ($U=\{1,2,3\}$, $A=\{1,2\}$, $B=\{2,3\}$; columns are the points $1,2,3$; all $64$ pairs of subsets of a $3$-set are checked by the program):
\begin{center}\small\begin{tabular}{lccc}\toprule &1&2&3\\\midrule
'''+SETT+r'''
\bottomrule\end{tabular}\end{center}
The areas of regions in a Venn diagram are not defined here and are not identified with values of the algebra.

\textbf{Sign rule of $\mathsf{AND}$ (Remark 1 of the paper).} For all integers, $4ab=(a+b)^2-(a-b)^2$, hence $\mathrm{sgn}(ab)=\mathrm{sgn}(|a+b|-|a-b|)$: the signs of a product are read off from whether the sum or the difference is the larger in size (equal signs make $|a+b|$ larger). This is the classical quarter-square identity. On $\D$ the sizes are $0,1,2$ and by Lemma 1 two different nonzero squares never occur, so the squares can be folded to signs first, which gives $\mathsf{AND}=\mathsf{DIF}(\mathsf{OR},\mathsf{XOR})$ (Theorem 1). On $\{-2,\dots,2\}$ the folded term fails (first failure $(a,b)=(2,1)$: $p^2-q^2=8>0$, but $\mathsf{DIF}(1,1)=0$), while the size comparison itself still holds. The program checks both statements. So the sign rules $-\cdot-=+$, $+\cdot-=-$ are consequences of the tables of $\mathsf{SYN}$, $\mathsf{DIF}$ and $\mathsf{MAG}$ on $\D$, but the statement ``multiplication is nothing but folded sums and differences'' is not claimed for other alphabets, and the wording about ``aligned axes'' is an interpretation.

\textbf{Status of the readings.}
\begin{center}\small\begin{longtable}{>{\raggedright\arraybackslash}p{6cm}>{\raggedright\arraybackslash}p{9.8cm}}\toprule statement&status\\\midrule\endhead
$\mathcal R(T_F)=T_{F\circ\rho}$, the table cycle, $\mathcal RZ=-iZ$&proved (Prop.~S14 and S15); finite instances checked by program.\\
$\mathsf{OR}=\mathsf{MAG}\circ\mathsf{SYN}$, $\mathsf{XOR}=\mathsf{MAG}\circ\mathsf{DIF}$; $\ker\mathsf{MAG}$ respected only by $\mathsf{AND},\mathsf{ANT}$; no nontrivial congruence&proved by evaluation at nine points and a short argument (Prop.~S16); checked by program.\\
Boolean restriction to $\{0,1\}$; set identities&classical facts expressed in this signature; checked by program.\\
$\mathrm{sgn}(ab)=\mathrm{sgn}(|a+b|-|a-b|)$; $\mathsf{AND}=\mathsf{DIF}(\mathsf{OR},\mathsf{XOR})$ on $\D$ (Remark 1, Thm.~1)&proved; checked by program.\\
Coin model: $-1,+1$ two faces, $0$ the edge; $\mathsf{ANT}$ turns the coin, $\mathsf{MAG}$ forgets the face; ``multiplication emerges from aligned axes''&interpretation, not a theorem. An analogy; a zero need not mean an edge in other settings.\\
The complex unit $i$ as the origin of the algebra; the algebra as a foundation of geometry, number theory, vectors or mass&\emph{not claimed.} Would need a defined space, distance, area and a proof; a research programme at most.\\
Sign-degree filtration, equational axiomatization, lattice of the clones generated by the 127 subsets&open.\\
\bottomrule\end{longtable}\end{center}
'''


# ---- S17: XY3 and G3 tables, generated from their definitions
STATES17=[(1,0),(0,0),(0,1)]                       # bit 0, 1, 2  <->  d = X-Y = +1, 0, -1
xy3f=lambda m,x,y:(m>>STATES17.index((x,y)))&1
XY3N={0:'0',1:'x',2:"x'y'",3:"y'",4:'y',5:r'x\vee y',6:"x'",7:'1'}
XY3F={0:lambda x,y:0,1:lambda x,y:x,2:lambda x,y:1-(x|y),3:lambda x,y:1-y,4:lambda x,y:y,5:lambda x,y:x|y,6:lambda x,y:1-x,7:lambda x,y:1}
for m_ in range(8):
    for (x_,y_) in STATES17: assert xy3f(m_,x_,y_)==XY3F[m_](x_,y_),(m_,x_,y_)
TRIN={0:'Kun',1:'Zhen',2:'Kan',3:'Dui',4:'Gen',5:'Li',6:'Xun',7:'Qian'}
def tri17(m_):
    """trigram of a mask: top line = bit 2, middle = bit 1, bottom = bit 0; solid (yang) = 1, broken (yin) = 0"""
    return '\\tri{%s}{%s}{%s}'%tuple('\\yang' if (m_>>b_)&1 else '\\yin' for b_ in (2,1,0))
def grid17(A_,a_):
    """Hexaogram (hexa-o-gram, o = orthogonal): the three bits of A are three HORIZONTAL lines (top to bottom: bit 2,1,0), the three bits of a
    are three VERTICAL lines (left to right: bit 2,1,0); bit 1 draws the line, bit 0 leaves it out.  Lines run from the first to the last
    of the nine crossing points and never beyond them; the crossings are always drawn as small dots, so that the position of an absent
    line stays readable.  Size: 0.41em x 0.40em per unit, line 0.6pt, dot radius 0.04em, the same baseline for all 27 symbols."""
    t='\\begin{tikzpicture}[x=0.41em,y=0.40em,baseline=0.25em]'
    for x_ in range(3):
        for y_ in range(3): t+='\\fill[black!75] (%d,%d) circle (0.04em);'%(x_,y_)
    for r_ in range(3):
        if (A_>>(2-r_))&1: t+='\\draw[line width=0.6pt,line cap=round] (0,%d)--(2,%d);'%(2-r_,2-r_)
    for c_ in range(3):
        if (a_>>(2-c_))&1: t+='\\draw[line width=0.6pt,line cap=round] (%d,0)--(%d,2);'%(c_,c_)
    return t+'\\end{tikzpicture}'
def bitri17(A_,a_):
    """Bi-trigram: the trigram of A and the trigram of a side by side"""
    return tri17(A_)+'\\ \\hspace{0.05em}'+tri17(a_)
XY3ROWS='\n'.join(r'$%d$&\texttt{0b%s}&%s \ %s&$%s$&%s&%s&%s\\'%(m_,format(m_,'03b'),tri17(m_),TRIN[m_],XY3N[m_],M(xy3f(m_,0,1)),M(xy3f(m_,0,0)),M(xy3f(m_,1,0))) for m_ in (0,1,2,3,4,5,6,7))
def trits17(g_):
    r_=[]
    for _ in range(3):
        t_=((g_+1)%3)-1; r_.append(t_); g_=(g_-t_)//3
    return tuple(reversed(r_))
def g3f(A_,a_):
    nA,na=XY3N[A_],XY3N[a_]
    if a_==0: return '$%s$'%nA
    if A_==0: return '$-%s$'%(('(%s)'%na) if '\\vee' in na else na)
    pa=('(%s)'%na) if '\\vee' in na else na
    return '$%s-%s$'%(nA,pa)
G3ROWS=[]
for gid_ in range(-13,14):
    t_=trits17(gid_); A_=sum(1<<(2-i) for i,v in enumerate(t_) if v==1); a_=sum(1<<(2-i) for i,v in enumerate(t_) if v==-1)
    gv=tuple(xy3f(A_,x,y)-xy3f(a_,x,y) for (x,y) in reversed(STATES17))     # g(-1),g(0),g(+1)
    assert A_&a_==0 and gv in catd, (gid_,gv)
    assert tuple({0:0,1:1,-1:2}[x] for x in gv)==tuple({0:0,1:1,-1:2}[x] for x in t_)
    size_,term_=catd[gv][0]
    for d_ in D: assert ev(term_,a=d_)==gv[d_+1],(gid_,d_)
    G3ROWS.append(r'$%d$&\texttt{%s}&(%d,%d)&%s&%s&%s&$[%s]$&%d&%s\\'%(gid_,esc(code(gv)),A_,a_,grid17(A_,a_),bitri17(A_,a_),g3f(A_,a_),','.join(sym[x] for x in gv),size_,tex(term_)))
assert len({r_ for r_ in G3ROWS})==27
G3T='\n'.join(G3ROWS)
SEC17=r"""
\section*{S17. The XY3 and G3 tables of the published divisional}
The eight XY3 masks are labeled by the eight trigrams (top line = bit 2, bottom line = bit 0, solid line = 1), following the labeling of the author's filings; each of the 27 G3 values is drawn below as a hexaogram and as a bi-trigram. The published divisional CN~122173050~A describes three configuration layers, XY3 (eight functions), G3 (27 functions) and C9 (a selector of three G3 branches). This section reproduces the first two from their definitions and places them in the algebra; nothing in the paper depends on it, and the notation $x,y$ is the pair of bits $(X,Y)$ of a BID.

\textbf{Definitions used.} The three BID states are $(X,Y)=(1,0),(0,0),(0,1)$, i.e.\ $d=X-Y=+1,0,-1$ (the state $(1,1)$ is sent to $d=0$ by the convention of S8). An XY3 function is a three-bit mask, bit $0$ for $(1,0)$, bit $1$ for $(0,0)$ and bit $2$ for $(0,1)$, read as a Boolean function of $x,y$ on these three states. A G3 function is the difference $\mathrm{XY3}(A)-\mathrm{XY3}(a)$ of two masks $A,a$ with no common one ($A\wedge a=0$); there are exactly $3^3=27$ such pairs, and they give the $27$ distinct unary functions of $\D$. The identifier $\mathrm{gid}\in\{-13,\dots,13\}$ is the balanced-ternary number whose three digits (most significant first; $0,1,2$ standing for $0,+1,-1$) are the values $g(-1),g(0),g(+1)$ at the states $d=-1,0,+1$. Hence the T3 code of $\mathrm{gid}$ is exactly the code of S5, and every G3 branch is a \chai{} term (last column, the minimal term over all seven operations of S5, variable $a=d$). The generator checks that every row below is computed from these definitions, that the $27$ functions are distinct, and that every term evaluates to its function at all three points; \texttt{verify\_direction\_rotation.py} checks the same statements independently.

\begin{center}\begin{tikzpicture}[x=9mm,y=9mm,>=stealth]
\foreach \x in {0,1,2}{\foreach \y in {0,1,2}{\fill[black!75] (\x,\y) circle (0.1);}}
\draw[line width=1.4pt,line cap=round] (0,0)--(2,0);
\draw[line width=1.4pt,line cap=round] (0,0)--(0,2);
\node[left] at (-0.6,2) {\scriptsize bit 2 ($d=-1$)};\node[left] at (-0.6,1) {\scriptsize bit 1 ($d=0$)};\node[left] at (-0.6,0) {\scriptsize bit 0 ($d=+1$)};
\node[above] at (0,2.6) {\scriptsize bit 2};\node[above] at (1,2.6) {\scriptsize bit 1};\node[above] at (2,2.6) {\scriptsize bit 0};
\node[right,align=left] at (3.3,1.7) {\small horizontal lines: the mask $A$ (here $A=1$, only bit 0)};
\node[right,align=left] at (3.3,1.0) {\small vertical lines: the mask $a$ (here $a=4$, only bit 2)};
\node[right,align=left] at (3.3,0.3) {\small dots: the nine crossings, always drawn};
\node[right,align=left] at (3.3,-0.4) {\small $g=\mathrm{XY3}(1)-\mathrm{XY3}(4)=x-y=[-,0,+]$, gid $-8$};
\end{tikzpicture}\\[2pt]{\small Figure S17. The hexaogram (hexa-o-gram, ``o'' for orthogonal) of a G3 function, here $(A,a)=(1,4)$: the bits of $A$ are horizontal lines, the bits of $a$ vertical lines.}\end{center}

\textbf{XY3.} Values are those of the mask at the states $(0,1),(0,0),(1,0)$ ($d=-1,0,+1$). The Boolean formulas $x,y,\dots$ agree with the masks on the three states only; at $(1,1)$ they take other values (for example $x\vee y=1$ there, whereas the zon state has the value $0$), which is why the domain is restricted to the three states or $(1,1)$ is sent to the zon state.
\begin{center}\small\renewcommand{\arraystretch}{1.35}\begin{tabular}{rlllccc}\toprule id&mask&trigram&function&$d=-1$&$d=0$&$d=+1$\\\midrule
"""+XY3ROWS+r"""
\bottomrule\end{tabular}\end{center}

\textbf{G3.} Pairs $(A,a)$ are the XY3 ids of the two masks; $g=[g(-1),g(0),g(+1)]$. The \emph{hexaogram} column draws the two masks superposed, as a $3\times3$ grid of three horizontal and three vertical lines, six lines in all (Figure S17): the three bits of $A$ are the horizontal lines, the three bits of $a$ the vertical lines, and a line is present exactly where its bit is $1$; the nine crossing points are always marked by small dots, so that the position of an absent line can be read. The \emph{bi-trigram} column writes the trigram of $A$ and the trigram of $a$ side by side (top line = bit 2). Both symbols are generated from the pair $(A,a)$ by the same program that builds the table, so symbol, code and function cannot disagree.
\begin{center}\scriptsize\renewcommand{\arraystretch}{1.5}\begin{longtable}{rllccccccp{4.1cm}}\toprule gid&T3 code&$(A,a)$&Hexaogram $(A,a)$&Bi-trigram $(A,a)$&$\mathrm{XY3}(A)-\mathrm{XY3}(a)$&$g$&size&\chai{} term\\\midrule\endhead
"""+G3T+r"""
\bottomrule\end{longtable}\end{center}
\textbf{Selectors.} On the three states the functions $x$, $x'y'$, $y$ (the masks $1$, $2$, $4$) are one-hot. At $(1,1)$ both $x$ and $y$ equal $1$ and $x'y'=0$, so a sum of the three products with these Boolean selectors would add two branches there; the total selector of the paper (S8) sends $(1,1)$ to the zon branch instead. The C9 layer is the three-branch selection of S9.
"""
import json as _json
_v1=_json.load(open(P('chaiv1_increments_results.json')))
assert _v1['status']=='ALL CHECKS PASSED'
SEC18=r'''
\section*{S18. Checks added for this version}
The standard-library program \texttt{chai\_verify\_v1\_increments.py} (%d checks, all passed; results in \texttt{chaiv1\_increments\_results.json}) evaluates, from the primitive definitions, the statements that are new in this version of the paper, and the explicit constructions of its proofs on \emph{every} binary function of the classes, not on samples:
\begin{itemize}\setlength\itemsep{0pt}
\item Section II: the three sections $\NOT=[+,+,0]$, $\ANT=[+,0,-]$, $\DIF(-1,\cdot)=[0,-,-]=-\NOT(-\cdot)$; that only the section at $0$ is an involution; $\NOT\circ\NOT=[0,0,+]$; the least $\DIF$-closed set containing $\{0,1\}$ is $\{-1,0,+1\}$ (iteration); $\DIF(-1,+1)=-1\ne-2$.
\item The printed Tables I, II and III against the definitions: all five binary tables and the three unary rows cell by cell; the five columns of Table II for the seven primitives; the numbers $3,81,12,3888,9,6561,27,19\,683$ of Table III.
\item Lemma 2 and Theorems 3--5: $\delta_p$ for the eight $p\ne0$; all $6561$ functions of $\Z_2$ from the five operations; all $19\,683$ binary functions as $\SYN\big(F',\AND(e_0,F(0))\big)$ (Theorem 3(c)); all $81$ odd binary functions from $\DIF$ alone with the terms $k_c$ and $s_p$; all $3888$ functions of $\mathcal B_2$ and $300$ random functions of $\mathcal B_3$ from $\SYN,\OR,\XOR,\AND,\NOT$.
\item Theorem 6: the criterion ``$\NOT\in\Sigma$ and ($\DIF\in\Sigma$ or $\{\SYN,\ANT\}\subseteq\Sigma$)'' agrees with all $127$ rows of \texttt{classification\_127.csv}, and $40$ rows are complete; the relation $\rho$ is preserved by $\ANT,\OR,\XOR,\AND,\NOT$ and not by $\SYN$ or $\DIF$; the witnesses of the necessity proof, $\DIF(0,+1)=-1\notin B$ and $(\DIF(+1,0),\DIF(0,+1))=(+1,-1)\notin\rho$; the five closures with constants of Table~IV of the paper against the shipped \texttt{results\_consts\_*.json} (size $19\,683$ and last height).
\item Sections III and VI: parity lemma and polarization at the nine pairs, the failure at $(2,1)$ on $\{-2,\dots,2\}$, $\Phi$ an involution with one fixed point, the congruence of the fold ($\AND$ and $\ANT$ only), the clockwise rotation of the five tables, $\MAX,\MIN$ as terms, the associativity failures $4,16,8,8$, the Shannon/Post expansion on $500$ random binary functions, the T9ID words of $\SYN$ and $\DIF$, the order-$8$ group and the Burnside count $2862$.
\end{itemize}
The large closures (the $6561$, $3888$ and $19\,683$ binary functions as closures of the projections) and the closures with nullary constants are not recomputed here; they are those of the earlier verified runs (\texttt{results.json}, \texttt{results\_consts\_*.json}) for the same seven operations.
'''%_v1['n_checks']

SEC19=r"""
\section*{S19. Structure of \chai{}: the diagram of the earlier version}
Figure S19 is the structure diagram of the earlier version of the paper, kept here because it shows the relations among the seven operations at a glance. Each arrow is a statement proved in the paper. Arrows from $(a,b)$ to $\SYN$ and $\DIF$: the pair $\Phi=(\SYN,\DIF)$ is an involution (Theorem~2). Vertical arrows labelled $\MAG$: folding the sum or the difference to its magnitude gives $\OR=\MAG\circ\SYN$ and $\XOR=\MAG\circ\DIF$. Arrows from $\OR$ and $\XOR$ to $\AND$: the polarization $\AND=\DIF(\OR,\XOR)$ (Theorem~1). Dashed double arrows: replacing the second argument $b$ by $\ANT(b)$ exchanges $\SYN\leftrightarrow\DIF$ and $\OR\leftrightarrow\XOR$, since $\SYN(a,b)=\DIF(a,\ANT(b))$ and $\OR(a,\ANT(b))=\XOR(a,b)$. The dotted rectangle is the zero-fixing clone $\Z$ ($f(0,\dots,0)=0$): the five binary operations stay inside it (Theorem~3a,b), while $\NOT$ is zero-breaking and, together with the five, gives all of $\Fone$ and $\Ftwo$ (Theorem~3c).
\begin{center}
\resizebox{0.8\textwidth}{!}{%
\begin{tikzpicture}[>=Stealth,font=\small,
box/.style={draw,rounded corners=2pt,inner sep=3pt,minimum width=1.0cm},
lab/.style={font=\scriptsize,inner sep=1.5pt}]
\node[box] (ab) at (2.2,4.1) {$(a,b)$};
\node[box] (syn) at (0,2.7) {$\SYN$};
\node[box] (dif) at (4.4,2.7) {$\DIF$};
\node[box] (or) at (0,1.2) {$\OR$};
\node[box] (xor) at (4.4,1.2) {$\XOR$};
\node[box] (and) at (2.2,0) {$\AND$};
\draw[->] (ab)--(syn); \draw[->] (ab)--(dif);
\node[lab] at (2.2,3.35) {$\Phi$ (involution)};
\draw[->] (syn)--node[lab,left]{$\MAG$}(or);
\draw[->] (dif)--node[lab,right]{$\MAG$}(xor);
\draw[->] (or)--(and); \draw[->] (xor)--(and);
\node[lab,align=center] at (2.2,-0.62) {$\AND=\DIF(\OR,\XOR)$ (polarization)};
\draw[<->,dashed] (syn)--node[lab,above]{$b\mapsto\ANT(b)$}(dif);
\draw[<->,dashed] (or)--node[lab,above]{$b\mapsto\ANT(b)$}(xor);
\draw[dotted,rounded corners=5pt] (-0.95,-0.95) rectangle (5.35,3.3);
\node[lab] at (-0.62,3.02) {$\Z$};
\node[box] (not) at (7.4,2.7) {$\NOT$};
\node[box,align=center] (all) at (7.4,0.6) {all of\\$\Fone,\Ftwo$};
\draw[->] (not)--node[lab,right,align=left]{zero-\\breaking}(all);
\draw[->,densely dotted] (5.3,1.2)--(all.west);
\end{tikzpicture}}\\[2pt]
{\small Figure S19. Structure of \chai{}. The five binary operations stay inside the zero-fixing clone $\Z$ (dotted; $f(0,\dots,0)=0$), and $\NOT$ leaves it. Dashed arrows: negation of the second argument.}
\end{center}
"""
env_=res['environment']
doc=r'''\documentclass[10pt,a4paper]{article}
\usepackage[T1]{fontenc}\usepackage{mathptmx}\usepackage[margin=2cm]{geometry}\usepackage{amsmath,amssymb,booktabs,longtable,array,pdflscape,hyperref}
\hypersetup{hidelinks}\setlength{\parindent}{0pt}\setlength{\parskip}{4pt}
\newcommand{\D}{D}\newcommand{\Z}{\mathcal Z}\newcommand{\SYN}{\mathsf{SYN}}\newcommand{\DIF}{\mathsf{DIF}}\newcommand{\OR}{\mathsf{OR}}\newcommand{\XOR}{\mathsf{XOR}}\newcommand{\AND}{\mathsf{AND}}\newcommand{\NOT}{\mathsf{NOT}}\newcommand{\ANT}{\mathsf{ANT}}\newcommand{\MAX}{\mathsf{MAX}}\newcommand{\MIN}{\mathsf{MIN}}\newcommand{\MAG}{\mathsf{MAG}}\newcommand{\chai}{Chai}\newcommand{\Fone}{\mathcal F_1}\newcommand{\Ftwo}{\mathcal F_2}\usepackage{tikz}\usetikzlibrary{arrows.meta,positioning}\newcommand{\yang}{\rule{0.8em}{0.18ex}}\newcommand{\yin}{\rule{0.25em}{0.18ex}\hspace{0.3em}\rule{0.25em}{0.18ex}}\newcommand{\tri}[3]{\raisebox{-0.01ex}{\shortstack[c]{#1\\[-0.3ex]#2\\[-0.3ex]#3}}}\newcommand{\iv}[1]{[\![#1]\!]}
\title{Chai: supplementary material\\[2pt]\large to ``Chai: A Signed Three-Valued Algebra from Parameterized Boolean Negation''}\author{Guangxian (Daway) Cai\\ Independent Researcher\\ ORCID: 0009-0003-2227-2162}\date{Supplement v1.0 for manuscript Chai2026v1}
\begin{document}\maketitle
This file is generated by \texttt{make\_supplement.py} from the results of \texttt{chai\_checker.py --full} and the other scripts of the package; it is not edited by hand, and the generator stops if any result file lacks data. Notation as in the paper; $-,0,+$ stand for $-1,0,+1$. The main paper contains the definitions, truth tables, and proofs of the arbitrary-arity generation theorems and of the completeness classification without primitive nullary constants. This supplement provides detailed finite closure data, additional term catalogues, and the computational evidence for the classification with nullary constants. Sections S0--S19 are cited from the paper as ``supplement Sn''; S14--S16 (quarter turn, complex table, folding, sign rule) and the constants case of S6 are checked by \texttt{verify\_direction\_rotation.py}; S17 reproduces the XY3 and G3 tables of the published divisional from their definitions, with hexaogram and bi-trigram symbols; S18 lists the checks added for this version (\texttt{chai\_verify\_v1\_increments.py}); S19 reproduces the structure diagram of the earlier version.

\section*{S0. What follows from classical theory and what is added here}
The table lists each result of the paper, the classical statement it rests on, and what the paper adds. The classical sources are Post (1921), S\l upecki (1939), Rosenberg (1970), Lau (2006) and, for strong Kleene logic, Robles and M\'endez (2022). The comparison is limited to the sources cited in the paper.
\begin{center}\small\begin{longtable}{>{\raggedright\arraybackslash}p{3.6cm}>{\raggedright\arraybackslash}p{5.6cm}>{\raggedright\arraybackslash}p{6.6cm}}
\toprule result&classical background&added here\\\midrule\endhead
Minimal carrier and the three sections (Prop.~1)&The least subset closed under an operation and containing a given set is a standard construction.&That $\{0,1\}$ closed under $\DIF(a,b)=\mathrm{sgn}(a-b)$ gives exactly $\{-1,0,+1\}$, relative to the chosen saturation; the three sections of $\DIF$ at the constants and the fact that only $\ANT$ is an involution. Elementary.\\
Boolean layer (Thm.~5)&The functions mapping $B^n$ into $B$ form a clone (an invariant-relation clone; Rosenberg, Lau).&$\SYN,\OR,\XOR,\AND,\NOT$ generate \emph{exactly} this clone, with explicit terms, and with arbitrary values at the points containing $-1$.\\
Zero-fixing barrier (Thm.~3)&$\Z=\mathrm{Pol}\{0\}$ is a maximal clone (Rosenberg; Lau). Hence any set inside $\Z$ is incomplete, and one operation outside $\Z$ can complete a set that generates $\Z$.&Explicit terms showing that the five binary operations generate \emph{all} of $\Z_n$ (not only a subclone), and that one zero-breaking operation then gives every operation (Thm.~3c reproves maximality with terms).\\
Odd operations (Thm.~4)&The functions commuting with $x\mapsto-x$ (self-dual, the centralizer of $\ANT$) form a clone; standard.&$\DIF$ alone generates this clone, with explicit terms $k_c$ and $s_p$ for every arity.\\
Classification, 40 of 127 (Thm.~6)&Post's and S\l upecki's completeness criteria; precomplete classes.&Exact condition for this seven-operation signature, minimal sets, counts with and without nullary constants, checked subset by subset by two methods (S6).\\
Polarization $\AND=\DIF(\OR,\XOR)$ (Thm.~1)&Elementary inclusion--exclusion on $\{0,1\}$-valued functions.&The identity at the sign level on $\D$, with the parity lemma, and its failure on $\{-m,\dots,m\}$, $m\ge2$.\\
Involution $\Phi=(\SYN,\DIF)$ (Thm.~2)&None specific.&Finite statement for this pair of operations; elementary.\\
Sign rule of $\AND$ (Remark~1, S16)&$ab=\tfrac14\big((a+b)^2-(a-b)^2\big)$, so $\mathrm{sgn}(ab)=\mathrm{sgn}(|a+b|-|a-b|)$ on all integers (quarter-square multiplication).&Only that on $\D$ the squares can be folded to signs before the comparison (Thm.~1); not a general principle (fails on $\{-2,\dots,2\}$).\\
Constants as nullary operations, 95 of 127 (Remark~3, S6)&The Post--S\l upecki criteria with constants; monotone operations and the relation $\varrho$ as invariants.&A short constructive proof for this signature (terms for $\NOT$ and $\ANT$ from $\SYN$ with one further operation) and the exact criterion; checked against the 127 subsets by a separate script.\\
Quarter turn of the truth table (Remark~2, S14)&Rotating a matrix is precomposition of the function with a rotation of the argument grid; multiplication by $i$ is the quarter turn of the complex plane.&That for \chai{} the turn carries $\SYN\to\DIF\to-\SYN$, $\OR\leftrightarrow\XOR$, $\AND\to-\AND$; the link with the order-8 symmetry group and the parity of functions. Elementary.\\
Complex table $\SYN+i\,\DIF$ (supplement only, S15)&$\mathrm{sgn}$ applied to both coordinates commutes with multiplication by $i$ (the dihedral symmetry of $\mathbb Z^2$).&The identification of the two tables with the saturated real and imaginary parts of $(1+i)\overline{(a+ib)}$, and with the sum and difference views of a BID. Elementary.\\
Folding is not a quotient (Prop.~5); simplicity and the rest of S16 (supplement only)&A functionally complete (primal) finite algebra with constants has no nontrivial congruence, subalgebra or automorphism; standard.&The exact list of the operations that respect $\ker\MAG$ (only $\AND,\ANT$), the folded complex table, and the Boolean/set readings. Elementary; the set statements are classical.\\
$K_3$, $L_3$ as terms (Prop.~3, S3)&Completeness implies that every function is a term.&Explicit minimal terms of familiar connectives, found by exhaustive search (S3, S5).\\
Difference lift, orbit count (Section VI, Prop.~6 and 7)&Elementary counting and Burnside's lemma.&The counts for this setting. The XY3/G3/C9 layers themselves are described in the published divisional CN~122173050~A and are background.\\
\bottomrule\end{longtable}\end{center}

\section*{S1. Algebraic laws of the seven binary operations}
Associativity failures are counted over the 27 ordered triples.
\begin{center}\small\begin{tabular}{lccccc}\toprule
operation & commutative & assoc.\ failures & idempotent & identity & absorbing\\\midrule
'''+LAWS+r'''
\bottomrule\end{tabular}\end{center}
Left-distributive pairs $(f,g)$, meaning $f(a,g(b,c))=g(f(a,b),f(a,c))$ for all $a,b,c$, among these seven operations: '''+LDT+r'''.

\section*{S2. The pseudocomplement $\NOT$}
On $\{0,+1\}$, $\NOT(x)=1-x$ is the Boolean complement. On $\D$, $\NOT(a)=\min\{x\in\{0,+1\}:\MAX(a,x)=+1\}$: it is the least Boolean element that joins $a$ to $+1$, so $\NOT(-1)=\NOT(0)=+1$ and $\NOT(+1)=0$. It is a \emph{relative} dual pseudocomplement (relative to the sublattice $\{0,+1\}$), not a pseudocomplement of $\D$, and $\NOT\circ\NOT=[0,0,+1]\neq\mathrm{id}$.

\section*{S3. The explicit terms with their complete tables}
\begin{align*}
\delta_{+1}(x)&=\SYN(x,\OR(x,x)), & \delta_{-1}(x)&=\XOR(x,\OR(x,x)),\\
h_{\pm1}(u,v)&=\AND(\delta_{+1}(u),\delta_{\pm1}(v)), & h_0(u,v)&=\DIF(\delta_{+1}(u),\AND(\delta_{+1}(u),\OR(v,v))),\\
k_{-1}(u,v)&=\DIF(u,\DIF(v,\DIF(v,u))), & k_0(u,v)&=\DIF(v,\DIF(v,\DIF(u,v))),\\
k_{+1}(u,v)&=\DIF(\DIF(u,v),\DIF(\DIF(u,v),v)), & \mathsf{MAG}(x)&=\NOT(\DIF(x,\DIF(x,\NOT(x)))).
\end{align*}
Unary terms ($\delta_{+1}=\iv{x=+1}$, $\delta_{-1}=\iv{x=-1}$, $\mathsf{MAG}(x)=x^2$):
\begin{center}'''.replace(r'\iv{x=+1}','[x{=}+1]').replace(r'\iv{x=-1}','[x{=}-1]')+unary_small+r'''\end{center}
Lemma~2: $h_{+1},h_{-1},h_0$ equal $[u{=}{+1}][v{=}{+1}]$, $[u{=}{+1}][v{=}{-1}]$, $[u{=}{+1}][v{=}0]$ (rows $u$, columns $v$):
\begin{center}'''+tabs_h+r'''\end{center}
Theorem~4: $k_c(u,v)=u$ if $u\neq0$ and $v=cu$, and $0$ otherwise (so $k_{-1}$ is $u$ exactly where $v=-u$, $k_0$ where $v=0$, $k_{+1}$ where $v=u$):
\begin{center}'''+tabs_k+r'''\end{center}
All of these tables are asserted by the generator against the intended functions at every point of $\D$ or $\D^2$.

\textbf{Familiar binary connectives as terms} (moved here from the paper; $K_3$ and $L_3$ coding $-1=F$, $0=U$, $+1=T$). Size: number of operator nodes of the term (tree size, repeated subterms counted each time), equal to the minimum over the seven operators found by exhaustive search up to size 5; each term below is evaluated by the generator at all nine points and compared with the stated function.
\begin{center}\small\begin{longtable}{llcp{8.2cm}}\toprule connective&meaning&size&term\\\midrule\endhead
'''+CONNT+r'''
\bottomrule\end{longtable}\end{center}
The standard connectives are terms, not primitives, and the primitives differ from them: $\NOT(0)=+1\neq\ANT(0)=0$, $\AND(-1,-1)=+1\neq\min(-1,-1)=-1$, and $\OR$ is not the lattice join. The $K_3$ and $L_3$ implications differ only at $(0,0)$ ($U\to U$ is $U$ in $K_3$ and $T$ in $L_3$).

\section*{S4. Generation heights}
Closure of the binary slice $\D^{\D^2}$: number of new functions at heights $0,1,2,\dots$ (the two projections, and in the last six rows also the three constants, are at height $0$). Produced by \texttt{chai\_checker.py --full} and \texttt{--full-consts}; the logs are \texttt{full\_run.log} and \texttt{consts\_*.log}.
\begin{center}\small\begin{tabular}{lrl}\toprule
generators & functions & new per height\\\midrule
'''+HEIGHTS+r'''
\bottomrule\end{tabular}\end{center}
The closure of $\DIF$ alone without constants is the set of $81=3^4$ odd binary functions. The first row has $6561=3^8$ functions, the binary slice of $\Z$; the two rows for $M_1$ and $M_2$ stay below $19\,683$, as the classification requires.

\begin{landscape}
\section*{S5. The 27 unary functions as terms}
Code: digits $0,1,2$ for $0,+1,-1$, listing $g(-1),g(0),g(+1)$, prefix \texttt{0t}. Size: minimum number of operator nodes (tree size, repeated subterms counted each time; exhaustive search by size); term: the lexicographically smallest minimal term. Three signatures: all seven operations (the table of v9, Table~IV; maximum size '''+str(mx(0))+r'''), $\{\DIF,\NOT\}$ (maximum '''+str(mx(1))+r''') and $\{\SYN,\ANT,\NOT\}$ (maximum '''+str(mx_san)+r'''). The three odd unary functions $0$, $x$, $\ANT(x)$ are generated by $\DIF$ alone: $\DIF(x,x)$, $x$, $\DIF(\DIF(x,x),x)$.
\begin{center}\scriptsize\begin{longtable}{lcc>{\raggedright\arraybackslash}p{4.6cm}c>{\raggedright\arraybackslash}p{4.6cm}c>{\raggedright\arraybackslash}p{4.4cm}}
\toprule code&$g$&size&all seven&size&$\{\DIF,\NOT\}$&size&$\{\SYN,\ANT,\NOT\}$\\\midrule\endhead
'''+UNARY+r'''
\bottomrule\end{longtable}\end{center}
\end{landscape}

\section*{S6. Classification of all 127 nonempty subsets}
For each subset: complete without nullary constants (Theorem~6) and its certificate of incompleteness (an invariant shared by all its operations), and complete when $0,+1,-1$ are nullary operations (Remark~3 of the paper) with the reason if not. The criterion of Proposition~6 (the subset contains $\DIF$, or contains $\SYN$ and one more operation) reproduces the last two columns on all $127$ subsets (checked by \texttt{verify\_direction\_rotation.py}); the proof certificates are: a monotone-operations argument for $\{\SYN\}$ and preservation of $\varrho$ (S7) for the subsets of $\{\OR,\XOR,\AND,\NOT,\ANT\}$. Both columns were computed by two methods that agree on every subset (\texttt{classification\allowbreak\_127.csv}, \texttt{classification\allowbreak\_table.py}), and by a third script (\texttt{independent\allowbreak\_free\allowbreak\_constants.py}).
\begin{center}\scriptsize\begin{longtable}{rlcp{3.9cm}cp{3.9cm}}
\toprule \#&subset&no const.&certificate&with const.&reason\\\midrule\endhead
'''+CLS+r'''
\bottomrule\end{longtable}\end{center}

\section*{S7. The invariants used for incompleteness}
\begin{center}\small\begin{tabular}{lccc}\toprule operation&fixes $0$&preserves $\{0,+1\}$&preserves $\varrho$\\\midrule
'''+INV+r'''
\bottomrule\end{tabular}\end{center}
Here $\varrho=\{(a,b)\in\D^2:ab\ge0\}$ has $7$ elements and an operation preserves it if it maps every $k$-tuple of pairs of $\varrho$ to a pair of $\varrho$ (checked over the $7^2$ or $7$ cases); the constants $c\in\D$ preserve it because $(c,c)\in\varrho$. Witnesses: $\DIF(0,+1)=-1\notin\{0,+1\}$; $(+1,0),(0,-1)\in\varrho$ but $(\SYN(+1,0),\SYN(0,-1))=(+1,-1)\notin\varrho$.

\section*{S8. The four BID states, involution, difference lift and orbits}
The four states of a BID $(X,Y)$ with the difference view $d=X-Y$, the sum view $c=X+Y$ and the total selector of Section~VI of the paper ($S_c$ is the one-hot branch selector; $(1,1)$ goes to $S_0$ by convention):
\begin{center}\small\begin{tabular}{cccc}\toprule $(X,Y)$&$d$&$c$&selector\\\midrule
'''+bidrows+r'''
\bottomrule\end{tabular}\end{center}
$\Phi(a,b)=(\SYN(a,b),\DIF(a,b))$ on $\D^2$ (asserted: $\Phi\circ\Phi=\mathrm{id}$):
\begin{center}\small'''+PHI+r'''\end{center}
Difference lift for $n=1$: the $27$ unary functions have '''+LIFT+r'''; the total is $64=8^2$ pairs of $\{0,1\}$-valued functions (the $8$ XY3 masks). For $n=2$: $512$ masks, $262\,144$ pairs (all enumerated by \texttt{check\_c9\_chain.py}).
Orbits of binary truth tables under the order-$8$ group generated by $\mathsf{SWAP}$, $\tau_A$, $\tau_B$: $'''+str(orb['burnside'])+r'''$ in total, namely '''+ORB+r'''; cycle counts of the $8$ group elements: '''+', '.join(map(str,orb['cycle_counts']))+r'''.

\section*{S9. A worked example of the three-slice C9 term}
For the ternary sum modulo $3$, $F(a,b)=((a+b+1)\bmod 3)-1$, the three branches $G_c=F(c,\cdot)$ are unary functions with the minimal terms below (variable $b$), and $P_c=\AND(E_c(a),G_c(b))$ with $E_c(x)=\NOT(\XOR(x,c))$, where the constants are the terms $0=\DIF(x,x)$, $+1=\NOT(0)$, $-1=\ANT(+1)$:
\begin{center}\small\begin{longtable}{clp{4.4cm}p{7cm}}\toprule $c$&$G_c$&term of $G_c$&$P_c$ as a term\\\midrule\endhead
'''+C9+r'''
\bottomrule\end{longtable}\end{center}
Then $F(a,b)=\SYN\big(P_{-1},\SYN(P_0,P_{+1})\big)$. The generator evaluates this expression at all nine points and asserts that it equals $F$. \texttt{check\_c9\_chain.py} does the same for all $19\,683$ binary functions.

\section*{S10. Nine-trit configuration words (T9ID)}
Row-major reading of the truth table (rows $a=-,0,+$, columns $b=-,0,+$), digits $0,1,2$ for $0,+1,-1$. The word of $\SYN$ is the word that the published divisional CN~122173050~A lists for its operation named OR; the words of $\DIF$ and $\AND$ agree with it. Words of the first six operations marked below match \texttt{results.json}.
\begin{center}\small\begin{tabular}{llc}\toprule operation&T9ID word&rows $a=-,0,+$ (columns $b=-,0,+$)\\\midrule
'''+T9R+r'''
\bottomrule\end{tabular}\end{center}

\section*{S11. Where each claim of the paper is checked}
\begin{center}\scriptsize\begin{longtable}{>{\raggedright\arraybackslash}p{4.6cm}>{\raggedright\arraybackslash}p{5.2cm}>{\raggedright\arraybackslash}p{6.2cm}}\toprule claim&script&output key\\\midrule\endhead
'''+TRR+r'''
\bottomrule\end{longtable}\end{center}

\section*{S12. What the computations establish, and what they do not}
\textbf{Exhaustive:} the parity lemma and polarization on all $9$ pairs; the involution on all $9$ points; every function of $\Z_2$ ($6561$) rebuilt from the five operations; all $19\,683$ binary and $27$ unary functions rebuilt from $\{\DIF,\NOT\}$; the closures of S4; the classification of all $127$ subsets by two methods; the difference lift for $n=1$ and $n=2$ (all $64$ and $262\,144$ pairs); all $19\,683$ binary functions as three-slice C9 terms.
\textbf{Random samples (fixed seed 7):} $25$ functions each of arity $3$ and $4$ rebuilt from $\{\DIF,\NOT\}$; $25$ functions each of $\Z_3$ and $\Z_4$ from the five operations; $25$ odd functions each of arity $1$, $3$, $4$ from $\DIF$. These are evidence, not proofs.
\textbf{Direction, rotation and folding (S14--S16):} the program \texttt{verify\_direction\_rotation.py} (standard library only, '''+str(DR_TESTS)+r''' tests, Python '''+DR_PY+r''') checks the table cycle, the lemma $\mathcal R(T_F)=T_{F\circ\rho}$ for all $19\,683$ functions, $\mathcal RZ=-iZ$ in exact Gaussian integers, the folding statements, the congruences, subalgebras and automorphisms by exhaustion, the Boolean and set identities, the $C_4$ and dihedral orbit counts, and the $17$ sign-linear functions; it includes negative controls (wrong sign, wrong direction, reflection for rotation) that must come out false.
\textbf{Proofs:} the statements for every arity $n$ (Theorems 3, 4 and 5) rest on the explicit constructions of Lemma~2 and Theorems~3, 4 and~5 in the paper, not on the sampled cases; the code checks the building blocks (S3) exhaustively, and S18 evaluates the constructions on every binary function of the classes.
Environment: Python '''+env_['python']+r''', numpy '''+env_['numpy']+r'''; the full run took about %.0f seconds on one core. Hashes of all files of the package are in \texttt{SHA256SUMS}.

\section*{S13. Use of generative AI tools}
Generative AI tools, including ChatGPT by OpenAI, Claude by Anthropic, and Gemini by Google, were used to assist with language editing, to discuss candidate proof formulations, and to develop and review the verification scripts described in this supplement. The author independently checked the mathematical arguments, the cited sources and the computational results, and accepts responsibility for the final content.
'''%runtime+SEC1416+SEC17+SEC18+SEC19+'\n\\end{document}\n'
doc=doc.replace('--full','-{}-full')
open(P('chai2026v1_supplement.tex'),'w',encoding='utf-8').write(doc)
print('wrote chai2026v1_supplement.tex')
