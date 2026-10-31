"""Audit of the chain polarization -> difference lift -> G3 -> C9 (three slices).
Run: python3 check_c9_chain.py   (pure Python, a few seconds). Every item is an assertion."""
import itertools, sys, time
from itertools import product
D=(-1,0,1)
sgn=lambda z:(z>0)-(z<0)
SYN=lambda a,b:sgn(a+b); DIF=lambda a,b:sgn(a-b); OR=lambda a,b:sgn((a+b)**2)
XOR=lambda a,b:sgn((a-b)**2); AND=lambda a,b:sgn(4*a*b)
NOT=lambda a:DIF(1,a); ANT=lambda a:DIF(0,a)
T0=time.time()
def say(k,v): print('[%5.1fs] %s: %s'%(time.time()-T0,k,v)); sys.stdout.flush()

# 0. the printed truth tables (rows a=-1,0,+1; columns b=-1,0,+1; '-','0','+')
sym={-1:'-',0:'0',1:'+'}
def table(f): return [''.join(sym[f(a,b)] for b in D) for a in D]
PRINTED={'SYN':['--0','-0+','0++'],'DIF':['0--','+0-','++0'],'OR':['++0','+0+','0++'],
         'XOR':['0++','+0+','++0'],'AND':['+0-','000','-0+']}
FUN={'SYN':SYN,'DIF':DIF,'OR':OR,'XOR':XOR,'AND':AND}
for k in PRINTED: assert table(FUN[k])==PRINTED[k],k
assert [sym[NOT(a)] for a in D]==list('++0') and [sym[ANT(a)] for a in D]==list('+0-')
say('truth tables of the paper','match')

# 1. parity lemma and polarization
for a,b in product(D,D):
    assert (abs(a+b),abs(a-b)) in {(0,0),(0,2),(2,0),(1,1)}
    p,q=a+b,a-b
    assert p*p-q*q==4*a*b
    assert AND(a,b)==sgn(p*p-q*q)==DIF(OR(a,b),XOR(a,b))
bad=[(a,b) for a,b in product(range(-2,3),repeat=2) if sgn((a+b)**2-(a-b)**2)!=sgn(sgn((a+b)**2)-sgn((a-b)**2))]
assert (2,1) in bad
say('polarization on D; fails on {-2..2} at',bad[:4])

# 2. difference lift, n=1 and n=2, by enumerating all pairs (U,V)
def lift(n):
    pts=3**n
    cnt={}
    for U in product((0,1),repeat=pts):
        for V in product((0,1),repeat=pts):
            F=tuple(u-v for u,v in zip(U,V))
            cnt[F]=cnt.get(F,0)+1
    return cnt
c1=lift(1)
assert len(c1)==27 and all(c1[F]==2**F.count(0) for F in c1) and sum(c1.values())==64==8**2
c2=lift(2)
assert len(c2)==3**9 and all(c2[F]==2**F.count(0) for F in c2) and sum(c2.values())==4**9==512**2
say('difference lift','n=1: 27 functions, 64 = 8^2 pairs; n=2: 19683 functions, %d = 512^2 pairs, counts 2^z exact'%sum(c2.values()))

# 3. G3: every unary function is a term of the seven operations
ops2=[SYN,DIF,OR,XOR,AND]; ops1=[NOT,ANT]
ident=tuple(D)
reach={ident}; frontier={ident}
while frontier:
    new=set()
    cur=list(reach)
    for f in frontier:
        for o in ops1: new.add(tuple(o(x) for x in f))
        for g in cur:
            for o in ops2:
                new.add(tuple(o(x,y) for x,y in zip(f,g))); new.add(tuple(o(y,x) for x,y in zip(f,g)))
    frontier=new-reach; reach|=new
assert len(reach)==27
say('G3: unary functions reachable from the identity','27 of 27')
# a unary function g equals U-V for 2^z mask pairs: already implied by step 2 (n=1), checked above

# 4. C9: three slices, one-hot selectors, term expansion, difference form
def slice_triple(F): return tuple(tuple(F[3*i:3*i+3]) for i in range(3))  # F row-major: a=-1,0,+1
assert len({slice_triple(F) for F in c2})==3**9
const={-1:ANT(1),0:0,1:1}
E=lambda c,x:NOT(XOR(x,const[c]))
for x in D:
    assert [E(c,x) for c in D].count(1)==1 and [E(c,x) for c in D].count(-1)==0 and E(0,x)+E(1,x)+E(-1,x)==1
    for c in D: assert E(c,x)==(1 if x==c else 0)
n=0
for F in c2:
    g=slice_triple(F)           # g[i][j]=F(a=D[i], b=D[j])
    for a,b in product(D,D):
        i,j=a+1,b+1
        P=[AND(E(c,a),g[c+1][j]) for c in D]
        assert SYN(P[0],SYN(P[1],P[2]))==F[3*i+j]
    # difference form with canonical masks u_c=[g_c=+1], v_c=[g_c=-1]
    for a,b in product(D,D):
        i,j=a+1,b+1
        U=sum((1 if a==c else 0)*(1 if g[c+1][j]==1 else 0) for c in D)
        V=sum((1 if a==c else 0)*(1 if g[c+1][j]==-1 else 0) for c in D)
        assert U in (0,1) and V in (0,1) and U-V==F[3*i+j]
    n+=1
say('C9 three-slice term and difference form','all %d binary functions'%n)

# 5. BID states: total selector equals [d=c]; one-hot; value of C9 on a BID pair
for X,Y in product((0,1),repeat=2):
    d=X-Y
    St=(1-X)*Y; Sr=(1-X)*(1-Y)+X*Y; Si=X*(1-Y)       # S_{-1},S_0,S_{+1} as 0/1
    assert [St,Sr,Si]==[1 if d==c else 0 for c in D] and St+Sr+Si==1
# canonical-domain signals of the form S_T=Y, S_R=(1-X)(1-Y), S_I=X add up to 2 at (1,1)
assert [(Y,(1-X)*(1-Y),X) for X,Y in product((0,1),repeat=2)]==[(0,1,0),(1,0,0),(0,0,1),(1,0,1)]
say('BID states','total selector = [d=c], one-hot on all four; the plain signals Y,(1-X)(1-Y),X are not one-hot at (1,1)')

# 6. configuration words (row-major, digits 0,1,2 for 0,+1,-1) and the BID reading of d and c
dig={0:'0',1:'1',-1:'2'}
def t9(f): return '0t'+'_'.join(''.join(dig[f(a,b)] for b in D) for a in D)
words={'SYN':t9(SYN),'DIF':t9(DIF),'OR':t9(OR),'XOR':t9(XOR),'AND':t9(AND),
       'MAX':t9(lambda a,b:max(a,b)),
       'SUM3':t9(lambda a,b:{0:0,1:1,2:-1}[(a+b)%3]),
       'CARRY':t9(lambda a,b:1 if (a,b)==(1,1) else (-1 if (a,b)==(-1,-1) else 0))}
say('T9ID words',words)
for X,Y in product((0,1),repeat=2):
    c,d=X+Y,X-Y
    assert X==(c+d)//2 and Y==(c-d)//2 and c*c-d*d==4*X*Y and d==c-2*Y and c==d+2*Y
assert len({(X+Y,X-Y) for X,Y in product((0,1),repeat=2)})==4
assert len({X-Y for X,Y in product((0,1),repeat=2)})==3 and len({X+Y for X,Y in product((0,1),repeat=2)})==3
say('c,d coordinates','(c,d) reversible, neither alone is; c^2-d^2=4XY')
print('ALL ASSERTIONS PASSED')
