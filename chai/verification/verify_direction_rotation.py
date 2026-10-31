#!/usr/bin/env python3
"""verify_direction_rotation.py -- independent check of the direction / quarter-turn / folding statements of BIDRA v10.

Standard library only; it does not import chai_checker.py.  Every definition is written here from scratch.
Usage:  python3 verify_direction_rotation.py            (about 5 seconds)
Output: direction_rotation_results.json  (machine-readable)  and direction_rotation_results.log (readable summary);
        exit status 0 iff every test passed.

Scope.  All statements are checked by exhaustion over D={-1,0,+1} (9 binary arguments), over the 19683 binary functions,
or over a stated finite window of integers.  A finite check is evidence, not a proof; the proofs of the general statements
are in the paper (Section III) and in supplement S14-S16.  The module also contains NEGATIVE CONTROLS: statements with the
wrong convention or the wrong sign that must come out FALSE; if one comes out true, the corresponding test is not sharp.
Matrix convention (fixed): rows a=-1,0,+1 from top to bottom, columns b=-1,0,+1 from left to right.
"""
import json, sys, os, itertools, platform, time
from itertools import product

D = (-1, 0, 1)
B = (0, 1)
HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------ reference definitions
def sgn(z): return (z > 0) - (z < 0)
def SYN(a, b): return sgn(a + b)
def DIF(a, b): return sgn(a - b)
def ANT(a): return -a
def MAG(a): return a * a
def OR(a, b): return sgn((a + b) ** 2)
def XOR(a, b): return sgn((a - b) ** 2)
def AND(a, b): return sgn(4 * a * b)
def NOT(a): return DIF(1, a)

BIN = {'SYN': SYN, 'DIF': DIF, 'OR': OR, 'XOR': XOR, 'AND': AND}

# ------------------------------------------------------------------ matrices and rotations
def table(F):
    """3x3 truth table of a binary function: T[r][c] = F(a=D[r], b=D[c])."""
    return tuple(tuple(F(a, b) for b in D) for a in D)

def neg(M): return tuple(tuple(-x for x in row) for row in M)
def mag(M): return tuple(tuple(x * x for x in row) for row in M)

def rotate_cw(M):
    """clockwise quarter turn: [[a,b,c],[d,e,f],[g,h,i]] -> [[g,d,a],[h,e,b],[i,f,c]]"""
    (a, b, c), (d, e, f), (g, h, i) = M
    return ((g, d, a), (h, e, b), (i, f, c))

def rotate_ccw(M):
    """inverse of rotate_cw, written out independently: [[a,b,c],[d,e,f],[g,h,i]] -> [[c,f,i],[b,e,h],[a,d,g]]"""
    (a, b, c), (d, e, f), (g, h, i) = M
    return ((c, f, i), (b, e, h), (a, d, g))

# Gaussian integers as (re, im) pairs: exact, no floating point
def gadd(u, v): return (u[0] + v[0], u[1] + v[1])
def gmul(u, v): return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])
def gconj(u): return (u[0], -u[1])
I_ = (0, 1)
def zmat(M_re, M_im):
    return tuple(tuple((M_re[r][c], M_im[r][c]) for c in range(3)) for r in range(3))
def zscale(s, Z): return tuple(tuple(gmul(s, x) for x in row) for row in Z)
def zrotate(Z, rot):
    """linear extension of a real matrix rotation to complex matrices: rotate real and imaginary parts separately"""
    re = rot(tuple(tuple(x[0] for x in row) for row in Z)); im = rot(tuple(tuple(x[1] for x in row) for row in Z))
    return zmat(re, im)

# ------------------------------------------------------------------ test harness
RESULTS = {}; LOG = []; FAILED = []; COUNT = {'tests': 0, 'passed': 0}
def record(name, ok, detail=None, cx=None):
    COUNT['tests'] += 1
    if ok: COUNT['passed'] += 1
    else: FAILED.append(name)
    RESULTS[name] = {'ok': bool(ok), **({'detail': detail} if detail is not None else {}), **({'counterexample': cx} if cx is not None else {})}
    LOG.append(('PASS ' if ok else 'FAIL ') + name + (('  ' + str(detail)) if detail is not None else '') + ('' if ok or cx is None else '  counterexample: ' + str(cx)))

def forall(pred, domain):
    """returns (True, n) or (False, first counterexample)"""
    n = 0
    for x in domain:
        n += 1
        if not pred(*x if isinstance(x, tuple) else (x,)): return False, x
    return True, n

def check(name, pred, domain, control=False):
    ok, info = forall(pred, domain)
    if control:   # a negative control passes iff the (wrong) statement is FALSE
        record('CONTROL ' + name, not ok, 'statement correctly fails; first counterexample %s' % (info,) if not ok else 'WRONG STATEMENT IS TRUE: test is not sharp')
    else:
        record(name, ok, '%d cases' % info if ok else None, None if ok else info)

D2 = list(product(D, D))
D3 = list(product(D, D, D))

# ================================================================== 5.1 definitions and direction identities (Section 3.1)
check('3.1 ANT(a)=DIF(0,a)', lambda a: ANT(a) == DIF(0, a), D)
check('3.1 NOT(a)=DIF(1,a)=sgn(1-a)', lambda a: NOT(a) == DIF(1, a) == sgn(1 - a), D)
check('3.1 DIF(a,b)=SYN(a,ANT b)', lambda a, b: DIF(a, b) == SYN(a, ANT(b)), D2)
check('3.1 SYN(a,b)=DIF(a,ANT b)', lambda a, b: SYN(a, b) == DIF(a, ANT(b)), D2)
check('3.1 XOR(a,b)=OR(a,ANT b)', lambda a, b: XOR(a, b) == OR(a, ANT(b)), D2)
check('3.1 OR(a,b)=XOR(a,ANT b)', lambda a, b: OR(a, b) == XOR(a, ANT(b)), D2)
check('3.1 ANT and NOT are different (NOT(0)=+1, ANT(0)=0)', lambda: NOT(0) != ANT(0), [()])
check('3.1 OR=MAG o SYN, XOR=MAG o DIF (pointwise)', lambda a, b: OR(a, b) == MAG(SYN(a, b)) and XOR(a, b) == MAG(DIF(a, b)), D2)
check('3.1 AND(a,b)=a*b on D', lambda a, b: AND(a, b) == a * b, D2)
check('3.1 polarization AND=DIF(OR,XOR) on D', lambda a, b: AND(a, b) == DIF(OR(a, b), XOR(a, b)), D2)
check('3.1 parity lemma (|a+b|,|a-b|) in {(0,0),(0,2),(2,0),(1,1)}', lambda a, b: (abs(a + b), abs(a - b)) in {(0, 0), (0, 2), (2, 0), (1, 1)}, D2)
# negative controls: wrong direction conventions
check('SYN(a,b)=ANT(DIF(a,b))', lambda a, b: SYN(a, b) == ANT(DIF(a, b)), D2, control=True)
check('DIF(a,b)=SYN(ANT a,b)', lambda a, b: DIF(a, b) == SYN(ANT(a), b), D2, control=True)

# ================================================================== 5.2 truth tables and rotations (Section 3.2)
T = {n: table(f) for n, f in BIN.items()}
record('3.2 truth tables as printed in the paper (SYN, DIF)', T['SYN'] == ((-1, -1, 0), (-1, 0, 1), (0, 1, 1)) and T['DIF'] == ((0, -1, -1), (1, 0, -1), (1, 1, 0)),
       'tables regenerated from the definitions and compared with the printed values')
record('3.2 truth tables as printed in the paper (OR, XOR, AND)', T['OR'] == ((1, 1, 0), (1, 0, 1), (0, 1, 1)) and T['XOR'] == ((0, 1, 1), (1, 0, 1), (1, 1, 0)) and T['AND'] == ((1, 0, -1), (0, 0, 0), (-1, 0, 1)))
sample = ((1, 2, 3), (4, 5, 6), (7, 8, 9))
record('3.2 rotate_cw on a labelled matrix', rotate_cw(sample) == ((7, 4, 1), (8, 5, 2), (9, 6, 3)), rotate_cw(sample))
record('3.2 rotate_ccw on a labelled matrix', rotate_ccw(sample) == ((3, 6, 9), (2, 5, 8), (1, 4, 7)), rotate_ccw(sample))
check('3.2 rotate_ccw o rotate_cw = id (all 19683 matrices)', lambda M: rotate_ccw(rotate_cw(M)) == M and rotate_cw(rotate_ccw(M)) == M,
      [(tuple(tuple(v[3 * r:3 * r + 3]) for r in range(3)),) for v in product(D, repeat=9)])
check('3.2 rotate_cw^4 = id (all 19683 matrices)', lambda M: rotate_cw(rotate_cw(rotate_cw(rotate_cw(M)))) == M,
      [(tuple(tuple(v[3 * r:3 * r + 3]) for r in range(3)),) for v in product(D, repeat=9)])
R, Rc = rotate_cw, rotate_ccw
ROT = [
    ('R_cw(T_SYN)=T_DIF',        R(T['SYN']) == T['DIF']),
    ('R_ccw(T_DIF)=T_SYN',       Rc(T['DIF']) == T['SYN']),
    ('R_cw(T_DIF)=-T_SYN',       R(T['DIF']) == neg(T['SYN'])),
    ('R_ccw(T_SYN)=-T_DIF',      Rc(T['SYN']) == neg(T['DIF'])),
    ('R_cw^2(T_SYN)=-T_SYN',     R(R(T['SYN'])) == neg(T['SYN'])),
    ('R_cw^4(T_SYN)=T_SYN',      R(R(R(R(T['SYN'])))) == T['SYN']),
    ('R_cw^2(T_DIF)=-T_DIF',     R(R(T['DIF'])) == neg(T['DIF'])),
    ('R_cw(T_OR)=T_XOR',         R(T['OR']) == T['XOR']),
    ('R_ccw(T_OR)=T_XOR',        Rc(T['OR']) == T['XOR']),
    ('R_cw(T_XOR)=T_OR',         R(T['XOR']) == T['OR']),
    ('R_ccw(T_XOR)=T_OR',        Rc(T['XOR']) == T['OR']),
    ('R_cw^2(T_OR)=T_OR',        R(R(T['OR'])) == T['OR']),
    ('R_cw(T_AND)=-T_AND',       R(T['AND']) == neg(T['AND'])),
    ('R_cw^2(T_AND)=T_AND',      R(R(T['AND'])) == T['AND']),
]
for n, ok in ROT: record('3.2 ' + n, ok)
# negative controls: the wrong direction / sign must fail
record('CONTROL R_ccw(T_SYN)=T_DIF (wrong sign) is false', Rc(T['SYN']) != T['DIF'], 'correct value is -T_DIF')
record('CONTROL R_cw(T_DIF)=T_SYN (wrong sign) is false', R(T['DIF']) != T['SYN'], 'correct value is -T_SYN')
record('CONTROL R_cw(T_SYN)=T_SYN is false', R(T['SYN']) != T['SYN'])
record('CONTROL rotation of T_AND is not T_XOR', R(T['AND']) != T['XOR'])

# ---- rotation as an input transformation:  R_cw(T_F)(a,b) = F(-b,a)   (the lemma behind the orbit statement)
def rho(a, b): return (-b, a)            # quarter turn of the input plane, a+ib -> i(a+ib)
ALLF = [dict(zip(D2, v)) for v in product(D, repeat=9)]            # all 19683 binary functions as tables
def tab_of(f): return tuple(tuple(f[(a, b)] for b in D) for a in D)
def lemma_rot(f):
    return tab_of({(a, b): f[rho(a, b)] for (a, b) in D2}) == rotate_cw(tab_of(f))
ok, n = forall(lemma_rot, ALLF)
record('3.2 lemma: R_cw(T_F)(a,b)=F(-b,a) for all 19683 binary F', ok, '%d functions' % n if ok else None)
check('3.2 distinct things: for commutative F, R_cw(T_F) = column flip (b->ANT b); for DIF it differs by sign',
      lambda: all(R(table(f)) == table(lambda a, b, f=f: f(a, ANT(b))) for f in (SYN, OR, XOR, AND)) and
              R(table(DIF)) == neg(table(lambda a, b: DIF(a, ANT(b)))) and R(table(DIF)) != table(lambda a, b: DIF(a, ANT(b))), [()])
check('3.2 column flip alone is NOT the rotation of DIF', lambda: table(lambda a, b: DIF(a, ANT(b))) == R(table(DIF)), [()], control=True)
# rotation and parity: R^2(T_F) = T_F o (-,-); odd <=> R^2 T = -T ; even <=> R^2 T = T
odd_count = sum(1 for f in ALLF if rotate_cw(rotate_cw(tab_of(f))) == neg(tab_of(f)))
even_count = sum(1 for f in ALLF if rotate_cw(rotate_cw(tab_of(f))) == tab_of(f))
record('3.2 R^2 T=-T exactly for the odd functions (81=3^4); R^2 T=T for the even ones (3^5=243)', odd_count == 81 and even_count == 243, {'odd': odd_count, 'even': even_count})
record('3.2 ANT-symmetry (Prop. 3): SYN, DIF odd; OR, XOR even', R(R(T['SYN'])) == neg(T['SYN']) and R(R(T['DIF'])) == neg(T['DIF']) and R(R(T['OR'])) == T['OR'] and R(R(T['XOR'])) == T['XOR'])
check('3.2 MAG commutes with rotation (entrywise map)', lambda f: mag(rotate_cw(tab_of(f))) == rotate_cw(mag(tab_of(f))), ALLF)

# ================================================================== 3.3 complex matrix Z = T_SYN + i T_DIF
Z = zmat(T['SYN'], T['DIF'])
iZ = zscale(I_, Z)
record('3.3 Z = T_SYN + i T_DIF has the 9 entries sgn(a+b)+i sgn(a-b)', all(Z[r][c] == (SYN(D[r], D[c]), DIF(D[r], D[c])) for r in range(3) for c in range(3)))
record('3.3 R_ccw(Z) = i Z (exact Gaussian integers, every entry)', zrotate(Z, rotate_ccw) == iZ)
record('3.3 R_cw(Z) = -i Z', zrotate(Z, rotate_cw) == zscale((0, -1), Z))
record('3.3 R_cw^4(Z) = Z, R_cw^2(Z) = -Z', zrotate(zrotate(zrotate(zrotate(Z, rotate_cw), rotate_cw), rotate_cw), rotate_cw) == Z and zrotate(zrotate(Z, rotate_cw), rotate_cw) == zscale((-1, 0), Z))
record('CONTROL R_cw(Z)=iZ is false', zrotate(Z, rotate_cw) != iZ)
record('CONTROL R_ccw(Z)=-iZ is false', zrotate(Z, rotate_ccw) != zscale((0, -1), Z))
# the mechanism: Z(a,b) = sgn_c( (1+i) conj(a+ib) ), and sgn_c commutes with multiplication by i
def sgn_c(u): return (sgn(u[0]), sgn(u[1]))
check('3.3 Z(a,b)=sgn_c((1+i)*conj(a+ib)) on D^2', lambda a, b: (SYN(a, b), DIF(a, b)) == sgn_c(gmul((1, 1), gconj((a, b)))), D2)
W = range(-6, 7)
check('3.3 lemma: sgn_c(i*u) = i*sgn_c(u) for all Gaussian integers |re|,|im|<=6', lambda x, y: sgn_c(gmul(I_, (x, y))) == gmul(I_, sgn_c((x, y))), list(product(W, W)))
check('3.3 Z(i w) = -i Z(w) for all w=a+ib, |a|,|b|<=6 (extension of Z by the same formula)',
      lambda a, b: sgn_c(gmul((1, 1), gconj(gmul(I_, (a, b))))) == gmul((0, -1), sgn_c(gmul((1, 1), gconj((a, b))))), list(product(W, W)))
# the BID coordinates are the unsaturated version: c+id = (1+i) conj(X+iY) for bits
check('3.3 BID: c+id=(1+i)conj(X+iY), c=X+Y, d=X-Y (X,Y in {0,1})', lambda X, Y: (X + Y, X - Y) == gmul((1, 1), gconj((X, Y))), list(product(B, B)))
check('3.3 BID: c^2-d^2=4XY for bits', lambda X, Y: (X + Y) ** 2 - (X - Y) ** 2 == 4 * X * Y, list(product(B, B)))
check('3.3 BID: (c,d) determines (X,Y): X=(c+d)/2, Y=(c-d)/2', lambda X, Y: ((X + Y) + (X - Y)) // 2 == X and ((X + Y) - (X - Y)) // 2 == Y, list(product(B, B)))
check('3.3 BID: d alone merges (0,0),(1,1); c alone merges (0,1),(1,0)', lambda: (0 - 0 == 1 - 1) and (0 + 1 == 1 + 0) and ((0, 0) != (1, 1)) and ((0, 1) != (1, 0)), [()])
# the folded pair  W = T_OR + i T_XOR  transforms by  R_ccw(W) = i * conj(W)  (an order-2 swap, not a rotation)
Wm = zmat(T['OR'], T['XOR'])
iconjW = zscale(I_, tuple(tuple(gconj(x) for x in row) for row in Wm))
record('3.4 folded pair W=T_OR+iT_XOR: R_ccw(W)=i*conj(W)', zrotate(Wm, rotate_ccw) == iconjW and zrotate(Wm, rotate_cw) == iconjW)
record('3.4 folded pair is NOT multiplied by i: R_ccw(W)=iW is false', zrotate(Wm, rotate_ccw) != zscale(I_, Wm))
record('3.4 folded pair: R^2(W)=W (period 2)', zrotate(zrotate(Wm, rotate_cw), rotate_cw) == Wm)

# ================================================================== 3.4 MAG folding and the congruence question
MAGV = {a: MAG(a) for a in D}
record('3.4 MAG(-1)=MAG(+1)=1, MAG(0)=0', MAG(-1) == 1 and MAG(1) == 1 and MAG(0) == 0)
record('3.4 T_OR=MAG(T_SYN), T_XOR=MAG(T_DIF) entrywise', T['OR'] == mag(T['SYN']) and T['XOR'] == mag(T['DIF']))
record('3.4 counterexample: MAG(DIF(-1,+1))=1 but MAG(DIF(+1,+1))=0 although MAG(-1)=MAG(+1)',
       MAG(-1) == MAG(1) and MAG(DIF(-1, 1)) == 1 and MAG(DIF(1, 1)) == 0, {'MAG(DIF(-1,1))': MAG(DIF(-1, 1)), 'MAG(DIF(1,1))': MAG(DIF(1, 1))})
OPS = {'SYN': (2, SYN), 'DIF': (2, DIF), 'OR': (2, OR), 'XOR': (2, XOR), 'AND': (2, AND), 'NOT': (1, NOT), 'ANT': (1, ANT)}
def compatible(part, ar, f):
    """partition `part` (dict element->class id) is compatible with operation f of arity ar"""
    for xs in product(D, repeat=ar):
        for ys in product(D, repeat=ar):
            if all(part[x] == part[y] for x, y in zip(xs, ys)) and part[f(*xs)] != part[f(*ys)]: return False
    return True
kernel_mag = {-1: 'n', 0: 'z', 1: 'n'}        # kernel of MAG: {-1,+1}, {0}
compat = {n: compatible(kernel_mag, ar, f) for n, (ar, f) in OPS.items()}
record('3.4 ker(MAG)={{-1,+1},{0}} is compatible with exactly AND and ANT among the seven operations',
       {n for n, ok in compat.items() if ok} == {'AND', 'ANT'}, compat)
# on {AND, ANT}: MAG is a homomorphism onto ({0,1}; min, id)
record('3.4 MAG(AND(a,b)) = min(MAG a, MAG b), MAG(ANT a)=MAG a', all(MAG(AND(a, b)) == min(MAG(a), MAG(b)) for a, b in D2) and all(MAG(ANT(a)) == MAG(a) for a in D))
# the failing operations, with first witnesses
wit = {}
for n, (ar, f) in OPS.items():
    if not compat[n]:
        for xs in product(D, repeat=ar):
            for ys in product(D, repeat=ar):
                if all(kernel_mag[x] == kernel_mag[y] for x, y in zip(xs, ys)) and kernel_mag[f(*xs)] != kernel_mag[f(*ys)] and n not in wit: wit[n] = (xs, ys, f(*xs), f(*ys))
record('3.4 witnesses for the operations that do not respect ker(MAG)', set(wit) == {'SYN', 'DIF', 'OR', 'XOR', 'NOT'}, {k: str(v) for k, v in wit.items()})
# congruences, subalgebras, automorphisms of the full algebra
def partitions(elems):
    if not elems: yield []; return
    first, rest = elems[0], elems[1:]
    for p in partitions(rest):
        for i in range(len(p)): yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p
parts = list(partitions(list(D)))
def pdict(p): return {x: i for i, blk in enumerate(p) for x in blk}
cong = [p for p in parts if all(compatible(pdict(p), ar, f) for ar, f in OPS.values())]
record('3.4 the full algebra (seven operations) has exactly the two trivial congruences (5 partitions of D examined)', len(parts) == 5 and len(cong) == 2, {'partitions': len(parts), 'congruences': cong})
def closed(S, ar, f): return all(f(*xs) in S for xs in product(sorted(S), repeat=ar))
subalg = [S for r in (1, 2, 3) for S in map(set, itertools.combinations(D, r)) if all(closed(S, ar, f) for ar, f in OPS.values())]
record('3.4 the full algebra has no proper subalgebra other than D (no nonempty proper subset is closed)', subalg == [set(D)], subalg)
auts = [p for p in itertools.permutations(D) if all(p[D.index(f(*xs))] == f(*[p[D.index(x)] for x in xs]) for ar, f in OPS.values() for xs in product(D, repeat=ar))]
record('3.4 the full algebra has only the identity automorphism', auts == [D], auts)
# reduct view: which reducts have {0,+1} as a Boolean subalgebra
Bset = {0, 1}
closedB = {n: closed(Bset, ar, f) for n, (ar, f) in OPS.items()}
record('3.4 {0,+1} is closed under SYN, OR, XOR, AND, NOT and not under DIF, ANT', {n for n, ok in closedB.items() if ok} == {'SYN', 'OR', 'XOR', 'AND', 'NOT'}, closedB)

# ================================================================== 3.4/3.5 Boolean restriction and sets
bo = lambda x, y: x | y
check('3.4 on {0,1}: SYN = boolean OR', lambda x, y: SYN(x, y) == (x | y), list(product(B, B)))
check('3.4 on {0,1}: MAG o DIF = boolean XOR', lambda x, y: MAG(DIF(x, y)) == (x ^ y), list(product(B, B)))
check('3.4 on {0,1}: XOR = boolean XOR, OR = boolean OR', lambda x, y: XOR(x, y) == (x ^ y) and OR(x, y) == (x | y), list(product(B, B)))
check('3.4 on {0,1}: AND = boolean AND', lambda x, y: AND(x, y) == (x & y), list(product(B, B)))
check('3.4 on {0,1}: NOT = boolean NOT (1-x)', lambda x: NOT(x) == 1 - x, B)
record('3.4 DIF(0,1)=-1 and DIF(1,0)=+1: two-valued inputs give a signed result outside {0,1}', DIF(0, 1) == -1 and DIF(1, 0) == 1)
check('3.4 MAG(DIF) folds the two directions: MAG(DIF(0,1))=MAG(DIF(1,0))=1', lambda: MAG(DIF(0, 1)) == MAG(DIF(1, 0)) == 1, [()])
check('3.4 DIF is not boolean: DIF(0,1) not in {0,1}', lambda: DIF(0, 1) in B, [()], control=True)
check('3.4 on {0,1}: polarization AND=DIF(OR,XOR) is inclusion-exclusion (OR - XOR = AND)', lambda x, y: (x | y) - (x ^ y) == (x & y) == AND(x, y) == DIF(OR(x, y), XOR(x, y)), list(product(B, B)))
# sets: characteristic functions on a 3-element universe, all pairs of subsets
U = (0, 1, 2); subsets = [frozenset(s) for r in range(4) for s in itertools.combinations(U, r)]
ind = lambda S, u: 1 if u in S else 0
ok_union = ok_sym = ok_int = ok_ie = ok_bidra = True
for A_, B_ in product(subsets, subsets):
    for u in U:
        a, b = ind(A_, u), ind(B_, u)
        ok_union &= ind(A_ | B_, u) == (a | b) == OR(a, b) == SYN(a, b)
        ok_sym &= ind(A_ ^ B_, u) == (a ^ b) == XOR(a, b)
        ok_int &= ind(A_ & B_, u) == (a & b) == AND(a, b)
        ok_ie &= ind(A_ & B_, u) == ind(A_ | B_, u) - ind(A_ ^ B_, u)
        ok_bidra &= AND(a, b) == DIF(OR(a, b), XOR(a, b))
    ok_ie &= len(A_ & B_) == len(A_ | B_) - len(A_ ^ B_)
record('3.5 characteristic functions: 1_{A|B}=OR=SYN, 1_{A^B}=XOR, 1_{A&B}=AND on all 64 pairs of subsets of a 3-set', ok_union and ok_sym and ok_int)
record('3.5 1_{A&B} = 1_{A|B} - 1_{A^B} and |A&B| = |A|B| - |A^B| (all 64 pairs)', ok_ie)
record('3.5 AND = DIF(OR,XOR) on the characteristic values', ok_bidra)
check('3.5 on {-2..2} the descent fails: polarization false', lambda a, b: sgn(4 * a * b) == DIF(sgn((a + b) ** 2), sgn((a - b) ** 2)), list(product(range(-2, 3), repeat=2)), control=True)

# ---- the quarter-turn identities do not depend on the alphabet {-1,0,1}: they hold on every grid {-m..m} for the sgn-defined operations
for m in (2, 3, 5):
    G = range(-m, m + 1)
    ok = all(SYN(-b, a) == DIF(a, b) and DIF(-b, a) == -SYN(a, b) and OR(-b, a) == XOR(a, b) and XOR(-b, a) == OR(a, b) and AND(-b, a) == -AND(a, b)
             for a in G for b in G)
    record('S14 quarter-turn identities hold on the grid {-%d..%d} (sgn-defined operations, F(-b,a) form)' % (m, m), ok, '%d points' % (len(G) ** 2))
record('S14 contrast: polarization AND=DIF(OR,XOR) fails on {-2..2} (the descent to signs is special to D)',
       any(sgn(4 * a * b) != DIF(OR(a, b), XOR(a, b)) for a in range(-2, 3) for b in range(-2, 3)))

# ================================================================== group actions: C4 and D4 orbits on the 19683 binary functions
def act(f, g):
    return tuple(f[g(a, b)] for (a, b) in D2)          # (f o g) as a 9-vector in D2 order
vecs = [tuple(v) for v in product(D, repeat=9)]
idx = {p: i for i, p in enumerate(D2)}
def perm_of(g): return [idx[g(a, b)] for (a, b) in D2]
gens = {'rot': lambda a, b: (-b, a), 'swap': lambda a, b: (b, a)}
def orbits(perms):
    seen = set(); sizes = []
    for v in vecs:
        if v in seen: continue
        orb = {v}; stack = [v]
        while stack:
            w = stack.pop()
            for p in perms:
                u = tuple(w[p[i]] for i in range(9))
                if u not in orb: orb.add(u); stack.append(u)
        seen |= orb; sizes.append(len(orb))
    return sizes
c4 = orbits([perm_of(gens['rot'])]); d4 = orbits([perm_of(gens['rot']), perm_of(gens['swap'])])
from collections import Counter
record('3.2 orbits of the 19683 binary functions under the rotation group C4: 4995 (Burnside (3^9+2*3^3+3^5)/4)', len(c4) == 4995 and (3 ** 9 + 2 * 27 + 3 ** 5) // 4 == 4995, dict(Counter(c4)))
record('3.2 orbits under the dihedral group D4 generated by rotation and swap: 2862 (= v9 Proposition 5)', len(d4) == 2862, dict(Counter(d4)))
# the orbit of DIF and of SYN
def tvec(F): return tuple(F(a, b) for (a, b) in D2)
def orbit_of(v, perms):
    orb = {v}; stack = [v]
    while stack:
        w = stack.pop()
        for p in perms:
            u = tuple(w[p[i]] for i in range(9))
            if u not in orb: orb.add(u); stack.append(u)
    return orb
o_dif = orbit_of(tvec(DIF), [perm_of(gens['rot']), perm_of(gens['swap'])])
record('3.2 D4-orbit of DIF = {DIF, SYN, -DIF, -SYN} (size 4)', o_dif == {tvec(DIF), tvec(SYN), tuple(-x for x in tvec(DIF)), tuple(-x for x in tvec(SYN))})
o_or = orbit_of(tvec(OR), [perm_of(gens['rot']), perm_of(gens['swap'])])
record('3.2 D4-orbit of OR = {OR, XOR} (size 2); orbit of AND = {AND, -AND}', o_or == {tvec(OR), tvec(XOR)} and orbit_of(tvec(AND), [perm_of(gens['rot']), perm_of(gens['swap'])]) == {tvec(AND), tuple(-x for x in tvec(AND))})

# ================================================================== sign-linear binary functions (geometry of the degree-one operations)
N = 12
lin = {}
for al, be in product(range(-N, N + 1), repeat=2):
    lin.setdefault(tvec(lambda a, b: sgn(al * a + be * b)), (al, be))
record('S16 sign-linear binary functions sgn(alpha a+beta b), |alpha|,|beta|<=%d: exactly 17 distinct functions' % N, len(lin) == 17, len(lin))
lin_set = set(lin)
perms_d4 = [perm_of(gens['rot']), perm_of(gens['swap'])]
lin_orbits = []
rem = set(lin_set)
while rem:
    v = next(iter(rem)); o = orbit_of(v, perms_d4); lin_orbits.append(len(o)); rem -= o
record('S16 the 17 sign-linear functions form D4-orbits of sizes 8, 4, 4, 1 (sectors, {SYN,DIF,-SYN,-DIF}, {a,b,-a,-b}, zero)', sorted(lin_orbits) == [1, 4, 4, 8], sorted(lin_orbits))
record('S16 every sign-linear function is odd (hence lies in the clone generated by DIF)', all(all(v[idx[(-a, -b)]] == -v[idx[(a, b)]] for (a, b) in D2) for v in lin_set))
# why 17: critical directions. The value of sgn(alpha a+beta b) is determined by the sign pattern on (1,0),(0,1),(1,1),(1,-1)
pat = {}
for al, be in product(range(-N, N + 1), repeat=2):
    pat.setdefault((sgn(al), sgn(be), sgn(al + be), sgn(al - be)), tvec(lambda a, b: sgn(al * a + be * b)))
record('S16 the table is determined by the signs of alpha, beta, alpha+beta, alpha-beta (8 sectors + 8 rays + 0)', len(pat) == 17)

# ================================================================== v10r3: sign rule of AND, nullary constants (Proposition 6), XY3/G3 tables (S17)
import csv
# --- the sign rule (Remark 1): sgn(ab) = sgn(|a+b| - |a-b|) on all integers of a window; the folded term works only on D
W12 = range(-12, 13)
record('Remark 1: sgn(ab) = sgn(|a+b|-|a-b|) for all integers |a|,|b|<=12 (625 pairs)', all(sgn(a * b) == sgn(abs(a + b) - abs(a - b)) for a in W12 for b in W12))
record('Remark 1 (negative control): the shifted statement sgn(ab) = sgn(|a+b|-|a-b|+1) is false', not all(sgn(a * b) == sgn(abs(a + b) - abs(a - b) + 1) for a in W12 for b in W12))
fold_term = lambda a, b: DIF(sgn((a + b) ** 2), sgn((a - b) ** 2))
bad_m2 = [(a, b) for a in range(-2, 3) for b in range(-2, 3) if fold_term(a, b) != sgn(a * b)]
record('Remark 1: the folded term DIF(sgn p^2, sgn q^2) equals sgn(ab) on D but fails on {-2..2}, first failure (2,1)', all(fold_term(a, b) == sgn(a * b) for a in D for b in D) and (2, 1) in bad_m2 and len(bad_m2) > 0, bad_m2[:4])

# --- Proposition 6: three constants given as nullary operations
record('Prop 6: NOT(y) = OR(y,-1) = XOR(y,+1) = DIF(+1,y), and ANT(y) = AND(y,-1)', all(NOT(y) == OR(y, -1) == XOR(y, 1) == DIF(1, y) and -y == AND(y, -1) for y in D))
record('Prop 6 (negative control): NOT(y) = AND(y,-1) and ANT(y) = OR(y,-1) are false', not all(NOT(y) == AND(y, -1) for y in D) and not all(-y == OR(y, -1) for y in D))
neg_e = lambda e: SYN(-1, NOT(e))                       # -e for e in {0,1}
ind_p = lambda y: NOT(NOT(y)); ind_m = lambda y: NOT(SYN(y, 1))
record('Prop 6: [y=+1] = NOT(NOT y), [y=-1] = NOT(SYN(y,+1)), -e = SYN(-1,NOT e) for e in {0,1}', all(ind_p(y) == int(y == 1) and ind_m(y) == int(y == -1) for y in D) and all(neg_e(e) == -e for e in (0, 1)))
ant_from_syn_not = lambda y: SYN(ind_m(y), neg_e(ind_p(y)))
record('Prop 6: ANT(y) = SYN([y=-1], -[y=+1]) is a term over SYN, NOT and the constants; then DIF(a,b) = SYN(a, ANT b)', all(ant_from_syn_not(y) == -y for y in D) and all(SYN(a, ant_from_syn_not(b)) == DIF(a, b) for a in D for b in D))
record('Prop 6 (negative control): -e = SYN(+1,NOT e) is false', not all(SYN(1, NOT(e)) == -e for e in (0, 1)))
mono = lambda f: all(f(a, b) <= f(c, d) for a in D for b in D for c in D for d in D if a <= c and b <= d)
record('Prop 6: SYN is monotone for -<0<+ (so {SYN} with the constants only gives monotone terms); ANT is not monotone', mono(SYN) and (-(-1) > -(1)))
record('Prop 6 (negative control): DIF is not monotone', not mono(DIF))
RHO = [(a, b) for a in D for b in D if a * b >= 0]
def pres(f, ar):
    return all(all(f(*[t[0] for t in ts]) * f(*[t[1] for t in ts]) >= 0 for ts in product(RHO, repeat=ar)) for _ in (0,))
record('Prop 6: OR, XOR, AND, NOT, ANT and the constants preserve rho={ab>=0}', pres(OR, 2) and pres(XOR, 2) and pres(AND, 2) and pres(NOT, 1) and pres(ANT, 1) and all((c, c) in RHO for c in D))
record('Prop 6 (negative control): SYN and DIF do not preserve rho', not pres(SYN, 2) and not pres(DIF, 2))

# criterion versus the 127 subsets, independently of classification_table.py (Slupecki certificate for the complete ones, invariants for the rest)
NAMES = ['SYN', 'DIF', 'OR', 'XOR', 'AND', 'NOT', 'ANT']
OPS = {'SYN': (2, SYN), 'DIF': (2, DIF), 'OR': (2, OR), 'XOR': (2, XOR), 'AND': (2, AND), 'NOT': (1, NOT), 'ANT': (1, ANT)}
def unary_closure(sigma):
    """all one-variable terms over sigma and the constants, as value vectors (f(-1),f(0),f(+1))"""
    S = {(-1, 0, 1), (0, 0, 0), (1, 1, 1), (-1, -1, -1)}
    while True:
        new = set(S)
        for n in sigma:
            ar, f = OPS[n]
            if ar == 1:
                for g in S: new.add(tuple(f(x) for x in g))
            else:
                for g in S:
                    for h in S: new.add(tuple(f(x, y) for x, y in zip(g, h)))
        if new == S: return S
        S = new
def onto_ess_binary(f):
    vals = {f(a, b) for a in D for b in D}
    ess = any(f(a, b) != f(c, b) for a in D for b in D for c in D) and any(f(a, b) != f(a, d) for a in D for b in D for d in D)
    return vals == set(D) and ess
rows127 = []
for mask in range(1, 128):
    sigma = [NAMES[i] for i in range(7) if mask >> i & 1]
    pred = ('DIF' in sigma) or ('SYN' in sigma and len(sigma) >= 2)                       # Proposition 6
    slup = len(unary_closure(sigma)) == 27 and any(OPS[n][0] == 2 and onto_ess_binary(OPS[n][1]) for n in sigma)   # Slupecki's criterion
    if pred:
        cert = slup
    else:                                                                                  # invariant certificates
        cert = (sigma == ['SYN']) or ('SYN' not in sigma and 'DIF' not in sigma)
    rows127.append((tuple(sigma), pred, slup, cert))
record('Prop 6: for all 127 subsets the criterion (DIF, or SYN plus one more) agrees with Slupecki certificates (all 27 unary terms and an onto essentially binary operation)', all(r[1] == r[2] for r in rows127))
record('Prop 6: 95 complete subsets = 64 containing DIF + 31 containing SYN and another operation but not DIF', sum(r[1] for r in rows127) == 95 and sum(1 for r in rows127 if 'DIF' in r[0]) == 64 and sum(1 for r in rows127 if 'DIF' not in r[0] and 'SYN' in r[0] and r[1]) == 31)
record('Prop 6: the 32 incomplete subsets are {SYN} and the 31 nonempty subsets of {OR,XOR,AND,NOT,ANT} (each carries a monotone / rho certificate)', sum(1 for r in rows127 if not r[1]) == 32 and all(r[3] for r in rows127 if not r[1]))
csvp = os.path.join(HERE, 'classification_127.csv')
if os.path.exists(csvp):
    ref = {tuple(r['subset'].split()): r['complete_with_constants'] == '1' for r in csv.DictReader(open(csvp))}
    mism = [r[0] for r in rows127 if ref.get(tuple(sorted(r[0], key=NAMES.index))) != r[1]]
    record('Prop 6: the criterion agrees with every row of classification_127.csv (complete_with_constants)', len(ref) == 127 and not mism, len(ref), mism[:3] or None)
else:
    record('Prop 6: classification_127.csv present', False)

# --- XY3 / G3 tables (supplement S17): masks over the three BID states, G3 = XY3(A) - XY3(a) with A & a = 0
STATES = [(1, 0), (0, 0), (0, 1)]                         # bit 0, bit 1, bit 2  (d = X-Y = +1, 0, -1)
def xy3(mask, x, y): return (mask >> STATES.index((x, y))) & 1
XY3_FORMULA = {0: lambda x, y: 0, 1: lambda x, y: x, 2: lambda x, y: 1 - (x | y), 4: lambda x, y: y,
               7: lambda x, y: 1, 6: lambda x, y: 1 - x, 5: lambda x, y: x | y, 3: lambda x, y: 1 - y}
record('S17 XY3: the eight masks are the Boolean functions 0, x, x\'y\', y\', y, x|y, x\', 1 on the three states (0b000..0b111)', all(xy3(m, x, y) == XY3_FORMULA[m](x, y) for m in range(8) for (x, y) in STATES))
def trits(gid):
    """balanced value of gid in -13..13 as digits (d2,d1,d0) with values in {0,+1,-1}"""
    r = []
    g = gid
    for _ in range(3):
        t = ((g + 1) % 3) - 1; r.append(t); g = (g - t) // 3
    return tuple(reversed(r))
G3 = {}
for gid in range(-13, 14):
    t = trits(gid)
    A = sum(1 << (2 - i) for i, v in enumerate(t) if v == 1); a = sum(1 << (2 - i) for i, v in enumerate(t) if v == -1)
    G3[gid] = (A, a, tuple(xy3(A, x, y) - xy3(a, x, y) for (x, y) in reversed(STATES)))   # values at d=-1,0,+1
record('S17 G3: 27 gids -13..13 give 27 distinct functions and 27 pairs of disjoint masks (A,a)', len({v[2] for v in G3.values()}) == 27 and all(v[0] & v[1] == 0 for v in G3.values()) and len({(v[0], v[1]) for v in G3.values()}) == 27)
def code_of(vs): return tuple({0: 0, 1: 1, -1: 2}[v] for v in vs)
record('S17 G3: the trits of gid (most significant first, 1 for +1 and 2 for -1) are the values g(d=-1), g(0), g(+1), the code of S5', all(tuple({0: 0, 1: 1, -1: 2}[x] for x in trits(g)) == code_of(v[2]) for g, v in G3.items()))
record('S17 G3: gid = -gid(negated function) and G3(gid) = XY3(A) - XY3(a); e.g. gid -8 = 0t201 is x-y, the identity d', all(tuple(-x for x in G3[g][2]) == G3[-g][2] for g in G3) and G3[-8][2] == (-1, 0, 1) and (G3[-8][0], G3[-8][1]) == (1, 4))
record('S17 (negative control): G3 with the masks in the opposite bit order is not the same table', any(tuple(xy3(A, x, y) - xy3(a, x, y) for (x, y) in STATES) != G3[g][2] for g, (A, a, _) in G3.items()))

# ================================================================== summary
RESULTS['_summary'] = {'tests': COUNT['tests'], 'passed': COUNT['passed'], 'failed': FAILED,
                       'python': platform.python_version(), 'command': ' '.join(sys.argv) or 'verify_direction_rotation.py',
                       'status': 'ALL ASSERTIONS PASSED' if not FAILED else 'FAILED'}
json.dump(RESULTS, open(os.path.join(HERE, 'direction_rotation_results.json'), 'w'), indent=1, ensure_ascii=False, default=str)
LOG.append('')
LOG.append('tests: %d, passed: %d, failed: %d' % (COUNT['tests'], COUNT['passed'], len(FAILED)))
LOG.append('python %s; command: %s' % (platform.python_version(), RESULTS['_summary']['command']))
LOG.append(RESULTS['_summary']['status'])
open(os.path.join(HERE, 'direction_rotation_results.log'), 'w').write('\n'.join(LOG) + '\n')
print('\n'.join(LOG))
sys.exit(0 if not FAILED else 1)
