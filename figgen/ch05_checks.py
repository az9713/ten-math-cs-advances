"""Ten exact checks for chapter 5 (permanent lower bounds).

Stdlib only; exact integer / Fraction / Z[zeta] arithmetic throughout.
This same listing is embedded in ch05 section 15.
"""
import itertools
import random
from fractions import Fraction

random.seed(5)
OK = []


def check(name, cond):
    OK.append((name, bool(cond)))
    print(('PASS' if cond else 'FAIL'), name)


def per(a):
    """Ryser's formula with Gray-code updates; exact."""
    n = len(a)
    if n == 0:
        return 1
    rs = [0]*n
    tot = 0
    gray_prev = 0
    for g in range(1, 1 << n):
        gray = g ^ (g >> 1)
        diff = gray ^ gray_prev
        b = diff.bit_length() - 1
        sgn = 1 if gray & diff else -1
        for i in range(n):
            rs[i] += sgn*a[i][b]
        p = 1
        for i in range(n):
            p *= rs[i]
        kbits = bin(gray).count('1')
        tot += -p if kbits % 2 else p
        gray_prev = gray
    return -tot if n % 2 else tot


def per_perm(a):
    n = len(a)
    t = 0
    for s in itertools.permutations(range(n)):
        p = 1
        for i in range(n):
            p *= a[i][s[i]]
        t += p
    return t


def det_perm(a):
    n = len(a)
    t = 0
    for s in itertools.permutations(range(n)):
        sg = sign_of(s)
        p = sg
        for i in range(n):
            p *= a[i][s[i]]
        t += p
    return t


def sign_of(s):
    sg, seen = 1, [False]*len(s)
    for i in range(len(s)):
        if not seen[i]:
            j, c = i, 0
            while not seen[j]:
                seen[j] = True
                j = s[j]
                c += 1
            if c % 2 == 0:
                sg = -sg
    return sg


# --- check 1: 4x4 permanent cancellation (symbolic, monomial dicts) ---
def sym_per4(mat):
    out = {}
    for s in itertools.permutations(range(4)):
        key, coef = [], 1
        for i in range(4):
            e = mat[i][s[i]]
            if isinstance(e, str):
                key.append(e)
            else:
                coef *= e
        key = tuple(sorted(key))
        out[key] = out.get(key, 0) + coef
    return {k: v for k, v in out.items() if v}


M1 = [['u', 'v', 1, 1], ['w', 'z', 1, 1],
      ['p', 'q', 2, -2], ['r', 's', 2, -2]]
want = {('u', 'z'): -8, ('v', 'w'): -8, ('p', 's'): 2, ('q', 'r'): 2}
check('C1  per 4x4 = -8(uz+vw)+2(ps+qr)', sym_per4(M1) == want)


# --- check 2: 4x4 determinant analogue collapses ---
def sym_det4(mat):
    out = {}
    for s in itertools.permutations(range(4)):
        key, coef = [], sign_of(s)
        for i in range(4):
            e = mat[i][s[i]]
            if isinstance(e, str):
                key.append(e)
            else:
                coef *= e
        key = tuple(sorted(key))
        out[key] = out.get(key, 0) + coef
    return {k: v for k, v in out.items() if v}


w2 = {('q', 'u'): 4, ('s', 'u'): -4, ('q', 'w'): -4, ('s', 'w'): 4,
      ('p', 'v'): -4, ('r', 'v'): 4, ('p', 'z'): 4, ('r', 'z'): -4}
check('C2  det 4x4 = 4((u-w)(q-s)-(v-z)(p-r))', sym_det4(M1) == w2)

# --- check 3: derivative identity (5.4) at t=4, s=5, d=3 ---
t, s_, d = 4, 5, 3
q = d - 1
X = [[random.randint(-9, 9) for _ in range(s_)] for _ in range(t)]


def m_tsd(x, tt, ss, dd):
    tot = 0
    for ii in itertools.combinations(range(tt), dd):
        for jj in itertools.combinations(range(ss), dd):
            tot += per_perm([[x[i][j] for j in jj] for i in ii])
    return tot


def partitions(seq):
    if not seq:
        yield []
        return
    first, rest = seq[0], seq[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i+1:]
        yield [[first]] + p


def fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


def q_i(i, a, x):
    pb = {}
    for bset in range(1, 1 << t):
        rows = [j for j in range(t) if bset >> j & 1]
        pb[bset] = sum(prod(x[j][c] for j in rows) for c in range(s_))
    tot = 0
    for ii in itertools.combinations([j for j in range(t) if j != i], q):
        for part in partitions(list(ii)):
            mu = 1
            for blk in part:
                mu *= (-1)**(len(blk) - 1)*fact(len(blk) - 1)
            term = mu
            for blk in part:
                bs = sum(1 << j for j in blk)
                term *= pb[bs] - prod(x[j][a] for j in blk)
            tot += term
    return tot


def prod(it):
    r = 1
    for v in it:
        r *= v
    return r


ok3 = True
for i in range(t):
    for a in range(s_):
        x1 = [row[:] for row in X]
        x1[i][a] = 1
        x0 = [row[:] for row in X]
        x0[i][a] = 0
        dM = m_tsd(x1, t, s_, d) - m_tsd(x0, t, s_, d)
        ok3 = ok3 and dM == q_i(i, a, X)
check('C3  (5.4): dM/dx_ia == Q_i(column; power sums), all i,a', ok3)

# --- check 4: symmetric identities (7.c) at t=6, q=3 ---
u6 = [random.randint(-9, 9) for _ in range(6)]


def esym(vals, j):
    tot = 0
    for c in itertools.combinations(range(len(vals)), j):
        tot += prod(vals[i] for i in c)
    return tot


ok4 = True
for i in range(6):
    hat = [u6[j] for j in range(6) if j != i]
    rhs = sum((-u6[i])**(3 - j)*esym(u6, j) for j in range(4))
    ok4 = ok4 and esym(hat, 3) == rhs
tot = sum(esym([u6[j] for j in range(6) if j != i], 3) for i in range(6))
ok4 = ok4 and tot == (6 - 3)*esym(u6, 3)
check('C4  (7.c): both elementary-symmetric identities', ok4)


# --- check 5: Lemma 6.1 at b=2, t=d=3, s=4 over Z[zeta] ---
def zmul(x, y):
    a, b = x
    c, e = y
    return (a*c - b*e, a*e + b*c - b*e)


def zper(mat):
    n = len(mat)
    tot = (0, 0)
    for sg in itertools.permutations(range(n)):
        p = (1, 0)
        for i in range(n):
            p = zmul(p, mat[i][sg[i]])
        tot = (tot[0] + p[0], tot[1] + p[1])
    return tot


zeta = [(1, 0), (0, 1), (-1, -1)]
XB = [[random.randint(-4, 4) for _ in range(4)] for _ in range(6)]
ucols = []
for j in range(3):
    col = [(1, 0)]*3 + [zmul((2, 0), zeta[j])]*3
    ucols.append(col)
B7 = []
for i in range(6):
    B7.append([(XB[i][a], 0) for a in range(4)]
              + [ucols[j][i] for j in range(3)])
B7.append([(1, 0)]*4 + [(0, 0)]*3)
lhs = zper(B7)
m1 = m_tsd([XB[i] for i in range(3)], 3, 4, 3)
m2 = m_tsd([XB[i] for i in range(3, 6)], 3, 4, 3)
rhs = 6*(8*m1 + 1*m2)
ok5 = lhs == (rhs, 0)
for ii in itertools.combinations(range(6), 3):
    ic = [j for j in range(6) if j not in ii]
    mv = zper([[ucols[j][i] for j in range(3)] for i in ic])
    if set(ii) == {0, 1, 2}:
        ok5 = ok5 and mv == (48, 0)
    elif set(ii) == {3, 4, 5}:
        ok5 = ok5 and mv == (6, 0)
    else:
        ok5 = ok5 and mv == (0, 0)
check('C5  Lem 6.1: per(B)=6(8*M1+M2); U-minors 48/6/0', ok5)

# --- check 6: 5x5 warm-up Jacobian ---
p3, q3 = [2, 3, 5], [7, 11, 13]


def b5(ye, yf, wij):
    return [[ye, 0] + p3,
            [0, yf, 1, 1, 1],
            [1, q3[0]] + wij[0],
            [1, q3[1]] + wij[1],
            [1, q3[2]] + wij[2]]


def coef5(se, sf, wij):
    rows = [i for i in range(5) if not (i == 0 and se or i == 1 and sf)]
    mat = [[b5(0, 0, wij)[i][j] for j in rows] for i in rows]
    return per_perm(mat)


def dcoef(se, sf, a, b):
    w1 = [[1]*3 for _ in range(3)]
    w1[a][b] = 2
    w0 = [[1]*3 for _ in range(3)]
    w0[a][b] = 0
    return (coef5(se, sf, w1) - coef5(se, sf, w0))//2


gb = [sum(p3) - p3[b] for b in range(3)]
ha = [sum(q3) - q3[a] for a in range(3)]
ok6 = True
for a in range(3):
    for b in range(3):
        ok6 = ok6 and dcoef(1, 1, a, b) == 2
        ok6 = ok6 and dcoef(0, 1, a, b) == 2*gb[b]
        ok6 = ok6 and dcoef(1, 0, a, b) == 2*ha[a]
        ok6 = ok6 and dcoef(0, 0, a, b) == gb[b]*ha[a]
jm = []
for se, sf in [(1, 1), (0, 1), (1, 0), (0, 0)]:
    jm.append([dcoef(se, sf, a, b)
               for a in range(2) for b in range(2)])
dj = det_perm(jm)
ok6 = ok6 and dj == 8*(p3[0] - p3[1])**2*(q3[0] - q3[1])**2
check('C6  5x5: entries 2,2Gb,2Ha,GbHa; detJ=8(p1-p2)^2(q1-q2)^2', ok6)

# --- check 7: Lemma 9.1 end-to-end at ell=3, k=6, m=7, n=13 ---
ell, kk, mm = 3, 6, 7
nn = kk + mm
pv = [2, 3, 5, 7, 11, 13, 17]
qv = [19, 23, 29, 31, 37, 41, 43]


def big_matrix():
    mat = [[0]*nn for _ in range(nn)]
    for u in range(ell):
        for b in range(mm):
            mat[u][kk + b] = pv[b]**(2**u)
    for v in range(ell):
        for b in range(mm):
            mat[ell + v][kk + b] = 1
    for a in range(mm):
        for u in range(ell):
            mat[kk + a][u] = 1
        for v in range(ell):
            mat[kk + a][ell + v] = qv[a]**(2**v)
    for a in range(mm):
        for b in range(mm):
            mat[kk + a][kk + b] = 1
    return mat


def jac_entry(al, be, a, b):
    tset = ([u for u in range(ell) if al >> u & 1]
            + [ell + v for v in range(ell) if be >> v & 1])
    keep = tset + [kk + i for i in range(mm) if i != a]
    keepc = tset + [kk + j for j in range(mm) if j != b]
    mat = big_matrix()
    sub = [[mat[i][j] for j in keepc] for i in keep]
    return per(sub)


def inj_sum(ws, nodes):
    tot = 0
    for pr in itertools.permutations(nodes, len(ws)):
        tot += prod(pr[i]**ws[i] for i in range(len(ws)))
    return tot


def gpoly(ws, nodes):
    deg = sum(ws)
    if not ws:
        return [1]
    h = inj_sum(ws, nodes)
    co = [0]*(deg + 1)
    co[0] = h
    for i in range(len(ws)):
        rest = ws[:i] + ws[i+1:]
        sub = gpoly(rest, nodes)
        for j, cj in enumerate(sub):
            co[j + ws[i]] -= cj
    return co


def falling(mv, sv):
    r = 1
    for i in range(sv):
        r *= mv - i
    return r


J49 = [[jac_entry(al, be, a, b)
        for a in range(mm) for b in range(mm)]
       for al in range(mm) for be in range(mm)]
ok7 = True
for al in range(mm):
    for be in range(mm):
        pw = [2**u for u in range(ell) if al >> u & 1]
        qw = [2**v for v in range(ell) if be >> v & 1]
        gp = gpoly(pw, pv)
        hq = gpoly(qw, qv)
        lpq = (fact(mm - 1 - len(pw) - len(qw))
               * falling(mm - 1 - len(pw), len(qw))
               * falling(mm - 1 - len(qw), len(pw)))
        for a in range(mm):
            for b in range(mm):
                gv = sum(c*pv[b]**j for j, c in enumerate(gp))
                hv = sum(c*qv[a]**j for j, c in enumerate(hq))
                jv = J49[al*mm + be][a*mm + b]
                ok7 = ok7 and jv == lpq*gv*hv
check('C7a (9.9): all 49x49 Jacobian entries factor as L*g*h', ok7)


def bareiss_det(mat):
    m2 = [[Fraction(v) for v in row] for row in mat]
    n = len(m2)
    sg = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if m2[r][c]), None)
        if piv is None:
            return 0
        if piv != c:
            m2[c], m2[piv] = m2[piv], m2[c]
            sg = -sg
        for r in range(c + 1, n):
            f = m2[r][c]/m2[c][c]
            for j in range(c, n):
                m2[r][j] -= f*m2[c][j]
    r = Fraction(sg)
    for i in range(n):
        r *= m2[i][i]
    return r


check('C7b det of the 49x49 Jacobian is nonzero', bareiss_det(J49) != 0)

# --- check 8: packing + Theorem 1.2 numbers at n=32 ---
n8, l8 = 32, 5
k8, m8 = 2*l8, 32 - 2*l8
cells = {}
ok8 = True
for tau in range(n8):
    for j in range(n8//k8):
        for r in range(k8):
            cell = (j*k8 + r, (j*k8 + r + tau) % n8)
            ok8 = ok8 and cell not in cells
            cells[cell] = (tau, j)
nu = n8*(n8//k8)
ok8 = ok8 and nu == 96 and 2*nu*k8 >= n8*n8
lb = nu*m8*m8//4
ok8 = ok8 and lb == 11616 and lb*128*l8 >= n8**4
check('C8  packing n=32: disjoint, nu=96, sum >= n^4/(128 log n)', ok8)

# --- check 9: Theorem 1.1 chain at n=2^16 ---
n9 = 2**16
ell9 = 16
d9 = ell9//4
t9 = 4*d9
b9 = n9//(2*t9)
r9 = b9*t9
s9 = n9 - r9 + d9
m9 = r9*s9
k9 = m9//4
ok9 = d9 >= 3 and n9 >= 6*t9 and 4*(2**t9 - 1) <= t9*s9
ok9 = ok9 and 3*r9 >= n9 and 2*r9 <= n9 and 2*s9 >= n9
ok9 = ok9 and 6*m9 >= n9*n9 and 48*k9 >= n9*n9
big_d = b9*(2**t9 - 1 + s9*(d9 - 2))
ok9 = ok9 and k9 + big_d < m9 and 2*big_d < m9
bound = n9*n9*1//144
ok9 = ok9 and bound == 29826161
check('C9  Thm 1.1 at n=2^16: (7.2) holds, k+D<m, bound=29826161', ok9)

# --- check 10: Schur coefficient identity, k=3 matching in 7x7 ---
n10, k10 = 7, 3
A10 = [[Fraction(random.randint(-9, 9)) for _ in range(n10)]
       for _ in range(n10)]
for i in range(k10):
    A10[i][i] = Fraction(0)
W10 = [[A10[k10 + i][k10 + j] for j in range(n10 - k10)]
       for i in range(n10 - k10)]
delta = bareiss_det(W10)


def mat_inv(mat):
    n = len(mat)
    aug = [row[:] + [Fraction(int(i == j)) for j in range(n)]
           for i, row in enumerate(mat)]
    for c in range(n):
        piv = next(r for r in range(c, n) if aug[r][c])
        aug[c], aug[piv] = aug[piv], aug[c]
        f = aug[c][c]
        aug[c] = [v/f for v in aug[c]]
        for r in range(n):
            if r != c and aug[r][c]:
                g = aug[r][c]
                aug[r] = [aug[r][j] - g*aug[c][j] for j in range(2*n)]
    return [row[n:] for row in aug]


Wi = mat_inv(W10)
S10 = [[A10[i][j] - sum(A10[i][k10 + x]*Wi[x][y]*A10[k10 + y][j]
                        for x in range(n10 - k10)
                        for y in range(n10 - k10))
        for j in range(k10)] for i in range(k10)]
ok10 = True
for tt in range(8):
    tset = [i for i in range(k10) if tt >> i & 1]
    rows = [i for i in range(n10) if i not in tset]
    mat = [[A10[i][j] for j in rows] for i in rows]
    lhs = bareiss_det(mat)
    comp = [i for i in range(k10) if i not in tset]
    smat = [[S10[i][j] for j in comp] for i in comp]
    ok10 = ok10 and lhs == delta*bareiss_det(smat)
check('C10 det Schur: [prod_T y]detX = delta*det S_comp, all T', ok10)

print('---')
print('%d/%d checks pass' % (sum(c for _, c in OK), len(OK)))
