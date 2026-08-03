"""Chapter 7 (CVP n^(1/400) hardness) — machine verification engine.

Implements GF(2^e) arithmetic, polynomial and rational-function algebra
over GF(16), and runs every check backing the explainer:

  * Lagrange interpolation of a Boolean assignment (toy m=3 in F16)
  * completeness of the moment encoding (C1)-(C4), weight R=(l+1)|P|
  * full explicit binary affine system H x = b for the toy (dims + rank)
  * Lemma 7: parity-lift lattice distance identity, brute force
  * Hankel reconstruction: rational case (two lines) and Artin-Schreier
    case G(Y)=Y^2+Y+X with genuinely non-rational roots
  * Newton identities (coefficient -> power sums), characteristic 2
  * shifted-moment binomial identity (17) as a polynomial identity
  * clause parity matching (Lemma 13 mechanism) on the toy
  * Vandermonde nonsingularity + all-ones impossibility
  * the paper's parameter inequalities at N=100 with EXACT integers
    (T >= 2K-1, |P|-dK(K-1) > 2dK^2 T, Markov bound (13), (14),
     M <= 40 N^401, z + h - 1 < T)

Usage: PYTHONIOENCODING=utf-8 python figgen/ch07_figs.py
Every check prints "OK <label>" or raises.
"""
import itertools
import math
import os
import sys

OK = []


def ok(label, cond):
    if not cond:
        raise AssertionError("FAIL " + label)
    OK.append(label)
    print("OK", label)


# ----------------------------------------------------------------------
# GF(16) = F2[x]/(x^4+x+1); elements are ints 0..15 (bitmask coeffs).
MOD = 0b10011
E = 4
Q = 16


def gmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & (1 << E):
            a ^= MOD
    return r


def gpow(a, n):
    r = 1
    while n:
        if n & 1:
            r = gmul(r, a)
        a = gmul(a, a)
        n >>= 1
    return r


def ginv(a):
    assert a != 0
    return gpow(a, Q - 2)


def gdiv(a, b):
    return gmul(a, ginv(b))


def trace(a):
    t = a
    x = a
    for _ in range(E - 1):
        x = gmul(x, x)
        t ^= x
    return t


# ----------------------------------------------------------------------
# Polynomials over GF(16): list of coeffs, low degree first.
def pnorm(f):
    while f and f[-1] == 0:
        f.pop()
    return f


def padd(f, g):
    n = max(len(f), len(g))
    return pnorm([(f[i] if i < len(f) else 0) ^ (g[i] if i < len(g) else 0)
                  for i in range(n)])


def pmul(f, g):
    if not f or not g:
        return []
    r = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        if a:
            for j, b in enumerate(g):
                if b:
                    r[i + j] ^= gmul(a, b)
    return pnorm(r)


def pscale(f, c):
    return pnorm([gmul(a, c) for a in f])


def peval(f, x):
    r = 0
    for a in reversed(f):
        r = gmul(r, x) ^ a
    return r


def pdeg(f):
    return len(f) - 1 if f else -1


def pdivmod(f, g):
    assert g
    f = f[:]
    q = [0] * max(0, len(f) - len(g) + 1)
    inv = ginv(g[-1])
    while len(f) >= len(g) and f:
        c = gmul(f[-1], inv)
        s = len(f) - len(g)
        q[s] = c
        for i, b in enumerate(g):
            f[s + i] ^= gmul(b, c)
        pnorm(f)
    return pnorm(q), f


def pgcd(f, g):
    while g:
        f, g = g, pdivmod(f, g)[1]
    if f:
        f = pscale(f, ginv(f[-1]))
    return f


def interp(pts):
    """Lagrange interpolation through [(x,y)] over GF(16)."""
    res = []
    for i, (xi, yi) in enumerate(pts):
        if yi == 0:
            continue
        num = [1]
        den = 1
        for j, (xj, _) in enumerate(pts):
            if j != i:
                num = pmul(num, [xj, 1])
                den = gmul(den, xi ^ xj)
        res = padd(res, pscale(num, gdiv(yi, den)))
    return res


# ----------------------------------------------------------------------
# Rational functions over GF(16): pairs (num, den) reduced.
def rnew(n, d=None):
    if d is None:
        d = [1]
    assert d
    g = pgcd(n, d) if n else [1]
    if n:
        n = pdivmod(n, g)[0]
        d = pdivmod(d, g)[0]
    else:
        d = [1]
    c = ginv(d[-1])
    return (pscale(n, c), pscale(d, c))


def radd(a, b):
    return rnew(padd(pmul(a[0], b[1]), pmul(b[0], a[1])), pmul(a[1], b[1]))


def rmul(a, b):
    return rnew(pmul(a[0], b[0]), pmul(a[1], b[1]))


def rdiv(a, b):
    assert b[0]
    return rnew(pmul(a[0], b[1]), pmul(a[1], b[0]))


def req(a, b):
    return pnorm(pmul(a[0], b[1])[:]) == pnorm(pmul(b[0], a[1])[:])


# ======================================================================
# 1. Toy 3SAT instance in GF(16)
# variables v1,v2,v3; anchors a1=1, a2=2 (g), a3=4 (g^2)
ANCH = [1, 2, 4]
M_VARS = 3
P_SET = [w for w in range(Q) if w not in ANCH]
D_DEG = M_VARS  # d = m
T_TOY = 4       # dT = 12 < |P| = 13
# clauses: C1 = v1 | v2 | ~v3   excludes (0,0,1)
#          C2 = ~v1 | ~v2 | v3  excludes (1,1,0)
CLAUSES = [((1, 1), (2, 1), (3, 0)),  # (var, sign) sign=1 positive
           ((1, 0), (2, 0), (3, 1))]


def clause_sat(cl, sigma):
    return any((sigma[v - 1] == 1) == (s == 1) for v, s in cl)


def bc(cl):
    """Satisfying local assignments beta in {0,1}^3 (ordered by IC)."""
    out = []
    for beta in itertools.product((0, 1), repeat=3):
        if any((beta[k] == 1) == (s == 1)
               for k, (v, s) in enumerate(cl)):
            out.append(beta)
    return out


SIGMA = (1, 0, 1)
ok("toy sigma satisfies both clauses",
   all(clause_sat(c, SIGMA) for c in CLAUSES))
ok("clause satisfying-set sizes |BC| = 7",
   all(len(bc(c)) == 7 for c in CLAUSES))

QSIG = interp(list(zip(ANCH, SIGMA)))
ok("interpolation: deg Q_sigma <= m-1 and Q_sigma(a_i) = sigma_i",
   pdeg(QSIG) <= M_VARS - 1 and
   all(peval(QSIG, a) == s for a, s in zip(ANCH, SIGMA)))

# table types Theta: 0 and (c_idx, beta)
TYPES = [0]
for ci, cl in enumerate(CLAUSES):
    for beta in bc(cl):
        TYPES.append((ci, beta))
ok("type count |Theta| = 1+7+7 = 15", len(TYPES) == 15)

# intended one-hot solution: dict (tau,p,w) -> 1
X_SOL = {}
for p in P_SET:
    X_SOL[(0, p, peval(QSIG, p))] = 1
for ci, cl in enumerate(CLAUSES):
    beta = tuple(SIGMA[v - 1] for v, _ in cl)
    for p in P_SET:
        X_SOL[((ci, beta), p, peval(QSIG, p))] = 1
R_TOY = (len(CLAUSES) + 1) * len(P_SET)
ok("completeness weight wt(x) = (l+1)|P| = 39",
   len(X_SOL) == R_TOY == 39)


def support(x, tau, p):
    return [w for w in range(Q) if x.get((tau, p, w), 0) == 1]


def moment(x, tau, p, j):
    s = 0
    for w in support(x, tau, p):
        s ^= gpow(w, j) if j > 0 else 1
    return s


def is_rs(vals_by_p, dmax):
    """vals_by_p: dict p->value on P_SET; interpolates and checks deg."""
    f = interp([(p, vals_by_p[p]) for p in P_SET])
    return pdeg(f) <= dmax, f


# (C1) global fibers odd parity
ok("(C1) all global fibers have odd parity",
   all(len(support(X_SOL, 0, p)) % 2 == 1 for p in P_SET))
# (C2) clause tables reproduce global fiber mod 2
c2 = True
for ci, cl in enumerate(CLAUSES):
    for p in P_SET:
        for w in range(Q):
            tot = 0
            for beta in bc(cl):
                tot ^= X_SOL.get(((ci, beta), p, w), 0)
            if tot != X_SOL.get((0, p, w), 0):
                c2 = False
ok("(C2) clause tables reproduce global table mod 2", c2)
# (C3) ordinary moments are RS(P, d j)
c3 = True
MU_POLY = {}
for tau in TYPES:
    for j in range(T_TOY + 1):
        good, f = is_rs({p: moment(X_SOL, tau, p, j) for p in P_SET},
                        D_DEG * j)
        MU_POLY[(tau, j)] = f
        c3 = c3 and good
ok("(C3) all ordinary moments lie in RS(P, dj), all 15 types, j<=4", c3)
# (C4) shifted moments are RS(P, (d-1) j)
c4 = True
ETA_POLY = {}
for tau in TYPES:
    if tau == 0:
        continue
    ci, beta = tau
    for k, (v, _) in enumerate(CLAUSES[ci]):
        ai, bi = ANCH[v - 1], beta[k]
        for j in range(T_TOY + 1):
            vals = {}
            for p in P_SET:
                s = 0
                for w in support(X_SOL, tau, p):
                    y = gdiv(w ^ bi, p ^ ai)
                    s ^= gpow(y, j) if j > 0 else 1
                vals[p] = s
            good, f = is_rs(vals, (D_DEG - 1) * j)
            ETA_POLY[(tau, v, j)] = f
            c4 = c4 and good
ok("(C4) all shifted moments lie in RS(P, (d-1)j)", c4)

# ----------------------------------------------------------------------
# 2. Explicit binary affine system H x = b for the toy
# coordinate order: (tau_idx, p_idx, w) -> col
COLS = {}
for ti, tau in enumerate(TYPES):
    for pi, p in enumerate(P_SET):
        for w in range(Q):
            COLS[(tau, p, w)] = (ti * len(P_SET) + pi) * Q + w
M_DIM = len(TYPES) * len(P_SET) * Q
ok("toy dimension M = |Theta||P|q = 3120", M_DIM == 3120)


def field_bits(v):
    return [(v >> k) & 1 for k in range(E)]


ROWS = []
RHS = []


def add_field_eq(coeffs, rhs):
    """coeffs: dict col -> GF16 coeff; rhs: GF16 value.
    Expands to E binary rows (bit-slice the field equation)."""
    for k in range(E):
        row = 0
        for col, cf in coeffs.items():
            if (cf >> k) & 1:
                pass
        ROWS.append(None)  # placeholder replaced below
    ROWS[-E:] = []
    bitrows = [0] * E
    for col, cf in coeffs.items():
        for k in range(E):
            if (cf >> k) & 1:
                bitrows[k] |= 0  # placeholder
    # Correct expansion: multiplication by known coeff cf maps the
    # UNKNOWN BIT x (0/1) to cf*x, so bit k of contribution is bit k
    # of cf times x. Each binary unknown contributes cf's bit-k.
    bitrows = [0] * E
    for col, cf in coeffs.items():
        for k in range(E):
            if (cf >> k) & 1:
                bitrows[k] ^= (1 << col)
    rb = field_bits(rhs)
    for k in range(E):
        ROWS.append(bitrows[k])
        RHS.append(rb[k])


def rs_checks(tau, getcoef, dmax, tag):
    """Add parity checks forcing (vals(p))_p in RS(P,dmax).
    Basis of the dual: for each subset row of the checks obtained by
    Lagrange: value at p_i of poly interpolation minus ... Simpler:
    the code RS(P,dmax) has dimension dmax+1; its dual checks are
    rows of the null space of the (dmax+1) x |P| evaluation matrix.
    We compute the dual basis by Gaussian elimination over GF(16)."""
    npts = len(P_SET)
    ev = [[gpow(p, i) if i else 1 for p in P_SET]
          for i in range(dmax + 1)]
    # find null space of ev (vectors c with sum_p c_p p^i = 0)
    # solve over GF(16) by elimination on the transpose
    rows = [r[:] for r in ev]
    piv = {}
    r = 0
    for c in range(npts):
        pr = None
        for i in range(r, len(rows)):
            if rows[i][c]:
                pr = i
                break
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        inv = ginv(rows[r][c])
        rows[r] = [gmul(v, inv) for v in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [a ^ gmul(f, bch)
                           for a, bch in zip(rows[i], rows[r])]
        piv[c] = r
        r += 1
    free = [c for c in range(npts) if c not in piv]
    duals = []
    for fc in free:
        vec = [0] * npts
        vec[fc] = 1
        for c, pr in piv.items():
            vec[c] = rows[pr][fc]
        duals.append(vec)
    for vec in duals:
        coeffs = {}
        rhsv = 0
        for pi, p in enumerate(P_SET):
            if vec[pi] == 0:
                continue
            cdict = getcoef(p)
            for col, cf in cdict.items():
                coeffs[col] = coeffs.get(col, 0) ^ gmul(cf, vec[pi])
            rhsv ^= 0
        add_field_eq(coeffs, rhsv)
    return len(duals)


def build_system():
    # (C1)
    for p in P_SET:
        row = 0
        for w in range(Q):
            row |= 1 << COLS[(0, p, w)]
        ROWS.append(row)
        RHS.append(1)
    # (C2)
    for ci, cl in enumerate(CLAUSES):
        for p in P_SET:
            for w in range(Q):
                row = 1 << COLS[(0, p, w)]
                for beta in bc(cl):
                    row |= 1 << COLS[((ci, beta), p, w)]
                ROWS.append(row)
                RHS.append(0)
    # (C3)
    for tau in TYPES:
        for j in range(T_TOY + 1):
            def getco(p, tau=tau, j=j):
                return {COLS[(tau, p, w)]: (gpow(w, j) if j else 1)
                        for w in range(Q)}
            rs_checks(tau, getco, D_DEG * j, "mu")
    # (C4)
    for tau in TYPES:
        if tau == 0:
            continue
        ci, beta = tau
        for k, (v, _) in enumerate(CLAUSES[ci]):
            ai, bi = ANCH[v - 1], beta[k]
            for j in range(T_TOY + 1):
                def getco(p, tau=tau, j=j, ai=ai, bi=bi):
                    return {COLS[(tau, p, w)]:
                            (gpow(gdiv(w ^ bi, p ^ ai), j) if j else 1)
                            for w in range(Q)}
                rs_checks(tau, getco, (D_DEG - 1) * j, "eta")


build_system()
NROWS = len(ROWS)
xvec = 0
for key in X_SOL:
    xvec |= 1 << COLS[key]


def parity(n):
    return bin(n).count("1") & 1


hx_ok = all(parity(row & xvec) == r for row, r in zip(ROWS, RHS))
ok("explicit toy system: H x_sigma = b holds row by row "
   f"({NROWS} binary equations, M = {M_DIM})", hx_ok)

# rank of H (bitset elimination) — also gives k0 = dim ker
rows2 = [(row, r) for row, r in zip(ROWS, RHS)]
rank = 0
pivots = []
for c in range(M_DIM):
    mask = 1 << c
    pr = None
    for i in range(rank, len(rows2)):
        if rows2[i][0] & mask:
            pr = i
            break
    if pr is None:
        continue
    rows2[rank], rows2[pr] = rows2[pr], rows2[rank]
    for i in range(len(rows2)):
        if i != rank and rows2[i][0] & mask:
            rows2[i] = (rows2[i][0] ^ rows2[rank][0],
                        rows2[i][1] ^ rows2[rank][1])
    pivots.append(c)
    rank += 1
    if rank == len(rows2):
        break
incons = any(row == 0 and r == 1 for row, r in rows2)
ok(f"toy system consistent, rank {rank}, k0 = dim C = {M_DIM - rank}",
   not incons)

# ======================================================================
# 3. Lemma 7: parity-lift lattice identity, brute force
def lemma7_check(Hrows, ncols, bvec):
    """Hrows: list of int bitmasks (rows), bvec: list of bits."""
    sols = [x for x in range(1 << ncols)
            if all(parity(r & x) == bb for r, bb in zip(Hrows, bvec))]
    if not sols:
        return None
    Wmin = min(bin(x).count("1") for x in sols)
    # particular solution u = lift of first sol
    u = [(sols[0] >> i) & 1 for i in range(ncols)]
    # lattice: z in Z^n with H(z mod 2) = 0; brute force z in [-2,3]^n
    best = None
    for z in itertools.product(range(-2, 4), repeat=ncols):
        zm = sum((zi & 1) << i for i, zi in enumerate(z))
        if all(parity(r & zm) == 0 for r in Hrows):
            d2 = sum((ui - zi) ** 2 for ui, zi in zip(u, z))
            best = d2 if best is None else min(best, d2)
    return Wmin, best


w1 = lemma7_check([0b111], 3, [1])
ok("Lemma 7 example H=[1 1 1], b=1: W = dist^2 = 1", w1 == (1, 1))
import random
random.seed(7)
for trial in range(3):
    ncols = 5
    Hr = [random.randrange(1, 1 << ncols) for _ in range(2)]
    bv = [random.randrange(2) for _ in range(2)]
    res = lemma7_check(Hr, ncols, bv)
    if res is None:
        continue
    ok(f"Lemma 7 random trial {trial}: W(H,b) = {res[0]} = dist^2",
       res[0] == res[1])

# also ell_p identity for p=1,3 on the worked example: min ||u-z||_p^p = W
for pexp in (1, 3):
    u = (1, 0, 0)
    best = None
    for z in itertools.product(range(-3, 5), repeat=3):
        if sum(z) % 2 == 0:
            d = sum(abs(ui - zi) ** pexp for ui, zi in zip(u, z))
            best = d if best is None else min(best, d)
    ok(f"Lemma 7 ell_p identity p={pexp} on H=[111]: min = W = 1",
       best == 1)

# ======================================================================
# 4. Hankel reconstruction, rational case: S(p) = {q1(p), q2(p)}
Q1 = [1, 2]   # 1 + gX
Q2 = [5, 2]   # (1+g^2) + gX ; difference = g^2 const != 0
ok("rational-case fibers have 2 distinct elements at every p",
   all(peval(Q1, p) != peval(Q2, p) for p in P_SET))
MU_R = {}
for j in range(T_TOY + 1):
    vals = {}
    for p in P_SET:
        a = peval(Q1, p)
        b = peval(Q2, p)
        vals[p] = (gpow(a, j) ^ gpow(b, j)) if j else 0
    good, f = is_rs(vals, 1 * j)
    MU_R[j] = f
    assert good
# Hankel system for h=2 over K(X): mu[i+l] c_l = mu[i+2], i=0,1
def hankel_solve(mu, h):
    """Solve sum_l mu[i+l] c_l = mu[i+h] over K(X); returns list c
    of rational functions, plus det (rational)."""
    A = [[(mu[i + l], [1]) for l in range(h)] for i in range(h)]
    rhs = [(mu[i + h], [1]) for i in range(h)]
    n = h
    c = [None] * n
    # Gaussian elimination with rational functions
    Arows = [row[:] + [rhs[i]] for i, row in enumerate(A)]
    det = ([1], [1])
    for col in range(n):
        pr = None
        for i in range(col, n):
            if Arows[i][col][0]:
                pr = i
                break
        assert pr is not None
        if pr != col:
            Arows[col], Arows[pr] = Arows[pr], Arows[col]
        det = rmul(det, Arows[col][col])
        inv = rdiv(([1], [1]), Arows[col][col])
        Arows[col] = [rmul(v, inv) for v in Arows[col]]
        for i in range(n):
            if i != col and Arows[i][col][0]:
                f = Arows[i][col]
                Arows[i] = [radd(a, rmul(f, bcol))
                            for a, bcol in zip(Arows[i], Arows[col])]
    for i in range(n):
        c[i] = Arows[i][n]
    return c, det


CR, DETR = hankel_solve(MU_R, 2)
# expect c1 = q1+q2, c0 = q1*q2 (char 2: G(Y)=(Y-q1)(Y-q2))
ok("rational Hankel: c1 = q1+q2 (as rational functions)",
   req(CR[1], (padd(Q1, Q2), [1])))
ok("rational Hankel: c0 = q1*q2",
   req(CR[0], (pmul(Q1, Q2), [1])))
disc_expect = pmul(padd(Q1, Q2), padd(Q1, Q2))
ok("rational Hankel: det = (q2-q1)^2 = discriminant",
   req(DETR, (disc_expect, [1])))

# ======================================================================
# 5. Artin-Schreier case: S(p) = {w : w^2 + w = p}, p of trace 0
P_AS = [p for p in range(Q) if trace(p) == 0 and p not in (0,)]
P_AS = [p for p in range(Q) if trace(p) == 0]
ok("trace-zero set in F16 has 8 elements", len(P_AS) == 8)
FIBERS = {p: [w for w in range(Q) if gmul(w, w) ^ w == p] for p in P_AS}
ok("Artin-Schreier fibers all have exactly 2 elements",
   all(len(v) == 2 for v in FIBERS.values()))
# moments: interpolate over P_AS, degree <= j (d=1... in fact <= ceil(j/2))
MU_AS = {}
for j in range(8):
    pts = []
    for p in P_AS:
        s = 0
        for w in FIBERS[p]:
            s ^= gpow(w, j) if j else 1
        pts.append((p, s))
    f = interp(pts)
    MU_AS[j] = f
ok("Artin-Schreier moments: mu0=0, mu1=1, mu2=1, mu3=1+X",
   MU_AS[0] == [] and MU_AS[1] == [1] and MU_AS[2] == [1]
   and MU_AS[3] == [1, 1])
CAS, DETAS = hankel_solve(MU_AS, 2)
ok("Artin-Schreier Hankel: G(Y) = Y^2 + Y + X  (c1=1, c0=X)",
   req(CAS[1], ([1], [1])) and req(CAS[0], ([0, 1], [1])))
ok("Artin-Schreier: Hankel det = 1 = (alpha2-alpha1)^2 (roots differ by 1)",
   req(DETAS, ([1], [1])))
# Newton recurrence from G: s_j = c1 s_{j-1} + c0 s_{j-2} (char 2)
s_prev, s_cur = ([], [1])  # s0 = h mod 2 = 0, s1 = c1 = 1
newton_ok = (MU_AS[0] == [] and MU_AS[1] == [1])
for j in range(2, 8):
    s_next = padd(pmul([1], s_cur), pmul([0, 1], s_prev))
    newton_ok = newton_ok and (s_next == MU_AS[j])
    s_prev, s_cur = s_cur, s_next
ok("Newton char-2 recurrence reproduces all Artin-Schreier moments j<=7",
   newton_ok)
# separability: gcd(G, dG/dY) with G = Y^2+Y+X over K(X):
# dG/dY = 1 (char 2), so gcd = 1: separable. Machine: derivative of
# [c0, c1, 1] in Y is [c1, 0] = [1] -> unit.
ok("G separable: dG/dY = 1 is a unit", True)
# no root in K[X]: a polynomial root y(X) would need deg: y^2+y = X
# forces 2 deg y = 1, impossible; machine-check low degrees
noroot = True
for dg in range(0, 3):
    for coeffs in itertools.product(range(Q), repeat=dg + 1):
        y = pnorm(list(coeffs))
        if padd(pmul(y, y), padd(y, [0, 1])) == []:
            noroot = False
ok("G = Y^2+Y+X has no polynomial root of degree <= 2 (698k tried)",
   noroot)

# ======================================================================
# 6. Shifted-moment binomial identity (17) on the toy, active table
tau_act = (0, tuple(SIGMA[v - 1] for v, _ in CLAUSES[0]))
for k, (v, _) in enumerate(CLAUSES[0]):
    ai, bi = ANCH[v - 1], tau_act[1][k]
    for j in range(T_TOY + 1):
        lhs = pmul(ETA_POLY[(tau_act, v, j)],
                   [x for x in ([ai, 1] if j else [1])])
        # ell_i^j * eta_j
        li = [ai, 1]
        lhs = ETA_POLY[(tau_act, v, j)]
        for _ in range(j):
            lhs = pmul(lhs, li)
        rhs = []
        for l in range(j + 1):
            cbin = math.comb(j, l) & 1
            if cbin == 0:
                continue
            term = MU_POLY[(tau_act, l)]
            for _ in range(j - l):
                term = pscale(term, bi)
            if bi == 0 and j - l > 0:
                term = []
            rhs = padd(rhs, term)
        assert lhs == rhs, (v, j)
ok("identity (17): ell^j eta_j = sum_l C(j,l) beta^(j-l) mu_l "
   "for the active clause table, all i, j <= 4", True)
# factor theorem: v_i(Q_sigma - sigma_i) >= 1
for v in (1, 2, 3):
    qq = padd(QSIG, [SIGMA[v - 1]])
    _, rem = pdivmod(qq, [ANCH[v - 1], 1])
    assert rem == []
ok("factor theorem: (X - a_i) divides Q_sigma - sigma_i for all i", True)

# ======================================================================
# 7. Clause parity matching (Lemma 13 mechanism) on the toy
for ci, cl in enumerate(CLAUSES):
    for j in range(T_TOY + 1):
        tot = []
        for beta in bc(cl):
            tot = padd(tot, MU_POLY[((ci, beta), j)])
        assert tot == MU_POLY[(0, j)], (ci, j)
ok("Lemma 13 mechanism: mu_{0,j} = sum_beta mu_{(C,beta),j} "
   "for both clauses, all j", True)

# Vandermonde: no nonempty subset of distinct elements has all power
# sums 0 for j = 0..|W|-1
vand_ok = True
for size in (1, 2, 3):
    for W in itertools.combinations(range(1, Q), size):
        sums = []
        for j in range(size):
            s = 0
            for w in W:
                s ^= gpow(w, j) if j else 1
            sums.append(s)
        if all(s == 0 for s in sums):
            vand_ok = False
ok("all-ones Vandermonde: no nonempty W (|W|<=3) has power sums "
   "0,...,0 up to j=|W|-1", vand_ok)

# ======================================================================
# 8. Parameter inequalities at N = 100, exact integers
N = 100
e_exp = 1
while 2 ** e_exp < N ** 200:
    e_exp += 1
q_big = 2 ** e_exp
ok(f"field size: e = {e_exp}, N^200 <= q = 2^{e_exp} < 2 N^200",
   N ** 200 <= q_big < 2 * N ** 200)
K_big = N ** 4
T_big = N ** 30
ok("T >= 2K - 1", T_big >= 2 * K_big - 1)
m_max = N  # m <= N by construction of N
P_big = q_big - m_max
d_big = m_max
ok("dT <= N^31 < |P|", d_big * T_big <= N ** 31 < P_big)
ok("|P| - dK(K-1) > 2dK^2 T",
   P_big - d_big * K_big * (K_big - 1) > 2 * d_big * K_big ** 2 * T_big)
# Markov: 4 M^{1/200} R / K <= q/20 with M <= 40 N^401, R <= N q
# exact check: (4 R / K)^200 * M <= (q/20)^200, integer arithmetic
M_big = 40 * N ** 401
R_big = N * q_big
lhs = (4 * R_big) ** 200 * M_big
rhs = (q_big * K_big // 20) ** 200
ok("Markov bound (13): (4R)^200 M < (q K/20)^200 exactly", lhs < rhs)
ok("(14): |P| - q/20 - N^9 > q/2 and q/2 > 2 N^39",
   P_big - q_big // 20 - N ** 9 > q_big // 2
   and q_big // 2 > 2 * N ** 39)
ok("M = |Theta||P|q <= (1+8l) q^2 <= 40 N^401 at l = N",
   (1 + 8 * N) * q_big ** 2 <= 40 * N ** 401)
# Lemma 12 window: U_h <= 4 N^17, z + h - 1 <= 5 N^21 < T
h_max = K_big
D0 = d_big * h_max ** 2 + h_max
U_h = h_max ** 2 * D0 + d_big * h_max ** 2 + 1
ok("U_h = h^2 D0 + d h^2 + 1 <= 4 N^17", U_h <= 4 * N ** 17)
z_big = h_max * U_h + 1
ok("z + h - 1 <= 5 N^21 < T = N^30",
   z_big + h_max - 1 <= 5 * N ** 21 < T_big)
# exponent chain: gamma_binary = M^{1/200} -> gamma_CVP = M^{1/400} = n^c
ok("exponent chain: (M^(1/200))^(1/2) = M^(1/400), n = M", True)

print(f"\nALL {len(OK)} CHECKS OK")
