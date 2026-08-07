"""Nineteen validation checks for Chapter 2 (binary & spherical codes).

Standard library only.  Every check prints PASS/FAIL and a number.
Exact checks use fractions.Fraction; numeric checks state tolerances.
Run: python ch02_checks.py   (~1-2 min)
"""
import math
import random
from fractions import Fraction as Fr

random.seed(2)
OK = []


def report(name, ok, detail=""):
    OK.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}")


def H2(u):
    if u <= 0 or u >= 1:
        return 0.0
    return -u*math.log2(u) - (1-u)*math.log2(1-u)


def g(v):
    return H2((1 - math.sqrt(max(0.0, 1-v)))/2)


def Hsph(u):
    if u <= 0:
        return 0.0
    return (1+u)*math.log2(1+u) - u*math.log2(u)


def qfun(u):
    return math.sqrt(u*(1+u))/(1+2*u) if u > 0 else 0.0


# ---- check 1: MRRW objective endpoints ---------------------------------
def check01():
    worst = 0.0
    for d in (0.05, 0.1, 0.2, 0.3, 0.4):
        f1 = 1 + g((1-2*d)**2) - g((1-2*d)**2 + 2*d*(1-2*d) + 2*d)
        m1 = H2(0.5 - math.sqrt(d*(1-d)))
        f0 = 1 + g(0.0) - g(2*d)
        worst = max(worst, abs(f1-m1), abs(f0 - (1-g(2*d))))
    report("01 F_delta endpoints = M1 and 1-g(2delta)", worst < 1e-12,
           f"max err {worst:.2e}")


# ---- check 2: cube ladder algebra, exact -------------------------------
def subsets(n, k):
    out = []

    def rec(start, cur):
        if len(cur) == k:
            out.append(tuple(cur))
            return
        for r in range(start, n):
            cur.append(r)
            rec(r+1, cur)
            cur.pop()
    rec(0, [])
    return out


def rank_exact(mat):
    m = [row[:] for row in mat]
    rows = len(m)
    cols = len(m[0]) if m else 0
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        inv = Fr(1, 1)/m[r][c]
        m[r] = [v*inv for v in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f*b for a, b in zip(m[i], m[r])]
        r += 1
    return r


def check02():
    n = 8
    ok = True
    for k in (1, 2, 3):
        lo = subsets(n, k-1)
        hi = subsets(n, k)
        up = subsets(n, k+1)
        D = [[Fr(1) if set(a) < set(b) else Fr(0) for b in hi]
             for a in lo]
        rk = rank_exact(D)
        ok &= (len(hi) - rk == math.comb(n, k) - math.comb(n, k-1))
        for _ in range(3):
            v = [Fr(random.randint(-5, 5)) for _ in hi]
            Uv = [sum(v[j] for j, b in enumerate(hi) if set(b) < set(a))
                  for a in up]
            DUv = [sum(Uv[i] for i, a in enumerate(up) if set(b) < set(a))
                   for b in hi]
            Dv = [sum(v[j] for j, b in enumerate(hi) if set(a) < set(b))
                  for a in lo]
            UDv = [sum(Dv[i] for i, a in enumerate(lo) if set(a) < set(b))
                   for b in hi]
            ok &= all(DUv[i] - UDv[i] == Fr(n-2*k)*v[i]
                      for i in range(len(hi)))
    report("02 ladder: dim E_k exact + DU-UD=(n-2k)I  (n=8)", ok)


# ---- check 3: k=0 stochasticity + balance (19), exact ------------------
def check03():
    n = 30
    ok = True
    for k in (0, 1, 2, 3):
        for i in range(k, n-k):
            Q = Fr((i-k+1)*(n-i-k))
            p_up = Q/(n*(i+1))
            p_dn = Q/(n*(n-i))
            ok &= (math.comb(n, i)*p_up == math.comb(n, i+1)*p_dn)
    for i in range(1, 30):
        n = 30
        up = Fr((i+1)*(n-i), n*(i+1))
        dn = Fr(i*(n-i+1), n*(n-i+1))
        ok &= (up + dn == 1)
    report("03 cube weights: (19) exact; k=0 rows sum to 1", ok)


# ---- check 4: Lemma 2.2 convergence ------------------------------------
def lam_tri(diag, off):
    """Largest eigenvalue of a symmetric tridiagonal matrix, by
    Sturm-count bisection (robust: no bipartite oscillation)."""
    m = len(diag)

    def count_below(x):
        cnt, d = 0, 1.0
        for i in range(m):
            d = diag[i] - x - (off[i-1]**2/d if i > 0 else 0.0)
            if d == 0.0:
                d = -1e-300
            if d < 0:
                cnt += 1
        return cnt
    lo = min(diag) - 2*max((abs(x) for x in off), default=0.0) - 1.0
    hi = max(diag) + 2*max((abs(x) for x in off), default=0.0) + 1.0
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if count_below(mid) == m:
            hi = mid
        else:
            lo = mid
    return 0.5*(lo+hi)


def check04():
    a, b = 0.35, 0.10
    gam = 2*(a-b)*(1-a-b)/math.sqrt(a*(1-a))
    gaps = []
    for n in (500, 2000, 4000):
        k, L = int(b*n), int(a*n)
        off = [(i-k+1)*(n-i-k)/(n*math.sqrt((i+1)*(n-i)))
               for i in range(k, L)]
        lam = lam_tri([0.0]*(L-k+1), off)
        gaps.append(gam - lam)
    ok = gaps[0] > gaps[1] > gaps[2] > 0 and gaps[2] < 6e-3
    report("04 Lemma 2.2: lam_max -> Gamma_H from below",
           ok, "gaps " + " ".join("%.1e" % x for x in gaps))


# ---- check 5: Toy 5.1 numbers + Perron identity (22) -------------------
def cube_lam_perron(n, k, L):
    lvls = list(range(k, L+1))
    off = [(i-k+1)*(n-i-k)/(n*math.sqrt((i+1)*(n-i)))
           for i in range(k, L)]
    lam = lam_tri([0.0]*len(lvls), off)
    m = len(lvls)
    v = [1.0]*m
    for _ in range(4000):
        w = [v[i] + (off[i-1]*v[i-1] if i > 0 else 0.0)
             + (off[i]*v[i+1] if i < m-1 else 0.0) for i in range(m)]
        nrm = math.sqrt(sum(x*x for x in w))
        v = [x/nrm for x in w]
    wgt = [math.sqrt(math.comb(n, lvls[i]))*v[i] for i in range(m)]
    p = {}
    for i in range(k, L):
        Q = (i-k+1)*(n-i-k)
        p[(i, i+1)] = Q/(n*(i+1))
        p[(i+1, i)] = Q/(n*(n-i))
    err = 0.0
    for j in range(m):
        s = 0.0
        if j > 0:
            s += p[(lvls[j-1], lvls[j])]*wgt[j-1]
        if j < m-1:
            s += p[(lvls[j+1], lvls[j])]*wgt[j+1]
        err = max(err, abs(s - lam*wgt[j]))
    return lam, err


def check05():
    lam10, err22 = cube_lam_perron(10, 1, 4)
    s = 0.2
    dk = math.comb(10, 1) - 1
    Damb = sum(math.comb(10, i) for i in range(1, 5))
    bound = (1-s)/(dk*(lam10-s))*Damb
    ok = (abs(lam10 - 0.4622578) < 1e-6 and abs(bound - 130.49) < 0.01
          and err22 < 1e-9)
    report("05 Toy 5.1: lam=0.46226, bound 130.49; Perron id (22)",
           ok, f"lam={lam10:.7f} bound={bound:.2f} err={err22:.1e}")


# ---- checks 6-8: Hahn coefficients -------------------------------------
def hahn(n, w, p, q, j):
    N = n - w
    j1 = w/2 - p
    j2 = N/2 - q
    jj = n/2 - j
    m0 = n/2 - w
    Sg = j1 + j2
    Dj = j2 - j1
    mu = m0/2*(j2*(j2+1) - j1*(j1+1))/(jj*(jj+1)) if jj > 0 else 0.0
    rad = (jj*jj - m0*m0)*(jj*jj - Dj*Dj)*((Sg+1)**2 - jj*jj)
    nu = (math.sqrt(max(rad, 0.0))/(2*jj*math.sqrt((2*jj-1)*(2*jj+1)))
          if jj > 0.5 else 0.0)
    return (n*mu - m0*m0)/(w*N), n*nu/(w*N)


def check06():
    n, w = 25, 10
    DJ = [math.comb(n, j) - (math.comb(n, j-1) if j else 0)
          for j in range(w+2)]
    worst = 0.0
    for j in range(0, w+1):
        b0, _ = hahn(n, w, 0, 0, j)
        tot = b0 if j > 0 else 0.0
        if j < w:
            _, c = hahn(n, w, 0, 0, j)
            tot += c*math.sqrt(DJ[j+1]/DJ[j])
        if j > 0:
            _, c = hahn(n, w, 0, 0, j-1)
            tot += c*math.sqrt(DJ[j-1]/DJ[j])
        worst = max(worst, abs(tot - 1.0))
    report("06 Johnson p=q=0: directed weights sum to 1 on full path",
           worst < 1e-10, f"max err {worst:.1e}")


def check07():
    n, w, p, q = 40, 15, 3, 5
    N = n - w
    DJ = [math.comb(n, j) - (math.comb(n, j-1) if j else 0)
          for j in range(w+2)]
    jm, jp = p+q, min(w, w-p+q, N+p-q)
    worst = 0.0
    for j in range(jm, jp):
        bj, cj = hahn(n, w, p, q, j)
        b0, c0 = hahn(n, w, 0, 0, j)
        pij = cj*cj/c0*math.sqrt(DJ[j+1]/DJ[j])
        pji = cj*cj/c0*math.sqrt(DJ[j]/DJ[j+1])
        worst = max(worst, abs(DJ[j]*pij - DJ[j+1]*pji))
        if cj <= 0:
            worst = 1.0
    report("07 Johnson balance (39) + c^pq>0 inside window",
           worst < 1e-6, f"max err {worst:.1e}")


def check08():
    al, be, ga, u = 0.3, 0.06, 0.1, 0.2
    n = 4000
    w, p, q, L = int(al*n), int(be*n), int(ga*n), int(u*n)
    bj, cj = hahn(n, w, p, q, L)
    b0, c0 = hahn(n, w, 0, 0, L)
    z, m = 1-2*u, 1-2*al
    ze, xi = 1-2*be-2*ga, 1-2*al+2*be-2*ga
    B = (ze*xi - m*z*z)**2/(z*z*(1-m*m)*(1-z*z))
    C = (z*z-xi*xi)*(ze*ze-z*z)/(2*z*z*(1-m*m)*math.sqrt(1-z*z))
    e1 = abs(bj*bj/b0 - B)
    e2 = abs(cj*cj/c0 - C)
    worst48 = 0.0
    for _ in range(20):
        a2 = random.uniform(0.05, 0.45)
        u2 = random.uniform(0.01, a2-0.005)
        A2, U2 = a2*(1-a2), u2*(1-u2)
        z2, m2 = 1-2*u2, 1-2*a2
        B2 = (m2 - m2*z2*z2)**2/(z2*z2*(1-m2*m2)*(1-z2*z2))
        C2 = ((z2*z2-m2*m2)*(1-z2*z2)
              / (2*z2*z2*(1-m2*m2)*math.sqrt(1-z2*z2)))
        lam = B2 + 2*C2
        rhs = (A2-U2)/(A2*(1+2*math.sqrt(U2)))
        worst48 = max(worst48, abs(1-lam - rhs))
    ok = e1 < 2e-3 and e2 < 2e-3 and worst48 < 1e-12
    report("08 Lemma 3.5 limits at n=4000; identity (48)",
           ok, f"errB={e1:.1e} errC={e2:.1e} err48={worst48:.1e}")


# ---- checks 9-11: spherical dimensions and recurrence ------------------
def DS(m, n):
    if m == 0:
        return 1
    if m == 1:
        return n
    return math.comb(m+n-1, m) - math.comb(m+n-3, m-2)


def check09():
    ok = True
    for n in range(3, 61):
        for m in range(0, 201):
            rhs = (2*m+n-2)*math.comb(m+n-3, m)//(n-2)
            ok &= (DS(m, n) == rhs)
    report("09 (59): two dimension formulas agree exactly", ok)


def check10():
    ok = True
    for n in range(4, 40):
        for i in range(0, 30):
            ok &= (DS(i, n) == sum(DS(k, n-1) for k in range(i+1)))
    report("10 branching (62): D_i(n) = sum_k dim H_k(n-1) exactly", ok)


def check11():
    ok = True
    for n in range(4, 30):
        for i in range(0, 60):
            a0sq = Fr((i+1)*(i+n-2), (2*i+n-2)*(2*i+n))
            mu = Fr(i+1, 2*i+n)
            md = Fr(i+n-2, 2*i+n-2)
            ok &= (mu*md == a0sq)
            ok &= (mu*mu*Fr(DS(i+1, n), DS(i, n)) == a0sq)
    n, k = 6, 2
    eta = k + (n-2)/2

    def geg(j, t):
        c0, c1 = 1.0, 2*eta*t
        if j == 0:
            return c0
        for jj in range(1, j):
            c2 = (2*(jj+eta)*t*c1 - (jj+2*eta-1)*c0)/(jj+1)
            c0, c1 = c1, c2
        return c1
    M = 20000
    ts = [-1 + 2*(i+0.5)/M for i in range(M)]
    wt = [(1-t*t)**(eta-0.5)*2/M for t in ts]

    def ip(f1, f2):
        return sum(w*f1(t)*f2(t) for t, w in zip(ts, wt))
    worst = 0.0
    for i in range(k, k+4):
        j = i - k
        n1 = math.sqrt(ip(lambda t: geg(j, t), lambda t: geg(j, t)))
        n2 = math.sqrt(ip(lambda t: geg(j+1, t),
                          lambda t: geg(j+1, t)))
        alpha = ip(lambda t: t*geg(j, t),
                   lambda t: geg(j+1, t))/(n1*n2)
        pred = math.sqrt((i-k+1)*(i+k+n-2)/((2*i+n-2)*(2*i+n)))
        worst = max(worst, abs(alpha - pred))
    report("11 (66)/(67): exact rational ids; quadrature alpha^(k)",
           ok and worst < 1e-6, f"quad err {worst:.1e}")


# ---- check 12: Kravchuk weights (74) + Weyl (111), exact ---------------
def weyl(n, lam):
    r1 = len(lam)
    out = Fr(1)
    for i0 in range(r1):
        i = i0 + 1
        li = lam[i0]
        out *= Fr(2*li + n - 2*i, n - 2*i)
        num, den = Fr(1), Fr(1)
        for t in range(li):
            num *= (n - i - r1 - 1 + li - t)
            den *= (r1 - i + li - t)
        out *= num/den
    for i0 in range(r1):
        for j0 in range(i0+1, r1):
            i, j = i0+1, j0+1
            out *= Fr((lam[i0]-lam[j0]+j-i)*(lam[i0]+lam[j0]+n-i-j),
                      (j-i)*(n-i-j))
    return out


def pkrav(n, r, lam, mu, li, sgn):
    hl = [Fr(2*lam[i] + n - 2*(i+1), 2) for i in range(r+1)]
    hm = [Fr(2*mu[m] + n - 1 - 2*(m+1), 2) for m in range(r)]
    rho = Fr(n - 2*r - 2, 2)
    e = Fr(sgn, 2)
    num = (hl[li] + rho*sgn)
    for m in range(r):
        num *= (hl[li] + e)**2 - hm[m]**2
    den = 2*hl[li]
    for qq in range(r+1):
        if qq != li:
            den *= hl[li]**2 - hl[qq]**2
    return num/den


def interlace_ok(lam, mu):
    for i in range(len(mu)):
        if not (lam[i] >= mu[i] >= lam[i+1]):
            return False
    return all(lam[i] >= lam[i+1] for i in range(len(lam)-1))


def check12():
    ok = True
    tested = 0
    while tested < 40:
        r = random.choice([1, 2, 3])
        n = random.choice([2*r+4, 2*r+5, 2*r+8, 2*r+11])
        lam = sorted([random.randint(0, 9) for _ in range(r+1)],
                     reverse=True)
        mu = [random.randint(lam[i+1], lam[i]) for i in range(r)]
        if not interlace_ok(lam, mu):
            continue
        tested += 1
        tot = Fr(0)
        for li in range(r+1):
            for sgn in (1, -1):
                nl = lam[:]
                nl[li] += sgn
                valid = (all(nl[i] >= nl[i+1] for i in range(r))
                         and nl[-1] >= 0 and interlace_ok(nl, mu))
                if valid:
                    pv = pkrav(n, r, lam, mu, li, sgn)
                    tot += pv
                    if sgn == 1:
                        rev = pkrav(n, r, nl, mu, li, -1)
                        ok &= (weyl(n, lam)*pv == weyl(n, nl)*rev)
        ok &= (tot == 1)
    for n in range(8, 20):
        for m in range(0, 8):
            ok &= (weyl(n, [m]) == DS(m, n))
    n, k = 12, 4
    for i in range(k, k+4):
        for j in range(0, k+1):
            pip = pkrav(n, 1, [i, j], [k], 0, 1)
            expl = Fr((i-k+1)*(i+k+n-2)*(i+n-3),
                      (i-j+1)*(i+j+n-3)*(2*i+n-2))
            ok &= (pip == expl)
    report("12 (74)+(111): (75) sums to 1, (76) exact, (79) matches",
           ok)


# ---- check 13: level-1 degenerations + Lemma 7.1 -----------------------
def gamma_phi(av, bv):
    x = [a*(1+a) for a in av]
    y = [b*(1+b) for b in bv]
    r = len(bv)
    G = 0.0
    for li in range(r+1):
        num = 1.0
        for ym in y:
            num *= (x[li]-ym)
        den = 1.0
        for j in range(r+1):
            if j != li:
                den *= (x[li]-x[j])
        G += num/den*qfun(av[li])
    Phi = sum(Hsph(a) for a in av) - sum(Hsph(b) for b in bv)
    return G, Phi


def check13():
    worst = 0.0
    for _ in range(20):
        a = random.uniform(0.1, 2.0)
        b = random.uniform(0.01, a*0.9)
        G, _ = gamma_phi([a, 1e-12], [b])
        grow = (a-b)*(1+a+b)/((1+2*a)*math.sqrt(a*(1+a)))
        worst = max(worst, abs(G-grow))
    n = 100000
    av, bv = [0.5, 0.2, 0.05], [0.3, 0.1]
    lam = [int(a*n) for a in av]
    mu = [int(b*n) for b in bv]
    x = [a*(1+a) for a in av]
    y = [b*(1+b) for b in bv]
    worst2 = 0.0
    for li in range(3):
        num = 1.0
        for ym in y:
            num *= (x[li]-ym)
        den = 1.0
        for j in range(3):
            if j != li:
                den *= (x[li]-x[j])
        R = num/den
        for sgn in (1, -1):
            pv = float(pkrav(n, 2, lam, mu, li, sgn))
            pred = R*(2*av[li]+1+sgn)/(2*(1+2*av[li]))
            worst2 = max(worst2, abs(pv-pred))
    report("13 Gamma_1((a,0),b)=Gamma_row; Lemma 7.1 at n=1e5",
           worst < 1e-5 and worst2 < 1e-4,
           f"row err {worst:.1e}, 7.1 err {worst2:.1e}")


# ---- check 14: BKL(1/2) ------------------------------------------------
def check14():
    s = 0.5

    def a0(t):
        return 0.5*((1-t*t)**-0.5 - 1)

    def obj(t):
        return Hsph(a0(t)) + 0.5*math.log2((1-t)/(1-s))
    lo, hi = 1e-6, s
    for _ in range(200):
        m1 = lo + (hi-lo)/3
        m2 = hi - (hi-lo)/3
        if obj(m1) < obj(m2):
            hi = m2
        else:
            lo = m1
    val = obj((lo+hi)/2)
    ok = abs(val - 0.4009442) < 1e-6
    worst = max(abs(2*qfun(a0(t)) - t) for t in (0.2, 0.5, 0.8))
    report("14 BKL(1/2)=0.4009442; 2q(a0(s))=s", ok and worst < 1e-12,
           f"BKL={val:.7f}")


# ---- check 15: kissing certificates ------------------------------------
def scale_tuple(av, bv, c):
    def sc(u):
        x = c*u*(1+u)
        return (-1 + math.sqrt(1+4*x))/2
    return [sc(a) for a in av], [sc(b) for b in bv]


def check15():
    s = 0.5
    a1, b1 = scale_tuple([0.0977018, 0.0030061], [0.0141955],
                         1.00002)
    G1, P1 = gamma_phi(a1, b1)
    a2, b2 = scale_tuple(
        [9.05293507e-2, 5.65430077e-4, 2.46387119e-6],
        [6.93043943e-3, 4.36084190e-5], 1.00002)
    G2, P2 = gamma_phi(a2, b2)
    ok = (2*G1 > s and P1 < 0.400944 and 2*G2 > s and P2 < 0.39661)
    report("15 kissing: level-1 cert 0.3968, level-2 cert 0.3966",
           ok, f"2G2={2*G2:.7f} Phi2={P2:.5f}")


# ---- check 16: the two half-line integrals (15.a) ----------------------
def simpson(f, a, b, m=4000):
    h = (b-a)/m
    s = f(a) + f(b)
    for i in range(1, m):
        s += f(a+i*h)*(4 if i % 2 else 2)
    return s*h/3


def check16():
    c = 0.25
    worst = 0.0
    for x in (0.03, 0.4, 2.0, 11.0):
        i1 = simpson(lambda th: (2/math.pi)*x/(x+c*math.sin(th)**2),
                     0, math.pi/2)
        worst = max(worst, abs(i1 - math.sqrt(x/(x+c))))
        a = (-1 + math.sqrt(1+4*x))/2
        i2 = simpson(lambda v: 1/(x + c*(1-v*v)), 0, 1)
        pred = 2/(1+2*a)*math.log((1+a)/a)
        worst = max(worst, abs(i2 - pred))
    report("16 (15.a): arcsine and nu integrals vs closed forms",
           worst < 1e-9, f"max err {worst:.1e}")


# ---- check 17: Chebyshev construction (104)-(107) ----------------------
def check17():
    R = 6.0
    ZR = 2/math.pi*math.asin(1/math.sqrt(4*R+1))
    zs = []
    for N in (8, 40, 200):
        xs = [R/2*(1+math.cos((2*t-1)*math.pi/(2*N)))
              for t in range(1, N+1)]
        zs.append(1 - sum(math.sqrt(x/(x+0.25)) for x in xs)/N)
    ok = (abs(zs[-1]-ZR) < 1e-3 and abs(zs[0]-ZR) > abs(zs[-1]-ZR))
    lamstar = 0.5*math.log2(2*math.pi/math.e)
    vals = []
    for RR in (0.5, 2.0, 10.0):
        ZRR = 2/math.pi*math.asin(1/math.sqrt(4*RR+1))
        aR = (-1+math.sqrt(1+4*RR))/2
        hR = math.log(2)*Hsph(aR)
        vals.append((math.log(2/ZRR) - hR)/(2*math.log(2)))
    ok &= all(vals[i] < vals[i+1] < lamstar for i in range(2))
    ok &= lamstar - vals[-1] < 1e-4
    report("17 Chebyshev: Z_N -> Z_R; objective -> lambda*",
           ok, f"obj(R=10)={vals[-1]:.5f} vs {lamstar:.5f}")


# ---- check 18: Prop 8.1 entropy identity + Wallis ----------------------
def check18():
    av, bv = [0.6, 0.13], [0.31]
    G, P = gamma_phi(av, bv)
    Z = 1 - 2*G
    xs = [a*(1+a) for a in av]
    ys = [b*(1+b) for b in bv]
    c = 0.25

    def K(t):
        val = math.log((t+xs[-1])/t)
        for xm, ym in zip(xs[:-1], ys):
            val += math.log((t+xm)/(t+ym))
        return val

    def integrand(v):
        t = c*(1-v*v)
        return math.log(math.pi*math.sqrt(t)) + K(t) + math.log(Z)
    dkl = simpson(integrand, 1e-9, 1-1e-9, 20000)
    lamstar = 0.5*math.log2(2*math.pi/math.e)
    lhs = lamstar - (0.5*math.log2(2/Z) - P)
    ok = abs(lhs - dkl/(2*math.log(2))) < 2e-4

    def wallis_f(a):
        if a < 1e-9:
            return 1.0
        return math.exp(-2*a)*math.tanh(a)/a
    wall = simpson(wallis_f, 0.0, 40.0, 40000)
    ok &= abs(wall - math.log(math.pi/2)) < 1e-6
    report("18 Prop 8.1 identity (99) numerically; Wallis integral",
           ok, f"gap={lhs:.6f} DKL/2ln2={dkl/(2*math.log(2)):.6f}")


# ---- check 19: the cross-chapter constant ------------------------------
def check19():
    lamstar = 0.5*math.log2(2*math.pi/math.e)
    ok = abs(2**(-lamstar) - math.sqrt(math.e/(2*math.pi))) < 1e-15
    ok &= abs(2**(-lamstar) - 0.65774462) < 1e-8
    report("19 2^(-lambda*) = sqrt(e/2pi) = 0.65774462 (= Chapter 1)",
           ok, f"{2**(-lamstar):.8f}")


if __name__ == "__main__":
    for f in (check01, check02, check03, check04, check05, check06,
              check07, check08, check09, check10, check11, check12,
              check13, check14, check15, check16, check17, check18,
              check19):
        f()
    print(f"{sum(OK)}/{len(OK)} checks passed")
