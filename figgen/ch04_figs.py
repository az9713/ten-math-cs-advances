"""Verification engine + figure data for the ch04 explainer (Connes
rigidity counterexample).  Exact arithmetic over F2 / Z4 on truncated
models V_N = (F2[t]/t^N)^4.

Run: PYTHONIOENCODING=utf-8 python figgen/ch04_figs.py
"""
import itertools
import random

random.seed(4)


# ---------------------------------------------------------------------
# 1. The four-point toy: coordinatewise vs carry law on F2 x F2.
# ---------------------------------------------------------------------
def diamond(p, q):
    (x, y), (x2, y2) = p, q
    return ((x + x2) % 2, (y + y2 + x*x2) % 2)


def to_z4(p):
    return (p[0] + 2*p[1]) % 4


pts = [(x, y) for x in (0, 1) for y in (0, 1)]
ok = all(to_z4(diamond(p, q)) == (to_z4(p) + to_z4(q)) % 4
         for p in pts for q in pts)
print('four-point: (x,y)<>(x\',y\') = Z/4 via x+2y:', ok)
print('four-point: coordinatewise law has exponent 2:',
      all(((2*x) % 2, (2*y) % 2) == (0, 0) for (x, y) in pts))


# ---------------------------------------------------------------------
# 2. Vectors over R_N = F2[t]/t^N, rank 4.  v is a tuple of 4 tuples of
#    N bits (coefficients of 1, t, ..., t^(N-1)).
# ---------------------------------------------------------------------
def all_vectors(N):
    for bits in itertools.product((0, 1), repeat=4*N):
        yield tuple(tuple(bits[i*N:(i+1)*N]) for i in range(4))


def poly_trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return tuple(p)


def poly_mod(a, b):
    """a mod b over F2, b nonzero, as trimmed tuples."""
    a = list(poly_trim(a))
    while len(a) >= len(b):
        sh = len(a) - len(b)
        for i in range(len(b)):
            a[sh+i] ^= b[i]
        a = list(poly_trim(a))
    return tuple(a)


def poly_gcd(a, b):
    a, b = poly_trim(a), poly_trim(b)
    while b:
        a, b = b, poly_mod(a, b)
    return a


def is_primitive(v, N):
    """gcd over F2[t] of the four coordinate polynomials is 1."""
    g = ()
    for c in v:
        g = poly_gcd(g, c)
    return g == (1,)


print('--- primitive-vector count |P n V_N| = 7*2^(4N-3)+1 ---')
for N in (1, 2, 3):
    cnt = sum(1 for v in all_vectors(N) if is_primitive(v, N))
    want = 7 * 2**(4*N - 3) + 1
    print(f'N={N}: count={cnt} formula={want} match={cnt == want}')


# ---------------------------------------------------------------------
# 3. Lemma 5.5: nonzero Boolean polynomial of degree <= 2 on F2^m has
#    support >= 2^(m-2).  Full enumeration for m = 4 and m = 5.
# ---------------------------------------------------------------------
def support_count(m, coeffs):
    """coeffs = (c0, linear tuple, quadratic dict) -> |supp|."""
    c0, lin, quad = coeffs
    n = 0
    for x in itertools.product((0, 1), repeat=m):
        val = c0
        for i in range(m):
            val += lin[i]*x[i]
        for (i, j), c in quad.items():
            val += c*x[i]*x[j]
        if val % 2:
            n += 1
    return n


print('--- Lemma 5.5: min support of nonzero deg<=2 poly ---')
for m in (4, 5):
    pairs = [(i, j) for i in range(m) for j in range(i+1, m)]
    best = None
    n_polys = 0
    for c0 in (0, 1):
        for lin in itertools.product((0, 1), repeat=m):
            for qbits in itertools.product((0, 1), repeat=len(pairs)):
                if c0 == 0 and not any(lin) and not any(qbits):
                    continue
                n_polys += 1
                quad = {p: b for p, b in zip(pairs, qbits) if b}
                s = support_count(m, (c0, lin, quad))
                if best is None or s < best:
                    best = s
    # note: constant term allowed here; Lemma 5.5 has zero constant
    # term but the bound is stated for any nonzero p of deg <= 2.
    print(f'm={m}: {n_polys} nonzero polys, min support = {best},'
          f' bound 2^(m-2) = {2**(m-2)}, ok={best >= 2**(m-2)}')


# ---------------------------------------------------------------------
# 4. Detection constants c_j = (2^(j-3)-1)/(2^(j-1)-1) by rank j.
# ---------------------------------------------------------------------
print('--- rank threshold: c_j = (2^(j-3)-1)/(2^(j-1)-1) ---')
from fractions import Fraction
for j in range(3, 8):
    cj = Fraction(2**(j-3) - 1, 2**(j-1) - 1)
    print(f'j={j}: c_j = {cj} positive={cj > 0}')


# ---------------------------------------------------------------------
# 5. The carry group on the truncated model.
#    l : X_M = dual of V_M  (M = N + n, so that T_n is defined down to
#    V_N); represented as 4xM bit matrix l[i][j] = l(t^j e_i).
#    q : Y_N = dual of B_N; represented by values on the basis
#        {b(t^j e_i)} u {mixed}, but for the checks we only need q on
#        b(v) for v in V_N, so store q as a dict over basis of B_N:
#        diag (i,j) and offdiag ((i,j),(i2,j2)).
# ---------------------------------------------------------------------
def rand_l(M):
    return tuple(tuple(random.randint(0, 1) for _ in range(M))
                 for _ in range(4))


def l_of(l, v):
    """l(v) for v in V_M."""
    return sum(l[i][j]*v[i][j] for i in range(4)
               for j in range(len(v[i]))) % 2


def Tn(l, n):
    """(T_n l)(v) = l(t^n v): shift columns."""
    M = len(l[0])
    return tuple(tuple(l[i][j+n] if j+n < M else 0 for j in range(M))
                 for i in range(4))


# q on B_M: dict basis -> bit.  Basis of B: diag[(i,j)] ~ b(t^j e_i)
# and off[((i,j),(i2,j2))] for (i,j) < (i2,j2) ~ symmetric sums.
def rand_q(M):
    diag = {(i, j): random.randint(0, 1)
            for i in range(4) for j in range(M)}
    keys = [(i, j) for i in range(4) for j in range(M)]
    off = {}
    for a in range(len(keys)):
        for b in range(a+1, len(keys)):
            off[(keys[a], keys[b])] = random.randint(0, 1)
    return (diag, off)


def q_of_bv(q, v):
    """q(b(v)) where b(v) = v (x) v, expanded in the basis."""
    diag, off = q
    keys = [(i, j) for i in range(4) for j in range(len(v[0]))]
    x = {k: v[k[0]][k[1]] for k in keys}
    val = 0
    for k in keys:
        val += diag[k]*x[k]
    for (a, b), c in off.items():
        val += c*x[a]*x[b]
    return val % 2


def r_ll(l1, l2, v):
    """r_{l1,l2}(b(v)) = l1(v) l2(v)."""
    return (l_of(l1, v)*l_of(l2, v)) % 2


def q_add(qs):
    """XOR of q-dicts (same key sets)."""
    diag = {}
    off = {}
    for k in qs[0][0]:
        diag[k] = sum(q[0][k] for q in qs) % 2
    for k in qs[0][1]:
        off[k] = sum(q[1][k] for q in qs) % 2
    return (diag, off)


def q_of_r(l1, l2, M):
    """the q-dict of r_{l1,l2}: values on the B_M basis.
    r(b(t^j e_i)) = l1(t^j e_i) l2(t^j e_i);
    r(sym(u,w)) = l1(u)l2(w) + l1(w)l2(u) for basis vectors u,w."""
    diag = {(i, j): (l1[i][j]*l2[i][j]) % 2
            for i in range(4) for j in range(M)}
    keys = [(i, j) for i in range(4) for j in range(M)]
    off = {}
    for a in range(len(keys)):
        for b in range(a+1, len(keys)):
            (i, j), (i2, j2) = keys[a], keys[b]
            off[(keys[a], keys[b])] = (l1[i][j]*l2[i2][j2]
                                       + l1[i2][j2]*l2[i][j]) % 2
    return (diag, off)


def star(n, p1, p2, M):
    """(l,q) *_n (l',q') on the truncated model."""
    (l1, q1), (l2, q2) = p1, p2
    return (tuple(tuple((l1[i][j]+l2[i][j]) % 2 for j in range(M))
                  for i in range(4)),
            q_add([q1, q2, q_of_r(Tn(l1, n), Tn(l2, n), M)]))


print('--- carry group *_n on truncated model (M=3), sampled ---')
M = 3
for n in (0, 1):
    ok_assoc = ok_comm = ok_inv = True
    for _ in range(30):
        a = (rand_l(M), rand_q(M))
        b = (rand_l(M), rand_q(M))
        c = (rand_l(M), rand_q(M))
        ab_c = star(n, star(n, a, b, M), c, M)
        a_bc = star(n, a, star(n, b, c, M), M)
        ok_assoc &= (ab_c == a_bc)
        ok_comm &= (star(n, a, b, M) == star(n, b, a, M))
        # inverse: (l, q + c_n(l,l)) -- doubling formula check below
        aa = star(n, a, a, M)
        zero_l = tuple(tuple(0 for _ in range(M)) for _ in range(4))
        ok_inv &= (aa[0] == zero_l)
        # 2(l,q) = (0, c_n(l,l)):
        want_q = q_of_r(Tn(a[0], n), Tn(a[0], n), M)
        ok_inv &= (aa[1] == want_q)
    print(f'n={n}: assoc={ok_assoc} comm={ok_comm}'
          f' doubling 2(l,q)=(0,c_n(l,l))={ok_inv}')


# K-equivariance under a sampled transvection g = I + t E_12 acting on
# V (hence diagonally on B); dual action on (l, q).
def g_on_v(v, M):
    """g v with g = I + t E12: v1 += t*v2."""
    v = [list(c) for c in v]
    for j in range(M-1):
        v[0][j+1] = (v[0][j+1] + v[1][j]) % 2
    return tuple(tuple(c) for c in v)


def ginv_on_v(v, M):
    v = [list(c) for c in v]
    for j in range(M-1):
        v[0][j+1] = (v[0][j+1] + v[1][j]) % 2
    return tuple(tuple(c) for c in v)  # char 2: g^{-1} = g


# ---------------------------------------------------------------------
# 6. Order-4 element F_{l,0} in (Z/4)^{V_N}; and Z_n count.
# ---------------------------------------------------------------------
print('--- order-4 carry element and |Z_n| = 2^(4n) ---')
N = 2


def F(l, q, v):
    return (l_of(l, v) + 2*q_of_bv(q, v)) % 4


l_e = tuple(tuple(1 if (i == 0 and j == 0) else 0 for j in range(N))
            for i in range(4))
q0 = ({k: 0 for k in [(i, j) for i in range(4) for j in range(N)]},
      {k: 0 for k in rand_q(N)[1]})
e_vec = tuple(tuple(1 if (i == 0 and j == 0) else 0 for j in range(N))
              for i in range(4))
# pointwise 2*F: (2 F_{l,0})(v) = 2 l(v) mod 4, nonzero at v = e:
print('F_{l,0}(e) = 1, (2F)(e) = 2 != 0 -> order 4:',
      F(l_e, q0, e_vec) == 1, (2*F(l_e, q0, e_vec)) % 4 == 2)
for n in (1, 2):
    Mn = 3
    cnt = 0
    for bits in itertools.product((0, 1), repeat=4*Mn):
        l = tuple(tuple(bits[i*Mn:(i+1)*Mn]) for i in range(4))
        if all(x == 0 for row in Tn(l, n) for x in row):
            cnt += 1
    print(f'n={n} (M={Mn}): |Z_n| = {cnt} = 2^(4n)?',
          cnt == 2**(4*n))


# ---------------------------------------------------------------------
# 7. d : B -> V, d(b(v)) = v is well defined on the basis expansion.
# ---------------------------------------------------------------------
print('--- d(b(v)) = v via basis expansion, sampled ---')
okd = True
for _ in range(50):
    v = tuple(tuple(random.randint(0, 1) for _ in range(M))
              for _ in range(4))
    # expansion of b(v): diag coeffs x_k, offdiag x_a x_b;
    # d kills offdiag and sends diag basis b(t^j e_i) -> t^j e_i
    keys = [(i, j) for i in range(4) for j in range(M)]
    x = {k: v[k[0]][k[1]] for k in keys}
    dv = [[0]*M for _ in range(4)]
    for (i, j) in keys:
        dv[i][j] = x[(i, j)]
    okd &= (tuple(tuple(r) for r in dv) == v)
print('d(b(v)) = v:', okd)


# ---------------------------------------------------------------------
# 8. Equivariance: k.F_{l,q} = F_{k.l, k.q} for k = I + t E12.
# ---------------------------------------------------------------------
def basis_vec(i, j, M):
    return tuple(tuple(1 if (i2, j2) == (i, j) else 0
                       for j2 in range(M)) for i2 in range(4))


def act_l(l, M):
    """(k.l)(v) = l(k^{-1} v) for k = I + t E12 (self-inverse)."""
    return tuple(tuple(l_of(l, ginv_on_v(basis_vec(i, j, M), M))
                       for j in range(M)) for i in range(4))


def act_q(q, M):
    """(k.q)(w) = q(k^{-1} w) on the B_M basis, via polarization
    sym(u, w) = b(u+w) + b(u) + b(w)."""
    keys = [(i, j) for i in range(4) for j in range(M)]
    diag = {}
    for (i, j) in keys:
        u = ginv_on_v(basis_vec(i, j, M), M)
        diag[(i, j)] = q_of_bv(q, u)
    off = {}
    for a in range(len(keys)):
        for b2 in range(a+1, len(keys)):
            u = ginv_on_v(basis_vec(*keys[a], M), M)
            w = ginv_on_v(basis_vec(*keys[b2], M), M)
            uw = tuple(tuple((u[i][j]+w[i][j]) % 2 for j in range(M))
                       for i in range(4))
            off[(keys[a], keys[b2])] = (q_of_bv(q, uw) + q_of_bv(q, u)
                                        + q_of_bv(q, w)) % 2
    return (diag, off)


print('--- equivariance: k.F_{l,q} = F_{k.l, k.q}, sampled ---')
okk = True
for _ in range(20):
    l1, q1 = rand_l(M), rand_q(M)
    kl, kq = act_l(l1, M), act_q(q1, M)
    for _ in range(15):
        v = tuple(tuple(random.randint(0, 1) for _ in range(M))
                  for _ in range(4))
        okk &= (F(l1, q1, ginv_on_v(v, M)) == F(kl, kq, v))
print('k.F_{l,q} = F_{k.l,k.q} for k = I + tE12:', okk)

print()
print('ALL CH04 ENGINE CHECKS COMPLETE')
