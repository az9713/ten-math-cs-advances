"""Verification engine + figure data for the ch03 explainer (non-sofic
group).  Implements the binary Leavitt algebra R = L_F2(1,2) as operators
on depth-D binary words and machine-checks every identity quoted in the
chapter.

Run: PYTHONIOENCODING=utf-8 python figgen/ch03_figs.py
"""
import itertools


# ---------------------------------------------------------------------
# The algebra: elements are F2-sums of monomials s_a t_b, stored as a
# frozenset of pairs (a, b) of binary words (XOR semantics).
# ---------------------------------------------------------------------
def el(*pairs):
    out = set()
    for p in pairs:
        out ^= {p}
    return frozenset(out)


ONE = el(('', ''))
ZERO = el()


def s(a):
    return el((a, ''))


def t(b):
    return el(('', b))


def add(x, y):
    return frozenset(set(x) ^ set(y))


def mul_mono(m1, m2):
    (a, b), (c, d) = m1, m2
    if c.startswith(b):
        return (a + c[len(b):], d)
    if b.startswith(c):
        return (a, d + b[len(c):])
    return None


def mul(x, y):
    out = set()
    for m1 in x:
        for m2 in y:
            m = mul_mono(m1, m2)
            if m is not None:
                out ^= {m}
    return frozenset(out)


def codepth(x):
    return max((len(b) for (a, b) in x), default=0)


def act_word(x, w):
    """Image of the depth-|w| cylinder symbol [w]: XOR-set of words."""
    out = set()
    for (a, b) in x:
        if w.startswith(b):
            out ^= {a + w[len(b):]}
    return out


def is_zero(x, depth=None):
    """Exact zero test: act on all words of length co-depth (or given)."""
    d = codepth(x) if depth is None else depth
    for w in itertools.product('01', repeat=d):
        if act_word(x, ''.join(w)):
            return False
    return True


def eq(x, y, depth=None):
    return is_zero(add(x, y), depth)


def e(a):
    return mul(s(a), t(a))


print('--- Leavitt relations ---')
for i in '01':
    for j in '01':
        want = ONE if i == j else ZERO
        print(f't{i}s{j} = {"1" if i == j else "0"}:',
              eq(mul(t(i), s(j)), want))
print('s0t0+s1t1 = 1:', eq(add(e('0'), e('1')), ONE))
print('t0s0=1 but s0t0 proper idempotent:',
      eq(mul(t('0'), s('0')), ONE),
      not eq(e('0'), ONE),
      eq(mul(e('0'), e('0')), e('0')))

# prefix-code identities
print('--- prefix-code identities ---')
ok = True
for a in ('000', '001', '01'):
    for b in ('000', '001', '01'):
        want = ONE if a == b else ZERO
        ok &= eq(mul(t(a), s(b)), want)
print('ta_i s_a_j = delta_ij on alpha:', ok)
print('e_a = e_a0 + e_a1:', eq(e('01'), add(e('010'), e('011'))))

ALPHA = ('000', '001', '01')
BETA = ('1000', '1001', '101')
NU = ('1100', '1101', '111')
ZETA = ('100', '101', '11')
D9 = ALPHA + BETA + NU
eD = ZERO
for a in D9:
    eD = add(eD, e(a))
print('e_D = 1 (nine-leaf code complete):', eq(eD, ONE))
eA = ZERO
for a in ALPHA:
    eA = add(eA, e(a))
print('e_alpha = e_0:', eq(eA, e('0')))


# units from prefix replacement tables
def U(table):
    """table: dict source-word -> target-word (complete codes)."""
    out = ZERO
    for a, b in table.items():
        out = add(out, mul(s(b), t(a)))
    return out


u_tab = {}
v_tab = {}
for i in range(3):
    u_tab[ALPHA[i]] = ALPHA[i] + '0'
    u_tab[BETA[i]] = ALPHA[i] + '1'
    u_tab[NU[i]] = ZETA[i]
    v_tab[ALPHA[i]] = ALPHA[i] + '0'
    v_tab[BETA[i]] = ZETA[i]
    v_tab[NU[i]] = ALPHA[i] + '1'
uu = U(u_tab)
vv = U(v_tab)
uinv = U({b: a for a, b in u_tab.items()})
vinv = U({b: a for a, b in v_tab.items()})
print('--- the units u, v ---')
print('u unit:', eq(mul(uu, uinv), ONE), eq(mul(uinv, uu), ONE))
print('v unit:', eq(mul(vv, vinv), ONE), eq(mul(vinv, vv), ONE))


def conj(g, ginv, x):
    return mul(mul(g, x), ginv)


def elem_root(a, b, r):
    return add(ONE, mul(mul(s(a), r), t(b)))


print('--- conjugation u(1+s_ai r t_aj)u^-1 = 1+s_ai0 r t_aj0 ---')
ok = True
for r in (ONE, s('0'), t('1'), mul(s('1'), t('0'))):
    for (i, j) in ((0, 1), (1, 2), (2, 0)):
        lhs = conj(uu, uinv, elem_root(ALPHA[i], ALPHA[j], r))
        rhs = elem_root(ALPHA[i] + '0', ALPHA[j] + '0', r)
        ok &= eq(lhs, rhs)
        lhs = conj(vv, vinv, elem_root(ALPHA[i], ALPHA[j], r))
        ok &= eq(lhs, rhs)
print('holds for sample r in {1, s0, t1, s1t0}, both u and v:', ok)

# sai0 r taj0 = sai (s0 r t0) taj  -> conjugates lie in Gamma
x = mul(mul(s(ALPHA[0] + '0'), s('1')), t(ALPHA[1] + '0'))
y = mul(mul(s(ALPHA[0]), mul(mul(s('0'), s('1')), t('0'))), t(ALPHA[1]))
print('sai0 r taj0 = sai(s0 r t0)taj sample:', eq(x, y))

# cylinder swaps and the characteristic-2 identity
print('--- cylinder swaps ---')


def swap(rho, sig):
    return add(add(add(add(ONE, e(rho)), e(sig)),
                   mul(s(rho), t(sig))), mul(s(sig), t(rho)))


rho, sig = '0001', '01'
P = mul(s(rho), t(sig))
Q = mul(s(sig), t(rho))
lhs = mul(mul(add(ONE, P), add(ONE, Q)), add(ONE, P))
print('(1+P)(1+Q)(1+P) = swap:', eq(lhs, swap(rho, sig)))
tau = swap(rho, sig)
print('swap is an involution:', eq(mul(tau, tau), ONE))

# Steinberg-type identities
print('--- Steinberg identities ---')
a, b, c = '000', '001', '01'
r1, r2 = s('1'), t('0')
x1 = elem_root(a, c, r1)
x2 = elem_root(c, b, r2)
x1i = elem_root(a, c, r1)   # char 2: own inverse
comm = mul(mul(mul(x1, x2), x1i), elem_root(c, b, r2))
print('[1+sa r tc, 1+sc r2 tb] = 1+sa r r2 tb:',
      eq(comm, elem_root(a, b, mul(r1, r2))))
print('xij(r)xij(r2) = xij(r+r2):',
      eq(mul(elem_root(a, b, r1), elem_root(a, b, r2)),
         elem_root(a, b, add(r1, r2))))

# commuting supports: Gamma on [0], J on [1000]
print('--- orthogonal supports ---')
print('e0 e1000 = 0:', eq(mul(e('0'), e('1000')), ZERO))
g = elem_root(ALPHA[0], ALPHA[1], s('1'))          # in Gamma
jrep = swap('10000', '10001')                      # in J = V_1000
print('sample g in Gamma, j in J commute:',
      eq(mul(g, jrep), mul(jrep, g)))
print('g-1 in e0Re0:', is_zero(add(mul(mul(e('0'), add(g, ONE)),
                                       e('0')), add(g, ONE))))
print('j-1 in e1000Re1000:',
      is_zero(add(mul(mul(e('1000'), add(jrep, ONE)), e('1000')),
                  add(jrep, ONE))))

# u moves J into [0001] inside [0]
print('--- u J u^-1 = V_0001 ---')
jc = conj(uu, uinv, jrep)
print('u j u^-1 = swap(00010,00011):',
      eq(jc, swap('00010', '00011')))
print('u j u^-1 - 1 in e0001 R e0001:',
      is_zero(add(mul(mul(e('0001'), add(jc, ONE)), e('0001')),
                  add(jc, ONE))))

# the X_i r Y_j expansion used to recover G
print('--- block recovery identity ---')
Xi = add(mul(s(ALPHA[0]), t('0')), mul(s(BETA[0]), t('1')))
Yj = add(mul(s('0'), t(ALPHA[1])), mul(s('1'), t(BETA[1])))
B = {('0', '0'): ONE, ('0', '1'): s('1'),
     ('1', '0'): ZERO, ('1', '1'): t('0')}
r = ZERO
for (p, q), bpq in B.items():
    r = add(r, mul(mul(s(p), bpq), t(q)))
lhs = mul(mul(Xi, r), Yj)
rhs = ZERO
LL = {('0', '0'): (ALPHA[0], ALPHA[1]), ('0', '1'): (ALPHA[0], BETA[1]),
      ('1', '0'): (BETA[0], ALPHA[1]), ('1', '1'): (BETA[0], BETA[1])}
for (p, q), bpq in B.items():
    la, lb = LL[(p, q)]
    rhs = add(rhs, mul(mul(s(la), bpq), t(lb)))
print('Xi r Yj = sum s_l b t_l:', eq(lhs, rhs))
print("u^-1(1+s_a1 r t_a2)u = 1+X1 r Y2 sample:",
      eq(conj(uinv, uu, elem_root(ALPHA[0], ALPHA[1], s('1'))),
         add(ONE, mul(mul(Xi, s('1')), Yj))))

# U_g well-definedness: refining a table entry leaves U unchanged
print('--- table refinement invariance ---')
g1 = U({'0': '10', '10': '0', '11': '11'})
g2 = U({'00': '100', '01': '101', '10': '0', '11': '11'})
print('refined table gives same unit:', eq(g1, g2))
print('example swap unit is an involution:', eq(mul(g1, g1), ONE))
