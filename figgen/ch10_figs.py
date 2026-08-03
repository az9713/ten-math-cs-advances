"""Figure-data generator + numeric validation for the ch10 explainer.

Run: PYTHONIOENCODING=utf-8 python figgen/ch10_figs.py
Outputs: validation numbers quoted in the HTML, plus SVG snippet files in
.ignore/ch10_snips/ (bare <svg> content; captions live in the HTML).
"""
import io
import itertools
import math
import os

import numpy as np

OUT = '.ignore/ch10_snips'
os.makedirs(OUT, exist_ok=True)
LOG2 = math.log(2.0)


def h(x):
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)


# =====================================================================
# 1. The symplectic generalized quadrangle W(q), q = 2, 3
# =====================================================================
def canon(v, q):
    for c in v:
        if c % q:
            inv = pow(c, q - 2, q) if q > 2 else 1
            return tuple((x * inv) % q for x in v)
    return None


def symp(x, y, q):
    return (x[0] * y[1] - x[1] * y[0] + x[2] * y[3] - x[3] * y[2]) % q


def build_wq(q):
    vecs = [v for v in itertools.product(range(q), repeat=4)
            if any(v)]
    pts = sorted(set(canon(v, q) for v in vecs))
    # totally isotropic 2-spaces = lines
    lines = set()
    for i, p in enumerate(pts):
        for pp in pts[i + 1:]:
            if symp(p, pp, q) == 0:
                span = set()
                for a in range(q):
                    for b in range(q):
                        w = tuple((a * p[k] + b * pp[k]) % q
                                  for k in range(4))
                        if any(w):
                            span.add(canon(w, q))
                lines.add(frozenset(span))
    return pts, sorted(lines, key=sorted)


def girth(adj):
    n = len(adj)
    best = 10 ** 9
    for root in range(n):
        dist = {root: 0}
        par = {root: -1}
        queue = [root]
        while queue:
            nxt = []
            for u in queue:
                for v in adj[u]:
                    if v not in dist:
                        dist[v] = dist[u] + 1
                        par[v] = u
                        nxt.append(v)
                    elif v != par[u]:
                        best = min(best, dist[u] + dist[v] + 1)
            queue = nxt
    return best


def wq_checks(q):
    pts, lines = build_wq(q)
    npts, nlin = len(pts), len(lines)
    want = (q + 1) * (q * q + 1)
    n_v = npts + nlin
    adj = [set() for _ in range(n_v)]
    for j, ln in enumerate(lines):
        for p in ln:
            i = pts.index(p)
            adj[i].add(npts + j)
            adj[npts + j].add(i)
    degs_p = sorted(len(adj[i]) for i in range(npts))
    degs_l = sorted(len(adj[npts + j]) for j in range(nlin))
    g = girth(adj)
    e = sum(len(a) for a in adj) // 2
    nq = n_v
    print(f'W({q}): |P|={npts} |M|={nlin} (want {want} each), '
          f'point-deg {degs_p[0]}..{degs_p[-1]}, '
          f'line-size {degs_l[0]}..{degs_l[-1]} (want {q + 1})')
    print(f'      girth={g} (want 8), n_q={nq}, e_q={e} '
          f'(formula {(q + 1) ** 2 * (q * q + 1)}), '
          f'2^(-4/3) n^(4/3)={2 ** (-4 / 3) * nq ** (4 / 3):.1f}')
    # relations ~ on each side
    relP = [[(symp(pts[i], pts[j], q) == 0 and i != j)
             for j in range(npts)] for i in range(npts)]
    lin_l = [set(ln) for ln in lines]
    relM = [[(i != j and len(lin_l[i] & lin_l[j]) > 0)
             for j in range(nlin)] for i in range(nlin)]
    # cross-check ~ against common-neighbor definition in the graph
    for i in range(npts):
        for j in range(npts):
            cn = len(adj[i] & adj[j]) > 0 and i != j
            assert cn == relP[i][j], 'relP mismatch'
    for i in range(nlin):
        for j in range(nlin):
            cn = (len(adj[npts + i] & adj[npts + j]) > 0 and i != j)
            assert cn == relM[i][j], 'relM mismatch'
    return pts, lines, relP, relM


def side_patterns(rel, name):
    """Search J-patterns and S3-patterns on one side."""
    n = len(rel)
    idx = range(n)
    r_of = {}
    for t in itertools.combinations(idx, 3):
        a, b, c = t
        if rel[a][b] or rel[a][c] or rel[b][c]:
            continue
        r = sum(1 for w in idx
                if w not in t and rel[w][a] and rel[w][b] and rel[w][c])
        r_of[t] = r
    dist = {}
    for r in r_of.values():
        dist[r] = dist.get(r, 0) + 1
    s3 = sum(v for k, v in dist.items() if k >= 3)
    # J-pattern: triples {x,y,z},{x',y,z} with r>=2 each and x ~ x'
    jpat = 0
    good = [t for t, r in r_of.items() if r >= 2]
    by_pair = {}
    for t in good:
        for pair in itertools.combinations(t, 2):
            third = [v for v in t if v not in pair][0]
            by_pair.setdefault(pair, []).append(third)
    for pair, thirds in by_pair.items():
        for x, xp in itertools.combinations(thirds, 2):
            if rel[x][xp]:
                jpat += 1
    print(f'   side {name}: triad r-distribution {dist}  '
          f'S3-patterns={s3}  J-patterns={jpat}')
    return s3, jpat


print('--- generalized quadrangle witnesses ---')
ptsA, linA, relPA, relMA = wq_checks(2)
s3p, jp = side_patterns(relPA, 'P q=2')
s3m, jm = side_patterns(relMA, 'M q=2')
assert jp == 0 and jm == 0, 'W(2) must be J-pattern-free on both sides'
pts3, lin3, relP3, relM3 = wq_checks(3)
s3m3, jm3 = side_patterns(relM3, 'M q=3')
assert s3m3 == 0, 'W(3) dual side must be S3-pattern-free'

# =====================================================================
# 2. Doily drawing of W(2)  (duads/synthemes model)
# =====================================================================
duads = [frozenset(d) for d in itertools.combinations(range(1, 7), 2)]


def synthemes():
    out = []
    elems = set(range(1, 7))
    first = 1
    for a in range(2, 7):
        rest = sorted(elems - {first, a})
        b, c, d, e = rest
        for pair in ((b, c, d, e), (b, d, c, e), (b, e, c, d)):
            out.append(frozenset({frozenset({first, a}),
                                  frozenset({pair[0], pair[1]}),
                                  frozenset({pair[2], pair[3]})}))
    return sorted(set(out), key=lambda s: sorted(map(sorted, s)))


syns = synthemes()
assert len(syns) == 15
cnt = {d: 0 for d in duads}
for s in syns:
    for d in s:
        cnt[d] += 1
assert all(v == 3 for v in cnt.values()), 'each duad on 3 synthemes'
# incidence-graph girth of the doily model
d_ix = {d: i for i, d in enumerate(duads)}
adj_d = [set() for _ in range(30)]
for j, s in enumerate(syns):
    for d in s:
        adj_d[d_ix[d]].add(15 + j)
        adj_d[15 + j].add(d_ix[d])
print('doily model: girth =', girth(adj_d), '(want 8; GQ(2,2) is unique,'
      ' so this IS W(2))')

CX, CY = 210.0, 168.0


def ang(i):
    return math.radians(-90 + 72 * (i - 1))


def polar(d):
    d = set(d)
    if 6 in d:
        i = (d - {6}).pop()
        return 136.0, ang(i)
    a, b = sorted(d)
    if (b - a) in (1, 4):
        i = a if b - a == 1 else 5
        return 88.0, ang(i) + math.radians(36)
    i = a if b - a == 2 else b
    return 46.0, ang((i % 5) + 1) + math.radians(12)


def pos(d):
    r, th = polar(d)
    return (CX + r * math.cos(th), CY + r * math.sin(th))


def bez3(p1, p2, p3):
    """Quadratic Bezier p1->p3 passing through p2 at t=1/2."""
    cx = 2 * p2[0] - (p1[0] + p3[0]) / 2
    cy = 2 * p2[1] - (p1[1] + p3[1]) / 2
    return (f'M{p1[0]:.1f},{p1[1]:.1f} Q{cx:.1f},{cy:.1f} '
            f'{p3[0]:.1f},{p3[1]:.1f}')


def line_type(s):
    outer = [d for d in s if 6 in d]
    inner = [d for d in s if 6 not in d and
             (max(d) - min(d)) in (2, 3)]
    return len(inner)  # 1 -> L1, 0 -> L2, 2 -> L3


COLS = {1: '#c0392b', 0: '#2e6da4', 2: '#2a7d2a'}
svg = ['<svg class="setupfig" viewBox="0 0 420 336" width="86%" '
       'role="img" aria-label="The doily: the generalized quadrangle '
       'W(2), 15 points and 15 lines">']
for s in syns:
    t = line_type(s)
    ds = list(s)
    ps = [pos(d) for d in ds]
    # order: middle point = the one nearest the other two combined
    best, order = None, None
    for perm in itertools.permutations(range(3)):
        a, b, c = perm
        length = (math.dist(ps[a], ps[b]) + math.dist(ps[b], ps[c]))
        if best is None or length < best:
            best, order = length, perm
    p1, p2, p3 = (ps[k] for k in order)
    svg.append(f'  <path d="{bez3(p1, p2, p3)}" fill="none" '
               f'stroke="{COLS[t]}" stroke-width="1.6" opacity="0.85"/>')
for d in duads:
    x, y = pos(d)
    svg.append(f'  <circle cx="{x:.1f}" cy="{y:.1f}" r="5.2" '
               f'fill="#f6f4ef" stroke="#1a1a1a" stroke-width="1.4"/>')
for d in duads:
    r, th = polar(d)
    lab = ''.join(str(t) for t in sorted(d))
    rl = r + 15 if r > 50 else r - 17
    lx = CX + rl * math.cos(th)
    ly = CY + rl * math.sin(th) + 3.5
    svg.append(f'  <text x="{lx:.1f}" y="{ly:.1f}" font-size="10.5"'
               f' fill="#555" text-anchor="middle">{lab}</text>')
svg.append('</svg>')
io.open(f'{OUT}/doily.svg', 'w', encoding='utf-8').write('\n'.join(svg))
print('wrote doily.svg')

# =====================================================================
# 3. Projective planes PG(2,p): the C4-free witnesses for Section 1
# =====================================================================
def pg2(p):
    vecs = [v for v in itertools.product(range(p), repeat=3) if any(v)]
    seen = set()
    pts = []
    for v in vecs:
        for c in v:
            if c % p:
                inv = pow(c, p - 2, p)
                cv = tuple((x * inv) % p for x in v)
                break
        if cv not in seen:
            seen.add(cv)
            pts.append(cv)
    return pts


def pg_girth(p):
    pts = pg2(p)
    n = len(pts)
    adj = [set() for _ in range(2 * n)]
    for i, pt in enumerate(pts):
        for j, ln in enumerate(pts):
            if sum(a * b for a, b in zip(pt, ln)) % p == 0:
                adj[i].add(n + j)
                adj[n + j].add(i)
    return len(pts), girth(adj)


print('--- projective-plane witnesses (Section 1 toy) ---')
for p in (2, 3, 5):
    npt, g = pg_girth(p)
    print(f'PG(2,{p}): {npt} points (want {p * p + p + 1}), incidence '
          f'girth {g} (want 6 => C4-free)')

# exact ex(n, C4) for tiny n by brute force
def ex_c4(n):
    pairs = list(itertools.combinations(range(n), 2))
    best = 0
    for mask in range(1 << len(pairs)):
        if bin(mask).count('1') <= best:
            continue
        rows = [0] * n
        for b, (i, j) in enumerate(pairs):
            if mask >> b & 1:
                rows[i] |= 1 << j
                rows[j] |= 1 << i
        ok = True
        for i, j in pairs:
            if bin(rows[i] & rows[j]).count('1') >= 2:
                ok = False
                break
        if ok:
            best = bin(mask).count('1')
    return best


for n in (4, 5, 6):
    print(f'ex({n}, C4) = {ex_c4(n)} (brute force)')

# chart data: KST upper bound vs incidence-graph lower bound, log-log
X0, X1, Y0, Y1 = 60.0, 440.0, 210.0, 25.0
LX0, LX1 = 0.6, 3.6      # log10 n
LY0, LY1 = 0.3, 5.4      # log10 m


def cx(ln):
    return X0 + (ln - LX0) / (LX1 - LX0) * (X1 - X0)


def cy(lm):
    return Y0 + (lm - LY0) / (LY1 - LY0) * (Y1 - Y0)


pts_kst = []
for ln in np.linspace(LX0, LX1, 60):
    n = 10 ** ln
    m = n / 4 * (1 + math.sqrt(4 * n - 3))
    pts_kst.append(f'{cx(ln):.1f},{cy(math.log10(m)):.1f}')
print('KST_curve:', ' '.join(pts_kst))
lower = []
for p in (2, 3, 5, 7, 11, 13, 17, 19, 23):
    n = 2 * (p * p + p + 1)
    e = (p + 1) * (p * p + p + 1)
    lower.append((n, e))
    print(f'  PG(2,{p}): n={n} edges={e} '
          f'edge/KST-ratio={e / (n / 4 * (1 + math.sqrt(4 * n - 3))):.3f}')
print('PG_points:', ' '.join(
    f'{cx(math.log10(n)):.1f},{cy(math.log10(e)):.1f}' for n, e in lower))
ref = [f'{cx(ln):.1f},{cy(1.5 * ln - 0.85):.1f}'
       for ln in np.linspace(LX0, LX1, 2)]
print('slope32_ref:', ' '.join(ref))
for ln in (1, 2, 3):
    print(f'  xtick n=10^{ln}: x={cx(ln):.1f}', end='')
for lm in (1, 2, 3, 4, 5):
    print(f'  ytick m=10^{lm}: y={cy(lm):.1f}', end='')
print()

# =====================================================================
# 4. Entropy numerics: kappa, window, epsilon, sizes
# =====================================================================
print('--- entropy-side numerics ---')
kappa = 1.5 - 0.75 * math.log2(3)
r3 = math.sqrt(3)
tau_s = 1 / (1 + r3)


def A(t):
    return kappa + t * math.log2(3)


def C(t):
    return 2 * h(t) - 1


def f(t):
    return C(t) - A(t)


w_closed = 0.5 * math.log2(1 + (r3 - 1) ** 4 / (8 * r3 * (1 + r3 ** 2)))
print(f'kappa = {kappa:.6f}')
print(f'tau* = 1/(1+sqrt3) = {tau_s:.6f}')
print(f'A(tau*) = {A(tau_s):.6f}   C(tau*) = {C(tau_s):.6f}')
print(f'window width f(tau*) = {f(tau_s):.6f}  '
      f'(closed form {w_closed:.6f})')


def bisect(fun, lo, hi):
    for _ in range(200):
        mid = (lo + hi) / 2
        if fun(lo) * fun(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


t1 = bisect(f, 0.05, tau_s)
t2 = bisect(f, tau_s, 0.49)
print(f'window nonempty for tau in ({t1:.4f}, {t2:.4f})')
eps_grid = [(f(t) / (2 * (1 - A(t))), t)
            for t in np.linspace(t1 + 1e-6, t2 - 1e-6, 20001)]
eps_max, t_eps = max(eps_grid)
print(f'sup epsilon = {eps_max:.6f} at tau = {t_eps:.4f}  '
      f'=> exponent 3/2+eps < {1.5 + eps_max:.4f}')
beta = (A(tau_s) + C(tau_s)) / 2
delta = (beta - A(tau_s)) / 8
inc = 2 * (beta - A(tau_s) - 2 * delta)
s_min = math.floor(1 / inc) + 1
print(f'example: beta={beta:.6f} delta={delta:.2e} '
      f'increment 2(beta-A-2delta)={inc:.6f}  s_min={s_min}')


def eta(loc_l):
    return (1 + math.log2(3)) / loc_l + h(1 / loc_l) / 2


l0 = 4
while True:
    m_pairs = l0 * (l0 - 1) // 2
    ok1 = eta(l0) < delta
    ok2 = (l0 + 3 * math.log2(m_pairs + 1) - delta * m_pairs) < -1
    if ok1 and ok2:
        break
    l0 += 1
print(f'minimal L0 for (14): {l0}  (eta(L0)={eta(l0):.2e} < '
      f'delta={delta:.2e})')
a0 = math.log2(l0)
a_s = 1 + 2 ** s_min * (a0 - 1)      # log2 of top-layer size, approx
print(f'log2 L_s ~ 1 + 2^{s_min}(log2 L0 - 1) = 2^{s_min} * '
      f'{a0 - 1:.2f}; log10(log2|V(H)|) ~ '
      f'{s_min * math.log10(2) + math.log10(a0 - 1):.1f}')

# Hamming ball size check
for m in (200, 2000):
    k = math.floor(tau_s * m)
    lg = 0.0
    tot = 0.0
    logs = [sum(math.log2((m - i) / (i + 1)) for i in range(j))
            for j in range(k + 1)]
    mx = max(logs)
    tot = sum(2 ** (v - mx) for v in logs)
    lg = mx + math.log2(tot)
    print(f'm={m}: (1/m) log2 D_m = {lg / m:.5f} vs h(tau*)='
          f'{h(tau_s):.5f}')

# brute-force check of Lemma 5.1 over all binary conditional laws
grid = np.linspace(0.0, 1.0, 161)
r00 = grid[:, None, None]
r11 = grid[None, :, None]
r01 = grid[None, None, :]


def hv(x):
    x = np.clip(x, 1e-12, 1 - 1e-12)
    return -x * np.log2(x) - (1 - x) * np.log2(1 - x)


gmax, arg = -9.9, None
for qv in np.linspace(0.02, 0.98, 49):
    u = 2 * qv - 1
    a0_, a1_, ast = (1 - u) ** 2 / 4, (1 + u) ** 2 / 4, (1 - u ** 2) / 2
    hz = a0_ * hv(r00) + a1_ * hv(r11) + ast * hv(r01)
    v = a0_ * r00 + a1_ * r11 + ast * r01
    d = a0_ * r00 + a1_ * (1 - r11) + ast / 2
    gg = hz - math.log2(3) * d - (hv(v) - h(qv)) / 2
    i = np.unravel_index(np.argmax(gg), gg.shape)
    if gg[i] > gmax:
        gmax, arg = float(gg[i]), (qv, grid[i[0]], grid[i[1]],
                                   grid[i[2]])
print(f'brute-force max of H(Z|XY)-(log2 3)d-(h(v)-h(q))/2 = '
      f'{gmax:.6f} (Lemma 5.1 bound kappa={kappa:.6f}) at q={arg[0]:.3f}')
assert gmax <= kappa + 1e-9

# g(u) curve for the Lemma 5.1 finish + quadratic envelope
GX0, GX1, GY0, GY1 = 60.0, 440.0, 200.0, 30.0
GU0, GU1 = -1.0, 1.0
GV0, GV1 = 0.05, 0.34


def gx(u):
    return GX0 + (u - GU0) / (GU1 - GU0) * (GX1 - GX0)


def gy(v):
    return GY0 + (v - GV0) / (GV1 - GV0) * (GY1 - GY0)


def g_rhs(u):
    return (0.5 * h((1 + u) / 2) + 1 - 0.75 * math.log2(3)
            + u * u / 4 * math.log2(4 / 3)
            + math.log2(math.e) * (math.sqrt(1 + u * u / 4) - 1))


coef = -(math.log2(math.e) / 8 - math.log2(4 / 3) / 4)
print(f'quadratic envelope coefficient = {coef:.6f} (must be < 0)')
us = np.linspace(-1, 1, 121)
print('g_curve:', ' '.join(f'{gx(u):.1f},{gy(g_rhs(u)):.1f}' for u in us))
print('env_curve:', ' '.join(
    f'{gx(u):.1f},{gy(kappa + coef * u * u):.1f}' for u in us))
print('kappa_y:', f'{gy(kappa):.1f}')
print('g ticks: u=-1:', f'{gx(-1):.1f}', ' u=0:', f'{gx(0):.1f}',
      ' u=1:', f'{gx(1):.1f}', ' v=0.1:', f'{gy(0.1):.1f}',
      ' v=0.2:', f'{gy(0.2):.1f}', ' v=0.3:', f'{gy(0.3):.1f}')
print(f'g(-1)=g(1)={g_rhs(1.0):.6f}  g(0)={g_rhs(0.0):.6f}')

# window plot (main + zoom)
WX0, WX1, WY0, WY1 = 55.0, 240.0, 205.0, 30.0
WT0, WT1 = 0.02, 0.48
WV0, WV1 = -0.8, 1.1
ZX0, ZX1, ZY0, ZY1 = 285.0, 445.0, 205.0, 30.0
ZT0, ZT1 = 0.33, 0.40
ZV0, ZV1 = 0.82, 0.95


def wx(t):
    return WX0 + (t - WT0) / (WT1 - WT0) * (WX1 - WX0)


def wy(v):
    return WY0 + (v - WV0) / (WV1 - WV0) * (WY1 - WY0)


def zx(t):
    return ZX0 + (t - ZT0) / (ZT1 - ZT0) * (ZX1 - ZX0)


def zy(v):
    return ZY0 + (v - ZV0) / (ZV1 - ZV0) * (ZY1 - ZY0)


ts = np.linspace(WT0, WT1, 90)
print('A_main:', ' '.join(f'{wx(t):.1f},{wy(A(t)):.1f}' for t in ts))
print('C_main:', ' '.join(f'{wx(t):.1f},{wy(C(t)):.1f}' for t in ts))
tz = np.linspace(ZT0, ZT1, 90)
print('A_zoom:', ' '.join(f'{zx(t):.1f},{zy(A(t)):.1f}' for t in tz))
print('C_zoom:', ' '.join(f'{zx(t):.1f},{zy(C(t)):.1f}' for t in tz))
tw = np.linspace(t1, t2, 60)
shade = ([f'{zx(t):.1f},{zy(C(t)):.1f}' for t in tw]
         + [f'{zx(t):.1f},{zy(A(t)):.1f}' for t in tw[::-1]])
print('win_shade:', ' '.join(shade))
print('zoom ticks: t=0.34:', f'{zx(0.34):.1f}', 't=0.39:',
      f'{zx(0.39):.1f}', 'v=0.85:', f'{zy(0.85):.1f}', 'v=0.90:',
      f'{zy(0.90):.1f}', 'v=0.94:', f'{zy(0.94):.1f}',
      'tau*:', f'{zx(tau_s):.1f}',
      't1:', f'{zx(t1):.1f}', 't2:', f'{zx(t2):.1f}',
      'A(t1):', f'{zy(A(t1)):.1f}', 'A(t2):', f'{zy(A(t2)):.1f}')
print('main ticks: t=0.1:', f'{wx(0.1):.1f}', 't=0.3:', f'{wx(0.3):.1f}',
      't=0.45:', f'{wx(0.45):.1f}',
      'v=0:', f'{wy(0.0):.1f}', 'v=0.5:', f'{wy(0.5):.1f}',
      'v=1.0:', f'{wy(1.0):.1f}', 'v=-0.5:', f'{wy(-0.5):.1f}',
      'zoombox:', f'{wx(ZT0):.1f},{wy(ZV1):.1f} to '
      f'{wx(ZT1):.1f},{wy(ZV0):.1f}')

# =====================================================================
# 5. The layered graph H (drawing data, L0 = 4)
# =====================================================================
v0 = list(range(1, 5))
v1 = list(itertools.combinations(v0, 2))
v2 = list(itertools.combinations(v1, 2))
print(f'--- layered H drawing: |V0|={len(v0)} |V1|={len(v1)} '
      f'|V2|={len(v2)} ---')
W, HY0, HY1, HY2 = 520.0, 292.0, 180.0, 52.0


def spread(k, width, margin=30.0):
    if k == 1:
        return [width / 2]
    return [margin + i * (width - 2 * margin) / (k - 1)
            for i in range(k)]


x0 = spread(4, W, 130)
x1 = spread(6, W, 90)
x2 = spread(15, W, 26)
p0 = {a: (x0[i], HY0) for i, a in enumerate(v0)}
p1 = {a: (x1[i], HY1) for i, a in enumerate(v1)}
p2 = {a: (x2[i], HY2) for i, a in enumerate(v2)}
seg1 = []
for a in v1:
    for par in a:
        seg1.append((p1[a], p0[par]))
seg2 = []
for a in v2:
    for par in a:
        seg2.append((p2[a], p1[par]))
print('H_edges1:', ' '.join(
    f'{s[0][0]:.0f},{s[0][1]:.0f},{s[1][0]:.0f},{s[1][1]:.0f}'
    for s in seg1))
print('H_edges2:', ' '.join(
    f'{s[0][0]:.0f},{s[0][1]:.0f},{s[1][0]:.0f},{s[1][1]:.0f}'
    for s in seg2))
print('H_v0:', ' '.join(f'{p[0]:.0f},{p[1]:.0f}' for p in p0.values()))
print('H_v1:', ' '.join(f'{p[0]:.0f},{p[1]:.0f}' for p in p1.values()))
print('H_v2:', ' '.join(f'{p[0]:.0f},{p[1]:.0f}' for p in p2.values()))
print('H_v1_labels:', ' '.join(
    f'{"".join(map(str, a))}@{p1[a][0]:.0f}' for a in v1))

# =====================================================================
# 6. Compactness gap: n^{1/48} punchline numbers
# =====================================================================
print('--- compactness gap numbers ---')
print('21/16 =', 21 / 16, ' 4/3 - 1/48 =', 4 / 3 - 1 / 48)
for e10 in (6, 12, 24, 48):
    print(f'  n=10^{e10}: family loses factor n^(1/48) = '
          f'{10 ** (e10 / 48):.3f}')
print('staircase: increment per layer =', f'{inc:.5f}',
      ' layers to exceed 1 =', s_min)
