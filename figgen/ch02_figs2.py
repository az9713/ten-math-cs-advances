"""Figures 18-34 for chapter 2 (see ledger in ch02_figs.py)."""
import math

from ch02_figlib import (RED, BLUE, GREEN, GOLD, INK, MUT, LT, LG, LY,
                         ACC, svg, ar, txt, sub, sup, line, poly, circ,
                         rect, path, mapper, ortho)
from ch02_figs import (H2, Hsph, qfun, a0fun, BKL, GammaH, kappaH, M1,
                       M2, lam_tri, smooth_lower)


# ================= 18 fig_layers =================
def fig_layers():
    b = ''
    n = 30
    x0, y0 = 250, 320
    for w in range(n+1):
        wd = 3 + 300*math.comb(n, w)/math.comb(n, 15)
        y = y0 - w*9.5
        col = RED if w == 12 else MUT
        op = 0.9 if w == 12 else 0.35
        b += rect(x0-wd/2, y-3.4, wd, 6.8, fill=col, stroke="none",
                  rx=2, opacity=op)
    b += txt(x0, y0+18, 'the cube, sliced into weight layers '
             '(n = 30, widths &#8733; C(n,w))', 10.5, INK, "middle")
    b += txt(x0+180, y0-12*9.5+4, 'layer w: C(n,w) words', 10.5, RED)
    b += line(x0+170, y0-12*9.5, x0+120, y0-12*9.5, RED, 1.3,
              marker=ar('fig_layers', 'R'))
    b += txt(x0, y0-31*9.5, 'weight 30', 9, MUT, "middle")
    bx = 545
    b += rect(bx, 90, 190, 150, fill=LT, stroke=INK, sw=1.1, rx=6)
    b += txt(bx+95, 115, 'Bassalygo&#8211;Elias:', 11, INK, "middle",
             weight="bold")
    b += txt(bx+95, 140, 'average over 2' + sup('', 'n', 10) +
             ' translates;', 10.5, INK, "middle")
    b += txt(bx+95, 160, 'some translate keeps a', 10.5, INK, "middle")
    b += txt(bx+95, 180, 'C(n,w)/2' + sup('', 'n', 10) +
             ' share of the code', 10.5, INK, "middle")
    b += txt(bx+95, 214, 'exponent toll: 1 &#8722; ' +
             sub('H', '2', 10) + '(&#945;)', 11, ACC, "middle")
    svg('fig_layers', 760, 360, b,
        'Weight layers of the cube and the averaging reduction')


# ================= 19 fig_supp =================
def fig_supp():
    b = ''
    n, w = 30, 12
    x0, y0, cell = 60, 70, 21
    for i in range(n):
        col = RED if i < w else BLUE
        b += rect(x0+i*cell, y0, cell-2, 26, fill=col, stroke="none",
                  rx=3, opacity=0.55 if i < w else 0.4)
    b += txt(x0+w*cell/2, y0-12, 'support of x: w = 12 coordinates',
             10.5, RED, "middle")
    b += txt(x0+w*cell+(n-w)*cell/2, y0-12,
             'complement: N = 18 coordinates', 10.5, BLUE, "middle")
    b += txt(x0+w*cell/2, y0+45, 'harmonics ' + sub('E', 'p', 10) +
             '(x), degree p = 3', 10.5, RED, "middle")
    b += txt(x0+w*cell+(n-w)*cell/2, y0+45, 'harmonics ' +
             sub('E', 'q', 10) + '(x' + sup('', 'c', 9) +
             '), degree q = 2', 10.5, BLUE, "middle")
    b += txt(x0+n*cell/2, y0+72, 'fiber ' + sub('E', 'x', 10) +
             sup('', 'p,q', 9) + ' = ' + sub('E', 'p', 10) +
             '(x) &#8855; ' + sub('E', 'q', 10) + '(x' +
             sup('', 'c', 9) + '),  dim = ' + sub('d', 'p', 10) +
             '(w)&#183;' + sub('d', 'q', 10) + '(N)', 11, GOLD,
             "middle")
    # branching window
    y1 = 210
    X = mapper(0, 12, 80, 680)
    b += line(70, y1, 692, y1, INK, 1.2)
    jm, jp = 5, 11
    b += rect(X(jm), y1-13, X(jp)-X(jm), 26, fill=GOLD, stroke="none",
              rx=4, opacity=0.35)
    for j in range(13):
        b += line(X(j), y1-4, X(j), y1+4, INK, 1.1)
        b += txt(X(j), y1+20, str(j), 9.5, MUT, "middle")
    b += txt(X(8), y1-22, 'Johnson degrees containing the fiber: ' +
             sub('j', '&#8722;', 9) + ' = p+q = 5  &#8804;  j  '
             '&#8804;  ' + sub('j', '+', 9) + ' = 11', 10.5, GOLD,
             "middle")
    b += txt(X(6), y1+40, 'window (31) for (w, N, p, q) = '
             '(12, 18, 3, 2):  ' + sub('j', '+', 9) +
             ' = min{12, 12&#8722;3+2, 18+3&#8722;2} = 11', 10.5, MUT,
             "middle")
    svg('fig_supp', 760, 270, b,
        'Support-complement split and the branching window')


# ================= 20 fig_hahn =================
def fig_hahn():
    b = ''
    xs = [90 + i*110 for i in range(6)]
    labels = [sub('j', '&#8722;', 9), 'j&#8722;1', 'j', 'j+1',
              '&#8943;', 'L']
    y = 150
    for i, (x, lb) in enumerate(zip(xs, labels)):
        if lb == '&#8943;':
            b += txt(x, y+5, '&#8943;', 16, MUT, "middle")
            continue
        b += circ(x, y, 24, LT, INK, 1.3)
        b += txt(x, y+4, lb, 11, INK, "middle", cls="v")
        b += poly([(x-12, y+30), (x, y+24), (x+12, y+30), (x, y+36)],
                  GOLD, 1.1, fill=GOLD, opacity=0.5, close=True)
        # self loop
        b += path(f'M{x-9:.0f},{y-22:.0f} C{x-16:.0f},{y-48:.0f} '
                  f'{x+16:.0f},{y-48:.0f} {x+9:.0f},{y-22:.0f}',
                  BLUE, 1.5, marker=ar('fig_hahn', 'B'))
        b += txt(x, y-52, sub('b&#770;', 'j' if i == 2 else '', 9),
                 10, BLUE, "middle", cls="v")
    for i in range(5):
        if labels[i] == '&#8943;' or labels[i+1] == '&#8943;':
            x1 = xs[i] + (26 if labels[i] != '&#8943;' else 14)
            x2 = xs[i+1] - (26 if labels[i+1] != '&#8943;' else 14)
        else:
            x1, x2 = xs[i]+26, xs[i+1]-26
        b += line(x1, y-5, x2, y-5, INK, 1.4)
        b += line(x2, y+5, x1, y+5, INK, 1.4)
    b += txt(xs[2]+55, y-14, sub('c&#770;', 'j', 9), 10.5, RED,
             "middle", cls="v")
    b += txt(340, 40, 'the layer graph: a path with self-loops '
             '(diagonal entries allowed)', 11.5, INK, "middle")
    b += txt(340, 230, 'square-over-classical weights (41)&#8211;(42):'
             '  b&#770; = (b' + sup('', 'p,q', 9) + ')&#178;/b&#8304;,'
             '   c&#770; = (c' + sup('', 'p,q', 9) +
             ')&#178;/c&#8304;', 11, MUT, "middle")
    b += txt(340, 252, 'fiber at every vertex: ' + sub('E', 'x', 10) +
             sup('', 'p,q', 9) + ' (tensor pair)', 10.5, GOLD,
             "middle")
    svg('fig_hahn', 680, 270, b,
        'The Johnson transition structure: path with loops')


# ================= 21 fig_m2gap =================
def fig_m2gap():
    b = ''
    try:
        import numpy as np
        from scipy.optimize import minimize
        have = True
    except Exception:
        have = False
    ds = [0.05, 0.08, 0.11, 0.14, 0.17, 0.20, 0.24, 0.28, 0.32, 0.36,
          0.40, 0.44]
    khs2 = smooth_lower(ds, [kappaH(d)[0] for d in ds])
    m1g = [max(M1(d) - v, 1e-7) for d, v in zip(ds, khs2)]
    m2g = []
    if have:
        def Lam(al, be, ga, u):
            z, m = 1-2*u, 1-2*al
            ze, xi = 1-2*be-2*ga, 1-2*al+2*be-2*ga
            B = (ze*xi-m*z*z)**2/(z*z*(1-m*m)*(1-z*z))
            C = ((z*z-xi*xi)*(ze*ze-z*z)
                 / (2*z*z*(1-m*m)*math.sqrt(1-z*z)))
            return B+2*C

        def obj(x):
            al, be, ga, u = x
            return (1-H2(al)+H2(u)-al*H2(be/al)
                    - (1-al)*H2(ga/(1-al)))
        rng = np.random.default_rng(5)
        for d in ds:
            m2v, _ = M2(d)
            best = min(m2v, kappaH(d)[0] + 0)  # kbin <= kH
            for _ in range(150 if d < 0.09 else 60):
                al = rng.uniform(d/2+1e-3, 0.499)
                mx = min(al, 1-al)
                be = rng.uniform(0, al/2*0.9)
                ga = rng.uniform(0, (1-al)/2*0.9)
                hi = min(al, al-be+ga, 1-al+be-ga)
                if hi - (be+ga) < 2e-3:
                    continue
                u = rng.uniform(be+ga+1e-3, hi-1e-3)

                def con(x):
                    al, be, ga, u = x
                    return (Lam(al, be, ga, u)
                            - (1 - d/(2*al*(1-al))) - 1e-10)
                cons = [{'type': 'ineq', 'fun': con},
                        {'type': 'ineq', 'fun': lambda x:
                         min(x[0], x[0]-x[1]+x[2], 1-x[0]+x[1]-x[2])
                         - x[3] - 1e-9},
                        {'type': 'ineq',
                         'fun': lambda x: x[3]-x[1]-x[2]-1e-9},
                        {'type': 'ineq',
                         'fun': lambda x: x[0]/2-x[1]},
                        {'type': 'ineq',
                         'fun': lambda x: (1-x[0])/2-x[2]}]
                try:
                    res = minimize(obj, [al, be, ga, u],
                                   method='SLSQP', constraints=cons,
                                   bounds=[(d/2+1e-9, 0.5-1e-9),
                                           (0, 0.25), (0, 0.5),
                                           (1e-6, 0.5)],
                                   options={'maxiter': 200,
                                            'ftol': 1e-13})
                    if res.success and con(res.x) > -1e-7:
                        best = min(best, res.fun)
                except Exception:
                    pass
            m2g.append(best)
        m2g = smooth_lower(ds, m2g)
        m2g = [max(M2(d)[0] - v, 1e-7) for d, v in zip(ds, m2g)]
    else:
        m2g = [max(g1*0.1, 1e-7) for g1 in m1g]
    x0, y0, w, h = 80, 250, 560, 195
    X = mapper(0.02, 0.47, x0, x0+w)
    Y = mapper(-6.3, -1.4, y0, y0-h)
    b += line(x0, y0, x0+w, y0, INK, 1.2)
    b += line(x0, y0, x0, y0-h, INK, 1.2)
    for e in (-2, -3, -4, -5, -6):
        b += line(x0-4, Y(e), x0, Y(e), INK, 1)
        b += txt(x0-7, Y(e)+3, '10' + sup('', str(e), 9), 9.5, MUT,
                 "end")
    for v in (0.1, 0.2, 0.3, 0.4):
        b += line(X(v), y0, X(v), y0+4, INK, 1)
        b += txt(X(v), y0+16, f'{v:.1f}', 9.5, MUT, "middle")
    b += txt(x0+w/2, y0+32, '&#948;', 13, INK, "middle", cls="v")
    b += poly([(X(d), Y(math.log10(v))) for d, v in zip(ds, m1g)],
              BLUE, 1.8)
    for d, v in zip(ds, m1g):
        b += circ(X(d), Y(math.log10(v)), 3.2, BLUE)
    b += poly([(X(d), Y(math.log10(v))) for d, v in zip(ds, m2g)],
              RED, 1.8)
    for d, v in zip(ds, m2g):
        b += circ(X(d), Y(math.log10(v)), 3.2, RED)
    b += txt(X(0.30), Y(-2.2), sub('M', '1', 10) + ' &#8722; &#954;' +
             sub('', 'H', 9), 11, BLUE)
    b += txt(X(0.13), Y(-4.9), sub('M', '2', 10) + ' &#8722; &#954;' +
             sub('', 'bin', 9), 11, RED)
    b += txt(x0+w/2, y0-h+6, 'computed record-breaking margins '
             '(numerical optimization; upper bounds on the true '
             'infima)', 10, MUT, "middle")
    svg('fig_m2gap', 700, 290, b,
        'Computed gaps between classical and new binary exponents')


# ================= 22 fig_repgraph =================
def fig_repgraph():
    b = ''
    b += rect(40, 40, 470, 250, fill="none", stroke=MUT, sw=1.2,
              rx=10, dash="6,4")
    b += txt(60, 62, '&#937; (retained)', 11, MUT)
    pos = {(0, 0): (120, 230), (1, 0): (240, 230), (2, 0): (360, 230),
           (0, 1): (180, 130), (1, 1): (300, 130), (0, 2): (240, 70)}
    edges = [((0, 0), (1, 0)), ((1, 0), (2, 0)), ((0, 0), (0, 1)),
             ((1, 0), (1, 1)), ((0, 1), (1, 1)), ((0, 1), (0, 2))]
    for e in edges:
        (x1, y1), (x2, y2) = pos[e[0]], pos[e[1]]
        b += line(x1, y1, x2, y2, INK, 1.5)
    (x1, y1) = pos[(2, 0)]
    b += line(x1+21, y1-9, x1+90, y1-40, MUT, 1.3, dash="4,3")
    b += txt(x1+96, y1-46, 'discarded edge', 9.5, MUT)
    for k, (x, y) in pos.items():
        b += circ(x, y, 22, LT, INK, 1.3)
        b += circ(x, y+8, 5, GOLD)
        b += txt(x, y-4, sub('V', '&#955;', 10), 11, INK, "middle",
                 cls="v")
    b += txt(186, 196, sub('p', '&#955;&#957;', 9), 10, RED, "middle",
             cls="v")
    b += txt(150, 262, 'gold dot: the one copy of the fiber E', 10,
             GOLD)
    b += rect(545, 60, 195, 190, fill=LY, stroke=ACC, sw=1.2, rx=8)
    b += txt(642, 88, 'requirements', 11, ACC, "middle",
             weight="bold")
    b += txt(556, 115, '&#8226; dim Hom(E, ' + sub('V', '&#955;', 10)
             + ') = 1', 10.5, INK)
    b += txt(556, 140, '&#8226; scalar transitions (52)', 10.5, INK)
    b += txt(556, 165, '&#8226; ' + sub('D', '&#955;', 10) +
             sub('p', '&#955;&#957;', 10) + ' = ' +
             sub('D', '&#957;', 10) + sub('p', '&#957;&#955;', 10) +
             '  (53)', 10.5, INK)
    b += txt(556, 190, '&#8226; connected, finite &#937;', 10.5, INK)
    b += txt(556, 222, '&#8658; bound (56) with &#923; = ' +
             '&#955;' + sub('', 'max', 8) + '(J)', 11, ACC)
    svg('fig_repgraph', 760, 300, b,
        'The general representation graph and its requirements')


# ================= 23 fig_zonal =================
def fig_zonal():
    b = ''

    def legendre(m, t):
        p0, p1 = 1.0, t
        if m == 0:
            return p0
        for j in range(1, m):
            p2 = ((2*j+1)*t*p1 - j*p0)/(j+1)
            p0, p1 = p1, p2
        return p1
    x0, y0, w, h = 55, 260, 330, 200
    X = mapper(-1, 1, x0, x0+w)
    Y = mapper(-1.05, 1.05, y0, y0-h)
    b += line(x0, Y(0), x0+w, Y(0), MUT, 1.0)
    b += line(x0, y0, x0, y0-h, INK, 1.2)
    b += line(x0, y0, x0+w, y0, INK, 1.2)
    for v in (-1, 0, 1):
        b += line(X(v), y0, X(v), y0+4, INK, 1)
        b += txt(X(v), y0+16, str(v), 9.5, MUT, "middle")
    b += txt(x0+w/2, y0+32, 't = &#9001;x, u&#9002;', 11, INK,
             "middle", cls="v")
    cols = {1: MUT, 2: GREEN, 3: BLUE, 6: RED}
    lx = x0 + 30
    for m, col in cols.items():
        pts = [(X(t), Y(legendre(m, t)))
               for t in [i/200-1 for i in range(401)]]
        b += poly(pts, col, 1.7 if m == 6 else 1.2)
        b += txt(lx, 32, f'm = {m}', 10, col)
        lx += 62
    # right: zonal bands on sphere
    cx, cy, R = 560, 160, 105
    b += circ(cx, cy, R, "none", INK, 1.4)
    xdir = (0.35, 0.25, 0.9)
    nn = math.sqrt(sum(v*v for v in xdir))
    xdir = tuple(v/nn for v in xdir)
    # orthonormal basis of xperp
    e1 = (-xdir[1], xdir[0], 0)
    n1 = math.sqrt(sum(v*v for v in e1))
    e1 = tuple(v/n1 for v in e1)
    e2 = (xdir[1]*e1[2]-xdir[2]*e1[1], xdir[2]*e1[0]-xdir[0]*e1[2],
          xdir[0]*e1[1]-xdir[1]*e1[0])
    for tt in [i/12-1+1/24 for i in range(24)]:
        val = legendre(6, tt)
        col = RED if val > 0 else BLUE
        rr = math.sqrt(max(0.0, 1-tt*tt))
        visible = []
        hidden = []
        for j in range(121):
            ph = 2*math.pi*j/120
            p3 = tuple(tt*xdir[i] + rr*(math.cos(ph)*e1[i]
                       + math.sin(ph)*e2[i]) for i in range(3))
            assert abs(sum(v*v for v in p3) - 1) < 1e-9
            px, py, dep = ortho(p3, 0.0, 0.0)
            pt = (cx + R*px, cy - R*py)
            (visible if dep >= 0 else hidden).append((pt, j))
        for grp, op in ((visible, 0.75), (hidden, 0.0)):
            runs = []
            cur = []
            last = None
            for pt, j in grp:
                if last is not None and j != last+1:
                    runs.append(cur)
                    cur = []
                cur.append(pt)
                last = j
            if cur:
                runs.append(cur)
            if op > 0:
                for run in runs:
                    if len(run) > 1:
                        b += poly(run, col, 3.2,
                                  opacity=op*min(1, abs(val)+0.25))
    px, py, dep = ortho(xdir, 0.0, 0.0)
    b += circ(cx+R*px, cy-R*py, 4, INK)
    nn2 = math.sqrt(px*px+py*py)
    ox2, oy2 = px/nn2, py/nn2
    b += line(cx+R*px+4*ox2, cy-R*py-4*oy2, cx+1.14*R*ox2,
              cy-1.14*R*oy2, INK, 0.9)
    b += txt(cx+1.24*R*ox2, cy-1.24*R*oy2+4, 'x', 12, INK, "middle",
             cls="v")
    b += txt(560, 292, sub('P', '6', 10) + '(&#9001;x,&#183;&#9002;) '
             'painted on the sphere: zonal bands', 10.5, INK,
             "middle")
    svg('fig_zonal', 760, 320, b,
        'Gegenbauer profiles and a zonal harmonic on the sphere')


# ================= 24 fig_tangentsphere =================
def fig_tangentsphere():
    b = ''
    cx, cy, R = 300, 175, 130
    az, el = 0.45, 0.32
    b += circ(cx, cy, R, "none", INK, 1.5)
    xdir = (0.15, -0.25, 0.96)
    nn = math.sqrt(sum(v*v for v in xdir))
    xdir = tuple(v/nn for v in xdir)
    e1 = (-xdir[1], xdir[0], 0)
    n1 = math.sqrt(sum(v*v for v in e1))
    e1 = tuple(v/n1 for v in e1)
    e2 = (xdir[1]*e1[2]-xdir[2]*e1[1], xdir[2]*e1[0]-xdir[0]*e1[2],
          xdir[0]*e1[1]-xdir[1]*e1[0])

    def pr(p3):
        px, py, dep = ortho(p3, az, el)
        return (cx + R*px, cy - R*py), dep
    # great circle S(x-perp)
    vis, hid = [], []
    for j in range(241):
        ph = 2*math.pi*j/240
        p3 = tuple(math.cos(ph)*e1[i] + math.sin(ph)*e2[i]
                   for i in range(3))
        assert abs(sum(v*v for v in p3)-1) < 1e-9
        assert abs(sum(p3[i]*xdir[i] for i in range(3))) < 1e-9
        pt, dep = pr(p3)
        (vis if dep >= 0 else hid).append((pt, j))
    for grp, dsh in ((hid, "4,4"), (vis, None)):
        runs, cur, last = [], [], None
        for pt, j in grp:
            if last is not None and j != last+1:
                runs.append(cur)
                cur = []
            cur.append(pt)
            last = j
        if cur:
            runs.append(cur)
        for run in runs:
            if len(run) > 1:
                b += poly(run, GOLD, 2.2 if not dsh else 1.4,
                          dash=dsh)
    xpt, _ = pr(xdir)
    b += line(cx, cy, xpt[0], xpt[1], INK, 1.6,
              marker=ar('fig_tangentsphere'))
    b += txt(xpt[0]+10, xpt[1]-4, 'x (code point)', 11.5, INK, cls="v")
    t = 0.60
    xi3 = tuple(math.cos(0.9)*e1[i] + math.sin(0.9)*e2[i]
                for i in range(3))
    u3 = tuple(t*xdir[i] + math.sqrt(1-t*t)*xi3[i] for i in range(3))
    assert abs(sum(v*v for v in u3)-1) < 1e-9
    upt, _ = pr(u3)
    b += line(cx, cy, upt[0], upt[1], BLUE, 1.6,
              marker=ar('fig_tangentsphere', 'B'))
    b += txt(upt[0]+8, upt[1]+12, 'u', 12, BLUE, cls="v")
    txpt, _ = pr(tuple(t*v for v in xdir))
    b += line(cx, cy, txpt[0], txpt[1], RED, 2.4)
    b += circ(txpt[0], txpt[1], 3, RED)
    b += txt(txpt[0]-8, txpt[1]-8, 't&#183;x', 10.5, RED, "end",
             cls="v")
    b += line(txpt[0], txpt[1], upt[0], upt[1], GREEN, 1.6, dash="5,3",
              marker=ar('fig_tangentsphere', 'G'))
    xipt, _ = pr(xi3)
    b += circ(xipt[0], xipt[1], 4, GOLD)
    b += txt(xipt[0]-12, xipt[1]+24, '&#958; &#8712; S(x' +
             sup('', '&#8869;', 9) + ')', 11, GOLD, "end", cls="v")
    b += line(xipt[0]-16, xipt[1]+18, xipt[0]-3, xipt[1]+5, GOLD,
              0.9)
    b += txt(300, 330, 'u = t&#183;x + &#8730;(1&#8722;t&#178;)&#183;'
             '&#958;:  latitude t toward x, direction &#958; on the '
             'tangent sphere (gold)', 11, INK, "middle")
    b += rect(490, 60, 245, 120, fill=LT, stroke=INK, sw=1.0, rx=6)
    b += txt(612, 84, 'the fiber at x:', 10.5, INK, "middle",
             weight="bold")
    b += txt(612, 106, sub('E', 'x', 11) + ' = ' + 'harmonics of '
             'degree k', 10.5, INK, "middle")
    b += txt(612, 126, 'on S(x' + sup('', '&#8869;', 9) + ')', 10.5,
             INK, "middle")
    b += txt(612, 152, 'eigenfunction shape (63):', 10, MUT, "middle")
    b += txt(612, 170, 'p(t)&#183;(1&#8722;t&#178;)' +
             sup('', 'k/2', 9) + '&#183;Y(&#958;)', 10.5, ACC,
             "middle")
    svg('fig_tangentsphere', 760, 350, b,
        'The tangent sphere and the latitude decomposition, computed '
        'orthographic projection')


# ================= 25 fig_rowpath =================
def fig_rowpath():
    b = ''
    y = 70
    xs = [70 + i*90 for i in range(4)]
    labs = [sub('V', 'k', 10) + sup('', 'S', 8),
            sub('V', 'k+1', 10) + sup('', 'S', 8), '&#8943;',
            sub('V', 'L', 10) + sup('', 'S', 8)]
    for x, lb in zip(xs, labs):
        if lb == '&#8943;':
            b += txt(x, y+5, '&#8943;', 15, MUT, "middle")
            continue
        b += circ(x, y, 22, LT, INK, 1.2)
        b += txt(x, y+4, lb, 10.5, INK, "middle", cls="v")
        b += circ(x, y+13, 4, GOLD)
    for i in range(3):
        x1 = xs[i] + (24 if labs[i] != '&#8943;' else 12)
        x2 = xs[i+1] - (24 if labs[i+1] != '&#8943;' else 12)
        b += line(x1, y, x2, y, INK, 1.4)
    b += txt(205, 25, 'fiber ' + sub('H', 'k', 10) + '(x' +
             sup('', '&#8869;', 8) + ') at every level', 10.5, GOLD,
             "middle")
    # weight profiles
    x0, y0, w, h = 470, 250, 250, 195
    n = 600
    X = mapper(0, 0.7, x0, x0+w)
    Y = mapper(0, 0.5, y0, y0-h)
    b += line(x0, y0, x0+w, y0, INK, 1.2)
    b += line(x0, y0, x0, y0-h, INK, 1.2)
    b += txt(x0+w/2, y0+30, 'i/n   (n = 600)', 10.5, INK, "middle")
    b += txt(x0-8, y0-h-8, 'edge weight (71)', 10, INK, "end")
    for v in (0.2, 0.4, 0.6):
        b += line(X(v), y0, X(v), y0+4, INK, 1)
        b += txt(X(v), y0+15, f'{v:.1f}', 9, MUT, "middle")
    for v in (0.2, 0.4):
        b += line(x0-4, Y(v), x0, Y(v), INK, 1)
        b += txt(x0-7, Y(v)+3, f'{v:.1f}', 9, MUT, "end")
    leg = [0]
    for bb, col in ((0.0, RED), (0.06, BLUE), (0.12, GREEN)):
        k = int(bb*n)
        pts = []
        for i in range(k, int(0.7*n)):
            c = ((i-k+1)*(i+k+n-2)
                 / math.sqrt((i+1)*(i+n-2)*(2*i+n-2)*(2*i+n)))
            pts.append((X(i/n), Y(c)))
        b += poly(pts, col, 1.6)
        lim = []
        for i in range(max(k, 8), int(0.7*n), 6):
            u = i/n
            gr = ((u-bb)*(1+u+bb)/((1+2*u)*math.sqrt(u*(1+u))))
            lim.append((X(u), Y(gr)))
        b += poly(lim, col, 1.0, dash="3,3", opacity=0.7)
        b += line(x0+w-95, y0-80+18*leg[0], x0+w-75, y0-80+18*leg[0],
                  col, 1.8)
        b += txt(x0+w-69, y0-76+18*leg[0], f'b = {bb:.2f}', 9.5, col)
        leg[0] += 1
    b += txt(x0+w/2, y0-h-24, 'raising the fiber degree b lowers the '
             'spectral profile&#8230;', 10.5, INK, "middle")
    b += txt(220, 250, '&#8230;but subtracts ' + sub('H', 'sph', 10) +
             '(b) from the exponent:', 11, INK, "middle")
    b += txt(220, 275, '&#934;' + sub('', 'row', 9) + '(a,b) = ' +
             sub('H', 'sph', 10) + '(a) &#8722; ' +
             sub('H', 'sph', 10) + '(b)', 12, ACC, "middle")
    b += txt(220, 300, 'net effect for small b &gt; 0: strictly '
             'better than the 1978 record', 10.5, MUT, "middle")
    svg('fig_rowpath', 760, 310, b,
        'The one-row spherical path and its computed weight profiles')


# ================= 26 fig_young =================
def fig_young():
    b = ''
    cell = 24
    lam = [7, 4, 2]
    mu = [5, 3]
    x0, y0 = 60, 60
    for r, ln in enumerate(lam):
        for c in range(ln):
            b += rect(x0+c*cell, y0+r*2*cell, cell-2, cell-2,
                      fill=LT, stroke=INK, sw=1.0, rx=2)
    for r, ln in enumerate(mu):
        for c in range(ln):
            b += rect(x0+c*cell+cell*0.35, y0+r*2*cell+cell,
                      cell-2, cell-2, fill=GOLD, stroke=ACC, sw=1.0,
                      rx=2, opacity=0.55)
    b += txt(x0+7.6*cell, y0+14, '&#955; = (7, 4, 2) ambient rows',
             11, INK)
    b += txt(x0+7.6*cell, y0+cell+14, '&#956; = (5, 3) stabilizer '
             'rows (gold)', 11, ACC)
    b += txt(x0+120, y0+6*cell+16, 'interlacing (72):  7 &#8805; 5 '
             '&#8805; 4 &#8805; 3 &#8805; 2 &#8805; 0', 11.5, INK)
    # right: one-box moves
    x1 = 480
    b += txt(x1+90, y0-15, 'coordinate multiplication:', 11, INK,
             "middle")
    b += txt(x1+90, y0+3, 'add or remove ONE box', 11, INK, "middle")
    lam2 = [4, 2]
    for r, ln in enumerate(lam2):
        for c in range(ln):
            b += rect(x1+c*cell, y0+30+r*cell, cell-2, cell-2,
                      fill=LT, stroke=INK, sw=1.0, rx=2)
    b += rect(x1+4*cell, y0+30, cell-2, cell-2, fill=LG,
              stroke=GREEN, sw=1.4, rx=2, dash="3,2")
    b += txt(x1+4*cell+30, y0+30+16, '+' + sub('e', '1', 9) +
             ' (5,2) &#10003;', 10, GREEN)
    b += rect(x1+2*cell, y0+30+cell, cell-2, cell-2, fill=LG,
              stroke=GREEN, sw=1.4, rx=2, dash="3,2")
    b += txt(x1+3*cell+8, y0+30+cell+16, '+' + sub('e', '2', 9) +
             ' (4,3) &#10003;', 10, GREEN)
    b += rect(x1, y0+30+2*cell, cell-2, cell-2, fill="#f7dcdc",
              stroke=RED, sw=1.4, rx=2, dash="3,2")
    b += txt(x1+cell+8, y0+30+2*cell+16, '+' + sub('e', '3', 9) +
             ': a third row breaks (72) for one-row &#956;', 10, RED)
    b += txt(x1+90, y0+30+3*cell+24, 'moves must keep rows weakly '
             'decreasing', 9.5, MUT, "middle")
    b += txt(x1+90, y0+30+3*cell+40, 'and keep interlacing &#956;',
             9.5, MUT, "middle")
    svg('fig_young', 760, 260, b,
        'Young diagrams, interlacing, and one-box moves')


# ================= 27 fig_tworow =================
def fig_tworow():
    b = ''
    n, k = 40, 4
    i0, i1 = k, k+8
    x0, y0 = 90, 280
    dx, dy = 72, 52

    def P(i, j):
        return (x0 + (i-i0)*dx, y0 - j*dy)

    def p79(i, j):
        pi_p = ((i-k+1)*(i+k+n-2)*(i+n-3)
                / ((i-j+1)*(i+j+n-3)*(2*i+n-2)))
        pi_m = ((i-k)*(i+k+n-3)*(i+1)
                / ((i-j+1)*(i+j+n-3)*(2*i+n-2)))
        pj_p = ((k-j)*(j+k+n-3)*(j+n-4)
                / ((i-j+1)*(i+j+n-3)*(2*j+n-4)))
        pj_m = (j*(k-j+1)*(j+k+n-4)
                / ((i-j+1)*(i+j+n-3)*(2*j+n-4)))
        return pi_p, pi_m, pj_p, pj_m
    for j in range(0, k+1):
        for i in range(i0, i1+1):
            if i < i1:
                w1 = math.sqrt(p79(i, j)[0]*p79(i+1, j)[1])
                b += line(*P(i, j), *P(i+1, j), BLUE if j == 0 else
                          INK, 0.8+9*w1, opacity=0.75)
            if j < k:
                w2 = math.sqrt(p79(i, j)[2]*p79(i, j+1)[3])
                b += line(*P(i, j), *P(i, j+1), GOLD, 0.8+9*w2,
                          opacity=0.85)
    for j in range(0, k+1):
        for i in range(i0, i1+1):
            b += circ(*P(i, j), 4.5, INK)
    b += rect(x0+6.1*dx, y0-4.4*dy, 2.3*dx, 1.9*dy, fill="none",
              stroke=MUT, sw=1.2, rx=6, dash="5,4")
    b += txt(x0+7.25*dx, y0-2.28*dy, 'asymptotically', 9.5, MUT,
             "middle")
    b += txt(x0+7.25*dx, y0-2.28*dy+13, 'constant corner', 9.5, MUT,
             "middle")
    b += txt(x0-14, y0+4, 'j = 0', 10.5, BLUE, "end")
    b += txt(x0-14, y0-k*dy+4, 'j = k', 10.5, INK, "end")
    b += txt(x0+4*dx, y0+30, 'i = k &#8230; (ambient first row); '
             'blue floor = the one-row path of &#167;11', 10.5, INK,
             "middle")
    b += txt(x0+4*dx, y0-k*dy-28, 'vertical (gold) edges exist only '
             'because k &gt; 0: the doorway weight ' +
             sub('p', 'j,+', 9) + '(i,0) = k(k+n&#8722;3)/'
             '((i+1)(i+n&#8722;3))', 10, GOLD, "middle")
    b += txt(x0+4*dx, y0+52, 'edge thickness = computed symmetric '
             'weight (77) for n = 40, k = 4', 9.5, MUT, "middle")
    svg('fig_tworow', 760, 360, b,
        'The two-row transition graph with computed edge weights')


# ================= 28 fig_box =================
def fig_box():
    b = ''
    # left: region + corner box
    x0, y0 = 60, 250
    b += poly([(x0, y0), (x0+240, y0), (x0+240, y0-80), (x0+120, y0-80),
               (x0+120, y0-160), (x0, y0-160)], INK, 1.4, fill=LT,
              close=True)
    b += rect(x0+175, y0-78, 62, 62, fill=GOLD, stroke=ACC, sw=1.3,
              rx=3, opacity=0.45)
    b += txt(x0+206, y0-90, 'box &#8459;', 11, ACC, "middle")
    b += txt(x0+120, y0-180, 'ambient weights &#937;' +
             sub('', 'n', 9), 11, INK, "middle")
    b += txt(x0+120, y0+22, 'corner ' + sub('&#955;', '*', 10) +
             ': all rows at their caps', 10.5, MUT, "middle")
    b += line(x0+237, y0-17, x0+237, y0-17, INK, 1)
    b += circ(x0+237, y0-78+62, 0.1, INK)
    b += circ(x0+237, y0-16, 4.5, RED)
    b += txt(x0+248, y0-12, sub('&#955;', '*', 9), 11, RED)
    # right: product-sine heatmap
    m = 25
    hx, hy, cs = 420, 60, 9.2
    for a in range(m):
        for c in range(m):
            v = math.sin(math.pi*(a+1)/(m+1))*math.sin(
                math.pi*(c+1)/(m+1))
            g = int(255 - v*180)
            col = f'#{g:02x}{g:02x}ff'
            b += (f'<rect x="{hx+a*cs:.1f}" y="{hy+(m-1-c)*cs:.1f}" '
                  f'width="{cs:.1f}" height="{cs:.1f}" fill="{col}" '
                  f'stroke="none"/>')
    b += rect(hx, hy, m*cs, m*cs, fill="none", stroke=INK, sw=1.2,
              rx=0)
    b += txt(hx+m*cs/2, hy+m*cs+20, 'product-of-sines ground state '
             'on the box (m = 25)', 10.5, INK, "middle")
    b += txt(hx+m*cs/2, hy+m*cs+40, 'Rayleigh quotient 2(' +
             sub('w', 'i', 9) + '+' + sub('w', 'j', 9) +
             ')cos(&#960;/(m+1)) &#8594; 2&#915;' + sub('', 'r', 9),
             10.5, MUT, "middle")
    svg('fig_box', 720, 340, b,
        'The corner box and its sine-product ground state')


# ================= 29 fig_escape =================
def fig_escape():
    b = ''
    x0, w = 70, 620
    y1 = 110
    X = mapper(0, 10, x0, x0+w)
    b += line(x0, y1, x0+w, y1, INK, 1.3)
    pts1 = [(0.4, 'x&#8323;', RED), (1.1, 'y&#8322;', BLUE),
            (1.9, 'x&#8322;', RED), (6.2, 'y&#8321;', BLUE),
            (8.0, 'x&#8321;', RED)]
    for v, lb, col in pts1:
        b += circ(X(v), y1, 5 if col == RED else 4, col)
        b += txt(X(v), y1-12, lb, 11, col, "middle", cls="v")
    b += line(X(6.2), y1+16, X(9.4), y1+16, BLUE, 1.4,
              marker=ar('fig_escape', 'B'))
    b += line(X(8.0), y1+30, X(9.7), y1+30, RED, 1.4,
              marker=ar('fig_escape', 'R'))
    b += txt(X(8.6), y1+46, 'pair diverges, ratio ' +
             sub('y', '1', 9) + '/' + sub('x', '1', 9) +
             ' &#8594; &#945;&#8321;', 10.5, INK, "middle")
    y2 = 230
    b += line(x0, y2, x0+w, y2, INK, 1.3)
    for v, lb, col in pts1[:3]:
        b += circ(X(v), y2, 5 if col == RED else 4, col)
        b += txt(X(v), y2-12, lb, 11, col, "middle", cls="v")
    b += circ(X(9.6), y2, 9, "none", ACC, 1.6)
    b += txt(X(9.6), y2-16, 'atom at &#8734;, mass 1&#8722;c&#178;',
             10.5, ACC, "middle")
    b += txt(X(9.6), y2+22, '= cap cost &#189;' + 'log&#8322;'
             '((1&#8722;t)/(1&#8722;s))', 10, ACC, "middle")
    b += txt(x0+w/2, 290, 'limit: a lower-level tuple plus an atom at '
             'infinity &#8212; where the weight &#8730;(x/(x+&#188;)) '
             'saturates at 1', 10.5, MUT, "middle")
    b += txt(x0+w/2, 40, 'a minimizing sequence at level 3', 11, INK,
             "middle")
    svg('fig_escape', 760, 305, b,
        'Escape to infinity equals a spherical cap (Lemma 7.5)')


# ================= 30 fig_hemisphere =================
def fig_hemisphere():
    b = ''
    cx, cy, R = 330, 200, 120
    b += circ(cx, cy, R, "none", INK, 1.5)
    b += line(cx-R-24, cy, cx+R+24, cy, MUT, 1.0, dash="4,3")
    ty = cy - R - 46
    b += line(cx-260, ty, cx+260, ty, INK, 1.3)
    b += txt(cx-270, ty+4, 'packing in ' + '&#8477;' +
             sup('', 'n', 9), 11, INK, "end")
    centers = [-215, -160, -105, -45, 25, 95, 160, 220]
    for c in centers:
        b += circ(cx+c, ty, 3.5, RED)
        b += line(cx+c-24, ty-7, cx+c+24, ty-7, RED, 2.2, opacity=0.4)
    for c in centers:
        px, py = cx+c, ty
        vx, vy = px-cx, py-cy
        nn = math.sqrt(vx*vx+vy*vy)
        sx, sy = cx+R*vx/nn, cy+R*vy/nn
        b += line(px, py, sx, sy, MUT, 0.9, dash="2,3")
        b += circ(sx, sy, 4, BLUE)
    b += txt(cx, cy+34, 'radial lift', 10.5, MUT, "middle")
    b += txt(cx+R+30, cy-46, 'spherical code on the upper hemisphere,',
             10.5, BLUE)
    b += txt(cx+R+30, cy-30, 'inner products &#8804; s', 10.5, BLUE)
    b += txt(cx, cy+R+34, '&#916;' + sub('', 'n', 9) +
             ' &#8804; ((1&#8722;s)/2)' + sup('', 'n/2', 9) +
             ' &#183; A(n+1, s)   (89)', 12, ACC, "middle")
    svg('fig_hemisphere', 760, 370, b,
        'The Sidelnikov upper-hemisphere lift')


# ================= 31 fig_measures =================
def fig_measures():
    b = ''
    av, bv = [0.6, 0.13], [0.31]
    xs = [a*(1+a) for a in av]
    ys = [v*(1+v) for v in bv]
    c = 0.25

    def K(t):
        val = math.log((t+xs[-1])/t)
        for xm, ym in zip(xs[:-1], ys):
            val += math.log((t+xm)/(t+ym))
        return val

    def GammaPhi():
        x = xs
        y = ys
        G = 0.0
        for li in range(2):
            num = 1.0
            for ym in y:
                num *= (x[li]-ym)
            den = 1.0
            for j2 in range(2):
                if j2 != li:
                    den *= (x[li]-x[j2])
            G += num/den*qfun(av[li])
        Phi = sum(Hsph(a) for a in av) - sum(Hsph(v) for v in bv)
        return G, Phi
    G, Phi = GammaPhi()
    Z = 1 - 2*G
    lamstar = 0.5*math.log2(2*math.pi/math.e)
    gap = lamstar - (0.5*math.log2(2/Z) - Phi)
    # left: indicator + K
    x0, y0, w, h = 60, 150, 280, 90
    X = mapper(0, 1.1, x0, x0+w)
    b += line(x0, y0, x0+w, y0, INK, 1.2)
    b += rect(X(0), y0-18, X(xs[1])-X(0), 18, fill=GOLD, stroke="none",
              rx=0, opacity=0.5)
    b += rect(X(ys[0]), y0-18, X(xs[0])-X(ys[0]), 18, fill=GOLD,
              stroke="none", rx=0, opacity=0.5)
    for v, lb in ((xs[1], sub('x', '2', 9)), (ys[0], sub('y', '1', 9)),
                  (xs[0], sub('x', '1', 9))):
        b += line(X(v), y0-3, X(v), y0+3, INK, 1.2)
        b += txt(X(v), y0+16, lb, 10, INK, "middle", cls="v")
    b += txt(x0+w/2, y0-30, 'interlacing indicator &#958; '
             '(quadratic coordinates)', 10.5, INK, "middle")
    y0b, hb = 300, 90
    Xk = mapper(0.002, 0.25, x0, x0+w)
    Yk = mapper(0, K(0.002), y0b, y0b-hb)
    b += line(x0, y0b, x0+w, y0b, INK, 1.2)
    b += line(x0, y0b, x0, y0b-hb, INK, 1.2)
    pts = [(Xk(t), Yk(K(t)))
           for t in [0.002 + i*(0.248/200) for i in range(201)]]
    b += poly(pts, ACC, 1.8)
    b += txt(x0+w/2, y0b+18, 't &#8712; (0, &#188;)', 10.5, MUT,
             "middle")
    b += txt(x0+70, y0b-hb+12, 'potential K(t)  (96)', 10.5, ACC)
    # right: three densities vs v-coordinate t=c(1-v^2)
    x1, y1, w1, h1 = 420, 300, 300, 240

    def rho_d(t):
        return 1/(math.pi*math.sqrt(t*(c-t)))

    def nu_d(t):
        return 1/math.sqrt(c-t)

    def rhoK_d(t):
        return math.exp(-K(t))/Z * rho_d(t)
    X2 = mapper(0, 0.25, x1, x1+w1)
    Y2 = mapper(0, 9.0, y1, y1-h1)
    b += line(x1, y1, x1+w1, y1, INK, 1.2)
    b += line(x1, y1, x1, y1-h1, INK, 1.2)
    ts = [0.0008 + i*(0.2482/400) for i in range(401)]
    nu_pts = [(X2(t), Y2(min(9, nu_d(t)))) for t in ts]
    rk_pts = [(X2(t), Y2(min(9, rhoK_d(t)))) for t in ts]
    shade = nu_pts + rk_pts[::-1]
    b += poly(shade, "none", 0, fill=RED, close=True, opacity=0.12)
    b += poly([(X2(t), Y2(min(9, rho_d(t)))) for t in ts], MUT, 1.2,
              dash="4,3")
    b += poly(nu_pts, BLUE, 2.0)
    b += poly(rk_pts, RED, 2.0)
    b += txt(X2(0.05), Y2(rho_d(0.05))-8, '&#961; (arcsine)', 10, MUT)
    b += txt(X2(0.135), Y2(4.35), '&#957; (target)', 10.5, BLUE)
    b += txt(X2(0.105), Y2(rhoK_d(0.105))+16, '&#961;' +
             sub('', 'K', 9) + ' (tilt)', 10.5, RED)
    b += txt(x1+w1/2, y1+20, 't', 11, INK, "middle", cls="v")
    b += txt(x1+w1/2, y1-h1-16, 'gap to &#955;* = ' +
             'D(&#957;&#8214;&#961;' + sub('', 'K', 9) +
             ')/(2 log 2) = %.4f bits' % gap, 11, ACC, "middle")
    svg('fig_measures', 760, 350, b,
        'The interlacing indicator, its potential, and the three '
        'measures of Proposition 8.1, computed')
    return gap


# ================= 32 fig_cheby =================
def fig_cheby():
    b = ''
    N, R = 8, 6.0
    xs = [R/2*(1+math.cos((2*li-1)*math.pi/(2*N)))
          for li in range(1, N+1)]
    ys = [R/2*(1+math.cos(m*math.pi/N)) for m in range(1, N)]
    x0, w = 70, 620
    X = mapper(0, R, x0, x0+w)
    y1 = 70
    b += line(x0, y1, x0+w, y1, INK, 1.3)
    for v in xs:
        b += circ(X(v), y1, 5, RED)
    for v in ys:
        b += line(X(v), y1-8, X(v), y1+8, BLUE, 2)
    b += txt(x0+w/2, y1-22, 'roots ' + sub('x', '&#8467;', 9) +
             ' (red) and critical points ' + sub('y', 'm', 9) +
             ' (blue) of the shifted Chebyshev polynomial, N = 8, '
             'R = 6', 10.5, INK, "middle")
    b += txt(x0+w+8, y1+4, 'all ' + sub('R', '&#8467;', 9) +
             ' = 1/N', 10.5, ACC)
    # indicator
    y2 = 150
    b += line(x0, y2, x0+w, y2, INK, 1.2)
    segs = [(0, xs[-1])] + [(ys[m], xs[m]) for m in range(N-2, -1, -1)]
    for a, bb2 in segs:
        b += rect(X(a), y2-14, X(bb2)-X(a), 14, fill=GOLD,
                  stroke="none", rx=0, opacity=0.55)
    b += txt(x0+w/2, y2+18, 'the indicator &#958;' + sub('', 'N', 9) +
             ' fills half of each equal &#952;-interval &#8594; '
             'weak-* limit &#189;', 10.5, GOLD, "middle")
    # convergence chart
    x0c, y0c, wc, hc = 110, 420, 520, 190
    lamstar = 0.5*math.log2(2*math.pi/math.e)
    Xc = mapper(math.log10(0.1), math.log10(100), x0c, x0c+wc)
    Yc = mapper(0.50, 0.615, y0c, y0c-hc)
    b += line(x0c, y0c, x0c+wc, y0c, INK, 1.2)
    b += line(x0c, y0c, x0c, y0c-hc, INK, 1.2)
    for v in (0.1, 1, 10, 100):
        b += line(Xc(math.log10(v)), y0c, Xc(math.log10(v)), y0c+4,
                  INK, 1)
        b += txt(Xc(math.log10(v)), y0c+16, str(v), 9.5, MUT, "middle")
    for v in (0.52, 0.56, 0.60):
        b += line(x0c-4, Yc(v), x0c, Yc(v), INK, 1)
        b += txt(x0c-7, Yc(v)+3, f'{v:.2f}', 9.5, MUT, "end")
    b += txt(x0c+wc/2, y0c+32, 'R  (log scale)', 11, INK, "middle")
    b += line(x0c, Yc(lamstar), x0c+wc, Yc(lamstar), RED, 1.3,
              dash="6,3")
    b += txt(x0c+wc, Yc(lamstar)-8, '&#955;* = 0.60440', 10.5, RED,
             "end")
    pts = []
    for i in range(100):
        Rv = 10**(-1 + i*0.03)
        ZR = 2/math.pi*math.asin(1/math.sqrt(4*Rv+1))
        aR = (-1+math.sqrt(1+4*Rv))/2
        hR = math.log(2)*Hsph(aR)
        gv = (math.log(2/ZR) - hR)/(2*math.log(2))
        if gv > 0.50:
            pts.append((Xc(math.log10(Rv)), Yc(gv)))
    b += poly(pts, BLUE, 1.9)
    b += txt(x0c+180, Yc(0.575), 'N &#8594; &#8734; packing objective '
             '(107) at cutoff R', 10.5, BLUE)
    for NN, col in ((4, MUT), (8, GREEN)):
        Rv = 6.0
        xs2 = [Rv/2*(1+math.cos((2*li-1)*math.pi/(2*NN)))
               for li in range(1, NN+1)]
        ys2 = [Rv/2*(1+math.cos(m*math.pi/NN)) for m in range(1, NN)]
        ZN = 1 - sum(math.sqrt(x/(x+0.25)) for x in xs2)/NN
        Phi = (sum(Hsph((-1+math.sqrt(1+4*x))/2) for x in xs2)
               - sum(Hsph((-1+math.sqrt(1+4*y))/2) for y in ys2))
        gv = 0.5*math.log2(2/ZN) - Phi
        b += circ(Xc(math.log10(Rv)), Yc(gv), 4, col)
        b += txt(Xc(math.log10(Rv))+10, Yc(gv)+4, f'N = {NN}', 9.5,
                 col)
    svg('fig_cheby', 760, 465, b,
        'Chebyshev interlacing tuples and convergence of the packing '
        'objective to lambda-star')


# ================= 33 fig_sech =================
def fig_sech():
    b = ''
    x0, y0, w, h = 70, 270, 560, 210
    X = mapper(0, 3.0, x0, x0+w)
    Y = mapper(0, 1.05, y0, y0-h)
    b += line(x0, y0, x0+w, y0, INK, 1.2)
    b += line(x0, y0, x0, y0-h, INK, 1.2)
    for v in (1, 2, 3):
        b += line(X(v), y0, X(v), y0+4, INK, 1)
        b += txt(X(v), y0+16, str(v), 9.5, MUT, "middle")
    b += txt(x0+w/2, y0+32, 'hyperbolic coordinate a  (t = &#188;'
             'sech&#178;a)', 11, INK, "middle")
    aa = [i*0.015 for i in range(201)]
    sech2 = [(X(a), Y(1/math.cosh(a)**2)) for a in aa]
    b += poly(sech2, ACC, 2.6)
    b += txt(X(0.8), Y(1/math.cosh(0.8)**2)-12, 'sech&#178;a  '
             '(= Chapter 1 Mellin measure)', 11, ACC)
    leg2 = 0
    for T, col in ((0.1, MUT), (1.0, BLUE), (10.0, GREEN)):
        ZT = 2/math.pi*math.asin(1/math.sqrt(4*T+1))
        pts = [(X(a), Y(2/(math.pi*ZT)/math.cosh(a)
                        / math.sqrt(1+4*T*math.cosh(a)**2)))
               for a in aa]
        b += poly(pts, col, 1.4)
        b += line(X(2.15), Y(0.86-0.09*leg2), X(2.32),
                  Y(0.86-0.09*leg2), col, 1.6)
        b += txt(X(2.37), Y(0.855-0.09*leg2), 'T = %g' % T, 9.5, col)
        leg2 += 1
    b += txt(x0+w/2, 30, 'tilted spherical measures &#961;' +
             sub('', 'K', 9) + sub('', 'T', 8) +
             ' converge to the Euclidean limit measure', 11.5, INK,
             "middle")
    svg('fig_sech', 700, 310, b,
        'The hyperbolic-coordinate dictionary between the two '
        'chapters')


# ================= 34 fig_landscape =================
def fig_landscape():
    b = ''
    # binary track
    x0, y0, w = 80, 130, 600
    X = mapper(0.25, 0.50, x0, x0+w)
    b += line(x0, y0, x0+w, y0, INK, 1.3)
    gv, kb, m2v = 1-H2(0.2), 0.4597, 0.4614
    for v, lb, col, dy, anc in (
            (gv, 'GV 0.2781', GREEN, -10, "middle"),
            (kb, '&#954;bin 0.4597', ACC, -24, "end"),
            (m2v, sub('M', '2', 9) + ' 0.4614', RED, -10, "start")):
        b += line(X(v), y0-6, X(v), y0+6, col, 2.2)
        xoff = -6 if anc == "end" else (6 if anc == "start" else 0)
        b += txt(X(v)+xoff, y0+dy, lb, 10, col, anc)
    b += rect(X(gv), y0-4, X(0.36)-X(gv), 8, fill=MUT, stroke="none",
              rx=2, opacity=0.25)
    b += txt(X(0.315), y0+22, 'LP-optimum floors (Samorodnitsky '
             'et al.): no two-point certificate descends here', 9.5,
             MUT, "middle")
    b += txt(x0+w/2, y0-46, 'binary, at &#948; = 0.2: the unknown '
             'interval and the moved ceiling', 11, INK, "middle")
    b += line(X(m2v), y0+14, X(kb), y0+14, ACC, 1.6,
              marker=ar('fig_landscape'))
    # spherical track
    y1 = 280
    X2 = mapper(0.3952, 0.4022, x0, x0+w)
    b += line(x0, y1, x0+w, y1, INK, 1.3)
    vals = [(0.4009442, sub('&#954;&#773;', '0', 9) + ' = ' +
             sub('B', 'KL', 9), RED),
            (0.39731, sub('&#954;&#773;', 'row', 9), GOLD),
            (0.39674, sub('&#954;&#773;', '1', 9), BLUE),
            (0.39660, sub('&#954;&#773;', '2', 9), GREEN)]
    for i, (v, lb, col) in enumerate(vals):
        b += line(X2(v), y1-6, X2(v), y1+6, col, 2.2)
        b += txt(X2(v), y1-12 - (13 if i % 2 else 0), lb, 10, col,
                 "middle")
    b += txt(X2(0.3987), y1+24, '&#8230;the ladder continues; its '
             's&#8594;1 endpoint is pinned at &#955;* = the exact LP '
             'optimum (Ch. 1)', 9.5, MUT, "middle")
    b += txt(x0+w/2, y1-46, 'spherical, at s = &#189; (kissing '
             'exponent): the strict hierarchy', 11, INK, "middle")
    b += line(X2(0.4009442), y1+42, X2(0.3968), y1+42, ACC, 1.6,
              marker=ar('fig_landscape'))
    # other regimes boxes
    y2 = 330
    boxes = [('fixed small n: SDP hierarchies', 90),
             ('d = n/2 &#8722; &#920;(&#8730;n): [PMP23]', 300),
             ('dims 8, 24: LP exactly tight', 500)]
    for lb, xx in boxes:
        b += rect(xx, y2, 180, 32, fill=LT, stroke=MUT, sw=1.0, rx=6)
        b += txt(xx+90, y2+20, lb, 9.5, INK, "middle")
    svg('fig_landscape', 760, 385, b,
        'The landscape after this chapter with computed anchors')


def main():
    fig_layers()
    fig_supp()
    fig_hahn()
    fig_m2gap()
    fig_repgraph()
    fig_zonal()
    fig_tangentsphere()
    fig_rowpath()
    fig_young()
    fig_tworow()
    fig_box()
    fig_escape()
    fig_hemisphere()
    gap = fig_measures()
    print("fig_measures DKL gap =", gap)
    fig_cheby()
    fig_sech()
    fig_landscape()


if __name__ == "__main__":
    main()
