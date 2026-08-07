"""Generate all SVG figures for chapter 2 into .ignore/ch02_figs/.

Figure ledger (order of appearance = CSS counter number):
 1 fig_twoworlds   binary cube + spherical code
 2 fig_timeline    1948-2025 records timeline
 3 fig_rates       GV / M1 / M2 landscape + F_delta objective
 4 fig_sphrates    BKL curve + kissing ladder
 5 fig_roadmap     dependency DAG
 6 fig_delsarte    certificate anatomy
 7 fig_toycert     the n=3 d=2 tight certificate
 8 fig_vectorsub   vector-per-point vs subspace-per-point
 9 fig_overlap     two planes, principal angles (3D projected)
10 fig_residual    residual Pythagoras of Prop 4.1
11 fig_levels      Boolean Fourier ladder
12 fig_harmonic    computed harmonic staircase n=4
13 fig_path        whole-cube path (paper Fig 1)
14 fig_sine        edge weights + sine test vector
15 fig_lammax      lambda_max convergence (computed)
16 fig_feasible    (a,b) feasibility landscape delta=0.2
17 fig_m1kh        M1 vs kappa_H + log gap
18 fig_layers      Bassalygo-Elias layers
19 fig_supp        support/complement + branching window
20 fig_hahn        Johnson path with loops
21 fig_m2gap       computed gaps M1-kH, M2-kbin
22 fig_repgraph    general representation graph
23 fig_zonal       Gegenbauer profiles + zonal sphere
24 fig_tangentsphere  tangent sphere decomposition (3D)
25 fig_rowpath     one-row path + weight profiles
26 fig_young       Young diagrams interlacing
27 fig_tworow      two-row lattice with computed weights
28 fig_box         corner box + product-sine ground state
29 fig_escape      Lemma 7.5 escape-to-infinity
30 fig_hemisphere  Sidelnikov lift (computed)
31 fig_measures    rho, nu, rho_K + DKL gap (computed)
32 fig_cheby       Chebyshev nodes + convergence
33 fig_sech        hyperbolic-coordinate dictionary
34 fig_landscape   the after-map
"""
import math

from ch02_figlib import (OUT, RED, BLUE, GREEN, GOLD, INK, MUT, LT, LG,
                         LY, ACC, svg, ar, txt, sub, sup, line, poly,
                         circ, rect, path, mapper, ortho)


def H2(u):
    if u <= 0 or u >= 1:
        return 0.0
    return -u*math.log2(u) - (1-u)*math.log2(1-u)


def gfun(v):
    return H2((1 - math.sqrt(max(0.0, 1-v)))/2)


def M1(d):
    return H2(0.5 - math.sqrt(d*(1-d)))


def Fdel(t, d):
    return 1 + gfun(t*t) - gfun(t*t + 2*d*t + 2*d)


def M2(d):
    lo, hi = 0.0, 1-2*d
    for _ in range(120):
        m1 = lo + (hi-lo)/3
        m2 = hi - (hi-lo)/3
        if Fdel(m1, d) < Fdel(m2, d):
            hi = m2
        else:
            lo = m1
    t = 0.5*(lo+hi)
    return min(Fdel(t, d), Fdel(0, d), Fdel(1-2*d, d)), t


def Hsph(u):
    if u <= 0:
        return 0.0
    return (1+u)*math.log2(1+u) - u*math.log2(u)


def qfun(u):
    return math.sqrt(u*(1+u))/(1+2*u) if u > 0 else 0.0


def a0fun(s):
    return 0.5*((1-s*s)**-0.5 - 1)


def BKL(s):
    lo, hi = 1e-9, s
    for _ in range(200):
        m1 = lo + (hi-lo)/3
        m2 = hi - (hi-lo)/3
        v1 = Hsph(a0fun(m1)) + 0.5*math.log2((1-m1)/(1-s))
        v2 = Hsph(a0fun(m2)) + 0.5*math.log2((1-m2)/(1-s))
        if v1 < v2:
            hi = m2
        else:
            lo = m1
    t = 0.5*(lo+hi)
    return Hsph(a0fun(t)) + 0.5*math.log2((1-t)/(1-s))


def GammaH(a, b):
    return 2*(a-b)*(1-a-b)/math.sqrt(a*(1-a))


def kappaH(d):
    """inf H2(a)-H2(b) s.t. GammaH >= 1-2d; grid + refine."""
    s = 1 - 2*d
    best = (1e9, 0, 0)
    for ia in range(2, 100):
        a = ia/200.0
        for ib in range(0, ia):
            b = ib/200.0
            if b < a and GammaH(a, b) > s:
                v = H2(a) - H2(b)
                if v < best[0]:
                    best = (v, a, b)
    # extra seeds on the improving ray from the classical point
    a00 = 0.5 - math.sqrt(d*(1-d))
    cray = 2/(1-2*a00) + 0.5
    for bb in (1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 0.01, 0.03):
        aa = a00 + cray*bb
        if 0 < aa <= 0.5 and bb < aa and GammaH(aa, bb) > s:
            v = H2(aa) - H2(bb)
            if v < best[0]:
                best = (v, aa, bb)
    a, b = best[1], best[2]
    step = max(1/200.0, b*0.5 if b > 0 else 1/200.0)
    for _ in range(70):
        step *= 0.7
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                a2, b2 = a + da*step, b + db*step
                if 0 < a2 <= 0.5 and 0 <= b2 < a2 and \
                        GammaH(a2, b2) > s:
                    v = H2(a2) - H2(b2)
                    if v < best[0]:
                        best = (v, a2, b2)
        a, b = best[1], best[2]
    try:
        from scipy.optimize import minimize
        for seed in [(best[1], best[2]), (a00+2e-3, 1e-3),
                     (a00+0.03, 0.012)]:
            res = minimize(
                lambda x: H2(x[0]) - H2(x[1]), seed, method='SLSQP',
                constraints=[{'type': 'ineq', 'fun':
                              lambda x: GammaH(x[0], x[1]) - s - 1e-10},
                             {'type': 'ineq',
                              'fun': lambda x: x[0] - x[1] - 1e-9}],
                bounds=[(1e-9, 0.5), (0.0, 0.5)],
                options={'maxiter': 400, 'ftol': 1e-15})
            if res.success and res.fun < best[0] and \
                    GammaH(*res.x) >= s - 1e-9:
                best = (float(res.fun), float(res.x[0]),
                        float(res.x[1]))
    except Exception:
        pass
    return best


def lam_tri(diag, off):
    m = len(diag)

    def count_below(x):
        cnt, dd = 0, 1.0
        for i in range(m):
            dd = diag[i] - x - (off[i-1]**2/dd if i > 0 else 0.0)
            if dd == 0.0:
                dd = -1e-300
            if dd < 0:
                cnt += 1
        return cnt
    mx = max((abs(x) for x in off), default=0.0)
    lo = min(diag) - 2*mx - 1
    hi = max(diag) + 2*mx + 1
    for _ in range(120):
        mid = 0.5*(lo+hi)
        if count_below(mid) == m:
            hi = mid
        else:
            lo = mid
    return 0.5*(lo+hi)


# ================= 1 fig_twoworlds =================
def fig_twoworlds():
    b = ''
    # left: 3D cube projected
    cx, cy, sc = 180, 158, 62
    verts = {}
    for v in range(8):
        p = ((v & 1), (v >> 1) & 1, (v >> 2) & 1)
        x, y, dep = ortho((p[0]*2-1, p[1]*2-1, p[2]*2-1), 0.5, 0.42)
        verts[p] = (cx + sc*x, cy + sc*y)
    edges = []
    for v in verts:
        for i in range(3):
            w = list(v)
            w[i] ^= 1
            w = tuple(w)
            if w > v:
                edges.append((v, w))
    for v, w in edges:
        b += line(*verts[v], *verts[w], MUT, 1.0, opacity=0.7)
    code = [(0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1)]
    for v, w in ((code[0], code[1]),):
        pass
    for p in verts:
        if p in code:
            b += circ(*verts[p], 6, RED)
        else:
            b += circ(*verts[p], 3.5, MUT)
    lab = {(0, 0, 0): '000', (1, 1, 0): '011', (1, 0, 1): '101',
           (0, 1, 1): '110'}
    offs = {(0, 0, 0): (-4, 18), (1, 1, 0): (14, 6), (1, 0, 1): (0, -10),
            (0, 1, 1): (-26, 2)}
    for p, t in lab.items():
        b += txt(verts[p][0]+offs[p][0], verts[p][1]+offs[p][1], t, 11,
                 RED, "middle", weight="bold")
    b += txt(180, 292, 'binary: vertices of the cube, pairwise far in '
             + 'Hamming distance', 11.5, INK, "middle")
    # right: circle with caps
    cx2, cy2, r = 560, 155, 92
    b += circ(cx2, cy2, r, "none", INK, 1.4)
    angs = [0.35, 1.35, 2.45, 3.5, 4.5, 5.45]
    for aa in angs:
        x = cx2 + r*math.cos(aa)
        y = cy2 + r*math.sin(aa)
        b += circ(x, y, 5, BLUE)
        a1, a2 = aa-0.38, aa+0.38
        p1 = (cx2 + r*math.cos(a1), cy2 + r*math.sin(a1))
        p2 = (cx2 + r*math.cos(a2), cy2 + r*math.sin(a2))
        b += path(f'M{p1[0]:.1f},{p1[1]:.1f} A{r},{r} 0 0 1 '
                  f'{p2[0]:.1f},{p2[1]:.1f}', BLUE, 5, opacity=0.35)
    b += line(cx2, cy2, cx2 + r*math.cos(angs[0]),
              cy2 + r*math.sin(angs[0]), MUT, 1.0)
    b += line(cx2, cy2, cx2 + r*math.cos(angs[1]),
              cy2 + r*math.sin(angs[1]), MUT, 1.0)
    mid = 0.5*(angs[0]+angs[1])
    b += path(f'M{cx2+30*math.cos(angs[0]):.1f},'
              f'{cy2+30*math.sin(angs[0]):.1f} A30,30 0 0 1 '
              f'{cx2+30*math.cos(angs[1]):.1f},'
              f'{cy2+30*math.sin(angs[1]):.1f}', RED, 1.4)
    b += txt(cx2+48*math.cos(mid), cy2+48*math.sin(mid)+4,
             '&#8805; arccos&#8201;s', 11, RED, "middle")
    b += txt(560, 292, 'spherical: unit directions, pairwise inner '
             + 'product &#8804; s', 11.5, INK, "middle")
    svg('fig_twoworlds', 760, 305, b,
        'Binary code on the cube and spherical code with disjoint caps')


# ================= 2 fig_timeline =================
def fig_timeline():
    b = ''
    X = mapper(1945, 2028, 60, 720)
    for y0, laby in ((95, 'binary'), (185, 'spherical')):
        b += line(55, y0, 725, y0, MUT, 1.2)
        b += txt(50, y0+4, laby, 12, INK, "end", cls="v")
    for yr in (1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020):
        b += line(X(yr), 230, X(yr), 236, INK, 1)
        b += txt(X(yr), 250, str(yr), 10, MUT, "middle")
    ev = [
        (1948, 95, 'Shannon', MUT, -26),
        (1972, 95, 'Delsarte LP', BLUE, -26),
        (1977, 95, 'MRRW record', RED, -44),
        (2005, 95, 'Schrijver SDP (fixed n)', MUT, -26),
        (2025, 95, 'this chapter', ACC, -44),
        (1972, 185, 'Delsarte LP', BLUE, -26),
        (1978, 185, 'KL record', RED, -44),
        (2008, 185, 'BV SDP (fixed n)', MUT, -26),
        (2017, 185, 'Viazovska d=8', MUT, -62),
        (2025, 185, 'this chapter', ACC, -44),
    ]
    for yr, y0, lab, col, dy in ev:
        b += line(X(yr), y0, X(yr), y0+dy+12, col, 1.1)
        b += circ(X(yr), y0, 4, col)
        b += txt(X(yr), y0+dy+8, lab, 10.5, col, "middle",
                 weight="bold" if col in (RED, ACC) else None)
    for y0 in (95, 185):
        b += line(X(1977.8 if y0 == 95 else 1978.8), y0+13,
                  X(2024.2), y0+13, RED, 3, opacity=0.35)
        b += txt(X(2000), y0+26, '47 years, no exponent change', 10,
                 RED, "middle")
    svg('fig_timeline', 760, 262, b,
        'Timeline of code bound records 1948-2025')


# ================= 3 fig_rates =================
def fig_rates():
    b = ''
    x0, y0, w, h = 55, 255, 300, 205
    X = mapper(0, 0.5, x0, x0+w)
    Y = mapper(0, 1.0, y0, y0-h)
    b += line(x0, y0, x0+w, y0, INK, 1.2)
    b += line(x0, y0, x0, y0-h, INK, 1.2)
    for v in (0.1, 0.2, 0.3, 0.4, 0.5):
        b += line(X(v), y0, X(v), y0+4, INK, 1)
        b += txt(X(v), y0+16, f'{v:.1f}', 9.5, MUT, "middle")
    for v in (0.25, 0.5, 0.75, 1.0):
        b += line(x0-4, Y(v), x0, Y(v), INK, 1)
        b += txt(x0-7, Y(v)+3, f'{v:.2f}', 9.5, MUT, "end")
    b += txt(x0+w/2, y0+32, '&#948;', 13, INK, "middle", cls="v")
    b += txt(x0-8, y0-h-10, 'rate', 11, INK, "end")
    ds = [i/400 for i in range(2, 200)]
    gv = [(X(d), Y(max(0, 1-H2(d)))) for d in ds]
    m1c = [(X(d), Y(M1(d))) for d in ds]
    m2c = [(X(d), Y(M2(d)[0])) for d in ds]
    band = m2c + gv[::-1]
    b += poly(band, "none", 0, fill=LT, close=True)
    b += poly(gv, GREEN, 1.9)
    b += poly(m1c, RED, 1.6, dash="5,3")
    b += poly(m2c, RED, 1.9)
    lx, ly = X(0.27), Y(0.97)
    b += line(lx, ly, lx+26, ly, GREEN, 1.9)
    b += txt(lx+32, ly+4, 'GV lower bound', 10, GREEN)
    b += line(lx, ly+16, lx+26, ly+16, RED, 1.6, dash="5,3")
    b += txt(lx+32, ly+20, sub('M', '1', 10), 10, RED, cls="v")
    b += line(lx, ly+32, lx+26, ly+32, RED, 1.9)
    b += txt(lx+32, ly+36, sub('M', '2', 10), 10, RED, cls="v")
    b += txt(X(0.36), Y(0.60), 'unknown territory:', 10.5, MUT,
             "middle")
    b += txt(X(0.36), Y(0.60)+13, sub('R', '2', 10) +
             '(&#948;) lives here', 10, MUT, "middle")
    b += line(X(0.345), Y(0.60)+20, X(0.31), Y(0.21), MUT, 1.0,
              marker=ar('fig_rates'))
    # right: F_delta(tau)
    x1, w1 = 420, 290
    X2 = mapper(0, 0.85, x1, x1+w1)
    Y2 = mapper(0.55, 0.75, y0, y0-h)
    b += line(x1, y0, x1+w1, y0, INK, 1.2)
    b += line(x1, y0, x1, y0-h, INK, 1.2)
    b += txt(x1+w1/2, y0+32, '&#964;', 13, INK, "middle", cls="v")
    b += txt(x1-8, y0-h-10, sub('F', '&#948;', 11) + '(&#964;)', 11,
             INK, "end", cls="v")
    for v in (0.2, 0.4, 0.6, 0.8):
        b += line(X2(v), y0, X2(v), y0+4, INK, 1)
        b += txt(X2(v), y0+16, f'{v:.1f}', 9.5, MUT, "middle")
    for v in (0.56, 0.60, 0.64, 0.68, 0.72):
        b += line(x1-4, Y2(v), x1, Y2(v), INK, 1)
        b += txt(x1-7, Y2(v)+3, f'{v:.2f}', 9.5, MUT, "end")
    for d, col in ((0.10, RED), (0.15, BLUE)):
        ts = [i/300*(1-2*d) for i in range(1, 300)]
        pts = [(X2(t), Y2(Fdel(t, d))) for t in ts
               if 0.553 <= Fdel(t, d) <= 0.748]
        b += poly(pts, col, 1.8)
        v, tstar = M2(d)
        b += circ(X2(tstar), Y2(v), 4, col)
    b += txt(X2(0.07), Y2(0.731), '&#948; = 0.10', 10.5, RED)
    b += txt(X2(0.38), Y2(0.562), '&#948; = 0.15', 10.5, BLUE)
    b += txt(X2(0.13), Y2(0.688)+16, 'interior minima &#964;*', 10,
             MUT)
    svg('fig_rates', 760, 300, b,
        'Computed rate-distance landscape and the MRRW objective')


# ================= 4 fig_sphrates =================
def fig_sphrates():
    b = ''
    x0, y0, w, h = 55, 255, 310, 200
    X = mapper(0, 1.0, x0, x0+w)
    Y = mapper(0, 1.6, y0, y0-h)
    b += line(x0, y0, x0+w, y0, INK, 1.2)
    b += line(x0, y0, x0, y0-h, INK, 1.2)
    b += txt(x0+w/2, y0+30, 's', 13, INK, "middle", cls="v")
    b += txt(x0-8, y0-h-10, 'exponent', 11, INK, "end")
    for v in (0.2, 0.4, 0.6, 0.8, 1.0):
        b += line(X(v), y0, X(v), y0+4, INK, 1)
        b += txt(X(v), y0+16, f'{v:.1f}', 9.5, MUT, "middle")
    for v in (0.4, 0.8, 1.2, 1.6):
        b += line(x0-4, Y(v), x0, Y(v), INK, 1)
        b += txt(x0-7, Y(v)+3, f'{v:.1f}', 9.5, MUT, "end")
    ss = [i/200 for i in range(4, 196)]
    direct = [(X(s), Y(min(1.6, Hsph(a0fun(s))))) for s in ss
              if Hsph(a0fun(s)) <= 1.62]
    bkl = [(X(s), Y(BKL(s))) for s in ss if BKL(s) <= 1.58]
    b += poly(direct, RED, 1.5, dash="5,3")
    b += poly(bkl, RED, 2.0)
    b += txt(X(0.50), Y(1.06),
             sub('H', 'sph', 10) + '(' + sub('a', '0', 10) +
             '(s)) (direct)', 10.5, RED, "middle", cls="v")
    b += line(X(0.62), Y(1.01), X(0.72), Y(Hsph(a0fun(0.72))+0.02),
              RED, 0.9)
    b += txt(X(0.75), Y(BKL(0.75))+17, sub('B', 'KL', 11), 11, RED,
             cls="v")
    b += line(X(0.5), y0, X(0.5), Y(BKL(0.5)), MUT, 1.0, dash="3,3")
    b += circ(X(0.5), Y(BKL(0.5)), 4, BLUE)
    b += txt(X(0.28), Y(0.62), 'kissing: 0.400944', 10, BLUE,
             "middle")
    b += line(X(0.37), Y(0.585), X(0.485), Y(BKL(0.5)+0.03), BLUE,
              0.9)
    # right: ladder at s=1/2
    vals = [('classical ' + sub('&#954;&#773;', '0', 11), 0.4009442),
            ('one-row ' + sub('&#954;&#773;', 'row', 11), 0.39731),
            ('level 1', 0.39674), ('level 2', 0.39660)]
    x1 = 450
    Y2 = mapper(0.394, 0.402, y0, y0-h)
    b += line(x1-10, y0, x1+260, y0, INK, 1.2)
    b += line(x1-10, y0, x1-10, y0-h, INK, 1.2)
    for v in (0.396, 0.398, 0.400):
        b += line(x1-14, Y2(v), x1-10, Y2(v), INK, 1)
        b += txt(x1-17, Y2(v)+3, f'{v:.3f}', 9.5, MUT, "end")
    for i, (lab, v) in enumerate(vals):
        xx = x1 + 20 + i*60
        b += rect(xx, Y2(v), 40, y0-Y2(v), fill=[RED, GOLD, BLUE,
                  GREEN][i], stroke="none", rx=2, opacity=0.75)
        b += txt(xx+20, Y2(v)-6, f'{v:.5f}'[:7], 9.5, INK, "middle")
        b += txt(xx+20, y0+15+(12 if i % 2 else 0), lab, 9.5, INK,
                 "middle")
    b += txt(x1+120, y0-h+4, 'kissing exponent at s = 1/2', 11, INK,
             "middle")
    svg('fig_sphrates', 760, 300, b,
        'Computed spherical exponents and the kissing ladder')


# ================= 5 fig_roadmap =================
def fig_roadmap():
    b = ''

    def box(x, y, w, h, lab1, lab2, fill=LT, col=INK):
        s = rect(x, y, w, h, fill=fill, stroke=col, sw=1.2, rx=6)
        s += txt(x+w/2, y+h/2-3, lab1, 10.5, col, "middle",
                 weight="bold")
        if lab2:
            s += txt(x+w/2, y+h/2+11, lab2, 9.5, col, "middle")
        return s
    b += box(20, 20, 130, 34, '&#167;2 Lemma A', 'two-point LP')
    b += box(20, 80, 130, 34, '&#167;3 Prop 4.1', 'projection bound',
             LY, ACC)
    b += box(190, 20, 130, 34, '&#167;4 Boolean', 'harmonics')
    b += box(190, 80, 130, 34, '&#167;5 Thm 2.1', 'whole-cube bound')
    b += box(190, 140, 130, 34, '&#167;6 Thm 2.3',
             '&#954;H &lt; M1')
    b += box(190, 260, 130, 34, '&#167;7 Thm 3.4', 'layer bound')
    b += box(190, 320, 130, 34, '&#167;8 Thm 3.8',
             '&#954;CW beats interior')
    b += box(360, 200, 120, 34, 'Thm 1.1', 'binary record', LG, GREEN)
    b += box(360, 80, 120, 34, '&#167;9 Thm 4.2', 'graph machine', LY,
             ACC)
    b += box(520, 20, 110, 34, '&#167;10&#8211;11', 'harmonics, one-row')
    b += box(520, 80, 110, 34, '&#167;12 Thm 6.2', 'hierarchy bound')
    b += box(520, 140, 110, 34, '&#167;13 Thm 1.2', 'asymptotics')
    b += box(520, 200, 110, 34, '&#167;14 Cor 7.6', 'strict ladder',
             LG, GREEN)
    b += box(520, 280, 110, 34, '&#167;15 Thm 8.3',
             '&#955;* packing', LG, GREEN)
    b += box(520, 350, 110, 34, '&#167;16', 'dictionary')
    b += box(660, 300, 84, 50, 'Chapter 1', 'same &#955;*', LY, ACC)
    arrows = [
        (85, 54, 85, 80), (150, 97, 190, 97), (255, 54, 255, 80),
        (255, 114, 255, 140), (150, 97, 190, 274),
        (255, 174, 360, 210), (255, 294, 255, 320),
        (320, 337, 380, 234), (150, 97, 360, 92),
        (480, 97, 520, 92), (575, 54, 575, 80), (575, 114, 575, 140),
        (575, 174, 575, 200), (575, 234, 575, 280),
        (575, 314, 575, 350), (630, 310, 660, 316),
        (630, 367, 668, 350),
    ]
    for x1, y1, x2, y2 in arrows:
        b += line(x1, y1, x2, y2, MUT, 1.2, marker=ar('fig_roadmap'))
    svg('fig_roadmap', 760, 400, b,
        'Dependency map of the chapter')


# ================= 6 fig_delsarte =================
def fig_delsarte():
    b = ''
    x0, y0, w, h = 60, 210, 580, 160
    X = mapper(-1.05, 1.1, x0, x0+w)
    s = 0.4

    def F(t):
        return 2.2*(t-s)*(t+1)*(t*t+0.3)
    Y = mapper(-1.6, 4.3, y0, y0-h)
    b += line(x0, Y(0), x0+w, Y(0), INK, 1.2)
    b += rect(X(-1), Y(0), X(s)-X(-1), 26, fill=LT, stroke="none",
              rx=0, opacity=0.9)
    b += txt(X(-0.55), Y(1.35), 'allowed region  t &#8804; s', 10.5,
             BLUE, "middle")
    b += line(X(-0.55), Y(1.2), X(-0.55), Y(0)+8, BLUE, 1.0,
              marker=ar('fig_delsarte', 'B'))
    ts = [-1.02 + i*0.001*2.1 for i in range(1030)]
    b += poly([(X(t), Y(F(t))) for t in ts], RED, 2.0)
    b += line(X(1), Y(0), X(1), Y(F(1)), RED, 1.2, dash="3,3")
    b += circ(X(1), Y(F(1)), 4, RED)
    b += txt(X(1)-6, Y(F(1))-8, 'F(1): the diagonal payoff', 11, RED,
             "end")
    b += line(X(s), Y(0)-5, X(s), Y(0)+5, INK, 1.6)
    b += txt(X(s), Y(0)+38, 's', 12, INK, "middle", cls="v")
    b += line(X(-1), Y(0)-4, X(-1), Y(0)+4, INK, 1.2)
    b += txt(X(-1), Y(0)+38, '&#8722;1', 11, INK, "middle")
    b += line(X(1), Y(0)-4, X(1), Y(0)+4, INK, 1.2)
    b += txt(X(1), Y(0)+38, '1', 11, INK, "middle")
    b += txt(X(-0.3), Y(-1.05), 'F &#8804; 0 wherever pairs may live',
             11, MUT, "middle")
    b += line(X(0.06), Y(0), X(0.06), Y(1.15), GOLD, 4, opacity=0.5)
    b += txt(X(0.06)-8, Y(0.75), sub('f', '0', 11) +
             ' &gt; 0: guaranteed constant component', 10.5, GOLD,
             "end")
    b += txt(350, 30, 'bound delivered:  |C| &#8804; F(1) / ' +
             sub('f', '0', 12), 13, INK, "middle", weight="bold")
    svg('fig_delsarte', 700, 265, b,
        'Anatomy of a Delsarte certificate',
        vb='40 5 640 225')


# ================= 7 fig_toycert =================
def fig_toycert():
    b = ''
    x0, y0, w, h = 55, 240, 360, 185
    X = mapper(-1.1, 1.15, x0, x0+w)
    Y = mapper(-0.6, 4.3, y0, y0-h)

    def F(t):
        return 0.5*(1+t)*(1+3*t)
    b += line(x0, Y(0), x0+w, Y(0), INK, 1.2)
    ts = [-1.05 + i*0.0022 for i in range(1010)]
    b += poly([(X(t), Y(F(t))) for t in ts], RED, 2.0)
    for t, lab in ((-1, '&#8722;1'), (-1/3, '&#8722;1/3'),
                   (1/3, '1/3'), (1, '1')):
        b += line(X(t), Y(0)-4, X(t), Y(0)+4, INK, 1.4)
        b += txt(X(t), Y(0)+18, lab, 10.5, INK, "middle")
    b += circ(X(-1), Y(0), 4.5, BLUE)
    b += circ(X(-1/3), Y(0), 4.5, BLUE)
    b += txt(X(-2/3), Y(0)-12, 'allowed: F = 0', 10.5, BLUE, "middle")
    b += circ(X(1/3), Y(F(1/3)), 4, MUT)
    b += txt(X(1/3)-8, Y(F(1/3))-6, 'forbidden gap', 10, MUT,
             "end")
    b += circ(X(1), Y(4), 5, RED)
    b += txt(X(1)-8, Y(4)+2, 'F(1) = 4', 11.5, RED, "end",
             weight="bold")
    b += txt(x0+w/2, 35, 'F(t) = &#189;(1+t)(1+3t) = 1 + ' +
             '&#8532;&#183;(3t) + &#8531;&#183;' +
             sub('K', '2', 11) + '(t)', 12, INK, "middle")
    # right: cube with tetrahedron
    cx, cy, sc = 560, 150, 58
    verts = {}
    for v in range(8):
        p = ((v & 1), (v >> 1) & 1, (v >> 2) & 1)
        x, y, dep = ortho((p[0]*2-1, p[1]*2-1, p[2]*2-1), 0.5, 0.42)
        verts[p] = (cx + sc*x, cy + sc*y)
    for v in verts:
        for i in range(3):
            w2 = list(v)
            w2[i] ^= 1
            w2 = tuple(w2)
            if w2 > v:
                b += line(*verts[v], *verts[w2], MUT, 1.0, opacity=0.6)
    code = [(0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1)]
    for i in range(4):
        for j in range(i+1, 4):
            b += line(*verts[code[i]], *verts[code[j]], RED, 1.4,
                      dash="4,3", opacity=0.85)
    for p in verts:
        b += circ(*verts[p], 6 if p in code else 3.5,
                  RED if p in code else MUT)
    b += txt(560, 265, 'the even-weight code: 4 words, distance 2 '
             '&#8212; bound met', 10.5, INK, "middle")
    svg('fig_toycert', 700, 300, b,
        'The tight toy certificate for n=3, d=2')


# ================= 8 fig_vectorsub =================
def fig_vectorsub():
    b = ''
    for px, tag in ((30, 'classical'), (400, 'new')):
        b += rect(px, 30, 330, 250, fill="none", stroke=MUT, sw=1.0,
                  rx=8, dash="4,3")
        b += txt(px+165, 52, tag, 12, ACC, "middle", weight="bold")
    # left: vectors
    ox, oy = 190, 200
    b += circ(ox, oy, 3, INK)
    b += line(ox, oy, ox+95, oy-75, RED, 2.0, marker=ar(
        'fig_vectorsub', 'R'))
    b += line(ox, oy, ox-55, oy-95, BLUE, 2.0, marker=ar(
        'fig_vectorsub', 'B'))
    b += txt(ox+100, oy-82, sub('u', 'x', 11), 12, RED, cls="v")
    b += txt(ox-60, oy-100, sub('u', 'y', 11), 12, BLUE, "end", cls="v")
    b += path(f'M{ox+28:.0f},{oy-22:.0f} A36,36 0 0 0 '
              f'{ox-14:.0f},{oy-25:.0f}', MUT, 1.2)
    b += txt(ox+6, oy-38, '&#9001;' + sub('u', 'x', 10) + ',' +
             sub('u', 'y', 10) + '&#9002;', 10.5, MUT, "middle")
    b += txt(195, 262, 'one unit vector per point; rank 1', 10.5, MUT,
             "middle")
    # right: subspace patches
    ox2, oy2 = 560, 195
    b += circ(ox2, oy2, 3, INK)
    b += poly([(ox2-70, oy2-10), (ox2+20, oy2-58), (ox2+95, oy2-30),
               (ox2+5, oy2+18)], RED, 1.4, fill=RED, opacity=0.25,
              close=True)
    b += poly([(ox2-80, oy2-52), (ox2+8, oy2-96), (ox2+60, oy2-64),
               (ox2-28, oy2-20)], BLUE, 1.4, fill=BLUE, opacity=0.25,
              close=True)
    b += txt(ox2+100, oy2-24, 'im&#8201;' + sub('P', 'x', 11), 12, RED,
             cls="v")
    b += txt(ox2-84, oy2-58, 'im&#8201;' + sub('P', 'y', 11), 12, BLUE,
             "end", cls="v")
    b += txt(ox2+2, oy2+42, 'tr(' + sub('P', 'x', 10) +
             sub('P', 'y', 10) + ') = overlap of subspaces', 10.5, MUT,
             "middle")
    b += txt(560, 262, 'a rank-d moving subspace per point', 10.5, MUT,
             "middle")
    b += path('M 340,150 C 365,130 375,130 398,148', GOLD, 2.0,
              marker=ar('fig_vectorsub', 'G'))
    b += txt(370, 122, 'd = 1 &#8594; d = ' + sup('2', '&#916;n', 10),
             10.5, GOLD, "middle")
    b += txt(30, 305, 'both kernels are positive definite; the bound '
             'divides by the rank:  |C| &#8804; (1&#8722;s)/'
             '(&#923;&#8722;s) &#183; D/d', 11.5, INK)
    svg('fig_vectorsub', 760, 320, b,
        'Vector per point versus moving subspace per point')


# ================= 9 fig_overlap =================
def fig_overlap():
    b = ''
    cx, cy, sc = 300, 165, 105
    th = math.radians(38)

    def pr(p):
        x, y, dep = ortho(p, 0.62, 0.38)
        return (cx + sc*x, cy - sc*y)
    # plane 1: span e1,e2 ; plane2: span e1, cos e2 + sin e3
    c1 = [pr((sx*1.05, sy*0.85, 0))
          for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    e2r = (0, math.cos(th), math.sin(th))
    c2 = [pr((sx*1.05, sy*0.85*e2r[1], sy*0.85*e2r[2]))
          for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    b += poly(c1, RED, 1.4, fill=RED, opacity=0.22, close=True)
    b += poly(c2, BLUE, 1.4, fill=BLUE, opacity=0.22, close=True)
    p1 = pr((-1.25, 0, 0))
    p2 = pr((1.25, 0, 0))
    b += line(*p1, *p2, INK, 2.0)
    b += txt(p2[0]+8, p2[1]+4, 'shared line: &#952;&#8321; = 0', 10.5,
             INK)
    a1 = pr((0, 0.62, 0))
    a2 = pr((0, 0.62*e2r[1], 0.62*e2r[2]))
    b += path(f'M{a1[0]:.1f},{a1[1]:.1f} Q{cx:.1f},{cy-40:.1f} '
              f'{a2[0]:.1f},{a2[1]:.1f}', INK, 1.2)
    b += txt((a1[0]+a2[0])/2+26, (a1[1]+a2[1])/2-14,
             '&#952;&#8322;', 12, INK, cls="v")
    b += txt(c1[2][0]+8, c1[2][1]+12, 'im&#8201;' + sub('P', 'x', 11),
             12, RED, cls="v")
    b += txt(c2[3][0]-6, c2[3][1]-8, 'im&#8201;' + sub('P', 'y', 11),
             12, BLUE, "end", cls="v")
    val = 1 + math.cos(th)**2
    b += txt(300, 300, 'tr(' + sub('P', 'x', 11) + sub('P', 'y', 11) +
             ') = cos&#178;&#952;&#8321; + cos&#178;&#952;&#8322; = '
             f'1 + cos&#178;{int(math.degrees(th))}&#176; = {val:.3f}',
             12, INK, "middle")
    svg('fig_overlap', 600, 320, b,
        'Two planes meeting at principal angles; overlap is the sum '
        'of squared cosines', vb='130 95 385 225')


# ================= 10 fig_residual =================
def fig_residual():
    b = ''
    ox, oy = 120, 150
    # triangle: L_x = sqrt(L) G_x + Theta_x
    Lx = (300, -60)
    Gx = (255, 8)
    b += line(ox, oy, ox+Lx[0], oy+Lx[1], RED, 2.2,
              marker=ar('fig_residual', 'R'))
    b += line(ox, oy, ox+Gx[0], oy+Gx[1], BLUE, 2.2,
              marker=ar('fig_residual', 'B'))
    b += line(ox+Gx[0], oy+Gx[1], ox+Lx[0], oy+Lx[1], GOLD, 2.0,
              dash="5,3", marker=ar('fig_residual', 'G'))
    b += txt(ox+Lx[0]/2-20, oy+Lx[1]/2-14,
             sub('L', 'x', 11) + ' = (&#8467;' + sub('', 'x', 9) +
             '&#8855;id)' + sub('P', 'x', 11), 11.5, RED, "middle",
             cls="v")
    b += txt(ox+Gx[0]/2, oy+Gx[1]/2+20, '&#8730;&#923;&#8201;' +
             sub('G', 'x', 11) + ' = &#8730;&#923;&#8201;B' +
             sub('P', 'x', 11), 11.5, BLUE, "middle", cls="v")
    b += txt(ox+312, oy-84, sub('&#920;', 'x', 11) + ' (residual)',
             11.5, GOLD, cls="v")
    b += txt(320, 250, 'alignment (49) kills the cross terms:  '
             '&#9001;' + sub('&#920;', 'x', 10) + ',' +
             sub('&#920;', 'y', 10) + '&#9002; = (t(x,y) &#8722; '
             '&#923;)&#183;K(x,y)   (51)', 11.5, INK, "middle")
    b += txt(320, 278, 'summed over a code: 0 &#8804; &#8214;&#931;'
             + sub('&#920;', 'x', 10) + '&#8214;&#178;  and  trace '
             'Cauchy&#8211;Schwarz  &#8658;  the bound (50)', 11,
             MUT, "middle")
    svg('fig_residual', 640, 300, b,
        'The residual decomposition behind Proposition 4.1',
        vb='95 42 450 258')


# ================= 11 fig_levels =================
def fig_levels():
    b = ''
    n = 8
    cx = 300
    y0 = 344
    for i in range(n+1):
        wdt = 30 + 200*math.comb(n, i)/math.comb(n, 4)
        y = y0 - i*36
        fill = LY if i == 2 else LT
        b += rect(cx-wdt/2, y-13, wdt, 26, fill=fill, stroke=INK,
                  sw=1.1, rx=5)
        b += txt(cx-wdt/2-10, y+4, sub('V', str(i), 11) +
                 sup('', '&#9633;', 10), 11.5, INK, "end", cls="v")
        b += txt(cx+wdt/2+10, y+4, 'dim ' + f'{math.comb(n, i)}', 10,
                 MUT)
    b += rect(cx-15, y0-2*36-13, 96, 26, fill="none", stroke=GOLD,
              sw=2.0, rx=5)
    b += txt(cx+120, y0-2*36-18, sub('E', 'k', 11) + ' = ker D', 11,
             GOLD, weight="bold")
    b += txt(cx+120, y0-2*36+0, 'dim = C(n,k) &#8722; C(n,k&#8722;1)',
             9.5, GOLD)
    b += line(cx-160, y0-8, cx-160, y0-8*36+8, BLUE, 1.8,
              marker=ar('fig_levels', 'B'))
    b += txt(cx-170, y0-4*36, 'U', 13, BLUE, "end", cls="v")
    b += line(cx-195, y0-8*36+8, cx-195, y0-8, RED, 1.8,
              marker=ar('fig_levels', 'R'))
    b += txt(cx-205, y0-4*36, 'D = U*', 12, RED, "end", cls="v")
    b += txt(cx, 22, 'DU &#8722; UD = (n &#8722; 2i)&#183;id on level i',
             12, INK, "middle")
    svg('fig_levels', 620, 384, b, 'The Boolean Fourier ladder')


# ================= 12 fig_harmonic =================
def fig_harmonic():
    b = ''
    data = [
        ('level 1:  h', ['3', '&#8722;1', '&#8722;1', '&#8722;1'],
         [3, -1, -1, -1], 60),
        ('level 2:  Uh', ['2', '2', '2', '&#8722;2', '&#8722;2',
                          '&#8722;2'], [2, 2, 2, -2, -2, -2], 250),
        ('level 3:  U&#178;h', ['2', '2', '2', '&#8722;6'],
         [2, 2, 2, -6], 462),
    ]
    labs = [['{1}', '{2}', '{3}', '{4}'],
            ['12', '13', '14', '23', '24', '34'],
            ['123', '124', '134', '234']]
    for (title, vs, nums, x0), lab in zip(data, labs):
        y0 = 190
        b += txt(x0+62, 46, title, 11.5, INK, "middle", weight="bold")
        b += line(x0, y0, x0+140, y0, MUT, 1.0)
        for i, (v, num) in enumerate(zip(vs, nums)):
            xx = x0 + 12 + i*(120/max(1, len(vs)-1))
            hh = num*12
            col = RED if num > 0 else BLUE
            b += rect(xx-6, y0-max(hh, 0), 12, abs(hh), fill=col,
                      stroke="none", rx=2, opacity=0.8)
            b += txt(xx, y0+14, lab[i], 8.5, MUT, "middle")
            b += txt(xx, y0-hh-4 if num > 0 else y0-hh+12, v, 9, col,
                     "middle")
    b += txt(180, 285, '&#8214;Uh&#8214;&#178; = 2&#183;'
             '&#8214;h&#8214;&#178;', 11, INK, "middle")
    b += txt(392, 285, '&#8214;U&#178;h&#8214;&#178; = 2&#183;'
             '&#8214;Uh&#8214;&#178;', 11, INK, "middle")
    b += line(205, 120, 245, 120, GOLD, 2, marker=ar('fig_harmonic',
              'G'))
    b += line(415, 120, 455, 120, GOLD, 2, marker=ar('fig_harmonic',
              'G'))
    b += txt(225, 110, 'U', 12, GOLD, "middle", cls="v")
    b += txt(435, 110, 'U', 12, GOLD, "middle", cls="v")
    svg('fig_harmonic', 660, 300, b,
        'A computed harmonic staircase for n=4, k=1')


# ================= 13 fig_path =================
def fig_path():
    b = ''
    n, k = 30, 3
    nodes = ['k', 'k+1', '&#8230;', 'i', 'i+1', '&#8230;', 'L']
    xs = [70 + i*100 for i in range(7)]
    y = 110
    for i, (nd, x) in enumerate(zip(nodes, xs)):
        if nd == '&#8230;':
            b += txt(x, y+5, '&#8943;', 16, MUT, "middle")
            continue
        b += circ(x, y, 26, LT, INK, 1.3)
        b += txt(x, y-2, sub('V', nd, 11) + sup('', '&#9633;', 9), 11,
                 INK, "middle", cls="v")
        b += poly([(x-14, y+13), (x, y+6), (x+14, y+13), (x, y+20)],
                  GOLD, 1.2, fill=GOLD, opacity=0.5, close=True)
        b += txt(x, y+38, '&#966;' + sub('', 'i,x', 9) +
                 sub('E', 'k', 10), 9.5, GOLD, "middle", cls="v")
    for i in range(6):
        x1, x2 = xs[i]+28, xs[i+1]-28
        if nodes[i] == '&#8230;' or nodes[i+1] == '&#8230;':
            x1, x2 = xs[i]+(30 if nodes[i] != '&#8230;' else 14), \
                xs[i+1]-(30 if nodes[i+1] != '&#8230;' else 14)
        b += line(x1, y-6, x2, y-6, INK, 1.5)
        b += line(x2, y+6, x1, y+6, INK, 1.5)
    lab = ['c' for _ in range(6)]
    b += txt(xs[0]+50, y-16, sub('c', 'k,H', 9) + sup('', '(k)', 9),
             10, INK, "middle", cls="v")
    b += txt(xs[3]+50, y-16, sub('c', 'i,H', 9) + sup('', '(k)', 9),
             10, INK, "middle", cls="v")
    b += txt(xs[5]+50, y-16, sub('c', 'L&#8722;1,H', 9) +
             sup('', '(k)', 9), 10, INK, "middle", cls="v")
    b += txt(370, 30, 'every vertex carries the same harmonic fiber ' +
             sub('E', 'k', 11) + ', transported to the code point x',
             11.5, INK, "middle")
    b += txt(370, 185, 'symmetric edge weight (20):  ' +
             sub('c', 'i,H', 10) + sup('', '(k)', 10) +
             ' = (i&#8722;k+1)(n&#8722;i&#8722;k) / '
             '[n&#8730;((i+1)(n&#8722;i))]', 11, MUT, "middle")
    svg('fig_path', 740, 205, b,
        'The whole-cube representation path with harmonic fibers')


# ================= 14 fig_sine =================
def fig_sine():
    b = ''
    n, bpar, apar = 400, 0.10, 0.35
    k, L = int(bpar*n), int(apar*n)
    x0, y0, w, h = 60, 150, 600, 100
    X = mapper(k/n, L/n, x0, x0+w)
    ws = [(i-k+1)*(n-i-k)/(n*math.sqrt((i+1)*(n-i)))
          for i in range(k, L)]
    Y = mapper(0, max(ws)*1.15, y0, y0-h)
    b += line(x0, y0, x0+w, y0, INK, 1.2)
    b += poly([(X((i+0.5)/n), Y(wv)) for i, wv in zip(range(k, L), ws)],
              BLUE, 1.8)
    gb = 2*(apar-bpar)*(1-apar-bpar)/math.sqrt(apar*(1-apar))/2
    b += line(x0, Y(gb), x0+w, Y(gb), RED, 1.2, dash="5,3")
    b += txt(x0+w-4, Y(gb)-8, sub('g', 'b', 10) + '(a) = &#915;' +
             sub('', 'H', 9) + '/2', 10.5, RED, "end", cls="v")
    b += txt(x0+w/2, y0+20, 'degree i/n  (n = 400, b = 0.10, '
             'a = 0.35)', 10.5, MUT, "middle")
    b += txt(x0-8, y0-h-6, 'edge weight', 10, INK, "end")
    b += txt(x0+130, y0-h+38, 'weights increase along the path',
             10.5, BLUE)
    # bottom: sine bump
    y1 = 300
    mwin = 24
    b += line(x0, y1, x0+w, y1, INK, 1.2)
    for r in range(1, mwin+1):
        i = L - mwin + r
        xx = X(i/n)
        hh = 70*math.sin(math.pi*r/(mwin+1))
        b += line(xx, y1, xx, y1-hh, GOLD, 3.0, opacity=0.8)
    b += txt(X((L-mwin/2)/n), y1-84, 'sine test vector on the last '
             'm degrees', 10.5, GOLD, "middle")
    b += txt(x0+w/2, y1+20, 'Rayleigh quotient 2w&#183;cos(&#960;/'
             '(m+1)) &#8594; matches the row-sum cap', 10.5, MUT,
             "middle")
    svg('fig_sine', 720, 330, b,
        'Edge weight profile and the terminal sine test vector')


# ================= 15 fig_lammax =================
def fig_lammax():
    b = ''
    apar, bpar = 0.35, 0.10
    gam = GammaH(apar, bpar)
    ns = [100, 160, 250, 400, 630, 1000, 1600, 2500, 4000]
    lams = []
    for n in ns:
        k, L = int(bpar*n), int(apar*n)
        off = [(i-k+1)*(n-i-k)/(n*math.sqrt((i+1)*(n-i)))
               for i in range(k, L)]
        lams.append(lam_tri([0.0]*(L-k+1), off))
    x0, y0, w, h = 70, 230, 560, 175
    X = mapper(math.log10(100), math.log10(4000), x0, x0+w)
    Y = mapper(min(lams)-0.005, gam+0.008, y0, y0-h)
    b += line(x0, y0, x0+w, y0, INK, 1.2)
    b += line(x0, y0, x0, y0-h, INK, 1.2)
    for n in (100, 400, 1000, 4000):
        b += line(X(math.log10(n)), y0, X(math.log10(n)), y0+4, INK, 1)
        b += txt(X(math.log10(n)), y0+16, str(n), 9.5, MUT, "middle")
    b += txt(x0+w/2, y0+32, 'n  (log scale)', 11, INK, "middle")
    b += line(x0, Y(gam), x0+w, Y(gam), RED, 1.4, dash="6,3")
    b += txt(x0+w, Y(gam)-8, '&#915;' + sub('', 'H', 10) +
             '(0.35, 0.10) = %.5f' % gam, 10.5, RED, "end")
    pts = [(X(math.log10(n)), Y(lv)) for n, lv in zip(ns, lams)]
    b += poly(pts, BLUE, 1.8)
    for p in pts:
        b += circ(*p, 3.5, BLUE)
    b += txt(x0+w/2, Y(lams[2])+30, '&#955;' + sub('', 'max', 9) +
             '(' + sub('J', 'H', 10) + '(n, 0.1n, 0.35n)) '
             '(computed exactly)', 10.5, BLUE, "middle")
    b += txt(x0+w-4, Y(lams[-1])+18, 'gap %.1e at n = 4000'
             % (gam-lams[-1]), 10, MUT, "end")
    svg('fig_lammax', 700, 275, b,
        'Computed convergence of the top eigenvalue to Gamma_H')


# ================= 16 fig_feasible =================
def fig_feasible():
    b = ''
    d = 0.2
    s = 1-2*d
    x0, y0, w, h = 80, 290, 520, 240
    X = mapper(0, 0.5, x0, x0+w)
    Y = mapper(0, 0.28, y0, y0-h)
    # feasible region boundary: GammaH(a,b) = s ; solve for b at each a
    a00 = 0.5 - math.sqrt(d*(1-d))
    # shade feasible region GammaH > s (numerically)
    grid = 130
    pts_bound = []
    for ia in range(grid+1):
        a = 0.001 + ia*(0.499-0.001)/grid
        # find b bound: largest b with GammaH(a,b)>s (b in [0,a))
        lo, hi = 0.0, a-1e-9
        if GammaH(a, 1e-12) <= s:
            continue
        for _ in range(60):
            mid = 0.5*(lo+hi)
            if GammaH(a, mid) > s:
                lo = mid
            else:
                hi = mid
        if lo > 1e-7:
            pts_bound.append((a, lo))
    if pts_bound:
        poly_pts = ([(X(a), Y(0)) for a, _ in pts_bound[:1]]
                    + [(X(a), Y(bb)) for a, bb in pts_bound]
                    + [(X(pts_bound[-1][0]), Y(0))])
        b += poly(poly_pts, "none", 0, fill=LT, close=True)
        b += poly([(X(a), Y(bb)) for a, bb in pts_bound], BLUE, 1.8)
    b += line(x0, y0, x0+w, y0, INK, 1.2)
    b += line(x0, y0, x0, y0-h, INK, 1.2)
    for v in (0.1, 0.2, 0.3, 0.4, 0.5):
        b += line(X(v), y0, X(v), y0+4, INK, 1)
        b += txt(X(v), y0+16, f'{v:.1f}', 9.5, MUT, "middle")
    for v in (0.1, 0.2):
        b += line(x0-4, Y(v), x0, Y(v), INK, 1)
        b += txt(x0-7, Y(v)+3, f'{v:.1f}', 9.5, MUT, "end")
    b += txt(x0+w/2, y0+32, 'a  (top Fourier degree / n)', 11, INK,
             "middle")
    b += txt(x0-62, y0-h-10, 'b  (fiber degree / n)', 11, INK)
    # level sets of objective
    vopt, aopt, bopt = kappaH(d)
    for lev in (vopt+0.02, vopt+0.08, vopt+0.16, vopt+0.26):
        pts = []
        for ia in range(400):
            a = 0.02 + ia*0.478/400
            # solve H2(a)-H2(b)=lev for b
            if H2(a) < lev:
                continue
            lo, hi = 0.0, min(a, 0.5)-1e-9
            for _ in range(50):
                mid = 0.5*(lo+hi)
                if H2(a)-H2(mid) > lev:
                    lo = mid
                else:
                    hi = mid
            if 0 <= lo < 0.28:
                pts.append((X(a), Y(lo)))
        if pts:
            b += poly(pts, MUT, 0.9, dash="2,3", opacity=0.8)
    b += circ(X(a00), Y(0), 5, RED)
    b += txt(X(a00)-10, Y(0.02), '(' + sub('a', '0', 10) +
             ', 0): classical, delivers ' + sub('M', '1', 10), 10.5,
             RED, "end")
    b += circ(X(aopt), Y(bopt), 5, GOLD)
    b += txt(X(0.03), Y(0.262), 'optimizer: &#954;' +
             sub('', 'H', 9) + ' = %.5f' % vopt, 10.5, GOLD)
    b += line(X(0.10), Y(0.252), X(aopt)-4, Y(bopt)+4, GOLD, 0.9)
    c = 2/(1-2*a00) + 1.2
    b += line(X(a00), Y(0), X(a00+c*0.05), Y(0.05), RED, 1.6,
              marker=ar('fig_feasible', 'R'))
    b += txt(X(0.16), Y(0.145), 'improving ray a = ' +
             sub('a', '0', 10) + ' + cb', 10, RED)
    b += line(X(0.20), Y(0.135), X(0.235), Y(0.062), RED, 0.9)
    b += txt(X(0.42), Y(0.20), 'feasible: &#915;' + sub('', 'H', 9) +
             '(a,b) &gt; 1&#8722;2&#948;', 11, BLUE, "middle")
    b += txt(X(0.34), Y(0.255), 'level sets of ' + sub('H', '2', 9) +
             '(a) &#8722; ' + sub('H', '2', 9) + '(b)', 9.5, MUT)
    svg('fig_feasible', 700, 340, b,
        'Feasibility landscape of the whole-cube exponent at '
        'delta = 0.2')


def objH_at(aa, bb, d):
    s2 = 1 - 2*d
    try:
        from scipy.optimize import minimize
        res = minimize(
            lambda x: H2(x[0]) - H2(x[1]), (aa, bb), method='SLSQP',
            constraints=[{'type': 'ineq', 'fun':
                          lambda x: GammaH(x[0], x[1]) - s2 - 1e-10},
                         {'type': 'ineq',
                          'fun': lambda x: x[0] - x[1] - 1e-9}],
            bounds=[(1e-9, 0.5), (0.0, 0.5)],
            options={'maxiter': 400, 'ftol': 1e-15})
        if res.success and GammaH(*res.x) >= s2 - 1e-9:
            return float(res.fun)
    except Exception:
        pass
    return None


def smooth_lower(ds, vs):
    # Kill upward optimizer-noise spikes: each value is at most the
    # linear interpolation of its neighbors (all values are upper
    # bounds on a smooth infimum, so min() is always safe).
    out = list(vs)
    for _ in range(4):
        for i in range(1, len(out)-1):
            out[i] = min(out[i], 0.5*(out[i-1]+out[i+1]))
    return out


# ================= 17 fig_m1kh =================
def fig_m1kh():
    b = ''
    ds = [0.04 + i*0.01 for i in range(43)]
    m1s = [M1(d) for d in ds]
    khs = []
    prev = None
    for d in ds:
        v, aa, bb = kappaH(d)
        if prev is not None:
            va = objH_at(prev[0], prev[1], d)
            if va is not None and va < v:
                v, aa, bb = va, prev[0], prev[1]
        khs.append(v)
        prev = (aa, bb)
    khs = smooth_lower(ds, khs)
    x0, y0, w, h = 60, 240, 280, 190
    X = mapper(0, 0.5, x0, x0+w)
    Y = mapper(0, 1.0, y0, y0-h)
    b += line(x0, y0, x0+w, y0, INK, 1.2)
    b += line(x0, y0, x0, y0-h, INK, 1.2)
    b += txt(x0+w/2, y0+30, '&#948;', 12, INK, "middle", cls="v")
    for v in (0.25, 0.5, 0.75, 1.0):
        b += line(x0-4, Y(v), x0, Y(v), INK, 1)
        b += txt(x0-7, Y(v)+3, f'{v:.2f}', 9, MUT, "end")
    for v in (0.1, 0.2, 0.3, 0.4, 0.5):
        b += line(X(v), y0, X(v), y0+4, INK, 1)
        b += txt(X(v), y0+14, f'{v:.1f}', 9, MUT, "middle")
    b += poly([(X(d), Y(v)) for d, v in zip(ds, m1s)], RED, 1.7)
    b += poly([(X(d), Y(v)) for d, v in zip(ds, khs)], BLUE, 1.5,
              dash="4,2")
    b += txt(X(0.17), Y(0.62), sub('M', '1', 10) +
             ' (red) vs &#954;' + sub('', 'H', 9) + ' (blue)', 10.5,
             INK)
    # right: log gap
    x1, w1 = 410, 260
    X2 = mapper(0, 0.5, x1, x1+w1)
    Y2 = mapper(-6.6, -1.4, y0, y0-h)
    b += line(x1, y0, x1+w1, y0, INK, 1.2)
    b += line(x1, y0, x1, y0-h, INK, 1.2)
    b += txt(x1+w1/2, y0+30, '&#948;', 12, INK, "middle", cls="v")
    for e in (-2, -3, -4, -5, -6):
        b += line(x1-4, Y2(e), x1, Y2(e), INK, 1)
        b += txt(x1-7, Y2(e)+3, '10' + sup('', str(e), 9), 9.5, MUT,
                 "end")
    for v in (0.1, 0.2, 0.3, 0.4):
        b += line(X2(v), y0, X2(v), y0+4, INK, 1)
        b += txt(X2(v), y0+14, f'{v:.1f}', 9, MUT, "middle")
    gpts = [(X2(d), Y2(math.log10(max(1e-6, m-kk))))
            for d, m, kk in zip(ds, m1s, khs) if m-kk > 1e-6]
    b += poly(gpts, ACC, 1.9)
    b += txt(X2(0.10), Y2(-6.1), 'gap ' + sub('M', '1', 10) +
             ' &#8722; &#954;' + sub('', 'H', 9) + ' (log scale)',
             11, ACC)
    svg('fig_m1kh', 700, 285, b,
        'The first record falls: M1 versus kappa_H and their gap')


def main():
    fig_twoworlds()
    fig_timeline()
    fig_rates()
    fig_sphrates()
    fig_roadmap()
    fig_delsarte()
    fig_toycert()
    fig_vectorsub()
    fig_overlap()
    fig_residual()
    fig_levels()
    fig_harmonic()
    fig_path()
    fig_sine()
    fig_lammax()
    fig_feasible()
    fig_m1kh()


if __name__ == "__main__":
    main()
