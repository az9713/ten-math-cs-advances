"""Chapter 7 figure generator -> .ignore/ch07_figs/*.svg

All ids/markers prefixed c7. House palette:
  accent red #c0392b, blue #2e6da4, purple #7d3c98, green #2a7d2a,
  gray #666, light fills #eef4fa / #fdf0ee / #f0f7f0.
Usage: PYTHONIOENCODING=utf-8 python figgen/ch07_svgs.py
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch07_figlib import Ax, fmt, svg_wrap, write_fig  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, ".ignore", "ch07_figs")

RED = "#c0392b"
BLUE = "#2e6da4"
PURPLE = "#7d3c98"
GREEN = "#2a7d2a"
GRAY = "#666"
DARK = "#333"
LBLUE = "#eef4fa"
LRED = "#fdf0ee"
LGREEN = "#f0f7f0"
LGRAY = "#f3f1ec"

FIGS = {}


def fig(name):
    def deco(fn):
        FIGS[name] = fn
        return fn
    return deco


def arrow_def(pid, color=DARK, scale=1.0):
    w = 7 * scale
    h = 5 * scale
    return (f'<defs><marker id="{pid}" viewBox="0 0 10 10" refX="9" '
            f'refY="5" markerWidth="{w}" markerHeight="{h}" '
            f'orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{color}"/>'
            f'</marker></defs>')


def marker2(pid, color):
    return (f'<marker id="{pid}" viewBox="0 0 10 10" refX="9" refY="5" '
            f'markerWidth="7" markerHeight="5" '
            f'orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>')


def line(x1, y1, x2, y2, color=DARK, w=1.5, dash=None, marker=None,
         opacity=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = f' marker-end="url(#{marker})"' if marker else ""
    o = f' opacity="{opacity}"' if opacity else ""
    return (f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" '
            f'y2="{fmt(y2)}" stroke="{color}" stroke-width="{w}"{d}{m}{o}/>')


def txt(x, y, s, fs=12, anchor="middle", color=DARK, cls=None, w=None):
    c = f' class="{cls}"' if cls else ""
    fw = f' font-weight="bold"' if w == "b" else ""
    return (f'<text x="{fmt(x)}" y="{fmt(y)}" font-size="{fs}" '
            f'text-anchor="{anchor}" fill="{color}"{c}{fw}>{s}</text>')


def rect(x, y, w, h, fill, stroke=DARK, sw=1.2, rx=6, dash=None,
         opacity=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' opacity="{opacity}"' if opacity else ""
    return (f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" '
            f'height="{fmt(h)}" rx="{rx}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{sw}"{d}{o}/>')


def circ(x, y, r, fill, stroke="none", sw=1.2, dash=None, opacity=None):
    s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke != "none" \
        else ' stroke="none"'
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' opacity="{opacity}"' if opacity else ""
    return (f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" '
            f'fill="{fill}"{s}{d}{o}/>')


def star(x, y, r, color=RED):
    pts = []
    for k in range(10):
        rr = r if k % 2 == 0 else r * 0.42
        a = -math.pi / 2 + k * math.pi / 5
        pts.append(f"{fmt(x + rr * math.cos(a))},"
                   f"{fmt(y + rr * math.sin(a))}")
    return f'<polygon points="{" ".join(pts)}" fill="{color}"/>'


# ======================================================================
# GF(16) minimal kit (self-contained copy from ch07_figs.py)
MODP = 0b10011
EE = 4
QF = 16


def gmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & (1 << EE):
            a ^= MODP
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
    return gpow(a, QF - 2)


def interp_gf(pts):
    res = [0] * len(pts)

    def padd(f, g):
        return [a ^ b for a, b in zip(f + [0] * (len(g) - len(f)),
                                      g + [0] * (len(f) - len(g)))]

    def pmul(f, g):
        r = [0] * (len(f) + len(g) - 1)
        for i, a in enumerate(f):
            for j, b in enumerate(g):
                r[i + j] ^= gmul(a, b)
        return r
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
        sc = gmul(yi, ginv(den))
        num = [gmul(c, sc) for c in num]
        res = padd(res, num) if res else num
    return res


def peval_gf(f, x):
    r = 0
    for a in reversed(f):
        r = gmul(r, x) ^ a
    return r


ANCH = [1, 2, 4]
P_SET = [w for w in range(QF) if w not in ANCH]
SIGMA = (1, 0, 1)
SIGMA2 = (1, 0, 0)
QS = interp_gf(list(zip(ANCH, SIGMA)))
QS2 = interp_gf(list(zip(ANCH, SIGMA2)))
QVALS = {p: peval_gf(QS, p) for p in P_SET}
QVALS2 = {p: peval_gf(QS2, p) for p in P_SET}


# ======================================================================
@fig("fig_cvp2d")
def f_cvp2d():
    W, H = 460, 300
    b1 = (52, 10)
    b2 = (18, 46)
    ox, oy = 88, 216
    pts = []
    for i in range(-4, 9):
        for j in range(-6, 7):
            x = ox + i * b1[0] + j * b2[0]
            y = oy - (i * b1[1] + j * b2[1])
            if 14 < x < W - 14 and 14 < y < H - 40:
                pts.append((x, y, i, j))
    tx, ty = 285, 120
    best = min(pts, key=lambda p: (p[0] - tx) ** 2 + (p[1] - ty) ** 2)
    r = math.hypot(best[0] - tx, best[1] - ty)
    o = [arrow_def("c7cvp_arr", DARK)]
    o.append(circ(tx, ty, r, LRED, RED, 1.4, dash="5 4", opacity=0.9))
    for (x, y, i, j) in pts:
        col = BLUE if not (i == 0 and j == 0) else DARK
        o.append(circ(x, y, 3.4, col))
    o.append(line(ox, oy, ox + b1[0] * 0.94, oy - b1[1] * 0.94, RED, 2.4,
                  marker="c7cvp_arr"))
    o.append(line(ox, oy, ox + b2[0] * 0.9, oy - b2[1] * 0.9, RED, 2.4,
                  marker="c7cvp_arr"))
    o.append(txt(ox + b1[0] + 4, oy - b1[1] + 14, "b&#8321;", 13, "start",
                 RED, cls="vb"))
    o.append(txt(ox + b2[0] - 14, oy - b2[1] - 6, "b&#8322;", 13, "end",
                 RED, cls="vb"))
    o.append(txt(ox - 10, oy + 16, "0", 12, "middle", DARK))
    o.append(star(tx, ty, 8, RED))
    o.append(txt(tx + 12, ty - 8, "t", 14, "start", RED, cls="vb"))
    o.append(line(tx, ty, best[0], best[1], GREEN, 2.2))
    o.append(circ(best[0], best[1], 5.4, "none", GREEN, 2.4))
    o.append(txt(best[0] - 8, best[1] + 18, "v&#9733;", 13, "middle",
                 GREEN, cls="v"))
    mx, my = (tx + best[0]) / 2, (ty + best[1]) / 2
    o.append(txt(mx + 10, my - 6, "dist&#8322;(t, L)", 12, "start", GREEN,
                 cls="v"))
    o.append(txt(16, H - 12,
                 "lattice L(B) = {i&#183;b&#8321; + j&#183;b&#8322; : "
                 "i, j &#8712; &#8484;}", 12, "start", BLUE))
    return svg_wrap("\n".join(o), W, H,
                    "Closest vector problem in the plane")


@fig("fig_landscape")
def f_landscape():
    W, H = 640, 262
    y0, y1 = 96, 148
    xa, xb = 30, W - 22
    segs = [
        (0.00, 0.34, LRED, RED, "NP-hard"),
        (0.34, 0.565, LGRAY, GRAY, "open"),
        (0.565, 0.80, LBLUE, BLUE, "NP&#8745;coNP"),
        (0.80, 1.00, LGREEN, GREEN, "poly-time"),
    ]
    o = [arrow_def("c7land_arr", DARK)]
    for (s, e, fill, col, lab) in segs:
        x = xa + s * (xb - xa)
        wdt = (e - s) * (xb - xa)
        o.append(rect(x, y0, wdt, y1 - y0, fill, col, 1.2, rx=0))
        o.append(txt(x + wdt / 2, y1 + 18, lab, 12, "middle", col, w="b"))
    o.append(line(xa, y1 + 34, xb + 8, y1 + 34, DARK, 1.4,
                  marker="c7land_arr"))
    o.append(txt(xb + 4, y1 + 52, "approximation factor &#947;(n)", 12,
                 "end", DARK))
    marks = [
        (0.03, "1", "exact: NP-hard", "vEB 1981", RED, -58),
        (0.115, "O(1)", "every constant", "ABSS 1997", RED, -30),
        (0.21, "n&#8201;<tspan baseline-shift=\"super\" "
         "font-size=\"75%\">a/log log n</tspan>", "almost-poly",
         "DKRS 2003", RED, -58),
        (0.315, "n^&#8201;1/400", "THIS CHAPTER", "deterministic, no PCP",
         RED, -30),
        (0.565, "&#8730;n&#8201;/&#8201;C", "barrier: in NP&#8745;coNP",
         "Aharonov&#8211;Regev 2005", BLUE, -58),
        (0.90, "2&#8201;<tspan baseline-shift=\"super\" "
         "font-size=\"75%\">n log log n/log n</tspan>",
         "achieved in poly time", "lattice reduction", GREEN, -30),
    ]
    for (px, top, lab1, lab2, col, dy) in marks:
        x = xa + px * (xb - xa)
        o.append(line(x, y0 + dy + 22, x, y0, col, 1.3, dash="3 3"))
        o.append(circ(x, y0 + 4, 3.2, col))
        o.append(txt(x, y0 + dy, top, 12, "middle", col, w="b"))
        o.append(txt(x, y0 + dy + 12, lab1, 9.5, "middle", col))
        o.append(txt(x, y0 + dy + 21.5, lab2, 9.5, "middle", GRAY))
    o.append(txt(xa + 0.44 * (xb - xa), y0 + 30,
                 "largest NP-hard exponent?", 10.5, "middle", GRAY))
    o.append(txt(xa + 0.44 * (xb - xa), y0 + 42,
                 "c &#8712; [1/400, 1/2)", 10.5, "middle", GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "The CVP approximation landscape")


@fig("fig_pipeline")
def f_pipeline():
    W, H = 640, 300
    o = ['<defs>' + marker2("c7pipe_a", DARK) +
         marker2("c7pipe_g", GREEN) + marker2("c7pipe_r", PURPLE) +
         '</defs>']
    boxes = [
        (18, 40, 130, 64, LGRAY, DARK, ["3SAT formula", "&#966;"],
         "m vars, &#8467; clauses"),
        (188, 40, 150, 64, LBLUE, BLUE,
         ["evaluation tables", "+ moment rules"], "&#167;6&#8211;7, over "
         "F&#8339; = F&#8322;&#8337;"),
        (378, 40, 130, 64, LBLUE, BLUE, ["binary affine", "system Hx = b"],
         "M &#8776; N&#8309;&#8304;&#8321; unknowns"),
        (18, 196, 210, 64, LRED, RED,
         ["parity-lift lattice &#923;&#8341;", "+ target u"],
         "Lemma 6: explicit basis"),
        (288, 196, 220, 64, LRED, RED,
         ["GapCVP instance", "(B&#8341;, u, r)"],
         "n = M, factor n^&#8201;1/400"),
    ]
    for (x, y, w, h, fill, col, lines_, sub) in boxes:
        o.append(rect(x, y, w, h, fill, col, 1.5))
        for k, s in enumerate(lines_):
            o.append(txt(x + w / 2, y + 26 + 18 * k, s, 12.5, "middle",
                         col, w="b"))
        o.append(txt(x + w / 2, y + h + 15, sub, 10.5, "middle", GRAY))
    o.append(line(148, 72, 184, 72, DARK, 1.6, marker="c7pipe_a"))
    o.append(txt(166, 62, "&#167;6", 10.5, "middle", GRAY))
    o.append(line(338, 72, 374, 72, DARK, 1.6, marker="c7pipe_a"))
    o.append(txt(356, 62, "&#167;7", 10.5, "middle", GRAY))
    o.append(line(420, 104, 190, 194, DARK, 1.6, marker="c7pipe_a"))
    o.append(txt(330, 148, "mod-2 lift (&#167;4)", 11, "middle", GRAY))
    o.append(line(228, 228, 284, 228, DARK, 1.6, marker="c7pipe_a"))
    o.append(txt(256, 218, "Lemma 7", 10.5, "middle", GRAY))
    o.append(line(80, 106, 80, 168, GREEN, 2.0, marker="c7pipe_g"))
    o.append(txt(88, 128, "completeness (Lemma 8):", 11, "start", GREEN))
    o.append(txt(88, 141, "&#966; satisfiable &#8658; weight &#8804; R "
                 "&#8658; dist &#8804; r", 11, "start", GREEN))
    xr = 560
    o.append(f'<path d="M{xr},196 C {xr + 46},120 {xr + 30},58 512,50" '
             f'fill="none" stroke="{PURPLE}" stroke-width="2" '
             f'stroke-dasharray="6 4" marker-end="url(#c7pipe_r)"/>')
    o.append(txt(xr + 12, 150, "soundness", 11, "start", PURPLE))
    o.append(txt(xr + 12, 163, "(Prop 14):", 11, "start", PURPLE))
    o.append(txt(514, 34, "any solution of weight &#8804; 4M^&#8201;1/200"
                 "&#8201;R &#8658; &#966; satisfiable", 10.5, "end",
                 PURPLE))
    return svg_wrap("\n".join(o), W, H,
                    "The reduction pipeline from 3SAT to CVP")


@fig("fig_basis")
def f_basis():
    W, H = 560, 240
    o = []
    for panel, (bb1, bb2, col, lab) in enumerate([
            (((44, 0)), (10, 38), GREEN, "basis 1: short, near-orthogonal"),
            (((44, 0)), (98, 38), RED,
             "basis 2: long, skewed &#8212; same lattice")]):
        ox = 80 + panel * 280
        oy = 168
        pts = set()
        for i in range(-6, 8):
            for j in range(-4, 6):
                x = ox + i * 44 + j * 10
                y = oy - (i * 0 + j * 38)
                if panel * 280 + 12 < x < panel * 280 + 268 and \
                        20 < y < H - 36:
                    pts.add((x, y))
        for (x, y) in sorted(pts):
            o.append(circ(x, y, 3.0, BLUE))
        pid = f"c7bas_arr{panel}"
        o.append(f'<defs>{marker2(pid, col)}</defs>')
        b1, b2 = bb1, bb2
        px = [(ox, oy), (ox + b1[0], oy - b1[1]),
              (ox + b1[0] + b2[0], oy - b1[1] - b2[1]),
              (ox + b2[0], oy - b2[1])]
        pp = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in px)
        o.append(f'<polygon points="{pp}" fill="{col}" opacity="0.15"/>')
        o.append(line(ox, oy, ox + b1[0] * 0.93, oy - b1[1] * 0.93, col,
                      2.2, marker=pid))
        o.append(line(ox, oy, ox + b2[0] * 0.93, oy - b2[1] * 0.93, col,
                      2.2, marker=pid))
        o.append(txt(panel * 280 + 140, H - 14, lab, 11.5, "middle", col))
        o.append(circ(ox, oy, 3.0, DARK))
    o.append(line(280, 18, 280, H - 30, GRAY, 1, dash="4 4"))
    return svg_wrap("\n".join(o), W, H,
                    "One lattice, two bases")


@fig("fig_gap")
def f_gap():
    W, H = 620, 300
    o = []
    for panel, kind in enumerate(["YES", "NO"]):
        ox = 20 + panel * 310
        cx, cy = ox + 145, 138
        r1, r2 = 34, 92
        o.append(circ(cx, cy, r2, LGRAY, GRAY, 1.3, dash="6 4",
                      opacity=0.85))
        o.append(circ(cx, cy, r1, LGREEN if kind == "YES" else LRED,
                      GREEN if kind == "YES" else RED, 1.6))
        b1 = (57, 8)
        b2 = (16, 49)
        if kind == "YES":
            latoff = (-12, 20)
        else:
            latoff = (-25, -32)
        for i in range(-4, 6):
            for j in range(-4, 6):
                x = cx + latoff[0] + i * b1[0] + j * b2[0]
                y = cy + latoff[1] - (i * b1[1] + j * b2[1])
                if ox + 8 < x < ox + 286 and 22 < y < 252:
                    d = math.hypot(x - cx, y - cy)
                    if kind == "NO" and d <= r2 + 6:
                        continue
                    o.append(circ(x, y, 3.0, BLUE))
        o.append(star(cx, cy, 7, DARK))
        o.append(txt(cx + 10, cy - 6, "t", 13, "start", DARK, cls="vb"))
        o.append(line(cx, cy, cx + r1 * 0.7071, cy - r1 * 0.7071,
                      GREEN if kind == "YES" else RED, 1.6))
        o.append(txt(cx + r1 * 0.75 + 4, cy - r1 * 0.75, "r", 12, "start",
                     GREEN if kind == "YES" else RED, cls="v"))
        o.append(line(cx, cy, cx - r2 * 0.5, cy + r2 * 0.866, GRAY, 1.4))
        o.append(txt(cx - r2 * 0.56 - 4, cy + r2 * 0.9, "&#947;(n)&#183;r",
                     12, "end", GRAY, cls="v"))
        lab = ("YES: some lattice point within r"
               if kind == "YES" else
               "NO: every lattice point beyond &#947;(n)&#183;r")
        o.append(txt(ox + 145, 278, lab, 12, "middle",
                     GREEN if kind == "YES" else RED, w="b"))
    o.append(txt(310, 24, "the annulus r &lt; d &#8804; &#947;r: "
                 "either answer allowed", 11, "middle", GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "The promise problem GapCVP")


@fig("fig_3sat")
def f_3sat():
    W, H = 620, 232
    o = []
    o.append(txt(18, 34, "&#966; =", 16, "start", DARK, cls="v"))
    lits = [[("v&#8321;", False), ("v&#8322;", False), ("v&#8323;", True)],
            [("v&#8321;", True), ("v&#8322;", True), ("v&#8323;", False)]]
    x = 62
    for ci, cl in enumerate(lits):
        x0 = x
        o.append(txt(x, 34, "(", 17, "start", DARK))
        x += 12
        for li, (v, neg) in enumerate(cl):
            s = ("&#172;" if neg else "") + v
            col = RED if neg else BLUE
            o.append(txt(x, 34, s, 14, "start", col, cls="v"))
            x += 34 if neg else 24
            if li < 2:
                o.append(txt(x, 34, "&#8744;", 13, "start", DARK))
                x += 22
        o.append(txt(x, 34, ")", 17, "start", DARK))
        x += 16
        o.append(rect(x0 - 4, 14, x - x0, 30, "none", GRAY, 1.0,
                      dash="3 3"))
        o.append(txt((x0 + x) / 2 - 2, 62, f"clause C{chr(0x2081 + ci)}",
                     11, "middle", GRAY))
        if ci == 0:
            o.append(txt(x + 2, 34, "&#8743;", 14, "start", DARK))
            x += 26
    o.append(txt(18, 96, "a clause is satisfied when at least one of its "
                 "three literals is true", 12, "start", DARK))
    rows = [("assignment &#963; = (1, 0, 1)", ""),
            ("C&#8321; = v&#8321;&#8744;v&#8322;&#8744;&#172;v&#8323;:",
             "v&#8321; = 1 &#10003;"),
            ("C&#8322; = &#172;v&#8321;&#8744;&#172;v&#8322;&#8744;"
             "v&#8323;:", "v&#8323; = 1 &#10003;")]
    o.append(rect(14, 112, 340, 92, LGREEN, GREEN, 1.2))
    for k, (a, b) in enumerate(rows):
        o.append(txt(26, 136 + k * 24, a, 12, "start",
                     GREEN if k == 0 else DARK,
                     w="b" if k == 0 else None))
        o.append(txt(230, 136 + k * 24, b, 12, "start", GREEN))
    o.append(rect(376, 112, 228, 92, LGRAY, GRAY, 1.2))
    o.append(txt(490, 136, "3SAT asks:", 12, "middle", DARK, w="b"))
    o.append(txt(490, 158, "does ANY of the 2^m assignments", 11.5,
                 "middle", DARK))
    o.append(txt(490, 176, "satisfy ALL &#8467; clauses?", 11.5, "middle",
                 DARK))
    o.append(txt(18, 224, "toy instance used throughout: m = 3, "
                 "&#8467; = 2; &#963; = (1,0,1) satisfies &#966;", 11.5,
                 "start", GRAY))
    return svg_wrap("\n".join(o), W, H, "Anatomy of a 3SAT formula")


@fig("fig_reduction")
def f_reduction():
    W, H = 620, 268
    o = ['<defs>' + marker2("c7red_a", DARK) + marker2("c7red_g", GREEN) +
         marker2("c7red_r", RED) + '</defs>']
    o.append(f'<ellipse cx="130" cy="130" rx="105" ry="96" fill="{LGRAY}"'
             f' stroke="{DARK}" stroke-width="1.3"/>')
    o.append(f'<ellipse cx="470" cy="130" rx="118" ry="96" '
             f'fill="{LBLUE}" stroke="{BLUE}" stroke-width="1.3"/>')
    o.append(line(60, 130, 200, 130, GRAY, 1.1, dash="4 3"))
    o.append(line(360, 130, 582, 130, GRAY, 1.1, dash="4 3"))
    o.append(txt(130, 48, "3SAT instances", 12.5, "middle", DARK, w="b"))
    o.append(txt(470, 48, "GapCVP instances", 12.5, "middle", BLUE,
                 w="b"))
    o.append(txt(130, 100, "satisfiable &#966;", 12, "middle", GREEN))
    o.append(txt(130, 172, "unsatisfiable &#966;", 12, "middle", RED))
    o.append(txt(470, 88, "YES: dist&#8322;(t, L) &#8804; r", 12,
                 "middle", GREEN))
    o.append(txt(470, 166, "NO: dist&#8322;(t, L) &gt; n^&#8201;1/400"
                 "&#8201;r", 12, "middle", RED))
    o.append(txt(470, 196, "(gap: nothing in between)", 10.5, "middle",
                 GRAY))
    o.append(line(214, 96, 366, 84, GREEN, 2.0, marker="c7red_g"))
    o.append(line(214, 168, 366, 162, RED, 2.0, marker="c7red_r"))
    o.append(txt(292, 74, "&#966; &#8614; (B, t, r)", 12, "middle", DARK,
                 cls="v"))
    o.append(txt(292, 116, "deterministic", 10.5, "middle", GRAY))
    o.append(txt(292, 128, "polynomial time", 10.5, "middle", GRAY))
    o.append(rect(96, 226, 430, 34, LRED, RED, 1.2))
    o.append(txt(311, 248, "an n^&#8201;1/400-approximator for CVP would "
                 "decide 3SAT &#8658; P = NP", 12, "middle", RED))
    return svg_wrap("\n".join(o), W, H,
                    "A gap-producing many-one reduction")


@fig("fig_gapamp")
def f_gapamp():
    W, H = 620, 190
    xa, xb = 60, 580
    y = 66
    scale = (xb - xa) / 10.0

    def X(v):
        return xa + v * scale
    o = [arrow_def("c7amp_arr", DARK)]
    o.append(line(xa - 10, y, xb + 14, y, DARK, 1.5, marker="c7amp_arr"))
    o.append(txt(xb + 6, y - 16, "distance d = dist&#8322;(t, L)", 11.5,
                 "end", DARK))
    for (v, lab) in [(1.6, "r"), (6.4, "&#947;r")]:
        o.append(line(X(v), y - 34, X(v), y + 6, GRAY, 1.2, dash="4 3"))
        o.append(txt(X(v), y - 40, lab, 13, "middle", DARK, cls="v"))
    o.append(rect(X(0), y - 8, X(1.6) - X(0), 16, LGREEN, GREEN, 1.0,
                  rx=3))
    o.append(rect(X(6.4), y - 8, X(10) - X(6.4), 16, LRED, RED, 1.0,
                  rx=3))
    o.append(txt(X(0.8), y + 26, "YES region", 11, "middle", GREEN))
    o.append(txt(X(8.2), y + 26, "NO region", 11, "middle", RED))
    o.append(txt(X(4.0), y + 26, "promise: never happens", 10.5, "middle",
                 GRAY))
    yy = 130
    o.append(txt(xa - 44, yy + 4, "output", 11, "start", PURPLE, w="b"))
    o.append(line(X(1.6), yy - 8, X(1.6), yy + 12, GRAY, 1, dash="3 3"))
    o.append(rect(X(0), yy - 6, X(6.4) - X(0), 12, "none", PURPLE, 1.4,
                  rx=3, dash="5 3"))
    o.append(txt(X(3.2), yy + 26,
                 "a &#947;-approximation d&#770; satisfies d &#8804; "
                 "d&#770; &#8804; &#947;d: on YES it lands &#8804; "
                 "&#947;r, on NO it lands &gt; &#947;r", 11, "middle",
                 PURPLE))
    o.append(txt(X(3.2), yy + 44,
                 "so comparing d&#770; with &#947;r decides the promise "
                 "problem &#8212; approximation &#8805; gap decision", 11,
                 "middle", DARK))
    return svg_wrap("\n".join(o), W, H,
                    "Why gap hardness means approximation hardness")


@fig("fig_exponents")
def f_exponents():
    W, H = 640, 258
    o = ['<defs>' + marker2("c7exp_a", DARK) + '</defs>']
    boxes = [
        (16, 30, 120, 56, "formula &#966;", "size s, N &#8776; s", LGRAY,
         DARK),
        (176, 30, 140, 56, "field F&#8339;", "q &#8776; N&#178;&#8304;"
         "&#8304;", LBLUE, BLUE),
        (356, 30, 150, 56, "dimension", "M &#8804; 40N&#8308;&#8304;"
         "&#185;", LBLUE, BLUE),
        (16, 156, 210, 56, "binary gap factor", "4&#8201;M^&#8201;1/200",
         LGREEN, GREEN),
        (286, 156, 200, 56, "Euclidean factor", "M^&#8201;1/400 = "
         "n^&#8201;1/400", LRED, RED),
    ]
    for (x, y, w, h, l1, l2, fill, col) in boxes:
        o.append(rect(x, y, w, h, fill, col, 1.4))
        o.append(txt(x + w / 2, y + 23, l1, 12.5, "middle", col, w="b"))
        o.append(txt(x + w / 2, y + 43, l2, 12.5, "middle", col))
    o.append(line(136, 58, 172, 58, DARK, 1.6, marker="c7exp_a"))
    o.append(txt(154, 48, "&#167;7", 10, "middle", GRAY))
    o.append(line(316, 58, 352, 58, DARK, 1.6, marker="c7exp_a"))
    o.append(txt(334, 48, "M &#8776; 9N&#183;q&#178;", 9.5, "middle",
                 GRAY))
    o.append(line(420, 86, 140, 152, DARK, 1.6, marker="c7exp_a"))
    o.append(txt(300, 118, "soundness threshold (Prop 14)", 10.5,
                 "middle", GRAY))
    o.append(line(226, 184, 282, 184, DARK, 1.6, marker="c7exp_a"))
    o.append(txt(254, 174, "&#8730; of Lemma 7", 10, "middle", GRAY))
    o.append(txt(20, 236, "q &#8776; N&#178;&#8304;&#8304; makes every "
                 "error term (Markov losses, discriminant zeros, valuation "
                 "windows) vanish against |P| &#8776; q", 10.5, "start",
                 GRAY))
    o.append(txt(20, 254, "M &#8776; N&#8308;&#8304;&#185; &#8776; "
                 "q&#178; &#8658; M^&#8201;1/200 &#8776; q^&#8201;1/100 "
                 "&#8776; N&#178;: a huge gap in absolute terms, a tiny "
                 "power of the dimension", 10.5, "start", GRAY))
    return svg_wrap("\n".join(o), W, H + 10,
                    "Where the exponent 1/400 comes from")


@fig("fig_hamming")
def f_hamming():
    W, H = 500, 300
    CODE = {0b000, 0b011, 0b101, 0b110}
    U = 0b100

    def proj(v):
        x = (v & 1)
        y = (v >> 1) & 1
        z = (v >> 2) & 1
        px = 120 + x * 170 + z * 78
        py = 236 - y * 150 - z * 62
        return px, py
    o = []
    edges = []
    for v in range(8):
        for k in range(3):
            u2 = v ^ (1 << k)
            if u2 > v:
                edges.append((v, u2))
    for (v, u2) in edges:
        x1, y1 = proj(v)
        x2, y2 = proj(u2)
        both = v in CODE and u2 in CODE
        o.append(line(x1, y1, x2, y2, GRAY, 1.1,
                      dash=None if True else None, opacity=0.7))
    for v in range(8):
        x, y = proj(v)
        incode = v in CODE
        col = GREEN if incode else "#bbb"
        if v == U:
            col = RED
        o.append(circ(x, y, 9 if incode or v == U else 6, col,
                      DARK if incode else "none", 1.2))
        lab = format(v, "03b")
        dy = -14 if y < 150 else 24
        o.append(txt(x, y + dy, lab, 12, "middle",
                     GREEN if incode else (RED if v == U else GRAY),
                     w="b" if incode or v == U else None))
    ux, uy = proj(U)
    for cw in (0b000, 0b101, 0b110):
        cx2, cy2 = proj(cw)
        o.append(line(ux, uy, cx2, cy2, RED, 1.6, dash="5 3",
                      opacity=0.85))
    o.append(txt(24, 286, "d&#8344;(u, C) = min wt(u &#8722; c) = "
                 "W(H, b) = 1  for H = [1 1 1], b = 1", 12, "start",
                 DARK))
    o.append(txt(24, 310, "code C = ker H = even-weight words (green)",
                 12, "start", GREEN))
    o.append(txt(24, 328, "received u = 100 (red): dashed = three "
                 "codewords at Hamming distance 1", 12, "start", RED))
    return svg_wrap("\n".join(o), W, H + 44,
                    "Nearest codeword in the 3-cube")


@fig("fig_paritylift")
def f_paritylift():
    W, H = 560, 330
    ox, oy = 120, 258
    s = 62

    def P(x, y):
        return ox + x * s, oy - y * s
    o = []
    o.append(txt(20, 26, "&#923;&#8341; = {z &#8712; &#8484;&#178; : "
                 "z&#8321; + z&#8322; even} = lift of C = {00, 11}", 12.5,
                 "start", DARK))
    for x in range(-1, 6):
        for y in range(-2, 3):
            px, py = P(x, y)
            if not (30 < px < W - 90 and 40 < py < H - 20):
                continue
            even = (x + y) % 2 == 0
            if even:
                base = (x % 2, y % 2)
                col = BLUE if base == (0, 0) else PURPLE
                o.append(circ(px, py, 6, col))
            else:
                o.append(circ(px, py, 4.5, "none", "#bbb", 1.3))
    ux, uy = P(1, 0)
    o.append(star(ux, uy, 9, RED))
    o.append(txt(ux + 3, uy - 13, "u = (1,0)", 12, "start", RED, cls="v"))
    r0 = 8
    o.append(f'<circle cx="{fmt(ux)}" cy="{fmt(uy)}" r="{r0}" '
             f'fill="none" stroke="{RED}" stroke-width="1.6" '
             f'stroke-dasharray="6 4">'
             f'<animate attributeName="r" values="6;{s};6" dur="4s" '
             f'repeatCount="indefinite"/>'
             f'<animate attributeName="opacity" values="0.9;0.55;0.9" '
             f'dur="4s" repeatCount="indefinite"/></circle>')
    o.append(circ(ux, uy, s, "none", RED, 1.0, dash="2 4", opacity=0.6))
    for (nx, ny) in [(0, 0), (2, 0), (1, 1), (1, -1)]:
        px, py = P(nx, ny)
        o.append(circ(px, py, 9.5, "none", GREEN, 2.0))
    o.append(txt(P(2, 1)[0] + 16, P(2, 1)[1],
                 "nearest lattice points:", 11.5, "start", GREEN))
    o.append(txt(P(2, 1)[0] + 16, P(2, 1)[1] + 15,
                 "dist&#8322;&#178; = 1 = W(H, b)", 11.5, "start", GREEN))
    o.append(txt(20, H - 34, "solid dots: lattice &#923;&#8341; "
                 "(blue = coset of 00, purple = coset of 11);  "
                 "open: odd parity, NOT in &#923;&#8341;", 11, "start",
                 GRAY))
    o.append(txt(20, H - 16, "animation: a circle grows from u until it "
                 "first touches the lattice &#8212; at radius 1", 11,
                 "start", GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "Parity-lift lattice for H = [1 1], b = 1")


@fig("fig_lemma7")
def f_lemma7():
    W, H = 620, 272
    coords = [0, -2, 1, 0, 3, -1, 0]
    xa = 90
    bw = 56
    y0 = 150
    sc = 26
    o = []
    o.append(txt(18, 28, "integer error vector u &#8722; &#955; "
                 "(coordinates):", 12.5, "start", DARK))
    for i, v in enumerate(coords):
        x = xa + i * (bw + 14)
        h = abs(v) * sc
        odd = v % 2 != 0
        col = RED if odd else "#9db6c9"
        if v != 0:
            o.append(rect(x, y0 - (h if v > 0 else 0), bw, h, col,
                          DARK, 1.0, rx=3, opacity=0.9))
        o.append(txt(x + bw / 2, y0 + 18 + (0 if v >= 0 else abs(v) * sc),
                     str(v), 12, "middle", DARK))
        o.append(txt(x + bw / 2, 62,
                     "odd" if odd else "even", 10.5, "middle",
                     RED if odd else GRAY))
    o.append(line(xa - 16, y0, xa + 7 * (bw + 14), y0, DARK, 1.3))
    o.append(line(xa - 16, y0 - sc, xa + 7 * (bw + 14), y0 - sc, RED,
                  1.1, dash="5 3"))
    o.append(txt(xa + 7 * (bw + 14) + 2, y0 - sc + 4, "1", 11, "start",
                 RED))
    o.append(txt(18, H - 34, "x = (u &#8722; &#955;) mod 2 selects the "
                 "odd coordinates (red). Each odd coordinate has "
                 "|u&#11388; &#8722; &#955;&#11388;| &#8805; 1,", 11.5,
                 "start", DARK))
    o.append(txt(18, H - 16, "so &#8721; |u&#11388; &#8722; "
                 "&#955;&#11388;|&#8348; &#8805; wt(x) &#8805; W(H, b) "
                 "&#8212; and the lift of a minimum-weight solution "
                 "achieves equality with &#177;1 entries.", 11.5, "start",
                 DARK))
    return svg_wrap("\n".join(o), W, H,
                    "Proof geometry of Lemma 7")


@fig("fig_h111")
def f_h111():
    W, H = 560, 320
    ex = (66, 0)
    ey = (26, -22)
    ez = (0, -60)
    ox, oy = 150, 220

    def proj(x, y, z):
        return (ox + x * ex[0] + y * ey[0] + z * ez[0],
                oy + x * ex[1] + y * ey[1] + z * ez[1])
    pts = []
    for x in range(-1, 4):
        for y in range(-1, 3):
            for z in range(-1, 3):
                if (x + y + z) % 2 == 0:
                    px, py = proj(x, y, z)
                    if 24 < px < W - 120 and 30 < py < H - 40:
                        pts.append((x, y, z, px, py))
    o = []
    for x in range(0, 2):
        for y in range(0, 2):
            p1 = proj(x, y, -1)
            p2 = proj(x, y, 2)
            o.append(line(*p1, *p2, "#ccc", 0.8, opacity=0.6))
    pts.sort(key=lambda t: -(t[1] + t[2]))
    for (x, y, z, px, py) in pts:
        depth = (y + z) / 4
        rr = 4.6 - 1.2 * depth
        o.append(circ(px, py, rr, BLUE, opacity=None))
    ux, uy = proj(1, 0, 0)
    o.append(star(ux, uy, 9, RED))
    o.append(txt(ux + 26, uy + 24, "u = (1,0,0)", 12, "start", RED,
                 cls="v"))
    for lam in [(0, 0, 0), (1, 1, 0), (1, 0, 1), (2, 0, 0), (1, -1, 0)]:
        lx, ly = proj(*lam)
        o.append(circ(lx, ly, 8.5, "none", GREEN, 2.0))
        o.append(line(ux, uy, lx, ly, GREEN, 1.4, dash="4 3",
                      opacity=0.8))
    o.append(txt(24, 30, "&#923;&#8341; = {z &#8712; &#8484;&#179; : "
                 "z&#8321;+z&#8322;+z&#8323; even},  H = [1 1 1], b = 1",
                 12.5, "start", DARK))
    o.append(txt(24, H - 44, "five of the lattice points at distance "
                 "exactly 1 from u (circled); no lattice point is closer:",
                 11.5, "start", GREEN))
    o.append(txt(24, H - 26, "dist&#8322;(u, &#923;&#8341;)&#178; = 1 = "
                 "min{wt(x) : x&#8321;+x&#8322;+x&#8323; = 1 over "
                 "F&#8322;}", 11.5, "start", GREEN))
    return svg_wrap("\n".join(o), W, H,
                    "The worked example in three dimensions")


# ======================================================================
# Part B: fields, encoding, constraints, reconstruction
@fig("fig_gf16")
def f_gf16():
    W, H = 640, 268
    o = []
    o.append(txt(18, 24, "F&#8321;&#8326; = F&#8322;[x]/(x&#8308;+x+1): "
                 "elements are 4-bit words; g = x generates the 15 "
                 "nonzero elements", 12, "start", DARK))
    logs = {1: 0}
    v = 1
    for k in range(1, 15):
        v = gmul(v, 2)
        logs[v] = k
    for w in range(16):
        cx = 24 + (w % 8) * 76
        cy = 44 + (w // 8) * 76
        o.append(rect(cx, cy, 66, 62, LBLUE if w else LGRAY, BLUE if w
                      else GRAY, 1.1, rx=5))
        o.append(txt(cx + 33, cy + 17, format(w, "04b"), 12, "middle",
                     DARK, w="b"))
        if w == 0:
            sub = "0"
        elif w == 1:
            sub = "1 = g&#8304;"
        else:
            sub = f"g^{logs[w]}"
        o.append(txt(cx + 33, cy + 34, sub, 11, "middle", BLUE))
        terms = []
        for k in range(3, -1, -1):
            if (w >> k) & 1:
                terms.append("x" + ("&#179;" if k == 3 else
                                    "&#178;" if k == 2 else
                                    "" if k == 1 else ""))
        poly = "+".join(terms).replace("x+x", "x + x") if terms else "0"
        if poly.endswith("x") and w & 1:
            poly = poly[:-1] + "1"
        o.append(txt(cx + 33, cy + 51, poly, 10, "middle", GRAY))
    o.append(txt(18, 216, "addition = bitwise XOR:  0110 &#8853; 0011 = "
                 "0101   (every element is its own negative: "
                 "char = 2)", 12, "start", DARK))
    o.append(txt(18, 236, "multiplication = polynomial product mod "
                 "x&#8308;+x+1:  g&#179;&#183;g&#178; = g&#8309; = "
                 "g&#178;+g = 0110", 12, "start", DARK))
    o.append(txt(18, 256, "the toy chapter field; the real construction "
                 "uses q = 2^&#8201;1329 &#8776; N&#178;&#8304;&#8304; "
                 "at N = 100", 11, "start", GRAY))
    return svg_wrap("\n".join(o), W, H, "The field with 16 elements")


@fig("fig_rs")
def f_rs():
    W, H = 640, 300
    ax = Ax(-0.7, 15.7, -0.8, 15.8, W=W, H=232, ml=46, mr=16, mt=30,
            mb=34)
    o = []
    o.append(ax.axes(xticks=(0, 5, 10, 15), yticks=(0, 5, 10, 15),
                     xlab="p (as integer 0&#8211;15)",
                     ylab="Q&#963;(p)"))
    for p in P_SET:
        o.append(ax.dot(p, QVALS[p], 4.0, BLUE))
    for a, s in zip(ANCH, SIGMA):
        o.append(ax.dot(a, s, 4.6, RED))
        o.append(ax.text(a, s, "anchor", 9.5, anchor="middle", dy=-9,
                         color=RED))
    o.append(txt(320, 20, "the same degree-2 polynomial, all 16 "
                 "evaluations &#8212; no visible smoothness", 11.5,
                 "middle", GRAY))
    y2 = 262
    o.append(txt(46, y2 - 8, "what “degree &#8804; D” still "
                 "buys: dimension D+1, so |P| &#8722; D &#8722; 1 "
                 "independent parity checks", 11.5, "start", DARK))
    xw = 540 / 13
    for i in range(13):
        col = LGREEN if i < 3 else LRED
        stroke = GREEN if i < 3 else RED
        o.append(rect(46 + i * xw, y2, xw - 3, 20, col, stroke, 1.0,
                      rx=3))
    o.append(txt(46 + 1.5 * xw, y2 + 34, "D+1 = 3 free", 10.5, "middle",
                 GREEN))
    o.append(txt(46 + 8 * xw, y2 + 34, "|P| &#8722; 3 = 10 parity checks "
                 "(each &#8594; e binary equations)", 10.5, "middle",
                 RED))
    return svg_wrap("\n".join(o), W, H,
                    "A Reed-Solomon codeword over F16")


@fig("fig_interp")
def f_interp():
    W, H = 640, 300
    o = []
    axl = Ax(0.4, 3.9, -0.45, 1.7, W=310, H=240, ml=40, mr=10, mt=26,
             mb=34)

    def qreal(x):
        return (x - 2) ** 2
    pts = [(0.72 + i * (3.28 - 0.72) / 80, qreal(0.72 + i * (3.28 - 0.72)
                                                 / 80)) for i in range(81)]
    o.append('<g>')
    o.append(axl.axes(xticks=(1, 2, 3), yticks=(0, 1),
                      xlab="X", ylab=""))
    o.append(axl.polyline(pts, BLUE, 2.0))
    for (a, s) in [(1, 1), (2, 0), (3, 1)]:
        o.append(axl.dot(a, s, 4.5, RED))
    o.append(axl.text(1, 1, "(a&#8321;, &#963;&#8321;)", 10, dy=-8,
                      color=RED, anchor="start"))
    o.append(axl.text(2, 0, "(a&#8322;, &#963;&#8322;)", 10, dy=16,
                      color=RED, anchor="middle"))
    o.append(axl.text(3, 1, "(a&#8323;, &#963;&#8323;)", 10, dy=-8,
                      color=RED, anchor="end"))
    o.append(txt(160, 16, "real-line analogy", 11, "middle", GRAY))
    o.append('</g>')
    xg = 340
    o.append(txt(xg + 4, 30, "the actual F&#8321;&#8326; toy: "
                 "Q&#963; = interpolation of &#963; = (1,0,1)", 11.5,
                 "start", DARK))
    o.append(txt(xg + 4, 48, "at anchors a&#8321;=1, a&#8322;=g, "
                 "a&#8323;=g&#178;", 11.5, "start", DARK))
    cells = [(a, s, True) for a, s in zip(ANCH, SIGMA)]
    cells += [(p, QVALS[p], False) for p in P_SET]
    cells.sort()
    cw = 34
    for k, (p, v, isanch) in enumerate(cells):
        x = xg + 4 + (k % 8) * cw
        y = 66 + (k // 8) * 62
        o.append(rect(x, y, cw - 4, 24, LRED if isanch else LBLUE,
                      RED if isanch else BLUE, 1.0, rx=3))
        o.append(txt(x + cw / 2 - 2, y + 16, str(p), 10.5, "middle",
                     DARK))
        o.append(rect(x, y + 24, cw - 4, 24, "none", GRAY, 0.8, rx=3))
        o.append(txt(x + cw / 2 - 2, y + 40, str(v), 10.5, "middle",
                     RED if isanch else BLUE, w="b" if isanch else None))
    o.append(txt(xg + 4, 196, "top: point p; bottom: Q&#963;(p).", 10.5,
                 "start", GRAY))
    o.append(txt(xg + 4, 214, "red anchors carry the Boolean values;",
                 10.5, "start", GRAY))
    o.append(txt(xg + 4, 232, "the 13 blue columns (p &#8712; P) are "
                 "the codeword", 10.5, "start", GRAY))
    o.append(txt(xg + 4, 262, "deg Q&#963; &#8804; m &#8722; 1 = 2, "
                 "checked by machine (&#167;17)", 10.5, "start", GREEN))
    return svg_wrap("\n".join(o), W, H,
                    "A Boolean assignment as one polynomial")


@fig("fig_tables")
def f_tables():
    W, H = 640, 468
    gx, gy = 58, 46
    cw, chh = 30, 17
    o = []
    o.append(txt(18, 24, "the global evaluation table x&#8320;: one row "
                 "per field value w, one column per point p &#8712; P",
                 12, "start", DARK))
    for wv in range(16):
        y = gy + wv * chh
        o.append(txt(gx - 8, y + chh - 5, str(wv), 8.5, "end", GRAY))
    for k, p in enumerate(P_SET):
        x = gx + k * cw
        o.append(txt(x + cw / 2, gy - 6, str(p), 8.5, "middle", GRAY))
    o.append(rect(gx, gy, cw * 13, chh * 16, "none", DARK, 1.2, rx=0))
    for wv in range(1, 16):
        o.append(line(gx, gy + wv * chh, gx + cw * 13, gy + wv * chh,
                      "#ddd", 0.5))
    for k in range(1, 13):
        o.append(line(gx + k * cw, gy, gx + k * cw, gy + chh * 16,
                      "#ddd", 0.5))
    o.append(txt(gx - 34, gy + chh * 8, "w", 12, "middle", DARK,
                 cls="v"))
    o.append(txt(gx + cw * 6.5, gy + chh * 16 + 16, "p", 12, "middle",
                 DARK, cls="v"))
    for k, p in enumerate(P_SET):
        x = gx + k * cw + cw / 2
        y1 = gy + QVALS[p] * chh + chh / 2
        y2 = gy + QVALS2[p] * chh + chh / 2
        o.append(f'<circle cx="{fmt(x)}" cy="{fmt(y1)}" r="5.2" '
                 f'fill="{RED}"><animate attributeName="opacity" '
                 f'values="1;1;0;0;1" keyTimes="0;0.45;0.5;0.95;1" '
                 f'dur="8s" repeatCount="indefinite"/></circle>')
        o.append(f'<circle cx="{fmt(x)}" cy="{fmt(y2)}" r="5.2" '
                 f'fill="{BLUE}" opacity="0">'
                 f'<animate attributeName="opacity" '
                 f'values="0;0;1;1;0" keyTimes="0;0.45;0.5;0.95;1" '
                 f'dur="8s" repeatCount="indefinite"/></circle>')
    lx = gx + cw * 13 + 18
    o.append(f'<g><text x="{lx}" y="{gy + 12}" font-size="11.5" '
             f'fill="{RED}">x&#8320;,&#8202;p,&#8202;w = 1 &#8660; '
             f'w = Q&#963;(p)</text>'
             f'<animate attributeName="opacity" '
             f'values="1;1;0;0;1" keyTimes="0;0.45;0.5;0.95;1" dur="8s" '
             f'repeatCount="indefinite"/></g>')
    o.append(f'<g opacity="0"><text x="{lx}" y="{gy + 12}" '
             f'font-size="11.5" fill="{BLUE}">now &#963; &#8594; '
             f'(1,0,0)</text>'
             f'<animate attributeName="opacity" values="0;0;1;1;0" '
             f'keyTimes="0;0.45;0.5;0.95;1" dur="8s" '
             f'repeatCount="indefinite"/></g>')
    o.append(txt(lx, gy + 44, "each column is a", 10.5, "start", GRAY))
    o.append(txt(lx, gy + 58, "“fiber”: the one-hot", 10.5,
                 "start", GRAY))
    o.append(txt(lx, gy + 72, "indicator of Q&#963;(p)", 10.5, "start",
                 GRAY))
    o.append(txt(lx, gy + 100, "animation: flipping", 10.5, "start",
                 GRAY))
    o.append(txt(lx, gy + 114, "&#963;&#8323; from 1 to 0", 10.5,
                 "start", GRAY))
    o.append(txt(lx, gy + 128, "moves every dot", 10.5, "start", GRAY))
    y3 = gy + chh * 16 + 40
    o.append(txt(18, y3, "clause tables: one copy per satisfying local "
                 "assignment &#946; &#8712; B&#8201;C; only the copy "
                 "&#946; = &#963;|&#8201;I&#8201;C is populated", 11.5,
                 "start", DARK))
    for kk in range(3):
        bx = 30 + kk * 205
        by = y3 + 14
        active = kk == 0
        o.append(rect(bx, by, 150, 64, LGREEN if active else LGRAY,
                      GREEN if active else GRAY, 1.1))
        for k in range(0, 13, 1):
            x = bx + 6 + k * 11
            if active:
                yv = by + 8 + (QVALS[P_SET[k]] % 8) * 6
                o.append(circ(x, yv, 2.0, GREEN))
        lab = ("(C&#8321;, &#946; = &#963;|I) &#8212; active" if kk == 0
               else ("(C&#8321;, &#946;&#8242;) &#8212; empty" if kk == 1
                     else "(C&#8321;, &#946;&#8243;) &#8212; empty"))
        o.append(txt(bx + 75, by + 78, lab, 10.5, "middle",
                     GREEN if active else GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "Evaluation tables and their fibers")


@fig("fig_types")
def f_types():
    W, H = 640, 240
    o = ['<defs>' + marker2("c7typ_a", DARK) + '</defs>']
    o.append(rect(250, 16, 140, 40, LGRAY, DARK, 1.3))
    o.append(txt(320, 41, "table types &#920;", 12.5, "middle", DARK,
                 w="b"))
    o.append(rect(30, 100, 110, 46, LGREEN, GREEN, 1.5))
    o.append(txt(85, 120, "&#964; = 0", 12, "middle", GREEN, w="b"))
    o.append(txt(85, 137, "global table", 10.5, "middle", GREEN))
    o.append(line(285, 58, 100, 96, DARK, 1.3, marker="c7typ_a"))
    o.append(line(320, 58, 320, 96, DARK, 1.3, marker="c7typ_a"))
    o.append(line(360, 58, 520, 96, DARK, 1.3, marker="c7typ_a"))
    for g, (cx, lab, act) in enumerate([
            (320, "clause C&#8321;: 7 subtypes (C&#8321;, &#946;)", 2),
            (520, "clause C&#8322;: 7 subtypes (C&#8322;, &#946;)", 4)]):
        for k in range(7):
            x = cx - 98 + k * 28
            active = k == act
            o.append(rect(x, 100, 24, 46, LGREEN if active else "#eee",
                          GREEN if active else "#aaa", 1.2, rx=4))
            if active:
                o.append(txt(x + 12, 128, "&#9679;", 11, "middle",
                             GREEN))
        o.append(txt(cx, 166, lab, 10.5, "middle", DARK))
    o.append(txt(320, 196, "|&#920;| = 1 + &#8721;|B&#8201;C| &#8804; "
                 "1 + 8&#8467;    (toy: 1 + 7 + 7 = 15)", 12, "middle",
                 DARK))
    o.append(txt(320, 220, "green = the tables a satisfying &#963; "
                 "populates: the global one and one subtype per clause "
                 "(&#946; = &#963;|I&#8201;C)", 10.5, "middle", GREEN))
    return svg_wrap("\n".join(o), W, H, "The table-type index set")


@fig("fig_constraints")
def f_constraints():
    W, H = 640, 350
    o = []
    gx, gy = 40, 60
    cw, chh = 22, 13
    o.append(txt(18, 24, "the four constraint families, drawn on the "
                 "tables they constrain", 12.5, "start", DARK))
    o.append(txt(gx + 66, gy - 16, "global table 0", 11, "middle", DARK))
    o.append(rect(gx, gy, cw * 6, chh * 10, "none", DARK, 1.1))
    for k in range(6):
        yv = gy + ((k * 3 + 1) % 10) * chh + chh / 2
        o.append(circ(gx + k * cw + cw / 2, yv, 3.4, RED))
    o.append(rect(gx + 2 * cw + 2, gy + 1, cw - 4, chh * 10 - 2, "none",
                  PURPLE, 1.6, dash="4 2"))
    o.append(txt(gx + 66, gy + chh * 10 + 16, "(C1): every global fiber "
                 "has an odd", 10.5, "middle", PURPLE))
    o.append(txt(gx + 66, gy + chh * 10 + 30, "number of ones "
                 "(&#8721;&#8202;x = 1 in F&#8322;)", 10.5, "middle",
                 PURPLE))
    bx = 250
    o.append(txt(bx + 90, gy - 16, "clause subtypes of C", 11, "middle",
                 DARK))
    for k in range(3):
        x = bx + k * 62
        o.append(rect(x, gy, 50, chh * 10, "none",
                      GREEN if k == 0 else "#aaa", 1.1))
        if k == 0:
            for kk in range(2):
                yv = gy + ((kk * 3 + 1) % 10) * chh + chh / 2
                o.append(circ(x + 12 + kk * 22, yv, 3.2, GREEN))
    o.append(txt(bx + 90, gy + chh * 10 + 16, "(C2): fiber-by-fiber, "
                 "the subtype tables", 10.5, "middle", GREEN))
    o.append(txt(bx + 90, gy + chh * 10 + 30, "XOR to the global fiber",
                 10.5, "middle", GREEN))
    y2 = 240
    o.append(txt(18, y2, "(C3): for each table and each j &#8804; T, "
                 "the j-th power sums of the fibers,  p &#8614; "
                 "&#8721;&#8202;w&#8201;x&#964;,p,w&#8201;w&#8202;^j,",
                 10.8, "start", BLUE))
    o.append(txt(18, y2 + 15, "form a Reed&#8211;Solomon codeword of "
                 "degree &#8804; dj", 10.8, "start", BLUE))
    o.append(txt(18, y2 + 36, "(C4): for each clause subtype and each "
                 "variable i of the clause, the shifted sums with "
                 "(w &#8722; &#946;&#7522;)/(p &#8722; a&#7522;)",
                 10.8, "start", PURPLE))
    o.append(txt(18, y2 + 51, "form a codeword of degree &#8804; "
                 "(d&#8722;1)j", 10.8, "start", PURPLE))
    o.append(txt(18, y2 + 72, "all four families are LINEAR over "
                 "F&#8322; in the indicators x&#964;,p,w &#8212; "
                 "field powers w^j are known coefficients, never "
                 "unknowns", 10.4, "start", DARK, w="b"))
    o.append(txt(18, y2 + 92, "together: one explicit binary affine "
                 "system Hx = b   (equation (3); toy: 8949 equations, "
                 "M = 3120 unknowns)", 10.4, "start", DARK))
    return svg_wrap("\n".join(o), W, H,
                    "The constraint system at a glance")


@fig("fig_secant")
def f_secant():
    W, H = 640, 330
    ax = Ax(0.3, 4.2, -0.5, 2.3, W=340, H=250, ml=42, mr=10, mt=24,
            mb=34)

    def q(x):
        return (x - 2) ** 2

    o = []
    o.append(ax.axes(xticks=(1, 2, 3, 4), yticks=(0, 1, 2), xlab="X",
                     ylab="Q(X)"))
    o.append(ax.polyline([(0.55 + i * 2.95 / 90,
                           q(0.55 + i * 2.95 / 90)) for i in range(91)],
                         BLUE, 2.0))
    a_i, b_i = 2, 0
    o.append(ax.dot(a_i, b_i, 5.0, RED))
    o.append(ax.text(a_i, b_i, "anchor (a&#7522;, &#946;&#7522;)", 10,
                     dy=17, color=RED, anchor="middle"))
    for p in (0.8, 1.3, 3.0, 3.45):
        o.append(ax.polyline([(a_i, b_i), (p, q(p))], "#d99", 1.0))
    xs = [0.7 + i * (3.45 - 0.7) / 40 for i in range(41)]
    xpix = [fmt(ax.X(x)) for x in xs]
    ypix = [fmt(ax.Y(q(x))) for x in xs]
    xv = ";".join(xpix + xpix[::-1])
    yv = ";".join(ypix + ypix[::-1])
    o.append(f'<line x1="{fmt(ax.X(a_i))}" y1="{fmt(ax.Y(b_i))}" '
             f'x2="{xpix[0]}" y2="{ypix[0]}" stroke="{PURPLE}" '
             f'stroke-width="2">'
             f'<animate attributeName="x2" values="{xv}" dur="10s" '
             f'repeatCount="indefinite"/>'
             f'<animate attributeName="y2" values="{yv}" dur="10s" '
             f'repeatCount="indefinite"/></line>')
    o.append(f'<circle cx="{xpix[0]}" cy="{ypix[0]}" r="4.6" '
             f'fill="{PURPLE}">'
             f'<animate attributeName="cx" values="{xv}" dur="10s" '
             f'repeatCount="indefinite"/>'
             f'<animate attributeName="cy" values="{yv}" dur="10s" '
             f'repeatCount="indefinite"/></circle>')
    ax2 = Ax(0.3, 4.2, -1.6, 2.1, W=250, H=250, ml=48, mr=10, mt=24,
             mb=34)
    o.append(f'<g transform="translate(370,0)">')
    o.append(ax2.axes(xticks=(1, 2, 3, 4), yticks=(-1, 0, 1, 2),
                      xlab="p", ylab="slope"))
    o.append(ax2.polyline([(x, x - 2) for x in (0.55, 3.95)], PURPLE,
                          2.0))
    o.append(ax2.text(2.4, 1.55, "(Q(p) &#8722; &#946;&#7522;)/"
                      "(p &#8722; a&#7522;)", 10.5, color=PURPLE,
                      anchor="middle"))
    o.append(ax2.text(2.4, 1.1, "= p &#8722; 2: degree dropped", 10.5,
                      color=PURPLE, anchor="middle"))
    o.append('</g>')
    o.append(txt(190, 16, "secants from the anchor (animated)", 11,
                 "middle", GRAY))
    o.append(txt(500, 16, "their slopes vs p", 11, "middle", GRAY))
    o.append(txt(18, H - 26, "real-line analogy of (C4): if the curve "
                 "passes through the anchor, the secant slope "
                 "(Q(X)&#8722;&#946;&#7522;)/(X&#8722;a&#7522;) is a "
                 "POLYNOMIAL of degree deg&#8202;Q &#8722; 1;", 11,
                 "start", DARK))
    o.append(txt(18, H - 10, "if not, it has a pole at a&#7522;. Over "
                 "F&#8339; the same dichotomy is the factor theorem "
                 "&#8212; and soundness will detect it with valuations "
                 "(&#167;10, &#167;12).", 11, "start", DARK))
    return svg_wrap("\n".join(o), W, H,
                    "Shifted moments are secant slopes")


@fig("fig_weight")
def f_weight():
    W, H = 620, 240
    o = []
    gx = 120
    cw = 26
    labels = ["global table 0", "table (C&#8321;, &#963;|I&#8321;)",
              "table (C&#8322;, &#963;|I&#8322;)"]
    for r in range(3):
        y = 46 + r * 44
        o.append(txt(gx - 10, y + 16, labels[r], 11, "end",
                     GREEN if r else RED))
        for k in range(13):
            o.append(rect(gx + k * cw, y, cw - 4, 24,
                          LRED if r == 0 else LGREEN,
                          RED if r == 0 else GREEN, 1.0, rx=3))
            o.append(txt(gx + k * cw + (cw - 4) / 2, y + 16, "1", 10,
                         "middle", RED if r == 0 else GREEN))
    o.append(txt(gx + 6.5 * cw, 30, "|P| = 13 ones per populated table",
                 11, "middle", GRAY))
    o.append(f'<path d="M {gx + 13 * cw + 8} 46 L {gx + 13 * cw + 22} '
             f'46 L {gx + 13 * cw + 22} 158 L {gx + 13 * cw + 8} 158" '
             f'fill="none" stroke="{DARK}" stroke-width="1.2"/>')
    o.append(txt(gx + 13 * cw + 30, 100, "&#8467; + 1 = 3", 11.5,
                 "start", DARK))
    o.append(txt(gx + 13 * cw + 30, 116, "tables", 11.5, "start", DARK))
    o.append(txt(18, 204, "wt(x&#963;) = (&#8467; + 1)&#8202;|P| = R "
                 "&#8212; toy: 3 &#215; 13 = 39;  real scale: R = "
                 "(&#8467;+1)(q &#8722; m) &#8776; N&#8202;q &#8776; "
                 "N&#178;&#8304;&#185;", 12, "start", DARK))
    o.append(txt(18, 226, "every other coordinate of x&#963; is zero: "
                 "13 of the 15 toy tables are empty", 11, "start", GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "The completeness weight budget")


@fig("fig_params")
def f_params():
    W, H = 640, 330
    o = ['<defs>' + marker2("c7par_a", GRAY) + '</defs>']
    rows = [
        ("N", "100 + s + m + &#8467;", "master size parameter", DARK),
        ("q", "least 2^e &#8805; N&#178;&#8304;&#8304; (e = 1329 at "
         "N = 100)", "field size: the big room", BLUE),
        ("K", "N&#8308;", "fiber-size cap in soundness", PURPLE),
        ("T", "N&#179;&#8304;", "moment budget (j &#8804; T)", GREEN),
        ("M", "|&#920;|&#8202;|P|&#8202;q &#8804; 40N&#8308;&#8304;"
         "&#185;", "number of binary unknowns = lattice rank n", RED),
    ]
    y = 34
    for (sym, val, role, col) in rows:
        o.append(rect(20, y, 44, 34, LGRAY, col, 1.3))
        o.append(txt(42, y + 22, sym, 14, "middle", col, w="b"))
        o.append(txt(76, y + 15, val, 12, "start", col))
        o.append(txt(76, y + 30, role, 10.5, "start", GRAY))
        y += 44
    xr = 400
    ineqs = [
        ("dT &#8804; N&#179;&#185; &lt; |P|", "moments are determined "
         "(RS uniqueness)"),
        ("T &#8805; 2K &#8722; 1", "Hankel matrix entries all exist"),
        ("|P| &#8722; dK(K&#8722;1) &gt; 2dK&#178;T", "zero-counting "
         "endgame of Lemma 10"),
        ("z + h &#8722; 1 &#8804; 5N&#178;&#185; &lt; T", "valuation "
         "amplification in Lemma 12"),
        ("discards &lt; q/20", "Markov: few oversized fibers"),
    ]
    yy = 34
    for (iq, why) in ineqs:
        o.append(rect(xr, yy, 224, 40, LBLUE, BLUE, 1.1))
        o.append(txt(xr + 112, yy + 17, iq, 11, "middle", BLUE, w="b"))
        o.append(txt(xr + 112, yy + 32, why, 9, "middle", GRAY))
        yy += 50
    for (sy, ty) in [(1, 0), (3, 1), (2, 2), (3, 2), (2, 3), (1, 4)]:
        o.append(line(258, 51 + sy * 44, 396, 54 + ty * 50, GRAY, 1.0,
                      dash="3 3", marker="c7par_a"))
    o.append(txt(20, 282, "every inequality was re-verified with exact "
                 "integer arithmetic at N = 100 (machine checks, "
                 "&#167;17)", 11.5, "start", GREEN))
    o.append(txt(20, 302, "the theme: q &#8776; N&#178;&#8304;&#8304; "
                 "is chosen so large that every loss the soundness "
                 "argument suffers is a vanishing fraction of |P| "
                 "&#8776; q", 11.5, "start", DARK))
    return svg_wrap("\n".join(o), W, H, "The parameter dashboard")


@fig("fig_family")
def f_family():
    W, H = 640, 300
    o = []

    def q1(x):
        return 0.55 + 0.32 * math.sin(1.35 * x + 0.4) + 0.09 * x

    def q2(x):
        return 1.75 - 0.22 * x + 0.13 * math.sin(2.1 * x)

    def q3(x):
        return 1.1 + 0.18 * math.cos(1.7 * x)
    axl = Ax(0, 4.3, 0, 2.3, W=310, H=240, ml=30, mr=8, mt=30, mb=30)
    o.append(axl.axes(xlab="p", ylab=""))
    for f, col in [(q1, BLUE), (q2, RED), (q3, GREEN)]:
        o.append(axl.polyline([(0.1 + i * 4 / 80,
                                f(0.1 + i * 4 / 80))
                               for i in range(81)], col, 1.8))
    ps = [0.5 + k * 0.5 for k in range(8)]
    for p in ps:
        for f, col in [(q1, BLUE), (q2, RED), (q3, GREEN)]:
            o.append(axl.dot(p, f(p), 3.4, col))
    o.append(txt(165, 18, "labeled: h global functions q&#8348;(X)", 11.5,
                 "middle", DARK))
    axr = Ax(0, 4.3, 0, 2.3, W=310, H=240, ml=30, mr=8, mt=30, mb=30)
    o.append('<g transform="translate(330,0)">')
    o.append(axr.axes(xlab="p", ylab=""))
    for p in ps:
        for f in (q1, q2, q3):
            o.append(axr.dot(p, f(p), 3.4, GRAY))
    o.append('</g>')
    o.append(txt(495, 18, "what the moments see: unordered fibers S(p)",
                 11.5, "middle", DARK))
    o.append(txt(320, H - 10, "same data, labels erased: the "
                 "reconstruction problem is to recover a global "
                 "description anyway &#8212; from power sums alone", 11.5,
                 "middle", GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "The reconstruction problem")


@fig("fig_hankel")
def f_hankel():
    W, H = 640, 300
    o = []
    cs = 52
    gx, gy = 60, 60
    hues = ["#fde8e6", "#fdf3dc", "#e8f2e2", "#e4eefa", "#f0e6f4",
            "#f7e6dc", "#e6f4f2"]
    strokes = [RED, "#b8860b", GREEN, BLUE, PURPLE, "#a0522d",
               "#2a7d7d"]
    for i in range(4):
        for l in range(4):
            k = i + l
            o.append(rect(gx + l * cs, gy + i * cs, cs - 4, cs - 4,
                          hues[k], strokes[k], 1.2, rx=4))
            o.append(txt(gx + l * cs + cs / 2 - 2,
                         gy + i * cs + cs / 2 + 3,
                         f"&#956;{chr(0x2080 + k)}", 13, "middle",
                         strokes[k]))
    o.append(txt(gx + 2 * cs - 2, gy - 18, "Hankel matrix "
                 "(&#956;&#7522;&#8330;&#8202;&#8467;(X)), h = 4", 12,
                 "middle", DARK))
    o.append(txt(gx + 2 * cs - 2, gy + 4 * cs + 16, "constant along "
                 "anti-diagonals: entry = &#956;&#7522;&#8330;&#8202;"
                 "&#8467;", 10.5, "middle", GRAY))
    xr = 350
    o.append(txt(xr, gy - 18, "why it is invertible at a full fiber "
                 "S(p) = {w&#8321;,&#8230;,w&#8344;}:", 12, "start",
                 DARK))
    o.append(txt(xr, gy + 8, "&#956;&#7522;&#8330;&#8202;&#8467;(p) = "
                 "&#8721;&#8348; w&#8348;^&#8202;i&#8202;&#183;&#8202;"
                 "w&#8348;^&#8202;&#8467;", 12.5, "start", DARK))
    o.append(txt(xr, gy + 36, "so  (&#956;&#7522;&#8330;&#8202;&#8467;"
                 "(p)) = V&#8202;V&#7511;,   V&#7522;&#8202;&#8348; = "
                 "w&#8348;^&#8202;i  (Vandermonde)", 12.5, "start",
                 DARK))
    o.append(txt(xr, gy + 66, "det = det(V)&#178; = &#8719;&#8202;"
                 "(w&#8348; &#8722; w&#8351;)&#178; &#8800; 0", 12.5,
                 "start", RED))
    o.append(txt(xr, gy + 96, "distinct elements &#8658; invertible;",
                 11.5, "start", GRAY))
    o.append(txt(xr, gy + 112, "a smaller fiber gives rank &lt; h "
                 "&#8658; det = 0", 11.5, "start", GRAY))
    o.append(txt(xr, gy + 142, "one system, X left symbolic:", 11.5,
                 "start", DARK))
    o.append(txt(xr, gy + 164, "&#8721;&#8202;&#8467; &#956;&#7522;"
                 "&#8330;&#8202;&#8467;&#8202;c&#8467; = &#8722;&#956;"
                 "&#7522;&#8330;&#8202;h   (0 &#8804; i &lt; h)", 12.5,
                 "start", BLUE))
    o.append(txt(xr, gy + 186, "solves ONCE over K(X); specializing at "
                 "any", 11.5, "start", GRAY))
    o.append(txt(xr, gy + 202, "full fiber recovers that fiber's "
                 "polynomial", 11.5, "start", GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "The Hankel system behind Lemma 10")


@fig("fig_sheets")
def f_sheets():
    W, H = 640, 330
    ax = Ax(-1.4, 2.6, -2.2, 1.4, W=360, H=280, ml=40, mr=10, mt=30,
            mb=36)
    o = []
    o.append(ax.axes(xticks=(0, 1, 2), yticks=(-2, -1, 0, 1), xlab="X",
                     ylab="Y"))
    up = []
    dn = []
    for i in range(121):
        x = -0.25 + i * (2.5 + 0.25) / 120
        s = math.sqrt(max(0.0, 1 + 4 * x))
        up.append((x, (-1 + s) / 2))
        dn.append((x, (-1 - s) / 2))
    o.append(ax.polyline(up, BLUE, 2.2))
    o.append(ax.polyline(dn, PURPLE, 2.2))
    o.append(ax.dot(-0.25, -0.5, 5, RED))
    o.append(ax.text(-0.25, -0.30, "branch point:", 10, color=RED,
                     anchor="end", dx=-9))
    o.append(ax.text(-0.25, -0.52, "disc = 0,", 10, color=RED,
                     anchor="end", dx=-9))
    o.append(ax.text(-0.25, -0.74, "fiber collapses", 10, color=RED,
                     anchor="end", dx=-9))
    for x in (0.6, 1.6):
        o.append(ax.vline(x, GRAY, 1.0, dash="3 3", y1=0.92))
        s = math.sqrt(1 + 4 * x)
        o.append(ax.dot(x, (-1 + s) / 2, 4.4, BLUE))
        o.append(ax.dot(x, (-1 - s) / 2, 4.4, PURPLE))
    o.append(ax.text(0.3, 1.2, "fiber S(p): the sheet values above p",
                     10.5, color=DARK, anchor="start"))
    o.append(txt(200, 18, "real-line analogy: the curve Y&#178; + Y = X",
                 11.5, "middle", GRAY))
    xr = 420
    o.append(txt(xr + 90, 30, "the F&#8321;&#8326; reality:", 11.5,
                 "middle", DARK))
    o.append(txt(xr + 90, 46, "fibers of Y&#178; + Y = p", 11.5,
                 "middle", DARK))
    fibs = []
    for p in range(16):
        ws = [w for w in range(16) if gmul(w, w) ^ w == p]
        if ws:
            fibs.append((p, ws))
    yv = 64
    o.append(txt(xr + 24, yv, "p", 11, "middle", GRAY))
    o.append(txt(xr + 120, yv, "S(p)", 11, "middle", GRAY))
    for (p, ws) in fibs:
        yv += 24
        o.append(rect(xr + 6, yv - 15, 40, 20, LBLUE, BLUE, 0.9, rx=3))
        o.append(txt(xr + 26, yv, str(p), 10.5, "middle", DARK))
        o.append(rect(xr + 66, yv - 15, 110, 20, LGRAY, GRAY, 0.9,
                      rx=3))
        o.append(txt(xr + 121, yv, f"{{{ws[0]}, {ws[1]}}}", 10.5,
                     "middle", DARK))
    o.append(txt(xr + 90, yv + 26, "G(Y) = Y&#178; + Y + X:", 11,
                 "middle", RED))
    o.append(txt(xr + 90, yv + 42, "roots &#8713; K(X) &#8212; a", 10.5,
                 "middle", GRAY))
    o.append(txt(xr + 90, yv + 56, "genuine extension needed", 10.5,
                 "middle", GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "A two-sheeted global description")


@fig("fig_lemma10flow")
def f_lemma10flow():
    W, H = 640, 344
    o = ['<defs>' + marker2("c7l10_a", DARK) + '</defs>']
    steps = [
        (20, 20, 280, 52, "1. factorize at a full fiber",
         "(&#956;(p)) = V&#8202;V&#7511; &#8658; &#916;&#8202;(p) "
         "&#8800; 0 &#8658; &#916; &#8800; 0", BLUE),
        (340, 20, 280, 52, "2. count zeros of &#916;",
         "deg &#916; &#8804; dh(h&#8722;1) &#8658; &#8805; |P| &#8722; "
         "dh(h&#8722;1) full fibers", BLUE),
        (20, 108, 280, 52, "3. solve the Hankel system over K(X)",
         "unique c&#8467; = n&#8467;/&#916;, deg n&#8467; &#8804; "
         "dh&#178; &#8658; G(Y)", PURPLE),
        (340, 108, 280, 52, "4. specialize at full fibers",
         "ev&#8346;(G) = &#8719;(Y &#8722; w), w &#8712; S(p)", PURPLE),
        (20, 196, 280, 52, "5. Newton: power sums of roots",
         "&#916;^j(&#956;&#8202;&#11388; &#8722; P&#8202;&#11388;) has "
         "deg &#8804; 2dK&#178;T, &gt; 2dK&#178;T zeros", RED),
        (340, 196, 280, 52, "6. conclude the moment identities",
         "&#956;&#11388;(X) = &#8721;&#8202;&#945;^&#8202;j over the "
         "root set R", RED),
        (180, 284, 300, 48, "7. separability for free",
         "&#916; = &#8719;(&#945;&#8348; &#8722; &#945;&#8351;)&#178; "
         "&#8800; 0 &#8658; roots distinct", GREEN),
    ]
    for (x, y, w, h, l1, l2, col) in steps:
        o.append(rect(x, y, w, h, LGRAY, col, 1.4))
        o.append(txt(x + w / 2, y + 20, l1, 11.5, "middle", col, w="b"))
        o.append(txt(x + w / 2, y + 38, l2, 10.5, "middle", DARK))
    o.append(line(300, 46, 336, 46, DARK, 1.4, marker="c7l10_a"))
    o.append(line(480, 72, 200, 104, DARK, 1.4, marker="c7l10_a"))
    o.append(line(300, 134, 336, 134, DARK, 1.4, marker="c7l10_a"))
    o.append(line(480, 160, 200, 192, DARK, 1.4, marker="c7l10_a"))
    o.append(line(300, 222, 336, 222, DARK, 1.4, marker="c7l10_a"))
    o.append(line(480, 248, 380, 280, DARK, 1.4, marker="c7l10_a"))
    return svg_wrap("\n".join(o), W, H,
                    "Proof plan of the reconstruction lemma")


@fig("fig_zerocount")
def f_zerocount():
    W, H = 620, 210
    o = []
    xa, xb = 60, 560
    y1, y2 = 70, 130
    o.append(txt(18, 30, "the counting hammer used twice in Lemma 10 "
                 "(and once in Lemma 13):", 12.5, "start", DARK))
    o.append(rect(xa, y1 - 16, (xb - xa) * 0.55, 26, LRED, RED, 1.2))
    o.append(txt(xa + (xb - xa) * 0.275, y1, "max #zeros = deg &#8804; "
                 "2dK&#178;T", 11.5, "middle", RED))
    o.append(rect(xa, y2 - 16, (xb - xa) * 0.9, 26, LGREEN, GREEN, 1.2))
    o.append(txt(xa + (xb - xa) * 0.45, y2, "#points where it provably "
                 "vanishes &#8805; |P| &#8722; dK(K&#8722;1) &gt; "
                 "2dK&#178;T", 11.5, "middle", GREEN))
    o.append(line(xa + (xb - xa) * 0.55, y1 + 14, xa + (xb - xa) * 0.55,
                  y2 - 18, GRAY, 1.2, dash="4 3"))
    o.append(txt(310, 172, "a nonzero polynomial cannot vanish more "
                 "often than its degree &#8658; the polynomial is "
                 "identically zero", 12, "middle", DARK, w="b"))
    o.append(txt(310, 194, "this is where the moment identities become "
                 "TRUE IDENTITIES in K(X), not just facts about sampled "
                 "points", 11, "middle", GRAY))
    return svg_wrap("\n".join(o), W, H, "Degree versus zero count")


# ======================================================================
# Part C: valuations, soundness, transfer, consequences
@fig("fig_ord")
def f_ord():
    W, H = 640, 300
    o = []
    a = 1.0
    ax = Ax(-0.1, 2.4, -2.4, 2.6, W=420, H=270, ml=42, mr=10, mt=22,
            mb=32)
    o.append(ax.axes(xticks=(0, 1, 2), yticks=(-2, -1, 0, 1, 2),
                     xlab="X", ylab=""))
    o.append(ax.vline(a, RED, 1.2, dash="5 3"))

    def clipped(f, lo=-2.3, hi=2.5, x0=0.0, x1=2.35, n=160, avoid=None):
        pts = []
        run = []
        for i in range(n + 1):
            x = x0 + (x1 - x0) * i / n
            if avoid and abs(x - avoid) < 0.02:
                if run:
                    pts.append(run)
                    run = []
                continue
            y = f(x)
            if lo <= y <= hi:
                run.append((x, y))
            else:
                if run:
                    pts.append(run)
                    run = []
        if run:
            pts.append(run)
        return pts

    curves = [
        (lambda x: 1.3 * (x - a), GREEN, "ord = 1"),
        (lambda x: 2.6 * (x - a) ** 3, BLUE, "ord = 3"),
        (lambda x: 0.35 / (x - a), PURPLE, "ord = &#8722;1 (pole)"),
        (lambda x: 1.55 + 0 * x, "#b8860b", "ord = 0 (unit)"),
    ]
    for f, col, lab in curves:
        for seg in clipped(f, avoid=a):
            o.append(ax.polyline(seg, col, 1.9))
    o.append(ax.text(0.06, 1.72, "ord = 0 (unit)", 10.5,
                     color="#b8860b"))
    o.append(ax.text(1.85, 0.72, "ord = 1", 10.5, color=GREEN))
    o.append(ax.text(1.9, 2.28, "ord = 3", 10.5, color=BLUE,
                     anchor="end"))
    o.append(ax.text(1.12, -1.55, "ord = &#8722;1", 10.5, color=PURPLE))
    o.append(ax.text(a, 2.32, "a", 12, color=RED, dx=4, cls="v"))
    xr = 452
    o.append(txt(xr, 40, "ord&#8202;&#8467;&#8342;(g) = how many", 11.5,
                 "start", DARK))
    o.append(txt(xr, 56, "factors of &#8467;&#8342; = X &#8722; a", 11.5,
                 "start", DARK))
    o.append(txt(xr, 72, "divide g", 11.5, "start", DARK))
    o.append(txt(xr, 100, "&gt; 0: zero at a", 11.5, "start", GREEN))
    o.append(txt(xr, 118, "= 0: finite, nonzero", 11.5, "start",
                 "#b8860b"))
    o.append(txt(xr, 136, "&lt; 0: pole at a", 11.5, "start", PURPLE))
    o.append(txt(xr, 168, "v(gh) = v(g) + v(h)", 11.5, "start", DARK))
    o.append(txt(xr, 186, "v(g+h) &#8805; min(v(g), v(h))", 11.5,
                 "start", DARK))
    o.append(txt(xr, 214, "the real plots are an analogy;", 10.5,
                 "start", GRAY))
    o.append(txt(xr, 230, "the definition is pure algebra", 10.5,
                 "start", GRAY))
    o.append(txt(xr, 246, "and works over any field", 10.5, "start",
                 GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "Orders of vanishing at a point")


@fig("fig_ultrametric")
def f_ultrametric():
    W, H = 640, 260
    o = []
    for panel in range(2):
        ox = 30 + panel * 320
        o.append(txt(ox + 130, 30,
                     "unique minimum" if panel == 0 else
                     "tied minimum", 12, "middle", DARK, w="b"))
        y0, y1 = 200, 56
        o.append(line(ox, y0 + 12, ox, y1 - 8, DARK, 1.2))
        for lv in range(-1, 3):
            y = y0 - (lv + 1) * 36
            o.append(line(ox - 4, y, ox + 250, y, "#ccc", 0.8,
                          dash="3 3"))
            o.append(txt(ox - 10, y + 4, str(lv), 10.5, "end", GRAY))
        vals = [(-1, RED), (1, BLUE), (2, GREEN)] if panel == 0 else \
               [(-1, RED), (-1, BLUE), (2, GREEN)]
        for k, (lv, col) in enumerate(vals):
            x = ox + 40 + k * 55
            y = y0 - (lv + 1) * 36
            o.append(circ(x, y, 6, col))
            o.append(txt(x, y - 12, f"x{chr(0x2081 + k)}", 11, "middle",
                         col, cls="v"))
        xs = ox + 210
        if panel == 0:
            y = y0 - 0 * 36
            o.append(circ(xs, y, 7, "none", DARK, 2))
            o.append(txt(xs, y - 14, "&#8721;x&#7522;", 11, "middle",
                         DARK))
            o.append(txt(ox + 130, 234, "v(&#8721;) = min v(x&#7522;) "
                         "exactly &#8212; the sum CANNOT vanish", 10.5,
                         "middle", GREEN))
        else:
            for lv2 in (0, 1, 3):
                y = y0 - (lv2 + 1) * 36
                o.append(circ(xs, y, 5, "none", GRAY, 1.4, dash="2 2"))
            o.append(txt(ox + 130, 234, "v(&#8721;) &#8805; min only: "
                         "the two lows can cancel", 10.5, "middle",
                         RED))
        o.append(txt(ox - 24, 40, "v", 12, "middle", DARK, cls="v"))
    return svg_wrap("\n".join(o), W, H,
                    "The minimum-valuation principle (Lemma 2)")


@fig("fig_ramify")
def f_ramify():
    W, H = 560, 270
    o = []
    x1, x2 = 150, 380
    y0 = 220
    sc = 60
    o.append(txt(x1, 36, "values on F = K(X)", 11.5, "middle", DARK))
    o.append(txt(x2, 36, "values on E &#8817; F, e = 2", 11.5, "middle",
                 DARK))
    for lv in (-1, 0, 1, 2):
        y = y0 - (lv + 1) * sc
        o.append(line(x1 - 40, y, x1 + 40, y, BLUE, 2))
        o.append(txt(x1 - 50, y + 4, str(lv), 11, "end", BLUE))
    for k in range(-2, 5):
        lv = k / 2
        y = y0 - (lv + 1) * sc
        full = k % 2 == 0
        o.append(line(x2 - 40, y, x2 + 40, y, PURPLE if not full else
                      BLUE, 2 if full else 1.2,
                      dash=None if full else "5 3"))
        if not full:
            o.append(txt(x2 + 50, y + 4, f"{k}/2", 11, "start", PURPLE))
    o.append(circ(x1, y0 - 2 * sc, 6, BLUE))
    o.append(txt(x1 + 12, y0 - 2 * sc + 4, "&#8467;&#8342;", 12, "start",
                 BLUE, cls="v"))
    o.append(circ(x2, y0 - 1.5 * sc, 6, PURPLE))
    o.append(txt(x2, y0 - 1.5 * sc + 22,
                 "u with u&#178; = &#8467;&#8342;", 11.5, "middle",
                 PURPLE))
    o.append(txt(280, 250, "2&#8202;v(u) = v(&#8467;&#8342;) = 1 "
                 "&#8658; v(u) = 1/2: new, fractional rungs &#8212; "
                 "ramification", 11.5, "middle", DARK))
    o.append(txt(280, 62, "extending ord&#8202;&#8467;&#8342; to a "
                 "finite separable E may refine the value ladder to "
                 "(1/e)&#8202;&#8484; (Lemma 4)", 10.5, "middle", GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "Ramification refines the valuation ladder")


@fig("fig_markov")
def f_markov():
    W, H = 640, 300
    import random as rnd
    rnd.seed(11)
    npts = 44
    sizes = []
    for k in range(npts):
        base = rnd.choice([0, 1, 1, 1, 2, 2, 3])
        sizes.append(base)
    for k in (7, 19, 33):
        sizes[k] = rnd.choice([9, 11, 13])
    KCAP = 6
    ax_x0, bw = 46, 12.4
    y0 = 210
    sc = 13
    o = []
    o.append(txt(18, 26, "fiber sizes |S&#964;(p)| of one table of a "
                 "low-weight solution (schematic data)", 12, "start",
                 DARK))
    for k, s in enumerate(sizes):
        x = ax_x0 + k * bw
        big = s > KCAP
        if s > 0:
            o.append(rect(x, y0 - s * sc, bw - 3, s * sc,
                          LRED if big else LBLUE, RED if big else BLUE,
                          0.9, rx=2))
        if big:
            o.append(line(x - 1, y0 - s * sc - 6, x + bw - 2, y0 + 6,
                          RED, 1.6))
            o.append(line(x - 1, y0 + 6, x + bw - 2, y0 - s * sc - 6,
                          RED, 1.6))
    o.append(line(ax_x0 - 8, y0, ax_x0 + npts * bw + 4, y0, DARK, 1.3))
    o.append(line(ax_x0 - 8, y0 - KCAP * sc, ax_x0 + npts * bw + 4,
                  y0 - KCAP * sc, PURPLE, 1.6, dash="6 4"))
    o.append(txt(ax_x0 + npts * bw + 2, y0 - KCAP * sc - 6,
                 "cap K = N&#8308;", 11.5, "end", PURPLE))
    o.append(txt(ax_x0 + npts * bw / 2, y0 + 18,
                 "evaluation points p &#8712; P", 11, "middle", GRAY))
    o.append(txt(18, 250, "total area = wt(x) &#8804; 4M^&#8201;1/200"
                 "&#8202;R, so at most wt(x)/K columns can exceed the "
                 "cap (Markov):", 12, "start", DARK))
    o.append(txt(18, 272, "|P &#8735; P&#964;| &#8804; 4M^&#8201;1/200"
                 "&#8202;R/N&#8308; &lt; q/20 &#8212; discard them "
                 "(&#10007;) and keep P&#964;, still &#8805; |P| &#8722; "
                 "q/20 points, separately for each table", 10.8, "start",
                 DARK))
    return svg_wrap("\n".join(o), W, H,
                    "Markov discard of oversized fibers")


@fig("fig_newton12")
def f_newton12():
    W, H = 640, 330
    o = []
    xa, xb = 70, 590
    y1 = 90

    def X(v):
        return xa + (v + 5.4) / 5.9 * (xb - xa)
    o.append(txt(18, 26, "where can t&#8348; = v&#7522;(y&#8348;) lie? "
                 "(y&#8348; = (&#945;&#8348; &#8722; &#946;&#7522;)/"
                 "&#8467;&#7522;, root of a monic equation with "
                 "coefficient valuations &#8805; &#8722;D&#8320;)", 12,
                 "start", DARK))
    o.append(line(xa - 12, y1, xb + 16, y1, DARK, 1.4))
    for (v, lab) in [(-5.0, "&#8722;D&#8320;"), (-1.0, "&#8722;1/h"),
                     (0.0, "0")]:
        o.append(line(X(v), y1 - 8, X(v), y1 + 8, DARK, 1.4))
        o.append(txt(X(v), y1 + 24, lab, 12, "middle", DARK))
    o.append(rect(X(-5.4), y1 - 7, X(-5.0) - X(-5.4), 14, LRED, RED,
                  1.0, rx=2))
    o.append(txt(X(-5.2), y1 - 16, "impossible: y^h would dominate "
                 "(Lemma 2)", 10, "middle", RED))
    o.append(rect(X(-1.0), y1 - 7, X(0) - X(-1.0), 14, LRED, RED, 1.0,
                  rx=2))
    o.append(txt(X(-0.5), y1 - 16, "impossible: two-term tie forces "
                 "t&#8348; &#8804; &#8722;1/h", 10, "middle", RED))
    o.append(rect(X(0), y1 - 7, X(0.5) - X(0), 14, LGREEN, GREEN, 1.0,
                  rx=2))
    o.append(txt(X(0.28), y1 - 16, "goal: t&#8348; &#8805; 0", 10,
                 "middle", GREEN))
    o.append(rect(X(-5.0), y1 - 7, X(-1.0) - X(-5.0), 14, "#fdf3dc",
                  "#b8860b", 1.0, rx=2))
    o.append(txt(X(-3.0), y1 - 16, "still allowed after step 1: "
                 "[&#8722;D&#8320;, &#8722;1/h]", 10, "middle",
                 "#b8860b"))
    y2 = 208
    o.append(txt(18, 148, "step 2 &#8212; amplification: take the z-th "
                 "moment window with z = hU&#8344; + 1 (moments up to "
                 "z + h &#8722; 1 &#8804; 5N&#178;&#185; &lt; T exist):",
                 11, "start", DARK))
    o.append(line(xa - 12, y2, xb + 16, y2, DARK, 1.4))
    for (v, lab) in [(-3.9, "&#8722;z/h"), (-1.9, "&#8722;U&#8344;"),
                     (0.0, "0")]:
        o.append(line(X(v), y2 - 8, X(v), y2 + 8, DARK, 1.4))
        o.append(txt(X(v), y2 + 24, lab, 12, "middle", DARK))
    o.append(line(X(-1.9), y2 - 36, X(-1.9), y2 - 8, BLUE, 1.4,
                  dash="4 3"))
    o.append(txt(X(-1.9), y2 - 40, "matrix bound: v&#7522;(y&#8348;^z) "
                 "&#8805; &#8722;U&#8344; (V&#8315;&#185; &#183; "
                 "nonneg moments)", 10.5, "middle", BLUE))
    o.append(line(X(-3.9), y2 + 36, X(-3.9), y2 + 8, RED, 1.4,
                  dash="4 3"))
    o.append(txt(X(-3.9), y2 + 50, "if t&#8348; &lt; 0: v&#7522;"
                 "(y&#8348;^z) = z&#8202;t&#8348; &#8804; &#8722;z/h "
                 "&lt; &#8722;U&#8344;", 10.5, "middle", RED))
    o.append(txt(310, 290, "the two bounds contradict &#8658; every "
                 "t&#8348; &#8805; 0 &#8658; v&#7522;(&#945; &#8722; "
                 "&#946;&#7522;) &#8805; 1: the root REMEMBERS the "
                 "clause bit", 12, "middle", DARK, w="b"))
    return svg_wrap("\n".join(o), W, H,
                    "The valuation window argument of Lemma 12")


@fig("fig_matching")
def f_matching():
    W, H = 640, 300
    o = []
    o.append(txt(18, 24, "clause C: the parity constraint (C2) forces "
                 "&#956;&#8320;,&#11388; = &#8721;&#946; &#956;"
                 "(C,&#946;),&#11388; &#8658; power sums of the "
                 "multiset difference vanish", 11.5, "start", DARK))
    o.append(rect(40, 60, 120, 170, LGRAY, DARK, 1.3))
    o.append(txt(100, 50, "global roots R&#8320;", 11.5, "middle",
                 DARK))
    roots = [(100, 100, RED), (100, 150, BLUE), (100, 200, GREEN)]
    for (x, y, col) in roots:
        o.append(circ(x, y, 7, col))
    boxes = [(260, 66, "R(C, &#946;&#8321;)", [RED]),
             (260, 146, "R(C, &#946;&#8322;)", [BLUE, GREEN]),
             (260, 226, "R(C, &#946;&#8323;)", [])]
    for (bx, by, lab, cols) in boxes:
        o.append(rect(bx, by - 26, 150, 52, LBLUE, BLUE, 1.2))
        o.append(txt(bx + 75, by - 34, lab, 11, "middle", BLUE))
        for k, col in enumerate(cols):
            o.append(circ(bx + 40 + k * 40, by, 7, col))
        if not cols:
            o.append(txt(bx + 75, by + 4, "(empty)", 10, "middle",
                         GRAY))
    o.append(line(107, 100, 296, 66, RED, 1.4, dash="4 3"))
    o.append(line(107, 150, 296, 146, BLUE, 1.4, dash="4 3"))
    o.append(line(107, 200, 336, 150, GREEN, 1.4, dash="4 3"))
    xr = 452
    o.append(txt(xr, 80, "if some &#945; &#8712; R&#8320; appeared in", 11,
                 "start", DARK))
    o.append(txt(xr, 96, "NO clause set, it would sit an odd", 11,
                 "start", DARK))
    o.append(txt(xr, 112, "number of times in the difference:", 11,
                 "start", DARK))
    o.append(txt(xr, 136, "&#8721;&#945;&#8712;W&#8202;&#945;^j = 0,",
                 12, "start", RED))
    o.append(txt(xr, 154, "j = 0, &#8230;, |W|&#8722;1", 11, "start",
                 RED))
    o.append(txt(xr, 178, "Vandermonde &#183; (1&#8230;1)&#7511; = 0", 11.5,
                 "start", DARK))
    o.append(txt(xr, 196, "with V invertible &#8212;", 11.5, "start",
                 DARK))
    o.append(txt(xr, 214, "impossible. So W = &#8709;:", 11.5, "start",
                 DARK))
    o.append(txt(xr, 240, "every global root matches", 11.5, "start",
                 GREEN))
    o.append(txt(xr, 258, "a satisfying local assignment", 11.5,
                 "start", GREEN))
    o.append(txt(18, 282, "|W| &#8804; 9K &lt; T + 1, so enough moment "
                 "equations exist to run the Vandermonde argument", 10.5,
                 "start", GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "Parity matching of clause roots (Lemma 13)")


@fig("fig_glue")
def f_glue():
    W, H = 640, 320
    o = ['<defs>' + marker2("c7glu_a", DARK) + '</defs>']
    o.append(circ(320, 96, 14, LRED, RED, 2))
    o.append(txt(320, 101, "&#945;", 13, "middle", RED, cls="v"))
    o.append(txt(320, 66, "one global root &#945; &#8712; R&#8320; "
                 "(Corollary 11: R&#8320; &#8800; &#8709;)", 11.5,
                 "middle", DARK))
    o.append(rect(80, 150, 190, 60, LBLUE, BLUE, 1.3))
    o.append(txt(175, 172, "clause C:  &#945; &#8712; R(C, &#946;)",
                 11.5, "middle", BLUE))
    o.append(txt(175, 192, "Lemma 12: v&#7522;(&#945; &#8722; "
                 "&#946;&#7522;) &#8805; 1", 11.5, "middle", BLUE))
    o.append(rect(370, 150, 190, 60, LBLUE, BLUE, 1.3))
    o.append(txt(465, 172, "clause C&#8242;:  &#945; &#8712; R(C&#8242;, "
                 "&#946;&#8242;)", 11.5, "middle", BLUE))
    o.append(txt(465, 192, "Lemma 12: v&#7522;(&#945; &#8722; &#946;"
                 "&#8242;&#7522;) &#8805; 1", 11.5, "middle", BLUE))
    o.append(line(308, 110, 200, 146, DARK, 1.4, marker="c7glu_a"))
    o.append(line(332, 110, 440, 146, DARK, 1.4, marker="c7glu_a"))
    o.append(txt(320, 136, "Lemma 13, for each clause", 10.5, "middle",
                 GRAY))
    o.append(rect(180, 238, 280, 40, LGREEN, GREEN, 1.4))
    o.append(txt(320, 255, "shared variable i: &#946;&#7522; = &#946;"
                 "&#8242;&#7522;", 12, "middle", GREEN, w="b"))
    o.append(txt(320, 271, "(Lemma 5: two constants close to &#945; "
                 "must be equal)", 10, "middle", GREEN))
    o.append(line(220, 214, 280, 234, DARK, 1.4, marker="c7glu_a"))
    o.append(line(420, 214, 360, 234, DARK, 1.4, marker="c7glu_a"))
    o.append(txt(320, 304, "so the local tuples glue into ONE global "
                 "assignment &#963; that satisfies every clause "
                 "&#8658; &#966; satisfiable (Prop 14)", 11.5, "middle",
                 DARK, w="b"))
    return svg_wrap("\n".join(o), W, H,
                    "Gluing local assignments through one root")


@fig("fig_branches")
def f_branches():
    W, H = 640, 330
    o = ['<defs>' + marker2("c7brn_a", DARK) + '</defs>']

    def box(x, y, w, h, l1, l2, fill, col):
        o.append(rect(x, y, w, h, fill, col, 1.3))
        o.append(txt(x + w / 2, y + 20, l1, 11.5, "middle", col, w="b"))
        if l2:
            o.append(txt(x + w / 2, y + 37, l2, 10.5, "middle", DARK))
    box(20, 20, 190, 48, "preprocess &#966;", "drop tautologies, "
        "repeats", LGRAY, DARK)
    box(20, 104, 190, 48, "no clauses left?", "&#8594; fixed YES "
        "instance", LGREEN, GREEN)
    box(20, 168, 190, 48, "empty clause?", "&#8594; ((2),(1),&#189;): "
        "NO instance", LRED, RED)
    box(250, 104, 170, 48, "build system (3)", "H, b in poly time",
        LBLUE, BLUE)
    box(250, 190, 170, 48, "Hx = b consistent?", "Gaussian elimination",
        LBLUE, BLUE)
    box(460, 190, 165, 62, "output", "(B&#8341;, u, &#8968;&#8730;R"
        "&#8969;)", LGRAY, DARK)
    box(250, 268, 170, 48, "inconsistent", "&#8594; fixed NO instance",
        LRED, RED)
    o.append(line(115, 68, 115, 100, DARK, 1.4, marker="c7brn_a"))
    o.append(line(210, 128, 246, 128, DARK, 1.4, marker="c7brn_a"))
    o.append(line(335, 152, 335, 186, DARK, 1.4, marker="c7brn_a"))
    o.append(line(335, 238, 335, 264, RED, 1.4, marker="c7brn_a"))
    o.append(txt(348, 254, "no", 10, "start", RED))
    o.append(line(420, 214, 456, 214, GREEN, 1.4, marker="c7brn_a"))
    o.append(txt(438, 206, "yes", 10, "middle", GREEN))
    o.append(txt(460, 280, "guarantees: &#966; sat &#8658; dist "
                 "&#8804; r;", 10.5, "start", DARK))
    o.append(txt(460, 296, "&#966; unsat &#8658; dist &gt; "
                 "n^&#8201;1/400&#8202;r", 10.5, "start", DARK))
    o.append(txt(460, 312, "(Lemma 8 &#8658; consistency on sat side)",
                 9.5, "start", GRAY))
    o.append(txt(24, 246, "trivial branches keep the map total:", 10.5,
                 "start", GRAY))
    o.append(txt(24, 262, "every &#966; gets a legal instance", 10.5,
                 "start", GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "The full reduction as an algorithm")


@fig("fig_lpnorms")
def f_lpnorms():
    W, H = 640, 280
    o = []
    cx, cy, s = 130, 130, 78
    o.append(rect(cx - s, cy - s, 2 * s, 2 * s, "none", GREEN, 1.8,
                  rx=0))
    o.append(circ(cx, cy, s, "none", BLUE, 1.8))
    pts = f"{cx},{cy - s} {cx + s},{cy} {cx},{cy + s} {cx - s},{cy}"
    o.append(f'<polygon points="{pts}" fill="none" stroke="{RED}" '
             f'stroke-width="1.8"/>')
    o.append(line(cx - s - 16, cy, cx + s + 16, cy, GRAY, 1))
    o.append(line(cx, cy - s - 16, cx, cy + s + 16, GRAY, 1))
    o.append(txt(cx + 30, cy - 32, "p = 1", 11, "middle", RED))
    o.append(txt(cx - s - 8, cy - 20, "p = 2", 11, "end", BLUE))
    o.append(txt(cx + s - 14, cy - s - 6, "p = &#8734;", 11, "middle",
                 GREEN))
    o.append(txt(cx, 250, "unit balls of &#8467;&#8346; norms", 11.5,
                 "middle", DARK))
    ax = Ax(0.8, 8.4, 0.0, 0.0062, W=340, H=240, ml=64, mr=12, mt=26,
            mb=40)
    o.append('<g transform="translate(290,6)">')
    o.append(ax.axes(xticks=(1, 2, 4, 8), yticks=(),
                     xlab="norm exponent p", ylab=""))
    for (yv, lab) in [(1 / 200, "1/200"), (1 / 400, "1/400"),
                      (1 / 800, "1/800")]:
        o.append(ax.hline(yv, "#bbb", 0.9, dash="3 3"))
        o.append(ax.text(0.86, yv, lab, 9.5, color=GRAY, dy=-3))
    o.append(ax.polyline([(p / 20, 1 / (200 * (p / 20)))
                          for p in range(20, 170)], PURPLE, 2.2))
    for p in (1, 2, 4, 8):
        o.append(ax.dot(p, 1 / (200 * p), 4.2, PURPLE))
    o.append(ax.text(2, 1 / 380, "hardness exponent 1/(200p)", 10.5,
                     color=PURPLE, dx=10, dy=-6))
    o.append('</g>')
    o.append(txt(460, 270, "Corollary 16: n&#8201;<tspan "
                 "baseline-shift=\"super\" font-size=\"75%\">1/(200p)"
                 "</tspan>-hardness for every fixed rational p &#8805; 1",
                 11, "middle", DARK))
    return svg_wrap("\n".join(o), W, H,
                    "Fixed finite norms inherit the gap")


@fig("fig_window")
def f_window():
    W, H = 640, 230
    xa, xb = 60, 590
    y = 110

    def X(c):
        return xa + c / 0.62 * (xb - xa)
    o = [arrow_def("c7win_arr", DARK)]
    o.append(line(xa - 10, y, xb + 18, y, DARK, 1.5,
                  marker="c7win_arr"))
    o.append(txt(xb + 10, y + 22, "exponent c in &#947; = n&#8202;^c",
                 11.5, "end", DARK))
    o.append(rect(X(0), y - 10, X(0.0025) - X(0), 20, RED, RED, 1.0,
                  rx=2))
    o.append(line(X(0.0025), y - 26, X(0.0025), y + 10, RED, 1.4))
    o.append(txt(X(0.0025) + 4, y - 32, "c = 1/400: NP-hard "
                 "(Theorem 1, this chapter)", 11.5, "start", RED))
    o.append(rect(X(0.0025), y - 10, X(0.5) - X(0.0025), 20, LGRAY,
                  GRAY, 1.0, rx=2))
    o.append(txt(X(0.25), y + 34, "open: is GapCVP(n^c) NP-hard "
                 "for any 1/400 &lt; c &lt; 1/2?", 11.5, "middle",
                 GRAY))
    o.append(line(X(0.5), y - 26, X(0.5), y + 26, BLUE, 2))
    o.append(txt(X(0.5), y - 34, "c = 1/2: GapCVP(C&#8730;n) &#8712; "
                 "NP &#8745; coNP", 11.5, "middle", BLUE))
    o.append(txt(X(0.5), y + 42, "NP-hardness here would give "
                 "NP = coNP", 10.5, "middle", BLUE))
    o.append(rect(X(0.5), y - 10, X(0.6) - X(0.5), 20, LBLUE, BLUE,
                  1.0, rx=2))
    o.append(txt(310, 200, "before this chapter, NO fixed c &gt; 0 was "
                 "known; the previous record factor n^&#8202;a/log log n "
                 "has exponent &#8594; 0", 11.5, "middle", DARK))
    return svg_wrap("\n".join(o), W, H,
                    "What is now known on the exponent axis")


@fig("fig_toysummary")
def f_toysummary():
    W, H = 640, 300
    o = ['<defs>' + marker2("c7toy_a", DARK) + '</defs>']
    steps = [
        (20, 30, 180, 74, DARK, LGRAY, "toy formula &#966;",
         ["m = 3 vars, &#8467; = 2 clauses", "&#963; = (1,0,1) "
          "satisfies"]),
        (240, 30, 180, 74, BLUE, LBLUE, "encoding field",
         ["q = 16, |P| = 13", "d = 3, T = 4, |&#920;| = 15"]),
        (460, 30, 164, 74, BLUE, LBLUE, "binary system",
         ["M = 3120 unknowns", "8949 equations"]),
        (20, 168, 180, 74, PURPLE, LGRAY, "linear algebra",
         ["rank H = 1839", "dim C = k&#8320; = 1281"]),
        (240, 168, 180, 74, GREEN, LGREEN, "completeness",
         ["H&#8202;x&#963; = b verified", "wt(x&#963;) = 39 = "
          "3&#183;13 = R"]),
        (460, 168, 164, 74, RED, LRED, "lattice instance",
         ["B&#8341;: 3120&#215;3120,", "det = 2&#185;&#8312;&#179;"
          "&#8313;, u &#8712; {0,1}&#8202;&#7481;"]),
    ]
    for (x, y, w, h, col, fill, l1, rest) in steps:
        o.append(rect(x, y, w, h, fill, col, 1.4))
        o.append(txt(x + w / 2, y + 20, l1, 12, "middle", col, w="b"))
        for k, s in enumerate(rest):
            o.append(txt(x + w / 2, y + 40 + k * 17, s, 10.5, "middle",
                         DARK))
    o.append(line(200, 67, 236, 67, DARK, 1.4, marker="c7toy_a"))
    o.append(line(420, 67, 456, 67, DARK, 1.4, marker="c7toy_a"))
    o.append(line(540, 104, 140, 164, DARK, 1.4, marker="c7toy_a"))
    o.append(line(200, 205, 236, 205, DARK, 1.4, marker="c7toy_a"))
    o.append(line(420, 205, 456, 205, DARK, 1.4, marker="c7toy_a"))
    o.append(txt(320, 272, "every arrow is executed by the chapter's "
                 "44-check battery; the in-page listing below "
                 "re-runs the core checks", 11, "middle", GRAY))
    return svg_wrap("\n".join(o), W, H,
                    "The toy instance end to end")


import re


def desup(svg):
    """Convert remaining literal ^exp runs in text content to tspan
    superscripts (check_svg forbids literal ^ in <text>)."""
    pat = re.compile(r"\^(?:&#820[12];)?((?:&#\d+;)|[0-9a-zA-Z/]+)")
    return pat.sub(lambda m: ('<tspan baseline-shift="super" '
                              'font-size="75%">' + m.group(1) +
                              '</tspan>'), svg)


def main():
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    for name, fn in FIGS.items():
        if only and name not in only:
            continue
        svg = desup(fn())
        assert "^" not in re.sub(r"<[^>]*>", "", svg), name
        write_fig(OUT, name, svg)
        print("wrote", name)


if __name__ == "__main__":
    main()
