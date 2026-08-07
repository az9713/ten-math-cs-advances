"""Generate all SVG figures for chapter 5 into .ignore/ch05_figs/.

Figure ledger (order of appearance = CSS counter number):
 1 fig_perdet      3x3 permanent as matchings
 2 fig_models     circuit vs formula vs division formula
 3 fig_results    the proved bounds, computed charts
 4 fig_roadmap    dependency map of both spines
 5 fig_agdim      cone / projectivization / avoidance count
 6 fig_powerfiber (d-1)^k fiber grid, computed roots
 7 fig_gateeq     gate equations + affine outputs
 8 fig_slice      Crit(x1x2x3) = axes + generic plane (3D proj)
 9 fig_project    direction map and kernel avoidance
10 fig_reverse    reverse-mode SMIL animation
11 fig_minorsum   M_{t,s,d} matchings in K_{4,5}
12 fig_partitions partition lattice |C|=3 with Mobius weights
13 fig_reduction  exponent staircase of (5.6)
14 fig_fourbyfour 4x4 cancellation cases
15 fig_blockB     anatomy of B(X)
16 fig_roots      root-of-unity filter, d=5
17 fig_params     parameter split of the n x n matrix
18 fig_boundchart computed circuit bound chart
19 fig_markedtree marked formula tree
20 fig_affinepath affine wrappers compose
21 fig_fivebyfive 5x5 warm-up anatomy
22 fig_blocks9    Lemma 9.1 block anatomy, l=3, m=7
23 fig_kron       the real 49x49 Jacobian heat map
24 fig_packing    cyclic matchings tile the 12x12 matrix
25 fig_flt        fractional linear wrappers
26 fig_schur      Schur compression funnel
27 fig_detchart   codimension asymmetry chart
28 fig_modelsmap  model hierarchy
29 fig_checks     validation surface grid
"""
import math
import os

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(
    __file__))), ".ignore", "ch05_figs")
os.makedirs(OUT, exist_ok=True)

RED = "#c0392b"
BLUE = "#2e6da4"
GREEN = "#2a7d2a"
GOLD = "#b8860b"
INK = "#333"
MUT = "#555"
BOX = "#f3f1ec"
ACC = "#7a1f1f"
LT = "#eef4fa"
LG = "#f1f7f1"


def svg(name, vb, body, aria):
    w, h = vb
    s = (f'<svg class="setupfig" viewBox="0 0 {w} {h}" width="100%" '
         f'role="img" aria-label="{aria}">\n' + body + '\n</svg>\n')
    with open(os.path.join(OUT, name + ".svg"), "w", encoding="utf-8") as f:
        f.write(s)
    print("wrote", name, len(s))


def txt(x, y, s, size=12, fill=INK, anchor="start", cls=None, weight=None):
    c = f' class="{cls}"' if cls else ""
    wgt = f' font-weight="{weight}"' if weight else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}"'
            f' fill="{fill}" text-anchor="{anchor}"{c}{wgt}>{s}</text>')


def sub(main, subscript, size=12):
    return (f'{main}<tspan baseline-shift="sub" '
            f'font-size="{size*0.72:.0f}">{subscript}</tspan>')


def sup(main, superscript, size=12):
    return (f'{main}<tspan baseline-shift="super" '
            f'font-size="{size*0.72:.0f}">{superscript}</tspan>')


def line(x1, y1, x2, y2, stroke=INK, w=1.4, dash=None, marker=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    mk = f' marker-end="url(#{marker})"' if marker else ""
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" '
            f'y2="{y2:.1f}" stroke="{stroke}" stroke-width="{w}"{d}{mk}/>')


def rect(x, y, w, h, fill=BOX, stroke=INK, sw=1.2, rx=0, opac=None):
    o = f' fill-opacity="{opac}"' if opac is not None else ""
    r = f' rx="{rx}"' if rx else ""
    st = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else \
        ' stroke="none"'
    return (f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" '
            f'height="{h:.1f}" fill="{fill}"{o}{st}{r}/>')


def circ(cx, cy, r, fill=INK, stroke=None, sw=1.2):
    st = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
    return f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r}" fill="{fill}"{st}/>'


def poly(pts, stroke=INK, w=1.6, fill="none", dash=None, close=False):
    p = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    tag = "polygon" if close else "polyline"
    return (f'<{tag} points="{p}" fill="{fill}" stroke="{stroke}" '
            f'stroke-width="{w}"{d}/>')


def arrow_defs(fid, colors):
    out = ["<defs>"]
    for i, col in enumerate(colors):
        out.append(
            f'<marker id="c5{fid}a{i}" markerWidth="9" markerHeight="9" '
            f'refX="7" refY="3" orient="auto">'
            f'<path d="M0,0 L7,3 L0,6 Z" fill="{col}"/></marker>')
    out.append("</defs>")
    return "\n".join(out)


# ---------------------------------------------------------------- fig 1
def fig_perdet():
    import itertools
    b = [arrow_defs("pd", [MUT])]
    b.append(txt(12, 20, "the 6 permutation patterns of a 3&#215;3 matrix "
                 "(one entry per row and column)", 12, MUT))
    perms = list(itertools.permutations(range(3)))
    labels = ["aei", "afh", "bdi", "bfg", "cdh", "ceg"]
    letters = [["a", "b", "c"], ["d", "e", "f"], ["g", "h", "i"]]
    x0 = 30
    for pi, pm in enumerate(perms):
        gx = x0 + pi*100
        gy = 36
        cell = 18
        for i in range(3):
            for j in range(3):
                chosen = (pm[i] == j)
                b.append(rect(gx + j*cell, gy + i*cell, cell, cell,
                              fill=(RED if chosen else "#fff"),
                              stroke="#999", sw=0.8))
                b.append(txt(gx + j*cell + 9, gy + i*cell + 13,
                             letters[i][j], 10.5,
                             "#fff" if chosen else "#888",
                             "middle",
                             cls="onfill v" if chosen else "v"))
        b.append(txt(gx + 27, gy + 70, "+" + labels[pi], 12, INK, "middle",
                     cls="v"))
    b.append(txt(12, 130, "the same six terms as perfect matchings "
                 "(rows &#8594; columns)", 12, MUT))
    for pi, pm in enumerate(perms):
        gx = x0 + pi*100
        gy = 142
        for i in range(3):
            b.append(circ(gx + 6, gy + i*16, 3.4, INK))
            b.append(circ(gx + 48, gy + i*16, 3.4, BLUE))
            b.append(line(gx + 10, gy + i*16, gx + 44, gy + pm[i]*16,
                          RED, 1.5))
    b.append(txt(12, 212, "permanent: add all six.  determinant: attach "
                 "sgn(&#963;) = &#177;1 first (afh, bdi, ceg flip sign).",
                 12, INK))
    svg("fig_perdet", (640, 224), "\n".join(b),
        "Permanent of a 3x3 matrix as six matchings")


# ---------------------------------------------------------------- fig 2
def fig_models():
    b = [arrow_defs("md", [MUT, RED])]

    def gate(x, y, lab, col=INK):
        return (circ(x, y, 11, "#fff", col, 1.6) +
                txt(x, y + 4.5, lab, 13, col, "middle"))

    def leaf(x, y, lab):
        return (rect(x - 10, y - 9, 20, 18, "#fff", "#999", 1, rx=4) +
                txt(x, y + 4, lab, 11.5, INK, "middle", cls="v"))

    # left: circuit
    b.append(txt(90, 18, "circuit: 4 gates", 12.5, ACC, "middle",
                 weight="bold"))
    xL, yv = 90, 0
    b.append(leaf(50, 210, "x"))
    b.append(leaf(110, 210, "y"))
    b.append(leaf(170, 210, "z"))
    b.append(gate(80, 160, "+"))     # u = x+y
    b.append(gate(80, 115, "&#215;"))   # v = u*u
    b.append(gate(125, 70, "+"))     # v+z
    b.append(gate(103, 30, "&#215;"))   # (v+z)*v
    for (x1, y1, x2, y2) in [(50, 201, 76, 170), (110, 201, 84, 170),
                             (80, 149, 80, 126), (80, 104, 76, 78),
                             (80, 104, 121, 79),
                             (170, 201, 129, 79), (125, 59, 108, 40),
                             (76, 108, 97, 36)]:
        b.append(line(x1, y1, x2, y2, MUT, 1.3, marker="c5mda0"))
    b.append(txt(48, 128, "u", 11, MUT, cls="v"))
    b.append(txt(60, 96, "v", 11, MUT, cls="v"))
    b.append(txt(30, 60, "v reused:", 11, RED))
    b.append(txt(30, 73, "2 out-edges", 11, RED))
    # middle: formula
    b.append(txt(320, 18, "formula: 7 gates, tree", 12.5, ACC, "middle",
                 weight="bold"))
    fx = 250
    b.append(gate(320, 30, "&#215;"))
    b.append(gate(275, 75, "+"))
    b.append(gate(370, 75, "&#215;"))
    b.append(gate(250, 120, "&#215;"))
    b.append(leaf(305, 120, "z"))
    b.append(gate(345, 120, "+"))
    b.append(gate(398, 120, "+"))
    b.append(gate(228, 165, "+"))
    b.append(gate(272, 165, "+"))
    b.append(leaf(220, 210, "x"))
    b.append(leaf(248, 210, "y"))
    b.append(leaf(282, 210, "x"))
    b.append(leaf(305, 210, "y"))
    b.append(leaf(333, 210, "x"))
    b.append(leaf(357, 210, "y"))
    b.append(leaf(390, 210, "x"))
    b.append(leaf(412, 210, "y"))
    for (x1, y1, x2, y2) in [(275, 64, 313, 38), (370, 64, 327, 38),
                             (250, 109, 268, 84), (305, 111, 281, 84),
                             (345, 109, 364, 84), (398, 109, 377, 84),
                             (228, 154, 244, 130), (272, 154, 256, 130),
                             (220, 201, 226, 176), (248, 201, 230, 176),
                             (282, 201, 268, 176), (305, 201, 276, 176),
                             (333, 201, 341, 131), (357, 201, 349, 131),
                             (390, 201, 396, 131), (412, 201, 400, 131)]:
        b.append(line(x1, y1, x2, y2, MUT, 1.2, marker="c5mda0"))
    b.append(txt(320, 232, "the shared (x+y)&#178; is recomputed",
                 11.5, RED, "middle"))
    # right: division formula
    b.append(txt(555, 18, "with division", 12.5, ACC, "middle",
                 weight="bold"))
    b.append(gate(555, 40, "&#247;"))
    b.append(leaf(520, 90, "1"))
    b.append(gate(590, 90, "+"))
    b.append(leaf(565, 140, "x"))
    b.append(leaf(615, 140, "y"))
    for (x1, y1, x2, y2) in [(520, 81, 548, 50), (590, 79, 562, 50),
                             (565, 131, 583, 100), (615, 131, 597, 100)]:
        b.append(line(x1, y1, x2, y2, MUT, 1.2, marker="c5mda0"))
    b.append(txt(555, 170, "1/(x+y): legal &#8212; the", 11.5, INK,
                 "middle"))
    b.append(txt(555, 184, "denominator x+y is not", 11.5, INK, "middle"))
    b.append(txt(555, 198, "the zero polynomial;", 11.5, INK, "middle"))
    b.append(txt(555, 212, "output must be polynomial", 11.5, INK,
                 "middle"))
    svg("fig_models", (640, 244), "\n".join(b),
        "Circuit versus formula versus division formula")


# ---------------------------------------------------------------- fig 3
def chart_axes(bx, by, bw, bh, xlab, ylab):
    out = [rect(bx, by, bw, bh, "#fff", "#bbb", 1)]
    out.append(txt(bx + bw/2, by + bh + 34, xlab, 11.5, MUT, "middle"))
    out.append(txt(max(4, bx - 62), by - 10, ylab, 11.5, MUT, "start"))
    return out


def fig_results():
    b = []
    # left panel: circuit advantage factor vs loglog n
    bx, by, bw, bh = 55, 40, 240, 170
    b += chart_axes(bx, by, bw, bh,
                    "log&#8322;log&#8322; n  (linear scale)",
                    "bound &#247; n&#178;")
    b.append(txt(175, 22, "circuits: advantage over the trivial n&#178;",
                 12, ACC, "middle", weight="bold"))
    xmin, xmax = 4.0, 160.0
    ymin, ymax = 0.0, 1.15

    def mx(v):
        return bx + (v - xmin)/(xmax - xmin)*bw

    def my(v):
        return by + bh - (v - ymin)/(ymax - ymin)*bh
    pts = []
    v = xmin
    while v <= xmax:
        pts.append((mx(v), my((v - 3.0)/144.0)))
        v += 2.0
    b.append(poly(pts, RED, 2.2))
    b.append(line(mx(xmin), my(1.0), mx(xmax), my(1.0), MUT, 1.3,
                  dash="6 4"))
    b.append(txt(mx(20), my(1.0) - 6, "trivial bound n&#178; &#8722; 1",
                 10.5, MUT))
    b.append(line(mx(147), by, mx(147), by + bh, "#999", 1, dash="2 3"))
    b.append(txt(mx(147), by + bh + 14, "147", 10, "#999", "middle"))
    for v, lab in [(4, "4"), (40, "40"), (80, "80"), (120, "120"),
                   (160, "160")]:
        b.append(txt(mx(v), by + bh + 14, lab, 10, MUT, "middle"))
    for v in [0.0, 0.5, 1.0]:
        b.append(txt(bx - 8, my(v) + 4, f"{v:.1f}", 10, MUT, "end"))
    b.append(txt(mx(60), my(0.28),
                 "(log&#8322;log&#8322;n &#8722; 3)/144", 11, RED))
    b.append(txt(mx(64), my(0.13), 'passes 1 only at n = 2'
                 '<tspan baseline-shift="super" font-size="8">2</tspan>'
                 '<tspan baseline-shift="super" font-size="6" dy="-4">147'
                 '</tspan>', 10.5, MUT))
    # right panel: formula bounds, log-log
    bx2, by2 = 360, 40
    b += chart_axes(bx2, by2, bw, bh, "n  (log scale)",
                    "leaf lower bound  (log scale)")
    b.append(txt(480, 22, "formulas: new vs classical", 12, ACC, "middle",
                 weight="bold"))
    lo, hi = math.log10(32), math.log10(1e8)

    def mx2(n):
        return bx2 + (math.log10(n) - lo)/(hi - lo)*bw
    yl, yh = 3.0, 31.0

    def my2(v):
        return by2 + bh - (math.log10(v) - yl)/(yh - yl)*bh
    p1, p2 = [], []
    n = 32.0
    while n <= 1e8:
        p1.append((mx2(n), my2(n**3)))
        p2.append((mx2(n), my2(n**4/(128*math.log2(n)))))
        n *= 1.25
    b.append(poly(p1, MUT, 1.8, dash="6 4"))
    b.append(poly(p2, RED, 2.2))
    ncross = 32.0
    while ncross**4/(128*math.log2(ncross)) < ncross**3:
        ncross += 1
    b.append(line(mx2(ncross), by2, mx2(ncross), by2 + bh, "#999", 1,
                  dash="2 3"))
    b.append(txt(mx2(ncross) + 3, by2 + 14, f"n = {int(ncross)}", 10,
                 "#999"))
    b.append(txt(mx2(1.2e4), my2(1e22), "n&#8308;/(128 log&#8322;n)",
                 11, RED, "middle"))
    b.append(txt(mx2(2.5e6), my2(1e12), "classical ~ n&#179;", 11, MUT,
                 "middle"))
    for n, lab in [(1e2, "10&#178;"), (1e4, "10&#8308;"),
                   (1e6, "10&#8310;"), (1e8, "10&#8312;")]:
        b.append(txt(mx2(n), by2 + bh + 14, lab, 10, MUT, "middle"))
    for e in [5, 15, 25]:
        b.append(txt(bx2 - 6, my2(10.0**e) + 4,
                     "10" + f'<tspan baseline-shift="super" '
                     f'font-size="7">{e}</tspan>', 10, MUT, "end"))
    svg("fig_results", (640, 258), "\n".join(b),
        "The proved bounds to scale")
    return int(ncross)


# ---------------------------------------------------------------- fig 4
def fig_roadmap():
    b = [arrow_defs("rm", [ACC])]

    def box(x, y, w, h, lines, col=BOX, tcol=INK):
        out = rect(x, y, w, h, col, "#999", 1.2, rx=7)
        for i, ln in enumerate(lines):
            out += txt(x + w/2, y + 16 + i*13.5, ln, 10.8, tcol, "middle")
        return out
    b.append(txt(160, 20, "circuit spine  (&#167;3&#8211;&#167;9)", 13,
                 ACC, "middle", weight="bold"))
    b.append(txt(490, 20, "formula spine  (&#167;10&#8211;&#167;12)", 13,
                 ACC, "middle", weight="bold"))
    b.append(box(20, 34, 130, 46, ["Lem 4.1 (&#167;4)", "B&#233;zout:",
                                   'e<tspan baseline-shift="super" '
                                   'font-size="7.5">k</tspan> &#8804; 2'
                                   '<tspan baseline-shift="super" '
                                   'font-size="7.5">q</tspan>']))
    b.append(box(170, 34, 130, 46, ["Lem 4.2 (&#167;5)", "slice + project",
                                    "the gradient"]))
    b.append(box(95, 106, 150, 46, ["Baur&#8211;Strassen (&#167;6)",
                                    "&#8711;P at 3&#215; cost"]))
    b.append(box(60, 178, 200, 40, ["Prop 4.3: codim k, degree d",
                                    "&#8658; 3L &#8805; k log&#8322;"
                                    "(d&#8722;1)"], LT))
    b.append(box(20, 244, 140, 46, ["Prop 5.1 + Cor 5.2 (&#167;7)",
                                    "critical loci of",
                                    "minor sums are small"]))
    b.append(box(180, 244, 140, 46, ["Lem 6.1 (&#167;8)",
                                     "block sum hides in",
                                     "one permanent"]))
    b.append(box(60, 316, 200, 42, ["&#167;9: d &#8776; &#188;log&#8322;n,"
                                    " k = &#920;(n&#178;)",
                                    "Thm 1.1: &#937;(n&#178;loglog n)"],
                 "#fdf7e3"))
    for (x1, y1, x2, y2) in [(85, 80, 140, 104), (235, 80, 190, 104),
                             (170, 152, 165, 176), (90, 290, 130, 314),
                             (250, 290, 210, 314), (160, 218, 160, 240)]:
        b.append(line(x1, y1, x2, y2, ACC, 1.4, marker="c5rma0"))
    b.append(box(350, 34, 130, 46, ["Lem 8.1 (&#167;10)", "tree occurrence",
                                    "bound: td &#8804; 4t"]))
    b.append(box(500, 34, 130, 46, ["Import 6 (&#167;3)", "Jacobian",
                                    "criterion"]))
    b.append(box(425, 106, 150, 46, ["Lem 9.1 (&#167;11)",
                                     "O(log n) matching has",
                                     "m&#178; independent coeffs"]))
    b.append(box(390, 178, 200, 40, ["&#167;12: pack &#957; &#8805; "
                                     "n&#178;/2k disjoint",
                                     "matchings; add the charges"], LT))
    b.append(box(350, 244, 130, 46, ["Lem 10.1 + 10.2", "division: td "
                                     "&#8804; 6t", "(&#167;12.2)"]))
    b.append(box(500, 244, 130, 46, ["Thms 1.2, 1.3",
                                     "&#937;(n&#8308;/log n)",
                                     "leaves"], "#fdf7e3"))
    for (x1, y1, x2, y2) in [(430, 80, 470, 104), (550, 80, 520, 104),
                             (495, 152, 490, 176), (455, 218, 425, 242),
                             (520, 218, 555, 242), (480, 267, 498, 267)]:
        b.append(line(x1, y1, x2, y2, ACC, 1.4, marker="c5rma0"))
    b.append(txt(480, 320, "division only relaxes 4t&#8722;2 to "
                 "6t&#8722;3", 10.8, MUT, "middle"))
    b.append(txt(480, 336, "(fractional linear wrappers + field "
                 "intersection)", 10.8, MUT, "middle"))
    svg("fig_roadmap", (650, 372), "\n".join(b),
        "Dependency map of both proofs")


# ---------------------------------------------------------------- fig 5
def fig_agdim():
    b = [arrow_defs("ag", [MUT])]
    # 3D cone, orthographic projection
    cx, cy = 105, 118
    az, el = 0.6, 0.32

    def proj(x, y, z):
        xr = x*math.cos(az) + y*math.sin(az)
        yr = -x*math.sin(az) + y*math.cos(az)
        return (cx + 52*xr, cy - 52*(z*math.cos(el) - yr*math.sin(el)))
    b.append(txt(105, 20, "a cone V &#8834; &#8450;&#179; (dim 2)", 12,
                 ACC, "middle", weight="bold"))
    rim = [proj(math.cos(t), math.sin(t), 1.0)
           for t in [i*2*math.pi/60 for i in range(61)]]
    rim2 = [proj(-math.cos(t), -math.sin(t), -1.0)
            for t in [i*2*math.pi/60 for i in range(61)]]
    b.append(poly(rim, BLUE, 1.6))
    b.append(poly(rim2, BLUE, 1.2, dash="4 3"))
    for t in [i*2*math.pi/8 for i in range(8)]:
        x, y = math.cos(t), math.sin(t)
        p1 = proj(x, y, 1.0)
        p2 = proj(-x, -y, -1.0)
        b.append(line(p1[0], p1[1], p2[0], p2[1], "#8fb3d4", 0.9))
    b.append(circ(*proj(0, 0, 0), 3, INK))
    b.append(txt(proj(0, 0, 0)[0] + 7, proj(0, 0, 0)[1] + 11, "0", 11,
                 INK))
    b.append(txt(105, 205, "each point rides a line through 0", 10.5,
                 MUT, "middle"))
    # middle: projectivization
    b.append(txt(320, 20, "&#8473;V &#8834; &#8473;&#178; (dim 1)", 12,
                 ACC, "middle", weight="bold"))
    b.append(circ(320, 118, 58, "#fff", "#999", 1.2))
    ppts = [(320 + 44*math.cos(t), 118 + 30*math.sin(t))
            for t in [i*2*math.pi/60 for i in range(61)]]
    b.append(poly(ppts, BLUE, 2))
    b.append(txt(320, 205, "directions only: one dimension spent", 10.5,
                 MUT, "middle"))
    b.append(line(180, 118, 245, 118, MUT, 1.2, marker="c5aga0"))
    b.append(txt(212, 106, "[&#183;]", 12, MUT, "middle"))
    # right: avoidance count
    b.append(txt(530, 20, "Lemma A: dodge if l + r &lt; N", 12, ACC,
                 "middle", weight="bold"))
    b.append(rect(430, 40, 200, 150, "#fff", "#999", 1.2, rx=6))
    spts = [(455 + 40*math.cos(t) + 18*math.sin(2.1*t),
             100 + 30*math.sin(t)) for t in
            [i*2*math.pi/40 for i in range(41)]]
    b.append(poly(spts, RED, 2))
    b.append(txt(457, 62, "S, dim r", 11, RED))
    b.append(line(560, 50, 610, 180, GREEN, 2.4))
    b.append(txt(600, 66, "generic", 10.5, GREEN, "end"))
    b.append(txt(608, 80, "l-plane", 10.5, GREEN, "end"))
    b.append(txt(530, 205, "planes through a point: codim N &#8722; l",
                 10.5, MUT, "middle"))
    svg("fig_agdim", (650, 218), "\n".join(b),
        "Cone, projectivization, avoidance principle")


# ---------------------------------------------------------------- fig 6
def fig_powerfiber():
    b = []
    d = 4

    def panel(cx, cy, a, lab):
        out = [circ(cx, cy, 2.5, INK)]
        out.append(line(cx - 62, cy, cx + 62, cy, "#bbb", 1))
        out.append(line(cx, cy - 62, cx, cy + 62, "#bbb", 1))
        rr = abs(a)**(1.0/3)
        th0 = math.atan2(a.imag, a.real)/3
        pts = []
        for j in range(3):
            th = th0 + j*2*math.pi/3
            x, y = cx + 26*rr*math.cos(th), cy - 26*rr*math.sin(th)
            pts.append((x, y))
            out.append(circ(x, y, 4, RED))
        out.append(poly(pts + [pts[0]], "#e8b7b0", 1, dash="3 3"))
        out.append(txt(cx, cy + 78, lab, 11.5, INK, "middle"))
        return out, pts
    p1, r1 = panel(95, 92, complex(1, 0),
                   "roots of z&#179; = a&#8321; = 1")
    p2, r2 = panel(255, 92, complex(0, 8),
                   "roots of z&#179; = a&#8322; = 8i")
    b += p1 + p2
    b.append(txt(480, 24, "the fiber: all (d&#8722;1)&#178; = 9 pairs",
                 12, ACC, "middle", weight="bold"))
    for i in range(3):
        for j in range(3):
            x = 405 + i*58
            y = 60 + j*52
            b.append(circ(x, y, 5.5, "#fff", RED, 1.8))
            b.append(circ(x - 2, y - 1.5, 1.8, RED))
            b.append(circ(x + 2.2, y + 1.8, 1.8, BLUE))
    b.append(rect(375, 36, 208, 152, "none", "#999", 1.1, rx=8))
    b.append(txt(480, 208, '9 isolated solutions of &#8711;H = a '
                 '&#8804; 2<tspan baseline-shift="super" '
                 'font-size="7.5">q</tspan> &#8658; q &#8805; '
                 'log&#8322;9', 11.5, INK, "middle"))
    svg("fig_powerfiber", (640, 222), "\n".join(b),
        "Computed fiber of the power-sum gradient")


# ---------------------------------------------------------------- fig 7
def fig_gateeq():
    b = [arrow_defs("ge", [MUT])]
    b.append(txt(130, 22, "circuit, q multiplications", 12.5, ACC,
                 "middle", weight="bold"))

    def cloud(x, y, w, lab):
        return (rect(x, y, w, 24, LT, BLUE, 1.1, rx=12) +
                txt(x + w/2, y + 16, lab, 10.5, BLUE, "middle"))

    def mgate(x, y, lab):
        return (circ(x, y, 12, "#fff", RED, 1.8) +
                txt(x, y + 4, "&#215;", 14, RED, "middle") +
                txt(x + 18, y + 4, lab, 11, RED, cls="v"))
    b.append(txt(60, 208, "u&#8321; &#8230; u&#8342;", 12, INK, "middle",
                 cls="v"))
    b.append(cloud(30, 158, 130, "affine (free)"))
    b.append(mgate(95, 122, "v&#8321;"))
    b.append(cloud(30, 82, 130, "affine (free)"))
    b.append(mgate(95, 48, "v&#8322;"))
    for (x1, y1, x2, y2) in [(60, 199, 80, 184), (95, 156, 95, 136),
                             (95, 110, 95, 108), (95, 106, 95, 108),
                             (95, 70, 95, 62)]:
        b.append(line(x1, y1, x2, y2, MUT, 1.2, marker="c5gea0"))
    b.append(line(120, 199, 140, 184, MUT, 1.2, marker="c5gea0"))
    b.append(txt(230, 60, "&#8658;", 22, INK, "middle"))
    b.append(rect(280, 34, 330, 168, BOX, "#999", 1.2, rx=8))
    b.append(txt(445, 58, "the fiber, described by the circuit:", 11.5,
                 INK, "middle"))
    b.append(txt(300, 86, "v&#8342; = A&#8342;(u, v&#8249;&#8342;) &#183; "
                 "B&#8342;(u, v&#8249;&#8342;)", 12.5, RED))
    b.append(txt(475, 86, "&#215; q,  degree &#8804; 2", 11, MUT))
    b.append(txt(300, 116, "H&#7522;(u, v) = a&#7522;", 12.5, BLUE))
    b.append(txt(475, 116, "&#215; k,  degree &#8804; 1", 11, MUT))
    b.append(txt(300, 152, 'affine B&#233;zout (Import 4): &#8804; 2'
                 '<tspan baseline-shift="super" font-size="7.5">q'
                 '</tspan> &#183; 1<tspan baseline-shift="super" '
                 'font-size="7.5">k</tspan> isolated solutions',
                 11.5, INK))
    b.append(txt(300, 178, "gate equations solve v from u uniquely "
                 "&#8658; bijection with fiber", 11, MUT))
    svg("fig_gateeq", (640, 220), "\n".join(b),
        "Gate equations cap the fiber size")


print("part 1 done")


# ---------------------------------------------------------------- fig 8
def fig_slice():
    b = []
    cx, cy = 200, 130
    az, el = 0.55, 0.38

    def proj(x, y, z):
        xr = x*math.cos(az) + y*math.sin(az)
        yr = -x*math.sin(az) + y*math.cos(az)
        return (cx + 62*xr, cy - 62*(z*math.cos(el) - yr*math.sin(el)))
    b.append(txt(200, 22, "Crit(x&#8321;x&#8322;x&#8323;) = the three "
                 "axes (a 1-dim cone)", 12.5, ACC, "middle",
                 weight="bold"))
    # generic plane spanned by v1=(1,.3,.55), v2=(-.25,1,.4)
    v1 = (1.0, 0.3, 0.55)
    v2 = (-0.25, 1.0, 0.4)
    corners = []
    for (sa, sb) in [(1.35, 1.35), (1.35, -1.35), (-1.35, -1.35),
                     (-1.35, 1.35)]:
        x = sa*v1[0] + sb*v2[0]
        y = sa*v1[1] + sb*v2[1]
        z = sa*v1[2] + sb*v2[2]
        corners.append(proj(x, y, z))
    b.append(poly(corners, GREEN, 1.4, fill="#dcecdc", close=True))
    b.append(f'<polygon points="' +
             " ".join(f"{x:.1f},{y:.1f}" for x, y in corners) +
             f'" fill="{GREEN}" fill-opacity="0.13" stroke="none"/>')
    # axes
    for vec, lab in [((1.8, 0, 0), "x&#8321;-axis"),
                     ((0, 1.8, 0), "x&#8322;-axis"),
                     ((0, 0, 1.55), "x&#8323;-axis")]:
        p1 = proj(*vec)
        p2 = proj(-vec[0], -vec[1], -vec[2])
        b.append(line(p2[0], p2[1], p1[0], p1[1], RED, 2.2))
        b.append(txt(p1[0] + 6, p1[1] - 4, lab, 10.5, RED))
    o = proj(0, 0, 0)
    b.append(circ(o[0], o[1], 3.5, INK))
    b.append(txt(o[0] + 8, o[1] + 13, "0 = only meeting point", 11, INK))
    p1 = proj(1.35*v1[0] + 1.35*v2[0], 1.35*v1[1] + 1.35*v2[1],
              1.35*v1[2] + 1.35*v2[2])
    b.append(txt(p1[0] - 4, p1[1] - 8, "W(&#8450;&#178;), a generic "
                 "2-dim slice", 11, GREEN, "end"))
    b.append(txt(200, 248, "k + &#948; = 2 + 1 &#8804; 3 = m: the plane "
                 "dodges the three critical lines", 11.5, MUT, "middle"))
    svg("fig_slice", (420, 262), "\n".join(b),
        "Generic plane meets the critical cone only at the origin")


# ---------------------------------------------------------------- fig 9
def fig_project():
    b = [arrow_defs("pj", [MUT, RED])]
    b.append(txt(95, 26, "&#8473;", 15, INK, "middle")
             + txt(104, 20, "k&#8722;1", 9, INK, "middle"))
    b.append(circ(95, 110, 55, "#fff", "#999", 1.3))
    b.append(txt(95, 190, "sliced input directions [u]", 10.5, MUT,
                 "middle"))
    b.append(line(160, 110, 235, 110, MUT, 1.3, marker="c5pja0"))
    b.append(txt(197, 98, "[&#8711;P(Wu)]", 11.5, MUT, "middle"))
    b.append(txt(390, 26, "&#8473;", 15, INK, "middle")
             + txt(400, 20, "m&#8722;1", 9, INK, "middle"))
    b.append(rect(250, 45, 280, 130, "#fff", "#999", 1.3, rx=10))
    im = [(300 + 90*t + 25*math.sin(6.28*t), 95 + 30*math.sin(3.1*t))
          for t in [i/40 for i in range(41)]]
    b.append(poly(im, RED, 2.2))
    b.append(txt(346, 66, "image, dim &#8804; k&#8722;1", 10.5, RED))
    b.append(line(268, 160, 512, 148, BLUE, 2))
    b.append(txt(500, 142, "&#8473;K, dim m&#8722;k&#8722;1", 10.5, BLUE,
                 "end"))
    b.append(txt(390, 196, "(m&#8722;k&#8722;1) + (k&#8722;1) = m&#8722;2 "
                 "&lt; m&#8722;1 &#8658; a generic K misses the image;",
                 10.8, INK, "middle"))
    b.append(txt(390, 212, "any rank-k map A with kernel K then kills no "
                 "gradient value", 10.8, INK, "middle"))
    svg("fig_project", (560, 226), "\n".join(b),
        "Projecting the gradient without creating zeros")


# --------------------------------------------------------------- fig 10
def fig_reverse():
    b = [arrow_defs("rv", [MUT, RED, BLUE])]
    DUR = 7.0

    def anim(vals, times):
        v = ";".join(str(x) for x in vals)
        t = ";".join(f"{x:.3f}" for x in times)
        return (f'<animate attributeName="opacity" values="{v}" '
                f'keyTimes="{t}" dur="{DUR}s" repeatCount="indefinite"/>')

    def gate(x, y, lab, name):
        return (circ(x, y, 13, "#fff", INK, 1.6) +
                txt(x, y + 4.5, lab, 13, INK, "middle") +
                txt(x + 19, y - 6, name, 11, MUT, cls="v"))

    def leaf(x, y, lab):
        return (rect(x - 11, y - 10, 22, 20, "#fff", "#999", 1.1, rx=4) +
                txt(x, y + 4, lab, 12, INK, "middle", cls="v"))
    # static skeleton
    b.append(txt(160, 20, "P = w&#183;v,  w = v&#183;x,  v = u&#183;u,  "
                 "u = x + y", 12.5, ACC, "middle", weight="bold"))
    nodes = {"x": (80, 272), "y": (170, 272), "u": (125, 217),
             "v": (125, 160), "w": (68, 104), "P": (100, 58)}
    edges = [("x", "u"), ("y", "u"), ("u", "v"), ("u", "v"),
             ("v", "w"), ("x", "w"), ("w", "P"), ("v", "P")]
    for (a, c) in edges:
        xa, ya = nodes[a]
        xc, yc = nodes[c]
        b.append(line(xa, ya - 12, xc, yc + 14, "#aaa", 1.2))
    b.append(leaf(*nodes["x"], "x"))
    b.append(leaf(*nodes["y"], "y"))
    b.append(gate(*nodes["u"], "+", "u"))
    b.append(gate(*nodes["v"], "&#215;", "v"))
    b.append(gate(*nodes["w"], "&#215;", "w"))
    b.append(gate(*nodes["P"], "&#215;", "P"))
    # phase 1: forward pulses on the three mult gates
    for i, nm in enumerate(["v", "w", "P"]):
        x, y = nodes[nm]
        t0 = 0.04 + i*0.09
        g = (f'<g opacity="0">' +
             circ(x, y, 17, "none", RED, 2.6) +
             anim([0, 0, 1, 1, 0, 0],
                  [0, t0, t0 + 0.02, t0 + 0.07, t0 + 0.1, 1]) + '</g>')
        b.append(g)
    fw = (f'<g>' + txt(258, 60, "phase 1: forward sweep,", 11.5, RED) +
          txt(258, 76, "store every product", 11.5, RED) +
          anim([1, 1, 0, 0], [0, 0.33, 0.37, 1]) + '</g>')
    b.append(fw)
    # phase 2: adjoints appear top-down
    adj = [("P", "1", 0.40), ("w", "v&#773;=v", 0.48),
           ("v", "v&#773;=w+u&#183;&#8230;", 0.56), ("u", "2u&#183;v&#773;",
                                                     0.64),
           ("x", "&#8706;P/&#8706;x", 0.74), ("y", "&#8706;P/&#8706;y",
                                              0.74)]
    for nm, lab, t0 in adj:
        x, y = nodes[nm]
        g = ('<g opacity="0">' +
             rect(x - 46, y - 8, 40, 17, "#eef4fa", BLUE, 1.1, rx=3) +
             txt(x - 26, y + 4.5, lab, 9.5, BLUE, "middle") +
             anim([0, 0, 1, 1, 1, 0], [0, t0, t0 + 0.03, 0.78, 0.97, 1])
             + '</g>')
        b.append(g)
    ph2 = ('<g opacity="0">' +
           txt(258, 108, "phase 2: reverse sweep &#8212;", 11.5, BLUE) +
           txt(258, 124, "a&#773; += v&#773;&#183;b,  b&#773; += "
               "v&#773;&#183;a", 11.5, BLUE) +
           txt(258, 140, "(&#8804; 2 mults per &#215; gate)", 11.5, BLUE) +
           anim([0, 0, 1, 1, 0, 0], [0, 0.38, 0.42, 0.72, 0.76, 1]) +
           '</g>')
    b.append(ph2)
    ph3 = ('<g opacity="0">' +
           txt(258, 176, "phase 3: every leaf holds its", 11.5, GREEN) +
           txt(258, 192, "derivative &#8212; the whole &#8711;P", 11.5,
               GREEN) +
           txt(258, 208, "for 3q = 9 multiplications", 11.5, GREEN) +
           anim([0, 0, 1, 1, 0], [0, 0.76, 0.8, 0.97, 1]) + '</g>')
    b.append(ph3)
    svg("fig_reverse", (450, 300), "\n".join(b),
        "Reverse-mode differentiation animation")


# --------------------------------------------------------------- fig 11
def fig_minorsum():
    b = []
    t, s = 4, 5
    matchings = [({0: 0, 1: 2, 3: 4}, RED), ({1: 1, 2: 3, 3: 0}, BLUE),
                 ({0: 3, 2: 2, 3: 1}, GREEN)]
    b.append(txt(160, 20, "three of the 240 matchings summed by "
                 "M&#8324;,&#8325;,&#8323;", 12, ACC, "middle",
                 weight="bold"))
    for pi, (mt, col) in enumerate(matchings):
        gx = 40 + pi*105
        for i in range(t):
            b.append(circ(gx, 44 + i*30, 4, INK))
        for j in range(s):
            b.append(circ(gx + 64, 38 + j*27, 4, BLUE))
        for i, j in mt.items():
            b.append(line(gx + 4, 44 + i*30, gx + 60, 38 + j*27, col, 2))
        b.append(txt(gx - 12, 36, "rows", 9.5, MUT, "middle"))
        b.append(txt(gx + 76, 30, "cols", 9.5, MUT, "middle"))
    b.append(txt(160, 178, "weight of a matching: &#8719; x&#8342;,"
                 "&#966;(j) over its d edges", 11, MUT, "middle"))
    # right: derivative
    gx = 430
    b.append(txt(505, 20, "&#8706;/&#8706;x&#8322;&#8322; freezes one "
                 "edge", 12, ACC, "middle", weight="bold"))
    for i in range(t):
        b.append(circ(gx, 44 + i*30, 4, INK if i != 1 else RED))
    for j in range(s):
        b.append(circ(gx + 64, 38 + j*27, 4, BLUE if j != 1 else RED))
    b.append(line(gx + 4, 74, gx + 60, 65, RED, 2.6))
    b.append(txt(gx + 30, 60, "(i,a) fixed", 9.5, RED, "middle"))
    for (i, j, col) in [(0, 0, "#bbb"), (2, 3, "#bbb"), (3, 4, "#bbb")]:
        b.append(line(gx + 4, 44 + i*30, gx + 60, 38 + j*27, col, 1.6,
                      dash="4 3"))
    b.append(txt(505, 178, "q = d&#8722;1 other rows, injected into "
                 "columns &#8800; a", 11, MUT, "middle"))
    svg("fig_minorsum", (640, 192), "\n".join(b),
        "Matchings summed by the minor-sum polynomial")


# --------------------------------------------------------------- fig 12
def fig_partitions():
    b = []
    parts = [([["1"], ["2"], ["3"]], 1, ["identity"], 1),
             ([["1", "2"], ["3"]], -1, ["(12)"], 1),
             ([["1", "3"], ["2"]], -1, ["(13)"], 1),
             ([["2", "3"], ["1"]], -1, ["(23)"], 1),
             ([["1", "2", "3"]], 2, ["(123)", "(132)"], -1)]
    b.append(txt(320, 20, "&#931;&#960; &#956;(&#960;) over partitions of "
                 "{1,2,3}  =  &#931;&#963; sgn(&#963;) over S&#8323;  =  "
                 "0", 12.5, ACC, "middle", weight="bold"))
    for pi, (blocks, mu, perms, sgn) in enumerate(parts):
        gx = 42 + pi*122
        b.append(rect(gx, 40, 104, 128, "#fff", "#999", 1.1, rx=7))
        y = 62
        for blk in blocks:
            wblk = 18*len(blk) + 10
            b.append(rect(gx + 12, y - 13, wblk, 22, LT, BLUE, 1, rx=10))
            for bi, el in enumerate(blk):
                b.append(txt(gx + 21 + bi*18, y + 3, el, 11.5, INK,
                             "middle"))
            y += 30
        b.append(txt(gx + 52, 138, f'&#956; = {mu:+d}', 12.5,
                     RED if mu < 0 else GREEN, "middle", weight="bold"))
        b.append(txt(gx + 52, 158, " ".join(perms), 10.5, MUT, "middle"))
    b.append(txt(320, 190, "weights +1 &#8722;1 &#8722;1 &#8722;1 +2 sum "
                 "to 0: any collision fiber of size &#8805; 2 cancels; "
                 "only injections survive", 11, MUT, "middle"))
    svg("fig_partitions", (640, 204), "\n".join(b),
        "Partition Mobius weights cancel non-injections")


# --------------------------------------------------------------- fig 13
def fig_reduction():
    b = [arrow_defs("rd", [RED])]
    q = 3
    bx, by, cell = 60, 40, 34
    b.append(txt(300, 20, "the rewriting (5.6): every exponent "
                 "&#8805; q collapses", 12.5, ACC, "middle",
                 weight="bold"))
    for e1 in range(6):
        for e2 in range(6):
            x = bx + e1*cell
            y = by + (5 - e2)*cell
            inbox = e1 < q and e2 < q
            b.append(rect(x, y, cell, cell,
                          LG if inbox else "#f7f0ef",
                          "#ccc", 0.7))
            b.append(circ(x + cell/2, y + cell/2, 2.2,
                          GREEN if inbox else "#c9a09a"))
    b.append(rect(bx, by + 3*cell, 3*cell, 3*cell, "none", GREEN, 2))
    for (f1, f2, t1, t2) in [(4, 1, 1, 1), (3, 3, 2, 3), (1, 4, 1, 1)]:
        x1 = bx + f1*cell + cell/2
        y1 = by + (5 - f2)*cell + cell/2
        x2 = bx + t1*cell + cell/2
        y2 = by + (5 - t2)*cell + cell/2
        b.append(f'<path d="M{x1:.0f},{y1:.0f} Q{(x1+x2)/2:.0f},'
                 f'{(y1+y2)/2 - 26:.0f} {x2:.0f},{y2:.0f}" fill="none" '
                 f'stroke="{RED}" stroke-width="1.8" '
                 f'marker-end="url(#c5rda0)"/>')
    b.append(txt(bx + 3*cell + 8, by + 3*cell - 8,
                 "u&#7522;-exponent &#8805; q: rewrite by (5.6)", 10.5,
                 RED))
    b.append(txt(bx + 1.5*cell, by + 6*cell + 18,
                 sub("exponent of u", "1", 11), 11, MUT, "middle"))
    b.append(txt(bx - 14, by + 1.5*cell,
                 sub("u", "2", 11), 11, MUT, "middle"))
    b.append(txt(bx + 1.5*cell, by + 4.4*cell,
                 "surviving box", 10, GREEN, "middle"))
    b.append(txt(bx + 1.5*cell, by + 4.4*cell + 13,
                 "{0,&#8230;,q&#8722;1}&#7511;", 10, GREEN, "middle"))
    b.append(txt(430, 110, "each arrow strictly lowers total", 11.5, INK))
    b.append(txt(430, 126, "column degree &#8658; termination;", 11.5,
                 INK))
    b.append(txt(430, 142, "what remains is a finite set of", 11.5, INK))
    b.append(txt(430, 158, "monomials &#8658; A is a finite", 11.5, INK))
    b.append(txt(430, 174, "&#8492;-module &#8658; dim A &#8804; dim "
                 "&#8492;", 11.5, INK))
    svg("fig_reduction", (620, 268), "\n".join(b),
        "Exponent staircase of the degree reduction")


# --------------------------------------------------------------- fig 14
def fig_fourbyfour():
    b = []
    mat = [["u", "v", "1", "1"], ["w", "z", "1", "1"],
           ["p", "q", "2", "&#8722;2"], ["r", "s", "2", "&#8722;2"]]
    cases = [([0, 1], "per", "&#8722;8", "upper&#8202;+&#8202;upper",
              "(uz+vw)&#183;(&#8722;8)"),
             ([2, 3], "per", "2", "lower&#8202;+&#8202;lower",
              "(ps+qr)&#183;2"),
             ([0, 2], "per", "0", "mixed: dies",
              "anything &#183; 0 = 0")]
    for ci, (rows, _, val, lab, prod) in enumerate(cases):
        gx = 30 + ci*205
        b.append(txt(gx + 72, 24, lab, 12, ACC, "middle", weight="bold"))
        cell = 26
        for i in range(4):
            for j in range(4):
                sel = (i in rows and j < 2)
                comp = (i not in rows and j >= 2)
                fill = RED if sel else (BLUE if comp else "#fff")
                b.append(rect(gx + j*cell, 36 + i*cell, cell, cell,
                              fill, "#999", 0.8,
                              opac=0.25 if (sel or comp) else None))
                b.append(txt(gx + j*cell + 13, 36 + i*cell + 17,
                             mat[i][j], 11,
                             INK, "middle",
                             cls="v" if j < 2 else None))
        b.append(txt(gx + 52, 168, "variable rows &#8594; cols 1,2",
                     9.5, RED, "middle"))
        b.append(txt(gx + 72, 186, "complementary constant minor: per = "
                     + val, 10, BLUE, "middle"))
        b.append(txt(gx + 72, 204, prod, 11, INK, "middle", cls="v"))
    svg("fig_fourbyfour", (650, 220), "\n".join(b),
        "The three cases of the 4x4 cancellation")


# --------------------------------------------------------------- fig 15
def fig_blockB():
    b = []
    # to visual scale: r=9 (b=3 blocks of t=3), s=8, d=2 => n=15
    t, bb, d, s = 3, 3, 2, 8
    r = t*bb
    u = 16
    x0, y0 = 40, 46
    cols = [GOLD, BLUE, GREEN]
    b.append(txt(x0 + (s + r - d)*u/2, 24,
                 "B(X), square of size r + s &#8722; d", 13, ACC,
                 "middle", weight="bold"))
    for h in range(bb):
        b.append(rect(x0, y0 + h*t*u, s*u, t*u, cols[h], "none",
                      opac=0.16))
    b.append(rect(x0, y0, s*u, r*u, "none", INK, 1.5))
    b.append(txt(x0 + s*u/2, y0 + r*u/2 + 4, "X  (r &#215; s variables)",
                 12.5, INK, "middle"))
    for h in range(bb):
        b.append(txt(x0 - 8, y0 + h*t*u + t*u/2 + 4, sub("R", str(h + 1),
                                                         11), 11,
                     cols[h], "end"))
    ux0 = x0 + s*u
    for h in range(bb):
        b.append(rect(ux0 + h*(t - d)*u, y0, (t - d)*u, r*u, "#fff",
                      "#999", 0.8))
        b.append(rect(ux0 + h*(t - d)*u, y0 + h*t*u, (t - d)*u, t*u,
                      cols[h], "none", opac=0.45))
    sat = bb*(t - d)
    for h in range(1, bb):
        gx = ux0 + sat*u + (h - 1)*d*u
        b.append(rect(gx, y0, d*u, r*u, "#fff", "#999", 0.8))
        b.append(rect(gx, y0, d*u, h*t*u, "#d7c4e8", "none", opac=0.8))
        b.append(rect(gx, y0 + h*t*u, d*u, t*u, "#9b59b6", "none",
                      opac=0.55))
    b.append(rect(ux0, y0, (r - d)*u, r*u, "none", INK, 1.5))
    b.append(txt(ux0 + (r - d)*u/2, y0 - 8,
                 "U: saturators, then root-of-unity filters", 10,
                 "#9b59b6", "middle"))
    yb = y0 + r*u
    b.append(rect(x0, yb, s*u, (s - d)*u, "#fbeee6", INK, 1.5))
    b.append(txt(x0 + s*u/2, yb + (s - d)*u/2 + 4, "all ones", 12, INK,
                 "middle"))
    b.append(rect(ux0, yb, (r - d)*u, (s - d)*u, "#eee", INK, 1.5))
    b.append(txt(ux0 + (r - d)*u/2, yb + (s - d)*u/2 + 4, "0", 13, INK,
                 "middle"))
    b.append(txt(x0 + s*u/2, yb + (s - d)*u + 16,
                 "s columns", 10.5, MUT, "middle"))
    b.append(txt(ux0 + (r - d)*u/2, yb + (s - d)*u + 16,
                 "r &#8722; d columns of U", 10.5, MUT, "middle"))
    b.append(txt((x0 + ux0 + (r - d)*u)/2, yb + (s - d)*u + 34,
                 "filter column of block h: 1 on R&#8321;&#8230;R"
                 "&#8342;&#8331;&#8321;, weight 2&#950;&#690; on "
                 "R&#8342;", 10, "#9b59b6", "middle"))
    svg("fig_blockB", (302, y0 + r*u + (s - d)*u + 48), "\n".join(b),
        "Anatomy of the bordered specialization matrix")


# --------------------------------------------------------------- fig 16
def fig_roots():
    b = []
    d = 5
    cx, cy, R = 120, 128, 62
    b.append(txt(cx, 22, "the d = 5 filter weights 2&#950;&#690;", 12.5,
                 ACC, "middle", weight="bold"))
    b.append(circ(cx, cy, R, "none", "#bbb", 1.1))
    b.append(line(cx - 85, cy, cx + 85, cy, "#ddd", 1))
    b.append(line(cx, cy - 85, cx, cy + 85, "#ddd", 1))
    sups = {0: "", 1: "", 2: "&#178;", 3: "&#179;", 4: "&#8308;"}
    for j in range(d):
        th = 2*math.pi*j/d
        x, y = cx + R*math.cos(th), cy - R*math.sin(th)
        b.append(line(cx, cy, x, y, "#e8b7b0", 1))
        b.append(circ(x, y, 4.5, RED))
        lab = "2" if j == 0 else "2&#950;" + sups[j]
        b.append(txt(cx + (R + 16)*math.cos(th),
                     cy - (R + 16)*math.sin(th) + 4, lab,
                     10.5, INK, "middle"))
    b.append(txt(cx, 232, "radius 2, angles 2&#960;j/5", 10.5, MUT,
                 "middle"))
    # right: elementary symmetric magnitudes
    ex0, ey0 = 290, 60
    b.append(txt(430, 22, "|2&#8305;&#8202;e&#7522;(&#950;&#8304;,"
                 "&#8230;,&#950;&#8308;)| in &#8719;(Y + 2&#950;&#690;"
                 "y&#8342;)", 12, ACC, "middle", weight="bold"))
    import cmath
    zs = [cmath.exp(2j*cmath.pi*j/5) for j in range(5)]
    import itertools as it
    vals = []
    for i in range(6):
        tot = 0
        for c in it.combinations(range(5), i):
            p = 1
            for j in c:
                p *= zs[j]
            tot += p
        vals.append(abs(2**i*tot))
    bw = 34
    for i, v in enumerate(vals):
        h = 4 + v*4.4
        x = ex0 + i*46
        b.append(rect(x, ey0 + 150 - h, bw, h,
                      RED if v > 0.5 else "#ccc", "none"))
        b.append(txt(x + bw/2, ey0 + 166, f"i={i}", 10, MUT, "middle"))
        b.append(txt(x + bw/2, ey0 + 150 - h - 6, f"{v:.0f}", 10.5,
                     INK, "middle"))
    b.append(txt(430, ey0 + 188, "only i = 0 (gives Y&#8309;) and i = 5 "
                 "(gives &#952;y&#8342;&#8309;, &#952; = 32)", 10.5, MUT,
                 "middle"))
    b.append(txt(430, ey0 + 202, "survive: every mixed term cancels "
                 "exactly", 10.5, MUT, "middle"))
    svg("fig_roots", (620, 274), "\n".join(b),
        "Root-of-unity filter kills mixed terms")


def abs_prod(c):
    p = 1
    for z in c:
        p *= z
    return p


print("part 2 done")


# --------------------------------------------------------------- fig 17
def fig_params():
    b = []
    n, d = 256, 3
    t = 4*d
    bb = n//(2*t)
    r = bb*t
    s = n - r + d
    sc = 340.0/n
    x0, y0 = 120, 40
    b.append(txt(x0 + n*sc/2, 24, f"n = {n}, d = {d}: t = {t}, b = {bb},"
                 f" r = {r}, s = {s}, m = rs = {r*s}", 12, ACC,
                 "middle", weight="bold"))
    cols = [GOLD, BLUE, GREEN]
    for h in range(bb):
        b.append(rect(x0, y0 + h*t*sc, s*sc, t*sc, cols[h % 3], "none",
                      opac=0.22))
    b.append(rect(x0, y0, s*sc, r*sc, "none", INK, 1.4))
    b.append(txt(x0 + s*sc/2, y0 + r*sc/2 + 5,
                 "X: r &#215; s variables", 12.5, INK, "middle"))
    b.append(rect(x0 + s*sc, y0, (r - d)*sc, r*sc, "#eee", INK, 1.2))
    b.append(txt(x0 + s*sc + (r - d)*sc/2, y0 + r*sc/2 + 5, "U", 13, INK,
                 "middle"))
    b.append(rect(x0, y0 + r*sc, s*sc, (s - d)*sc, "#fbeee6", INK, 1.2))
    b.append(txt(x0 + s*sc/2, y0 + r*sc + (s - d)*sc/2 + 5, "1", 13, INK,
                 "middle"))
    b.append(rect(x0 + s*sc, y0 + r*sc, (r - d)*sc, (s - d)*sc, "#f7f7f7",
                  INK, 1.2))
    b.append(txt(x0 + s*sc + (r - d)*sc/2, y0 + r*sc + (s - d)*sc/2 + 5,
                 "0", 13, INK, "middle"))
    b.append(txt(x0 - 8, y0 + r*sc/2, f"r = {r}", 11, MUT, "end"))
    b.append(txt(x0 - 8, y0 + r*sc + (s - d)*sc/2, f"s&#8722;d = {s-d}",
                 11, MUT, "end"))
    b.append(txt(x0 + s*sc/2, y0 + n*sc + 18, f"s = {s}", 11, MUT,
                 "middle"))
    b.append(txt(x0 + s*sc + (r - d)*sc/2, y0 + n*sc + 18,
                 f"r&#8722;d = {r-d}", 11, MUT, "middle"))
    b.append(txt(x0 + 44, y0 + 14, sub("R", "1", 10), 10.5, GOLD))
    b.append(txt(x0 + 44, y0 + 14 + t*sc, sub("R", "2", 10), 10.5, BLUE))
    b.append(txt(x0 + n*sc + 40, y0 + 40, "half the matrix", 11, MUT))
    b.append(txt(x0 + n*sc + 40, y0 + 56, "stays variable:", 11, MUT))
    b.append(txt(x0 + n*sc + 40, y0 + 72, "m &#8776; n&#178;/4", 11, MUT))
    svg("fig_params", (640, y0 + n*sc + 34), "\n".join(b),
        "Parameter split of the n by n permanent")


# --------------------------------------------------------------- fig 18
def fig_boundchart():
    b = []
    bx, by, bw, bh = 70, 40, 470, 190
    b += chart_axes(bx, by, bw, bh, "log&#8322; n  (log scale)",
                    "advantage factor log&#8322;log&#8322;n &#8722; 3")
    b.append(txt(bx + bw/2, 22, "Theorem 1.1's gain over n&#178;/144, "
                 "from threshold to astronomical n", 12.5, ACC, "middle",
                 weight="bold"))
    lo, hi = math.log10(16), math.log10(4096)

    def mx(l2n):
        return bx + (math.log10(l2n) - lo)/(hi - lo)*bw

    def my(v):
        return by + bh - v/10.0*bh
    pts = []
    l2 = 16.0
    while l2 <= 4096:
        pts.append((mx(l2), my(math.log2(l2) - 3)))
        l2 *= 1.05
    b.append(poly(pts, RED, 2.2))
    for l2 in [16, 32, 256, 4096]:
        fv = math.log2(l2) - 3
        b.append(circ(mx(l2), my(fv), 3.4, RED))
        b.append(txt(mx(l2) + 8, my(fv) + 13,
                     f"factor {fv:.0f} at n = 2"
                     f'<tspan baseline-shift="super" font-size="8">'
                     f"{int(l2)}</tspan>", 10, INK))
    b.append(txt(bx + 12, by + 16, "bound = (n&#178;/144)&#183;factor;  "
                 "it passes the trivial n&#178;&#8722;1 when the factor "
                 "reaches 144,", 10.2, MUT))
    b.append(txt(bx + 12, by + 31, "i.e. at n = 2"
                 '<tspan baseline-shift="super" font-size="8">2</tspan>'
                 '<tspan baseline-shift="super" font-size="6" dy="-4">'
                 "147</tspan>"
                 '<tspan dy="4"> &#8212; but it grows forever: '
                 "C(per&#8342;)/n&#178; &#8594; &#8734;</tspan>", 10.2,
                 MUT))
    for l2, lab in [(16, "16"), (64, "64"), (256, "256"),
                    (1024, "1024"), (4096, "4096")]:
        b.append(txt(mx(l2), by + bh + 14, lab, 10, MUT, "middle"))
    for v in [0, 2, 4, 6, 8, 10]:
        b.append(txt(bx - 6, my(v) + 4, str(v), 10, MUT, "end"))
    svg("fig_boundchart", (640, 288), "\n".join(b),
        "Computed circuit lower bound chart")


# --------------------------------------------------------------- fig 19
def fig_markedtree():
    b = []
    # tree: root at top; leaves at bottom. marked leaves: 3
    nodes = {
        "root": (300, 40, "m"), "g1": (180, 95, "m"), "g2": (420, 95, "u"),
        "g3": (120, 150, "m"), "g4": (250, 150, "u"),
        "g5": (370, 150, "u"), "g6": (470, 150, "u"),
        "g7": (70, 205, "m"), "g8": (170, 205, "m"),
        "l1": (40, 260, "Y"), "l2": (100, 260, "Z"),
        "l3": (140, 260, "Y"), "l4": (200, 260, "Z"),
        "l5": (225, 205, "Z"), "l6": (275, 205, "Y"),
        "l7": (345, 205, "Z"), "l8": (395, 205, "Z"),
        "l9": (445, 205, "Z"), "l10": (495, 205, "Z")}
    edges = [("root", "g1"), ("root", "g2"), ("g1", "g3"), ("g1", "g4"),
             ("g2", "g5"), ("g2", "g6"), ("g3", "g7"), ("g3", "g8"),
             ("g7", "l1"), ("g7", "l2"), ("g8", "l3"), ("g8", "l4"),
             ("g4", "l5"), ("g4", "l6"), ("g5", "l7"), ("g5", "l8"),
             ("g6", "l9"), ("g6", "l10")]
    marked_nodes = {"root", "g1", "g3", "g4", "g7", "g8", "l1", "l3",
                    "l6"}
    for (a, c) in edges:
        xa, ya = nodes[a][0], nodes[a][1]
        xc, yc = nodes[c][0], nodes[c][1]
        mk = a in marked_nodes and c in marked_nodes
        b.append(line(xa, ya, xc, yc, RED if mk else "#bbb",
                      2.2 if mk else 1.2))
    merge = {"g1", "g3"}
    for nm, (x, y, kind) in nodes.items():
        if nm.startswith("l"):
            isy = kind == "Y"
            b.append(rect(x - 11, y - 10, 22, 20,
                          RED if isy else "#fff",
                          RED if isy else "#999", 1.2, rx=4,
                          opac=0.25 if isy else None))
            b.append(txt(x, y + 4, "y" if isy else "z", 11.5,
                         RED if isy else MUT, "middle", cls="v"))
        elif nm in merge:
            b.append(rect(x - 9, y - 9, 18, 18, RED, RED, 1.4))
            b.append(txt(x, y + 4, "&#9632;", 1, RED, "middle"))
        else:
            mk = nm in marked_nodes
            b.append(circ(x, y, 9, "#fff", RED if mk else "#999",
                          1.8 if mk else 1.2))
    b.append(txt(300, 20, "t = 3 marked leaves &#8658; exactly t&#8722;1 "
                 "= 2 merge gates (&#9632;)", 12.5, ACC, "middle",
                 weight="bold"))
    b.append(txt(555, 100, "unmarked", 10.5, MUT, "middle"))
    b.append(txt(555, 114, "subtrees:", 10.5, MUT, "middle"))
    b.append(txt(555, 128, "free", 10.5, MUT, "middle"))
    b.append(txt(60, 165, "one-marked-child", 10, MUT))
    b.append(txt(60, 178, "chain: affine wrap", 10, MUT))
    svg("fig_markedtree", (620, 288), "\n".join(b),
        "Marked subtree of a formula")


# --------------------------------------------------------------- fig 20
def fig_affinepath():
    b = [arrow_defs("ap", [MUT])]
    x = 40
    steps = [("u", None), ("&#215; h&#8321;", "h&#8321;u"),
             ("+ h&#8322;", "h&#8321;u + h&#8322;"),
             ("h&#8323; &#8722; &#183;",
              "(h&#8323;&#8722;h&#8322;) &#8722; h&#8321;u")]
    for i, (op, res) in enumerate(steps):
        if i == 0:
            b.append(circ(x, 60, 14, "#fff", RED, 1.8))
            b.append(txt(x, 65, "u", 13, RED, "middle", cls="v"))
        else:
            b.append(rect(x - 26, 42, 52, 36, LT, BLUE, 1.3, rx=8))
            b.append(txt(x, 64, op, 12, BLUE, "middle"))
            b.append(txt(x, 100, res, 11, INK, "middle", cls="v"))
        if i < 3:
            b.append(line(x + 28, 60, x + 92, 60, MUT, 1.5,
                          marker="c5apa0"))
        x += 120
    b.append(rect(410, 34, 190, 60, "#fdf7e3", GOLD, 1.4, rx=8))
    b.append(txt(505, 56, "one affine map:", 11.5, INK, "middle"))
    b.append(txt(505, 76, "Au + B,  A = &#8722;h&#8321;, B = h&#8323;"
                 "&#8722;h&#8322;", 11.5, INK, "middle", cls="v"))
    b.append(txt(320, 130, "however long the chain, the marked value is "
                 "wrapped by just two elements of &#8450;[Z]", 11, MUT,
                 "middle"))
    body = ('<g transform="translate(0,-26)">' + "\n".join(b) + "</g>")
    svg("fig_affinepath", (640, 116), body,
        "Affine wrappers compose along a path")


# --------------------------------------------------------------- fig 21
def fig_fivebyfive():
    b = []
    labs = [["y&#8337;", "0", "p&#8321;", "p&#8322;", "p&#8323;"],
            ["0", "y&#8338;", "1", "1", "1"],
            ["1", "q&#8321;", "w&#8321;&#8321;", "w&#8321;&#8322;",
             "w&#8321;&#8323;"],
            ["1", "q&#8322;", "w&#8322;&#8321;", "w&#8322;&#8322;",
             "w&#8322;&#8323;"],
            ["1", "q&#8323;", "w&#8323;&#8321;", "w&#8323;&#8322;",
             "w&#8323;&#8323;"]]
    cell = 38
    x0, y0 = 150, 50
    for i in range(5):
        for j in range(5):
            if i < 2 and j < 2:
                fill = "#f3d9d5" if i == j else "#faf0ee"
            elif i < 2:
                fill = "#f9e8c8" if i == 0 else "#fdf6e8"
            elif j < 2:
                fill = "#d9e6f2" if j == 1 else "#ecf2f8"
            else:
                fill = LG
            b.append(rect(x0 + j*cell, y0 + i*cell, cell, cell, fill,
                          "#999", 0.8))
            b.append(txt(x0 + j*cell + cell/2, y0 + i*cell + cell/2 + 4,
                         labs[i][j], 12, INK, "middle", cls="v"))
    b.append(txt(x0 + cell, y0 - 20, "internal", 10.5, RED, "middle"))
    b.append(txt(x0 + 3.5*cell, y0 - 20, "external", 10.5, GREEN,
                 "middle"))
    b.append(txt(x0 - 12, y0 + cell + 4, "marked diagonal", 10, RED,
                 "end"))
    b.append(txt(x0 + 5*cell + 10, y0 + 0.5*cell + 4,
                 "p-signals exit row e", 10, GOLD))
    b.append(txt(x0 + 5*cell + 10, y0 + 1.5*cell + 4,
                 "ones exit row f", 10, MUT))
    b.append(txt(x0 - 12, y0 + 3.5*cell, "q-signals enter", 10, BLUE,
                 "end"))
    b.append(txt(x0 - 12, y0 + 3.5*cell + 13, "column f", 10, BLUE,
                 "end"))
    b.append(txt(x0 + 3.5*cell, y0 + 5*cell + 16,
                 "probe block W (&#8594; 1 after &#8706;/&#8706;w)", 10.5,
                 GREEN, "middle"))
    svg("fig_fivebyfive", (486, y0 + 5*cell + 32), "\n".join(b),
        "Anatomy of the 5x5 warm-up matrix")


# --------------------------------------------------------------- fig 22
def fig_blocks9():
    b = []
    ell, kk, mm = 3, 6, 7
    n = kk + mm
    cell = 22
    x0, y0 = 120, 60
    al, be = 5, 3   # sample: P(5)={e0,e2}, Q(3)={f0,f1}
    pset = [u for u in range(ell) if al >> u & 1]
    qset = [v for v in range(ell) if be >> v & 1]
    tset = pset + [ell + v for v in qset]
    for i in range(n):
        for j in range(n):
            if i < kk and j < kk:
                fill = "#f3d9d5" if i == j else "#faf0ee"
            elif i < kk:
                fill = "#f9e8c8" if i < ell else "#fbf3e2"
            elif j < kk:
                fill = "#d9e6f2" if j >= ell else "#ecf2f8"
            else:
                fill = LG
            b.append(rect(x0 + j*cell, y0 + i*cell, cell, cell, fill,
                          "#bbb", 0.5))
    for i in tset:
        b.append(rect(x0 + i*cell, y0 + i*cell, cell, cell, "none", RED,
                      2.2))
    for i in range(kk):
        if i not in tset:
            b.append(txt(x0 + i*cell + cell/2, y0 + i*cell + cell/2 + 4,
                         "S", 10, MUT, "middle"))
        else:
            b.append(txt(x0 + i*cell + cell/2, y0 + i*cell + cell/2 + 4,
                         "0", 10, RED, "middle"))
    for u in range(ell):
        b.append(txt(x0 - 8, y0 + u*cell + cell/2 + 4, sub("e", str(u),
                                                           10.5), 10.5,
                     GOLD, "end"))
        b.append(txt(x0 - 8, y0 + (ell + u)*cell + cell/2 + 4,
                     sub("f", str(u), 10.5), 10.5, BLUE, "end"))
    b.append(txt(x0 + kk*cell + mm*cell/2, y0 - 30,
                 "m = 7 external columns b", 10.5, GREEN, "middle"))
    b.append(txt(x0 + 1.5*cell, y0 - 30, "D: zeros +", 10, RED, "middle"))
    b.append(txt(x0 + 1.5*cell, y0 - 18, "marked diag", 10, RED,
                 "middle"))
    b.append(txt(x0 + kk*cell + mm*cell/2, y0 - 6,
                 "row e&#7524;: entries p&#7495;&#178;&#8202;"
                 "&#694;&#8202; &#183; row f&#7515;: ones", 9.5, GOLD,
                 "middle"))
    b.append(txt(x0 - 62, y0 + kk*cell + mm*cell/2, "V: col e ones,",
                 9.5, BLUE, "end"))
    b.append(txt(x0 - 62, y0 + kk*cell + mm*cell/2 + 13,
                 "col f&#7515;: q&#8342;&#178;&#8202;&#7515;", 9.5, BLUE,
                 "end"))
    b.append(txt(x0 + kk*cell + mm*cell/2,
                 y0 + kk*cell + mm*cell/2 + 4, "W &#8801; 1", 12, GREEN,
                 "middle"))
    b.append(txt(x0 + n*cell/2, y0 + n*cell + 20,
                 "sample (&#945;,&#946;) = (5,3): T = {e&#8320;,e&#8322;,"
                 "f&#8320;,f&#8321;} (red diagonal zeros stay in play; "
                 "S-cells are consumed)", 10, INK, "middle"))
    svg("fig_blocks9", (465, y0 + n*cell + 34), "\n".join(b),
        "Block anatomy of the Lemma 9.1 specialization")


# --------------------------------------------------------------- fig 23
def fig_kron():
    import itertools as it

    def fact(x):
        r = 1
        for i in range(2, x + 1):
            r *= i
        return r

    def falling(mv, sv):
        r = 1
        for i in range(sv):
            r *= mv - i
        return r

    def inj_sum(ws, nodes):
        tot = 0
        for pr in it.permutations(nodes, len(ws)):
            p = 1
            for i in range(len(ws)):
                p *= pr[i]**ws[i]
            tot += p
        return tot

    def gpoly(ws, nodes):
        if not ws:
            return [1]
        deg = sum(ws)
        co = [0]*(deg + 1)
        co[0] = inj_sum(ws, nodes)
        for i in range(len(ws)):
            rest = ws[:i] + ws[i+1:]
            sb = gpoly(rest, nodes)
            for j, cj in enumerate(sb):
                co[j + ws[i]] -= cj
        return co
    ell, mm = 3, 7
    pv = [2, 3, 5, 7, 11, 13, 17]
    qv = [19, 23, 29, 31, 37, 41, 43]
    J = [[0]*(mm*mm) for _ in range(mm*mm)]
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
                for bcol in range(mm):
                    gv = sum(c*pv[bcol]**j for j, c in enumerate(gp))
                    hv = sum(c*qv[a]**j for j, c in enumerate(hq))
                    J[al*mm + be][a*mm + bcol] = lpq*gv*hv
    mags = [[math.log10(abs(v)) if v else 0 for v in row] for row in J]
    mx = max(max(r) for r in mags)
    b = []
    cell = 7.2
    x0, y0 = 120, 52
    for i in range(49):
        for j in range(49):
            v = mags[i][j]/mx
            neg = J[i][j] < 0
            # warm scale for positive, blue for negative
            if neg:
                col = f"rgb({int(235-130*v)},{int(240-110*v)},245)"
            else:
                col = f"rgb(245,{int(235-150*v)},{int(225-190*v)})"
            b.append(f'<rect x="{x0+j*cell:.1f}" y="{y0+i*cell:.1f}" '
                     f'width="{cell:.1f}" height="{cell:.1f}" '
                     f'fill="{col}"/>')
    for g in range(8):
        b.append(line(x0 + g*7*cell, y0, x0 + g*7*cell, y0 + 49*cell,
                      "#fff", 1.4))
        b.append(line(x0, y0 + g*7*cell, x0 + 49*cell, y0 + g*7*cell,
                      "#fff", 1.4))
    b.append(rect(x0, y0, 49*cell, 49*cell, "none", INK, 1.2))
    # row scale strip: L_{P,Q}
    for al in range(mm):
        for be in range(mm):
            pw = bin(al).count("1")
            qw = bin(be).count("1")
            lpq = (fact(mm - 1 - pw - qw)*falling(mm - 1 - pw, qw)
                   * falling(mm - 1 - qw, pw))
            v = math.log10(lpq)/3.0
            col = f"rgb({int(240-140*v)},{int(230-60*v)},{int(200-60*v)})"
            b.append(f'<rect x="{x0-16}" y="{y0+(al*mm+be)*cell:.1f}" '
                     f'width="10" height="{cell:.1f}" fill="{col}"/>')
    b.append(txt(x0 - 22, y0 + 24.5*cell, sub("L", "P,Q", 10), 10.5, INK,
                 "end"))
    b.append(txt(x0 + 24.5*cell, y0 - 28, "the actual 49&#215;49 Jacobian"
                 " of the &#167;15 run (shaded by log-magnitude)",
                 11.5, ACC, "middle", weight="bold"))
    b.append(txt(x0 + 24.5*cell, y0 - 10, "columns (a,b) grouped by a "
                 "(7 blocks of 7); rows (&#945;,&#946;) grouped by "
                 "&#945;", 10.5, MUT, "middle"))
    b.append(txt(x0 + 24.5*cell, y0 + 49*cell + 18,
                 "coarse 7&#215;7 block pattern = A&#8318;&#7511; layer; "
                 "fine texture inside each block = A&#8346;&#7511; layer; "
                 "left strip = row scaling L", 10, MUT, "middle"))
    svg("fig_kron", (620, y0 + 49*cell + 34), "\n".join(b),
        "Computed Jacobian heat map showing Kronecker structure")


# --------------------------------------------------------------- fig 24
def fig_packing():
    b = []
    n, k = 12, 4
    cell = 21

    def hue(tau):
        h = tau*360.0/12
        return f"hsl({h:.0f},62%,62%)"
    x0, y0 = 60, 46
    b.append(txt(x0 + 2*cell, 20, "four sample matchings", 11.5, ACC,
                 "middle", weight="bold"))
    x1 = x0 + n*cell + 60
    b.append(txt(x1 + n*cell/2, 20, "all &#957; = 36 matchings: a "
                 "partition of the entries", 11.5, ACC, "middle",
                 weight="bold"))
    samples = [(0, 0), (3, 1), (7, 2), (10, 0)]
    for i in range(n):
        for j in range(n):
            b.append(rect(x0 + j*cell, y0 + i*cell, cell, cell, "#fff",
                          "#ccc", 0.6))
    for (tau, jb) in samples:
        for r in range(k):
            i = jb*k + r
            j = (i + tau) % n
            b.append(rect(x0 + j*cell + 1, y0 + i*cell + 1, cell - 2,
                          cell - 2, hue(tau), "none"))
    for (tau, jb) in samples:
        i = jb*k
        j = (i + tau) % n
        b.append(txt(x0 + j*cell + cell/2, y0 + i*cell - 3,
                     f"&#964;={tau}, j={jb}", 9, INK, "middle"))
    for i in range(n):
        for j in range(n):
            tau = (j - i) % n
            jb = i//k
            light = "62%" if jb == 1 else ("72%" if jb == 0 else "50%")
            h = tau*360.0/12
            b.append(rect(x1 + j*cell, y0 + i*cell, cell, cell,
                          f"hsl({h:.0f},58%,{light})", "#fff", 0.7))
    for jb in range(3):
        b.append(line(x1 - 4, y0 + jb*k*cell, x1 + n*cell + 4,
                      y0 + jb*k*cell, INK, 1.6))
        b.append(txt(x1 - 8, y0 + jb*k*cell + k*cell/2 + 4,
                     f"j={jb}", 10, INK, "end"))
    b.append(line(x1 - 4, y0 + n*cell, x1 + n*cell + 4, y0 + n*cell,
                  INK, 1.6))
    b.append(txt(x1 + n*cell/2, y0 + n*cell + 16, "hue = offset &#964;; "
                 "band = row block j; (&#964;,j) recoverable from any "
                 "cell", 9.5, MUT, "middle"))
    svg("fig_packing", (x1 + n*cell + 40, y0 + n*cell + 30), "\n".join(b),
        "Cyclic matchings tile the matrix")


# --------------------------------------------------------------- fig 25
def fig_flt():
    b = [arrow_defs("fl", [MUT])]
    ops = [("u + h", "1 h; 0 1"), ("u &#8722; h", "1 &#8722;h; 0 1"),
           ("h &#8722; u", "&#8722;1 h; 0 1"), ("h&#183;u", "h 0; 0 1"),
           ("u/h", "1 0; 0 h"), ("h/u", "0 h; 1 0")]
    b.append(txt(160, 20, "the six wrappers and their matrices", 12,
                 ACC, "middle", weight="bold"))
    for i, (op, mtx) in enumerate(ops):
        x = 40 + (i % 3)*105
        y = 42 + (i//3)*66
        b.append(rect(x, y, 92, 52, LT, BLUE, 1.2, rx=7))
        b.append(txt(x + 46, y + 20, op, 12, BLUE, "middle", cls="v"))
        a1, a2 = mtx.split("; ")
        b.append(txt(x + 46, y + 36, "(" + a1 + ")", 10, INK, "middle"))
        b.append(txt(x + 46, y + 47, "(" + a2 + ")", 10, INK, "middle"))
    b.append(txt(500, 20, "a valid path composes", 12, ACC, "middle",
                 weight="bold"))
    y = 48
    for i, lab in enumerate(["M&#8321;", "M&#8322;", "M&#8323;"]):
        b.append(rect(410 + i*62, y, 44, 30, "#fff", "#999", 1.1, rx=6))
        b.append(txt(432 + i*62, y + 19, lab, 12, INK, "middle"))
        if i < 2:
            b.append(line(456 + i*62, y + 15, 470 + i*62, y + 15, MUT,
                          1.4, marker="c5fla0"))
    b.append(txt(500, 104, "M = M&#8323;M&#8322;M&#8321; &#8800; 0:", 11.5,
                 INK, "middle"))
    b.append(txt(500, 121, "denominator factors as", 11, MUT, "middle"))
    b.append(txt(500, 137, "(c&#8342;u&#8320;+d&#8342;)(cu&#8342;+d)",
                 11, MUT, "middle", cls="v"))
    b.append(txt(500, 158, "normalize a nonzero entry to 1:", 11, MUT,
                 "middle"))
    b.append(txt(500, 174, "&#8804; 3 parameters per wrapper", 11.5, INK,
                 "middle"))
    svg("fig_flt", (620, 190), "\n".join(b),
        "Fractional linear wrappers under division")


# --------------------------------------------------------------- fig 26
def fig_schur():
    b = [arrow_defs("sc", [MUT])]
    b.append(txt(120, 22, "determinant: 2&#7511; coefficients", 12, ACC,
                 "middle", weight="bold"))
    for i in range(5):
        y = 40 + i*26
        b.append(rect(40, y, 120, 18, "#fff", "#999", 1, rx=4))
        lab = ["det &#931;", "det &#931;&#8321;&#8321;",
               "det &#931;&#8322;&#8322;",
               "det &#931;&#8321;&#8322;,&#8321;&#8322;", "&#8943;"][i]
        b.append(txt(100, y + 13, "&#948;&#183;" + lab if i < 4 else lab,
                     10.5, INK, "middle"))
    for i in range(5):
        b.append(line(162, 49 + i*26, 218, 100, MUT, 1.1,
                      marker="c5sca0"))
    b.append(rect(222, 72, 130, 56, "#fdf7e3", GOLD, 1.4, rx=8))
    b.append(txt(287, 94, "k&#178; entries of &#931;", 11.5, INK,
                 "middle"))
    b.append(txt(287, 112, "+ the scalar &#948;", 11.5, INK, "middle"))
    b.append(txt(287, 150, "td &#8804; k&#178; + 1", 13, RED, "middle",
                 weight="bold"))
    b.append(txt(287, 168, "k = &#920;(log n) &#8658; O(log&#178;n)", 11,
                 MUT, "middle"))
    # right: comparison bars
    b.append(txt(490, 22, "independent coefficients from a", 12, ACC,
                 "middle", weight="bold"))
    b.append(txt(490, 38, "k = 2&#8968;log&#8322;n&#8969; matching, "
                 "n = 1024", 11, ACC, "middle"))
    kv = 20
    detv = kv*kv + 1
    perv = (1024 - kv)**2
    bx = 420
    hdet = 12
    hper = 150.0
    b.append(rect(bx, 190 - hdet, 60, hdet, BLUE, "none"))
    b.append(txt(bx + 30, 205, "det:", 10.5, MUT, "middle"))
    b.append(txt(bx + 30, 218, f"{detv}", 10.5, BLUE, "middle"))
    b.append(rect(bx + 100, 190 - hper, 60, hper, RED, "none"))
    b.append(txt(bx + 130, 205, "per:", 10.5, MUT, "middle"))
    b.append(txt(bx + 130, 218, f"{perv:,}".replace(",", "&#8202;"),
                 10.5, RED, "middle"))
    b.append(txt(bx + 130, 190 - hper - 8, "m&#178;", 11, RED, "middle"))
    b.append(txt(490, 244, "(bars to scale would need a 2500&#215; "
                 "taller page)", 9.5, MUT, "middle"))
    svg("fig_schur", (640, 258), "\n".join(b),
        "Schur compression versus permanent independence")


# --------------------------------------------------------------- fig 27
def fig_detchart():
    b = []
    bx, by, bw, bh = 70, 40, 480, 190
    b += chart_axes(bx, by, bw, bh, "n  (log scale)",
                    "critical-locus codimension (log scale)")
    b.append(txt(bx + bw/2, 22, "why the geometric method is permanent-"
                 "specific", 12.5, ACC, "middle", weight="bold"))
    lo, hi = math.log10(2**16), math.log10(2**64)

    def mx(n):
        return bx + (math.log10(n) - lo)/(hi - lo)*bw

    def my(v):
        return by + bh - (math.log10(v) + 1.0)/41.0*bh
    pts_per, pts_det, pts_four, pts_full = [], [], [], []
    e = 16.0
    while e <= 64.0:
        n = 2**e
        d = max(3, int(e//4))
        t = 4*d
        bb = n/(2*t)
        r = bb*t
        s = n - r + d
        k = r*s/4
        pts_per.append((mx(n), my(k)))
        pts_det.append((mx(n), my(n - 1)))
        pts_four.append((mx(n), my(4)))
        pts_full.append((mx(n), my(2*n)))
        e += 1.0
    b.append(poly(pts_per, BLUE, 2.2))
    b.append(poly(pts_full, "#999", 1.4, dash="2 3"))
    b.append(poly(pts_det, RED, 2))
    b.append(poly(pts_four, RED, 1.3, dash="6 4"))
    # legend in the free upper-left corner
    lx, ly = bx + 14, by + 14
    entries = [
        (BLUE, None, 2.2, "permanent specialization: k = "
         "&#8970;rs/4&#8971; = &#920;(n&#178;)  (&#167;9)"),
        (RED, None, 2.0, "ceiling for EVERY det specialization: "
         "n&#8722;1  (Import 7)"),
        ("#999", "2 3", 1.4, "full per&#8342; itself: &#8804; 2n "
         "[BCMV25]"),
        (RED, "6 4", 1.3, "naive bordered det construction: 4")]
    for i, (col, dash, wd, lab) in enumerate(entries):
        yy = ly + i*18
        b.append(line(lx, yy, lx + 26, yy, col, wd, dash=dash))
        b.append(txt(lx + 33, yy + 4, lab, 10, INK))
    for e2 in [16, 32, 48, 64]:
        b.append(txt(mx(2**e2), by + bh + 14,
                     "2" + f'<tspan baseline-shift="super" '
                     f'font-size="7.5">{e2}</tspan>', 10, MUT, "middle"))
    for e2 in [0, 12, 24, 36]:
        b.append(txt(bx - 6, my(10.0**e2) + 4,
                     "10" + f'<tspan baseline-shift="super" '
                     f'font-size="7.5">{e2}</tspan>', 10, MUT, "end"))
    svg("fig_detchart", (640, 288), "\n".join(b),
        "Codimension comparison permanent versus determinant")


# --------------------------------------------------------------- fig 28
def fig_modelsmap():
    b = [arrow_defs("mm", [MUT])]

    def box(x, y, w, h, lines, col="#fff", stroke="#999", dash=None,
            tcol=INK):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        out = (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" '
               f'fill="{col}" stroke="{stroke}" stroke-width="1.3"{d}/>')
        for i, ln in enumerate(lines):
            out += txt(x + w/2, y + 17 + i*13.5, ln, 10.5, tcol, "middle")
        return out
    b.append(box(20, 40, 150, 50, ["formulas &#177; &#247;",
                                   "Thms 1.2&#8211;1.3:",
                                   "&#937;(n&#8308;/log n)"], BOX, ACC))
    b.append(box(210, 40, 170, 50, ["ABPs = weakly skew =",
                                    "det representations",
                                    "dc(per) &#8805; n&#178;/2 [MR04]"]))
    b.append(box(430, 40, 150, 50, ["general circuits",
                                    "Thm 1.1:",
                                    "&#937;(n&#178;loglog n)"], BOX, ACC))
    b.append(line(172, 65, 206, 65, MUT, 1.4, marker="c5mma0"))
    b.append(line(382, 65, 426, 65, MUT, 1.4, marker="c5mma0"))
    b.append(txt(400, 52, "quasi-poly", 8.5, MUT, "middle"))
    b.append(txt(190, 52, "poly", 9, MUT, "middle"))
    b.append(box(20, 130, 130, 46, ["monotone", "exponential [JS82]"],
                 "#fff", "#999", "5 4"))
    b.append(box(165, 130, 150, 46, ["synt. multilinear",
                                     'n<tspan baseline-shift="super" '
                                     'font-size="8">&#937;(log n)</tspan>'
                                     " [Raz09]"], "#fff", "#999", "5 4"))
    b.append(box(330, 130, 130, 46, ["homog. depth 3/4",
                                     "exp [NW97, GKKS14]"], "#fff",
                 "#999", "5 4"))
    b.append(box(475, 130, 140, 46, ["constant depth",
                                     "superpoly [LST25]"], "#fff", "#999",
                 "5 4"))
    b.append(txt(320, 200, "dashed = restricted models: stronger bounds, "
                 "no implication for the unrestricted boxes above", 10.5,
                 MUT, "middle"))
    b.append(txt(320, 216, "(a general circuit may be non-monotone, "
                 "non-multilinear, deep, and inhomogeneous at once)",
                 10.5, MUT, "middle"))
    svg("fig_modelsmap", (640, 230), "\n".join(b),
        "Model hierarchy and where the results sit")


# --------------------------------------------------------------- fig 29
def fig_checks():
    b = []
    cards = [
        ("C1", "per 4&#215;4 identity", "&#167;8.1", "symbolic, 24 terms"),
        ("C2", "det 4&#215;4 collapse", "&#167;13.2", "symbolic"),
        ("C3", "derivative id. (5.4)", "&#167;7.3",
         "t=4, s=5, d=3, exact"),
        ("C4", "symmetric ids. (7.c)", "&#167;7.5", "t=6, q=3, exact"),
        ("C5", "Lemma 6.1 end-to-end", "&#167;8",
         "&#8484;[&#950;], b=2, t=d=3, s=4"),
        ("C6", "5&#215;5 Jacobian", "&#167;11.1",
         "det J = 8&#183;36 = 288"),
        ("C7", "Lemma 9.1 end-to-end", "&#167;11",
         "n=13: 49&#178; entries + det"),
        ("C8", "packing + Thm 1.2", "&#167;12", "n = 32, &#957; = 96"),
        ("C9", "Thm 1.1 chain", "&#167;9", "n = 2&#8202;&#8321;&#8310;"),
        ("C10", "Schur identity", "&#167;13.1", "n=7, k=3, &#8474;")]
    for i, (cid, what, sec, par) in enumerate(cards):
        x = 24 + (i % 5)*122
        y = 20 + (i//5)*86
        b.append(rect(x, y, 110, 74, LG, GREEN, 1.3, rx=8))
        b.append(txt(x + 10, y + 18, cid, 12, GREEN, weight="bold"))
        b.append(txt(x + 100, y + 18, "PASS", 10, GREEN, "end",
                     weight="bold"))
        b.append(txt(x + 55, y + 36, what, 9.8, INK, "middle"))
        b.append(txt(x + 55, y + 51, sec, 9.5, ACC, "middle"))
        b.append(txt(x + 55, y + 65, par, 9, MUT, "middle"))
    svg("fig_checks", (640, 200), "\n".join(b),
        "Validation coverage map")


def main():
    fig_perdet()
    fig_models()
    nc = fig_results()
    print("formula crossing n =", nc)
    fig_roadmap()
    fig_agdim()
    fig_powerfiber()
    fig_gateeq()
    fig_slice()
    fig_project()
    fig_reverse()
    fig_minorsum()
    fig_partitions()
    fig_reduction()
    fig_fourbyfour()
    fig_blockB()
    fig_roots()
    fig_params()
    fig_boundchart()
    fig_markedtree()
    fig_affinepath()
    fig_fivebyfive()
    fig_blocks9()
    fig_kron()
    fig_packing()
    fig_flt()
    fig_schur()
    fig_detchart()
    fig_modelsmap()
    fig_checks()


if __name__ == "__main__":
    main()
