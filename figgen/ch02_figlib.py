"""Shared SVG helpers for ch02 figures. Marker ids prefixed c2_."""
import math
import os

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(
    __file__))), ".ignore", "ch02_figs")
os.makedirs(OUT, exist_ok=True)

RED = "#c0392b"
BLUE = "#2e6da4"
GREEN = "#2a7d2a"
GOLD = "#b8860b"
INK = "#333"
MUT = "#777"
LT = "#eef4fa"
LG = "#f1f7f1"
LY = "#fdf3d7"
ACC = "#7a1f1f"

_defs_used = set()


def svg(name, w, h, body, aria, vb=None):
    defs = (f'<defs>'
            f'<marker id="c2_ar_{name}" viewBox="0 0 10 10" refX="9" '
            f'refY="5" markerWidth="7" markerHeight="7" '
            f'orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{INK}"/></marker>'
            f'<marker id="c2_arR_{name}" viewBox="0 0 10 10" refX="9" '
            f'refY="5" markerWidth="7" markerHeight="7" '
            f'orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{RED}"/></marker>'
            f'<marker id="c2_arB_{name}" viewBox="0 0 10 10" refX="9" '
            f'refY="5" markerWidth="7" markerHeight="7" '
            f'orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{BLUE}"/></marker>'
            f'<marker id="c2_arG_{name}" viewBox="0 0 10 10" refX="9" '
            f'refY="5" markerWidth="7" markerHeight="7" '
            f'orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{GOLD}"/></marker>'
            f'</defs>')
    vbs = vb if vb else f"0 0 {w} {h}"
    s = (f'<svg class="setupfig" viewBox="{vbs}" width="100%" '
         f'role="img" aria-label="{aria}">\n{defs}\n{body}\n</svg>\n')
    with open(os.path.join(OUT, name + ".svg"), "w",
              encoding="utf-8") as f:
        f.write(s)
    print("wrote", name, len(s))


def ar(name, color=""):
    tag = {"": "c2_ar_", "R": "c2_arR_", "B": "c2_arB_",
           "G": "c2_arG_"}[color]
    return f'url(#{tag}{name})'


def txt(x, y, s, size=12, fill=INK, anchor="start", cls=None,
        weight=None):
    c = f' class="{cls}"' if cls else ""
    wgt = f' font-weight="{weight}"' if weight else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}"'
            f' fill="{fill}" text-anchor="{anchor}"{c}{wgt}>{s}</text>')


def sub(main, s, size=12):
    return (f'{main}<tspan baseline-shift="sub" '
            f'font-size="{size*0.72:.0f}">{s}</tspan>')


def sup(main, s, size=12):
    return (f'{main}<tspan baseline-shift="super" '
            f'font-size="{size*0.72:.0f}">{s}</tspan>')


def line(x1, y1, x2, y2, stroke=INK, w=1.4, dash=None, marker=None,
         opacity=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = f' marker-end="{marker}"' if marker else ""
    o = f' opacity="{opacity}"' if opacity is not None else ""
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" '
            f'y2="{y2:.1f}" stroke="{stroke}" stroke-width="{w}"{d}{m}{o}'
            f'/>')


def poly(pts, stroke=INK, w=1.6, fill="none", dash=None, opacity=None,
         close=False):
    p = " ".join(f"{x:.2f},{y:.2f}" for x, y in pts)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' opacity="{opacity}"' if opacity is not None else ""
    tag = "polygon" if close else "polyline"
    return (f'<{tag} points="{p}" fill="{fill}" stroke="{stroke}" '
            f'stroke-width="{w}"{d}{o}/>')


def circ(x, y, r, fill=INK, stroke="none", sw=1.0, opacity=None):
    o = f' opacity="{opacity}"' if opacity is not None else ""
    return (f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{sw}"{o}/>')


def rect(x, y, w, h, fill=LT, stroke=INK, sw=1.0, rx=4, opacity=None,
         dash=None):
    o = f' opacity="{opacity}"' if opacity is not None else ""
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" '
            f'height="{h:.1f}" rx="{rx}" fill="{fill}" stroke="{stroke}"'
            f' stroke-width="{sw}"{o}{d}/>')


def path(d, stroke=INK, w=1.6, fill="none", dash=None, opacity=None,
         marker=None):
    dd = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' opacity="{opacity}"' if opacity is not None else ""
    m = f' marker-end="{marker}"' if marker else ""
    return (f'<path d="{d}" fill="{fill}" stroke="{stroke}" '
            f'stroke-width="{w}"{dd}{o}{m}/>')


def axis_chart(x0, y0, w, h, xlab, ylab, name=None):
    """Simple L-shaped axes; returns svg string. y0 = bottom."""
    s = line(x0, y0, x0+w, y0, INK, 1.2)
    s += line(x0, y0, x0, y0-h, INK, 1.2)
    s += txt(x0+w/2, y0+28, xlab, 12, INK, "middle", cls="v")
    s += txt(x0-8, y0-h-8, ylab, 12, INK, "end", cls="v")
    return s


def chart_line(xs, ys, X, Y, color=RED, w=1.8, dash=None):
    pts = [(X(x), Y(y)) for x, y in zip(xs, ys)]
    return poly(pts, color, w, dash=dash)


def mapper(x0, x1, X0, X1):
    def f(v):
        return X0 + (v-x0)/(x1-x0)*(X1-X0)
    return f


def ortho(pt, az=0.6, el=0.35):
    """Orthographic projection of 3D point; returns (x, y, depth)."""
    ca, sa = math.cos(az), math.sin(az)
    ce, se = math.cos(el), math.sin(el)
    x, y, z = pt
    xr = ca*x + sa*y
    yr = -sa*x + ca*y
    X = xr
    Y = -z*ce + yr*se
    depth = yr*ce + z*se
    return X, Y, depth
