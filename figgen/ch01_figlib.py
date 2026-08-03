"""SVG plot helpers for ch01 figures (house style, light card bg).

Coordinate mapping: data (x,y) -> pixels inside a margin box.
All ids/markers must be prefixed 'c1' + figure key by the caller.
"""

SERIF = 'font-family="Georgia,serif"'  # CSS on page overrides to STIX


def fmt(v):
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return s if s != "-0" else "0"


class Ax:
    """Axes mapping data coords to pixel coords."""

    def __init__(self, x0, x1, y0, y1, W=460, H=230,
                 ml=46, mr=14, mt=12, mb=34):
        self.x0, self.x1, self.y0, self.y1 = x0, x1, y0, y1
        self.W, self.H = W, H
        self.ml, self.mr, self.mt, self.mb = ml, mr, mt, mb
        self.iw = W - ml - mr
        self.ih = H - mt - mb

    def X(self, x):
        return self.ml + (x - self.x0) / (self.x1 - self.x0) * self.iw

    def Y(self, y):
        return self.mt + (self.y1 - y) / (self.y1 - self.y0) * self.ih

    def P(self, x, y):
        return f"{fmt(self.X(x))},{fmt(self.Y(y))}"

    def polyline(self, pts, color, w=2.0, dash=None, opacity=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        o = f' opacity="{opacity}"' if opacity else ""
        p = " ".join(self.P(x, y) for x, y in pts)
        return (f'<polyline points="{p}" fill="none" stroke="{color}" '
                f'stroke-width="{w}"{d}{o}/>')

    def path_area(self, pts, color, opacity=0.18, ybase=0.0):
        first, last = pts[0], pts[-1]
        p = " ".join("L" + self.P(x, y) for x, y in pts)
        return (f'<path d="M{self.P(first[0], ybase)} {p} '
                f'L{self.P(last[0], ybase)} Z" fill="{color}" '
                f'opacity="{opacity}" stroke="none"/>')

    def hline(self, y, color="#999", w=1.0, dash="4 3", x0=None, x1=None):
        x0 = self.x0 if x0 is None else x0
        x1 = self.x1 if x1 is None else x1
        return (f'<line x1="{fmt(self.X(x0))}" y1="{fmt(self.Y(y))}" '
                f'x2="{fmt(self.X(x1))}" y2="{fmt(self.Y(y))}" '
                f'stroke="{color}" stroke-width="{w}" '
                f'stroke-dasharray="{dash}"/>')

    def vline(self, x, color="#999", w=1.0, dash="4 3", y0=None, y1=None):
        y0 = self.y0 if y0 is None else y0
        y1 = self.y1 if y1 is None else y1
        return (f'<line x1="{fmt(self.X(x))}" y1="{fmt(self.Y(y0))}" '
                f'x2="{fmt(self.X(x))}" y2="{fmt(self.Y(y1))}" '
                f'stroke="{color}" stroke-width="{w}" '
                f'stroke-dasharray="{dash}"/>')

    def axes(self, xticks=(), yticks=(), xlab="", ylab="", fs=11,
             xtickfmt=fmt, ytickfmt=fmt, y_axis_at=None, x_axis_at=None):
        """Draw axis lines, ticks, labels. Labels are plain text
        (caller supplies unicode; no _ or ^)."""
        xa = self.y0 if x_axis_at is None else x_axis_at
        ya = self.x0 if y_axis_at is None else y_axis_at
        out = []
        out.append(f'<line x1="{fmt(self.X(self.x0))}" '
                   f'y1="{fmt(self.Y(xa))}" x2="{fmt(self.X(self.x1))}" '
                   f'y2="{fmt(self.Y(xa))}" stroke="#333" '
                   f'stroke-width="1.2"/>')
        out.append(f'<line x1="{fmt(self.X(ya))}" '
                   f'y1="{fmt(self.Y(self.y0))}" x2="{fmt(self.X(ya))}" '
                   f'y2="{fmt(self.Y(self.y1))}" stroke="#333" '
                   f'stroke-width="1.2"/>')
        for t in xticks:
            xp = self.X(t)
            yp = self.Y(xa)
            out.append(f'<line x1="{fmt(xp)}" y1="{fmt(yp)}" '
                       f'x2="{fmt(xp)}" y2="{fmt(yp + 4)}" '
                       f'stroke="#333" stroke-width="1"/>')
            out.append(f'<text x="{fmt(xp)}" y="{fmt(yp + 16)}" '
                       f'font-size="{fs}" text-anchor="middle" '
                       f'fill="#333">{xtickfmt(t)}</text>')
        for t in yticks:
            xp = self.X(ya)
            yp = self.Y(t)
            out.append(f'<line x1="{fmt(xp - 4)}" y1="{fmt(yp)}" '
                       f'x2="{fmt(xp)}" y2="{fmt(yp)}" '
                       f'stroke="#333" stroke-width="1"/>')
            out.append(f'<text x="{fmt(xp - 7)}" y="{fmt(yp + 4)}" '
                       f'font-size="{fs}" text-anchor="end" '
                       f'fill="#333">{ytickfmt(t)}</text>')
        if xlab:
            out.append(f'<text x="{fmt(self.X(self.x1))}" '
                       f'y="{fmt(self.Y(xa) + 28)}" font-size="{fs + 1}" '
                       f'text-anchor="end" fill="#333">{xlab}</text>')
        if ylab:
            out.append(f'<text x="{fmt(self.X(ya) - 4)}" '
                       f'y="{fmt(self.Y(self.y1) + 7)}" '
                       f'font-size="{fs + 1}" text-anchor="start" '
                       f'fill="#333">{ylab}</text>')
        return "\n".join(out)

    def text(self, x, y, s, fs=12, anchor="start", color="#333", cls=None,
             dx=0, dy=0):
        c = f' class="{cls}"' if cls else ""
        return (f'<text x="{fmt(self.X(x) + dx)}" y="{fmt(self.Y(y) + dy)}"'
                f' font-size="{fs}" text-anchor="{anchor}" '
                f'fill="{color}"{c}>{s}</text>')

    def dot(self, x, y, r=3.5, color="#333"):
        return (f'<circle cx="{fmt(self.X(x))}" cy="{fmt(self.Y(y))}" '
                f'r="{r}" fill="{color}"/>')


def svg_wrap(body, W, H, label):
    return (f'<svg class="setupfig" viewBox="0 0 {W} {H}" width="100%" '
            f'role="img" aria-label="{label}">\n{body}\n</svg>')


def sample(f, x0, x1, n=240):
    pts = []
    for i in range(n + 1):
        x = x0 + (x1 - x0) * i / n
        pts.append((x, f(x)))
    return pts


def write_fig(outdir, name, svg):
    import os
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, name + ".svg"), "w",
              encoding="utf-8") as fh:
        fh.write(svg)
