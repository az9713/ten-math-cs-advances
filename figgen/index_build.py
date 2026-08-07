"""Build index.html: generate the two overview figures and splice them
into .ignore/index_parts/body.html at the <!--FIG:name--> markers.

Figures:
  atlas     -- the ten chapters on the PDF page axis, colored by
               outsider difficulty, annotated with pp + figure counts.
  machinery -- the shared-machinery map (primer arrows + the ch2->ch1
               s->1 back-edge).

Idempotent: rerun any time.  Usage:  python figgen/index_build.py
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARTS = os.path.join(ROOT, ".ignore", "index_parts", "body.html")
FIGD = os.path.join(ROOT, ".ignore", "index_figs")
OUT = os.path.join(ROOT, "index.html")

# ch, first page, last page, figs, difficulty key, short label
CH = [
    (1, 3, 28, 32, "med", "Sphere packing"),
    (2, 29, 78, 34, "high", "Binary + spherical codes"),
    (3, 79, 95, 24, "high", "Non-sofic group"),
    (4, 96, 113, 22, "highest", "Connes rigidity"),
    (5, 114, 153, 29, "med", "Permanent bounds"),
    (6, 154, 182, 29, "high", "Parallel repetition"),
    (7, 183, 218, 38, "medhigh", "CVP hardness"),
    (8, 219, 228, 13, "high", "Ehrhart volume"),
    (9, 229, 235, 15, "low", "Multicolor Ramsey"),
    (10, 236, 249, 20, "lowmed", "Compactness + degeneracy"),
]
DIFF = {  # key -> (fill color, legend label)
    "low": ("#2e7d32", "low"),
    "lowmed": ("#689023", "low–med"),
    "med": ("#b8860b", "medium"),
    "medhigh": ("#c05c20", "med–high"),
    "high": ("#a83232", "high"),
    "highest": ("#6d2a8a", "highest"),
}


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def atlas_svg():
    x0, x1 = 212.0, 648.0
    p0, p1 = 3.0, 249.0

    def X(p):
        return x0 + (p - p0) * (x1 - x0) / (p1 - p0)

    W, H = 740, 402
    ytop, rowh, barh = 52, 28, 16
    yax = ytop + 10 * rowh + 4
    s = []
    s.append('<svg viewBox="0 0 %d %d" width="%d" role="img" '
             'aria-label="The ten chapters on the PDF page axis">'
             % (W, H, W))
    # legend
    lx = 30.0
    s.append('<text x="%.1f" y="27" font-size="11.5" fill="#333" '
             'font-weight="bold">outsider difficulty:</text>' % lx)
    lx += 118
    for key in ["low", "lowmed", "med", "medhigh", "high", "highest"]:
        col, lab = DIFF[key]
        s.append('<rect x="%.1f" y="18" width="12" height="12" rx="2" '
                 'fill="%s"/>' % (lx, col))
        s.append('<text x="%.1f" y="28" font-size="11.5" fill="#333">%s'
                 '</text>' % (lx + 16, esc(lab)))
        lx += 16 + 6.2 * len(lab) + 18
    # gridlines + axis ticks
    for p in [3, 50, 100, 150, 200, 249]:
        x = X(p)
        s.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" '
                 'stroke="#e3ddd2" stroke-width="1"/>' % (x, ytop - 8, x, yax))
        s.append('<text x="%.1f" y="%d" font-size="11" fill="#666" '
                 'text-anchor="middle">%d</text>' % (x, yax + 16, p))
    s.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#333" '
             'stroke-width="1.2"/>' % (x0, yax, x1, yax))
    s.append('<text x="%.1f" y="%d" font-size="11.5" fill="#333" '
             'text-anchor="middle">PDF page</text>'
             % ((x0 + x1) / 2, yax + 34))
    # rows
    for i, (ch, a, b, nf, dk, lab) in enumerate(CH):
        y = ytop + i * rowh
        col = DIFF[dk][0]
        s.append('<text x="204" y="%.1f" font-size="12" fill="#333" '
                 'text-anchor="end"><tspan font-weight="bold">%d</tspan>'
                 '<tspan fill="#666" font-size="11"> · %s</tspan>'
                 '</text>' % (y + barh - 3.5, ch, esc(lab)))
        s.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%d" rx="2.5"'
                 ' fill="%s"/>' % (X(a), y, X(b) - X(a), barh, col))
        s.append('<text x="656" y="%.1f" font-size="11" fill="#666">'
                 '%d pp · %d figs</text>'
                 % (y + barh - 3.5, b - a + 1, nf))
    s.append('</svg>')
    return "\n".join(s)


BOXW, BOXH = 210, 54


def box(s, x, y, title, sub, primer):
    stroke = "#7a1f1f" if primer else "#8a8375"
    sw = 1.8 if primer else 1.1
    s.append('<rect x="%d" y="%d" width="%d" height="%d" rx="7" '
             'fill="#fffdf8" stroke="%s" stroke-width="%.1f"/>'
             % (x, y, BOXW, BOXH, stroke, sw))
    cx = x + BOXW / 2
    s.append('<text x="%.1f" y="%d" font-size="13" font-weight="bold" '
             'fill="#333" text-anchor="middle">%s</text>'
             % (cx, y + 22, esc(title)))
    s.append('<text x="%.1f" y="%d" font-size="10.5" fill="#666" '
             'text-anchor="middle">%s</text>' % (cx, y + 40, esc(sub)))


def harrow(s, x1, y, x2, label):
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#7a1f1f" '
             'stroke-width="1.6" marker-end="url(#ixah)"/>'
             % (x1, y, x2, y))
    s.append('<text x="%.1f" y="%d" font-size="11" fill="#7a1f1f" '
             'text-anchor="middle">%s</text>'
             % ((x1 + x2) / 2, y - 9, esc(label)))


def machinery_svg():
    W, H = 770, 448
    L, R = 45, 505
    y1, y2, y3 = 52, 162, 272
    s = []
    s.append('<svg viewBox="0 0 %d %d" width="%d" role="img" '
             'aria-label="Shared-machinery map between the chapters">'
             % (W, H, W))
    s.append('<defs>'
             '<marker id="ixah" viewBox="0 0 10 10" refX="9" refY="5" '
             'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
             '<path d="M0,0 L10,5 L0,10 z" fill="#7a1f1f"/></marker>'
             '<marker id="ixahd" viewBox="0 0 10 10" refX="9" refY="5" '
             'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
             '<path d="M0,0 L10,5 L0,10 z" fill="#34687f"/></marker>'
             '</defs>')
    # dashed back-edge ch2 -> ch1 over the top
    s.append('<path d="M %d %d C %d 4, %d 4, %d %d" fill="none" '
             'stroke="#34687f" stroke-width="1.5" stroke-dasharray="6 4" '
             'marker-end="url(#ixahd)"/>'
             % (R + BOXW // 2, y1 - 2, R + BOXW // 2, L + BOXW // 2,
                L + BOXW // 2, y1 - 2))
    s.append('<text x="380" y="34" font-size="11" fill="#34687f" '
             'text-anchor="middle">ch. 2 §15, s → 1 limit: '
             'recovers the ch. 1 packing exponent</text>')
    # boxes
    box(s, L, y1, "Ch. 1 — Sphere packing", "home of the LP-certificate spine", True)
    box(s, R, y1, "Ch. 2 — Codes", "Delsarte LP, second pass", True)
    box(s, L, y2, "Ch. 3 — Non-sofic group", "property (T) + expanders in §3", True)
    box(s, R, y2, "Ch. 4 — Connes rigidity", "reuses the §3 primer", False)
    box(s, L, y3, "Ch. 8 — Ehrhart volume", "lattice primer in §2", True)
    box(s, R, y3, "Ch. 7 — CVP hardness", "uses lattices + codes", False)
    # primer arrows
    harrow(s, L + BOXW, y1 + 27, R - 5, "LP certificates, reused")
    harrow(s, L + BOXW, y2 + 27, R - 5, "property (T) / expander primer")
    harrow(s, L + BOXW, y3 + 27, R - 5, "Blichfeldt + Minkowski, proved once")
    # codes edge ch2 -> ch7 around the right side
    s.append('<path d="M %d %d C 757 130, 757 240, %d %d" fill="none" '
             'stroke="#7a1f1f" stroke-width="1.6" marker-end="url(#ixah)"/>'
             % (R + BOXW, y1 + 40, R + BOXW + 5, y3 + 14))
    s.append('<text x="744" y="187" font-size="11" fill="#7a1f1f" '
             'text-anchor="middle" transform="rotate(-90 744 187)">'
             'codes primer (§§1–2)</text>')
    # self-contained strip
    s.append('<text x="380" y="372" font-size="11.5" fill="#666" '
             'text-anchor="middle" font-style="italic">self-contained '
             '— no shared primer:</text>')
    small = [("Ch. 5", "permanent bounds"), ("Ch. 6", "parallel repetition"),
             ("Ch. 9", "multicolor Ramsey"), ("Ch. 10", "compactness")]
    sw, gap = 140, 22
    x = (W - (len(small) * sw + (len(small) - 1) * gap)) / 2
    for t, sub in small:
        s.append('<rect x="%.1f" y="384" width="%d" height="40" rx="6" '
                 'fill="#fffdf8" stroke="#8a8375" stroke-width="1.1"/>'
                 % (x, sw))
        s.append('<text x="%.1f" y="401" font-size="12" font-weight="bold" '
                 'fill="#333" text-anchor="middle">%s</text>'
                 % (x + sw / 2, esc(t)))
        s.append('<text x="%.1f" y="416" font-size="10.5" fill="#666" '
                 'text-anchor="middle">%s</text>'
                 % (x + sw / 2, esc(sub)))
        x += sw + gap
    s.append('</svg>')
    return "\n".join(s)


def main():
    figs = {
        "atlas": ('<figure id="fig-atlas">\n%s\n<figcaption>The survey on '
                  'its own page axis: one bar per chapter, from first to '
                  'last PDF page, colored by outsider difficulty and '
                  'annotated with chapter length and the explainer’s '
                  'figure count. Pages 3–249; 256 figures in total.'
                  '</figcaption>\n</figure>') % atlas_svg(),
        "machinery": ('<figure id="fig-machinery">\n%s\n<figcaption>The '
                      'shared-machinery map. Boxes with a dark-red border '
                      'are primer homes; solid arrows are primers reused '
                      'forward (read the source section first, or accept '
                      'it as a black box); the dashed arrow is the one '
                      'backward flow of a <i>result</i>: chapter 2’s '
                      's → 1 limit independently recovers chapter '
                      '1’s packing exponent. The bottom row imports '
                      'nothing.</figcaption>\n</figure>') % machinery_svg(),
    }
    with open(PARTS, encoding="utf-8") as f:
        html = f.read()
    for name, block in figs.items():
        marker = "<!--FIG:%s-->" % name
        assert marker in html, "missing marker " + marker
        html = html.replace(marker, block)
    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write(html)
    os.makedirs(FIGD, exist_ok=True)
    for name in figs:
        fn = os.path.join(FIGD, "fig_%s.svg" % name)
        with open(fn, "w", encoding="utf-8", newline="\n") as f:
            f.write(figs[name])
    print("wrote", OUT, "(%d bytes)" % len(html.encode("utf-8")))


if __name__ == "__main__":
    main()
