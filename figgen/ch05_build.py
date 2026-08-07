"""Assemble ch05-permanent-bounds.html from part files + generated SVGs.

Parts: .ignore/ch05_parts/p*.html (sorted). Markers <!--FIG:name--> are
replaced by .ignore/ch05_figs/name.svg content; the marker <!--LISTING-->
is replaced by the HTML-escaped contents of figgen/ch05_checks.py wrapped
in the codewrap/copy-button block. Idempotent; rerun any time.
Usage: PYTHONIOENCODING=utf-8 python figgen/ch05_build.py
"""
import glob
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARTS = os.path.join(ROOT, ".ignore", "ch05_parts")
FIGS = os.path.join(ROOT, ".ignore", "ch05_figs")
LISTING = os.path.join(ROOT, "figgen", "ch05_checks.py")
OUT = os.path.join(ROOT, "ch05-permanent-bounds.html")

BTN = ('<button class="copybtn" type="button" onclick="copyCode(this)" '
       'aria-label="Copy code to clipboard"><svg viewBox="0 0 24 24" '
       'fill="none" stroke="currentColor" stroke-width="2" '
       'stroke-linecap="round" stroke-linejoin="round">'
       '<rect x="9" y="9" width="13" height="13" rx="2"/>'
       '<path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"'
       '/></svg><span>Copy</span></button>')


def main():
    parts = sorted(glob.glob(os.path.join(PARTS, "p*.html")))
    if not parts:
        sys.exit("no part files found")
    body = "\n".join(io.open(p, encoding="utf-8").read() for p in parts)

    code = io.open(LISTING, encoding="utf-8").read()
    esc = (code.replace("&", "&amp;").replace("<", "&lt;")
           .replace(">", "&gt;"))
    block = ('<div class="codewrap">' + BTN + "<pre><code>" + esc
             + "</code></pre></div>")
    body = body.replace("<!--LISTING-->", block)

    missing = []

    def splice(m):
        name = m.group(1)
        path = os.path.join(FIGS, name + ".svg")
        if not os.path.exists(path):
            missing.append(name)
            return m.group(0)
        return io.open(path, encoding="utf-8").read()

    body = re.sub(r"<!--FIG:([A-Za-z0-9_]+)-->", splice, body)
    io.open(OUT, "w", encoding="utf-8").write(body)
    nfig = len(re.findall(r"<figure", body))
    print(f"wrote {OUT}: {len(body)} bytes, {len(parts)} parts, "
          f"{nfig} figures")
    if missing:
        print("MISSING FIGS:", missing)
        sys.exit(1)


if __name__ == "__main__":
    main()
