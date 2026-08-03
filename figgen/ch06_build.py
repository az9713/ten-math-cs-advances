"""Assemble ch06-quantum-parallel-repetition.html from part files + generated SVGs.

Parts: .ignore/ch06_parts/p*.html (sorted). Markers <!--FIG:name--> are
replaced by .ignore/ch06_figs/name.svg content.
Usage: PYTHONIOENCODING=utf-8 python figgen/ch06_build.py
"""
import glob
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARTS = os.path.join(ROOT, ".ignore", "ch06_parts")
FIGS = os.path.join(ROOT, ".ignore", "ch06_figs")
OUT = os.path.join(ROOT, "ch06-quantum-parallel-repetition.html")


def main():
    parts = sorted(glob.glob(os.path.join(PARTS, "p*.html")))
    if not parts:
        sys.exit("no part files found")
    body = "\n".join(io.open(p, encoding="utf-8").read() for p in parts)

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
