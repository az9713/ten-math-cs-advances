"""Inject the HTML-escaped validation listing into p17.html.

Reads the tested listing from the scratchpad, escapes it, and replaces
the <!--LISTING--> placeholder (or a previously injected codewrap) so
the script is idempotent.
Usage: PYTHONIOENCODING=utf-8 python figgen/ch07_inject_listing.py <listing.py>
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PART = os.path.join(ROOT, ".ignore", "ch07_parts", "p17.html")

BTN = ('<div class="codewrap"><button class="copybtn" type="button" '
       'onclick="copyCode(this)" aria-label="Copy code to clipboard">'
       '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" '
       'stroke-width="2" stroke-linecap="round" stroke-linejoin="round">'
       '<rect x="9" y="9" width="13" height="13" rx="2"/>'
       '<path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/>'
       '</svg><span>Copy</span></button><pre><code>')


def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;")
            .replace("'", "&#x27;"))


def main():
    src = io.open(sys.argv[1], encoding="utf-8").read().rstrip("\n")
    block = BTN + esc(src) + "</code></pre></div>"
    html = io.open(PART, encoding="utf-8").read()
    new, n = re.subn(r"<!--LISTING-->", lambda m: block, html)
    if n == 0:
        new, n = re.subn(
            r'<div class="codewrap">.*?</code></pre></div>',
            lambda m: block, html, flags=re.S)
    if n != 1:
        sys.exit(f"expected 1 injection point, found {n}")
    io.open(PART, "w", encoding="utf-8").write(new)
    print(f"injected {len(src)} chars of listing into {PART}")


if __name__ == "__main__":
    main()
