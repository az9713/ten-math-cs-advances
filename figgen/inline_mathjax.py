"""Embed MathJax 3.2.2 inline so explainers render with zero network.

Usage: PYTHONIOENCODING=utf-8 python figgen/inline_mathjax.py FILE.html...
Replaces the CDN <script src=...tex-mml-svg.js ...> tag (with or without
SRI lines) by an inline <script> containing the full bundle. Idempotent:
files without a CDN tag are left untouched.
"""
import io
import re
import sys

BUNDLE = '.ignore/mathjax-3.2.2-tex-mml-svg.js'
CDN_RE = re.compile(
    r'<script src="https://cdn\.jsdelivr\.net/npm/mathjax@[^"]*"'
    r'[^>]*(?:\n[^>]*)*></script>')


def main():
    js = io.open(BUNDLE, encoding='utf-8').read()
    assert '</script' not in js, 'bundle would break inline embedding'
    inline = '<script>/* MathJax 3.2.2 tex-mml-svg, embedded for offline' \
             ' rendering */\n' + js + '\n</script>'
    for path in sys.argv[1:]:
        s = io.open(path, encoding='utf-8').read()
        s2, nsub = CDN_RE.subn(lambda _m: inline, s, count=1)
        if nsub:
            io.open(path, 'w', encoding='utf-8').write(s2)
        print(f'{path}: {"inlined" if nsub else "no CDN tag found"}')


if __name__ == '__main__':
    main()
