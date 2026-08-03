"""Make a checker-friendly copy of an explainer with the inline MathJax
bundle swapped back to a CDN tag (checkers would otherwise scan 2 MB of JS).

Usage: python figgen/strip_mathjax.py FILE.html OUT.html
"""
import io
import sys

MARK = '<script>/* MathJax 3.2.2'
CDN = ('<script src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/'
       'tex-mml-svg.js" integrity="sha384-msWocAZtTDh+149KSxjbRTGVfGbDji'
       'ff/j0JX1iu/A/GxxTWSu6ozNemlfkag/8d" crossorigin="anonymous"'
       ' async></script>')


def main():
    src, out = sys.argv[1], sys.argv[2]
    s = io.open(src, encoding='utf-8').read()
    i = s.find(MARK)
    if i < 0:
        io.open(out, 'w', encoding='utf-8').write(s)
        print('no inline bundle; copied as-is')
        return
    j = s.index('</script>', i) + len('</script>')
    s = s[:i] + CDN + s[j:]
    io.open(out, 'w', encoding='utf-8').write(s)
    print(f'stripped {j - i} bytes of inline MathJax -> {out}')


if __name__ == '__main__':
    main()
