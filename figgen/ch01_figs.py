"""Chapter 1 (Cohn-Elkies sphere packing LP rate): figure data + checks.

Emits SVG figures to .ignore/ch01_figs/ and prints a validation battery.
Run:  PYTHONIOENCODING=utf-8 python figgen/ch01_figs.py
"""
import math
import cmath
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch01_figlib import Ax, svg_wrap, sample, write_fig, fmt  # noqa: E402

import mpmath as mp  # noqa: E402

mp.mp.dps = 30

OUT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), ".ignore", "ch01_figs")

RED = "#c0392b"
BLUE = "#2e6da4"
GREEN = "#2a7d2a"
GRAY = "#888"
DARK = "#333"
MUT = "#555"
PUR = "#7d3c98"
ORA = "#b06a00"

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(f"CHECK {name}: {'OK' if ok else 'FAIL'}  {detail}")


def simpson(f, a, b, n=400):
    if n % 2:
        n += 1
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += f(a + i * h) * (4 if i % 2 else 2)
    return s * h / 3


# ----------------------------------------------------------------- constants
def constants():
    disp = mp.quad(lambda a: -mp.e**(-2 * a) / (2 * a**2 * mp.cosh(a))
                   * a * mp.sinh(a), [0, mp.inf])
    check("c01 displacement = -log(pi/2)/2",
          abs(disp + mp.log(mp.pi / 2) / 2) < 1e-20, f"{disp}")
    # Wallis partial products of (2k/(2k-1))(2k/(2k+1))
    P = mp.mpf(1)
    for k in range(1, 20001):
        P *= mp.mpf(4 * k * k) / ((2 * k - 1) * (2 * k + 1))
    check("c02 Wallis -> pi/2", abs(P - mp.pi / 2) < 1e-4, f"{P}")
    p = lambda u: mp.pi / 4 / mp.cosh(mp.pi * u / 2)**2  # noqa: E731
    for t in (0.6, 1.3):
        cf = mp.quad(lambda u: p(u) * mp.cos(t * u), [-mp.inf, mp.inf])
        check(f"c03 charfn t={t}", abs(cf - t / mp.sinh(t)) < 1e-18)
    for x in (0.0, 0.7, 2.0):
        lhs = mp.quad(lambda u: p(u) * mp.log(mp.sqrt(x * x + u * u)),
                      [-mp.inf, -1, 0, 1, mp.inf])
        rhs = mp.digamma((x + 1) / 2) + mp.log(2)
        check(f"c04 identity(22) x={x}", abs(lhs - rhs) < 1e-10)
    ipsi = mp.quad(lambda x: mp.digamma((x + 1) / 2), [0, 1])
    check("c05 int psi = -log pi", abs(ipsi + mp.log(mp.pi)) < 1e-18)
    for s in (0.9, 0.99):
        th = mp.pi * (1 + s) / 2
        Ps = lambda T: mp.sin(th) / (4 * (mp.cosh(mp.pi * T / 2)  # noqa
                                          - mp.cos(th)))
        Ms = (1 - s) / 2
        inner = lambda T: mp.quad(  # noqa: E731
            lambda x: mp.log(mp.sqrt(x * x + T * T / 4)), [0, 1])
        J = -mp.quad(lambda T: Ps(T) * inner(T),
                     [-mp.inf, 0, mp.inf]) / Ms
        check(f"c06 J(sigma={s})", abs(J - mp.log(mp.pi / 2)) < 0.02,
              f"J={J}")
    s = 0.3
    th = mp.pi * (1 + s) / 2
    M = mp.quad(lambda T: mp.sin(th) / (4 * (mp.cosh(mp.pi * T / 2)
                                             - mp.cos(th))),
                [-mp.inf, mp.inf])
    check("c07 mass M_sigma=(1-sigma)/2", abs(M - (1 - s) / 2) < 1e-18)
    a_star = mp.log(2 * mp.pi / mp.e) / (2 * mp.log(2))
    check("c23 2^(-alpha*) = sqrt(e/2pi)",
          abs(2**(-a_star) - mp.sqrt(mp.e / (2 * mp.pi))) < 1e-25,
          f"alpha*={a_star}")
    check("c24 1/8+1-log8 < 0", mp.mpf(1) / 8 + 1 - mp.log(8) < 0,
          f"{mp.mpf(1)/8+1-mp.log(8)}")
    v2 = mp.pi**1 / mp.gamma(2)
    v3 = mp.pi**1.5 / mp.gamma(2.5)
    check("c22 v2=pi, v3=4pi/3", abs(v2 - mp.pi) < 1e-20
          and abs(v3 - 4 * mp.pi / 3) < 1e-20)
    d = 4000
    vd = mp.pi**(d / 2) / mp.gamma(d / 2 + 1)
    check("c22b vd^(1/d)*sqrt(d)->sqrt(2 pi e)",
          abs(vd**(mp.mpf(1) / d) * mp.sqrt(d)
              - mp.sqrt(2 * mp.pi * mp.e)) < 0.01)
    # Poisson summation on Z for Gaussian e^{-pi n^2 s}
    ssum = lambda s: mp.nsum(lambda n: mp.e**(-mp.pi * n * n * s),  # noqa
                             [-mp.inf, mp.inf])
    check("c26 Poisson Z Gaussian s=0.7",
          abs(ssum(mp.mpf("0.7"))
              - ssum(1 / mp.mpf("0.7")) / mp.sqrt(mp.mpf("0.7"))) < 1e-25)
    B, T = 6.0, 2.2
    lhs = mp.quad(lambda a: 1 - mp.cos(a * T), [B, B + 1])
    rhs = 1 - mp.sin(T / 2) / (T / 2) * mp.cos((B + 0.5) * T)
    check("c28 resonance identity (58)", abs(lhs - rhs) < 1e-12)


# ------------------------------------------------- d=3 Hankel + eigenchecks
def d3_fourier(g, rho):
    """Radial Fourier transform in d=3: (2/rho) int g(r) r sin(2 pi r rho)."""
    pieces = [0, 1, 2, 3, 4, 6, 8]
    if rho == 0:
        return 4 * mp.pi * mp.quad(lambda r: g(r) * r * r, pieces)
    return 2 / rho * mp.quad(lambda r: g(r) * r * mp.sin(2 * mp.pi * r * rho),
                             pieces)


def mellin_checks():
    g = lambda r: r * r * mp.e**(-mp.pi * r * r) / (1 + r * r)  # noqa: E731
    d = 3
    lam = mp.mpf(d) / 2
    for z in (mp.mpf("1.3"), mp.mpf(2) + 0.7j):
        Mg = mp.quad(lambda r: g(r) * r**(d - z - 1), [0, 1, 2, 4, 8])
        Mhat = mp.quad(lambda rho: d3_fourier(g, rho) * rho**(z - 1),
                       [0, 1, 2, 3, 4, 6])
        rhs = mp.pi**(lam - z) * mp.gamma(z / 2) \
            / mp.gamma((d - z) / 2) * Mg
        rel = abs(Mhat - rhs) / abs(rhs)
        check(f"c08 Mellin-Hankel z={z}", rel < 1e-8, f"rel={rel}")
    # psi+ self-dual, psi- anti-self-dual in d=3
    psi_p = lambda r: mp.e**(-mp.pi * r * r)  # noqa: E731
    psi_m = lambda r: (r * r - mp.mpf(d) / (4 * mp.pi)) \
        * mp.e**(-mp.pi * r * r)  # noqa: E731
    errp = max(abs(d3_fourier(psi_p, rho) - psi_p(rho))
               for rho in (0.3, 0.9, 1.7))
    errm = max(abs(d3_fourier(psi_m, rho) + psi_m(rho))
               for rho in (0.3, 0.9, 1.7))
    check("c09 psi+ self-Fourier (d=3)", errp < 1e-12, f"{errp}")
    check("c09b psi- anti-self-Fourier (d=3)", errm < 1e-12, f"{errm}")
    # Laguerre eigenfunction k=2, alpha=d/2-1: hat = (+1)*phi2
    alpha = mp.mpf(d) / 2 - 1
    phi2 = lambda r: mp.laguerre(2, alpha, 2 * mp.pi * r * r) \
        * mp.e**(-mp.pi * r * r)  # noqa: E731
    err2 = max(abs(d3_fourier(phi2, rho) - phi2(rho))
               for rho in (0.3, 0.9, 1.7))
    check("c10 Laguerre k=2 eigen (+1)", err2 < 1e-15, f"{err2}")


# ------------------------------------------------------- mass concentration
def mass_fig():
    """g = phi2/phi2(0) - phi0, eigenvalue +1, g(0)=0. Radial mass density
    r^{d-1}|g(r)| vs r/sqrt(d), normalized to max 1."""
    curves = []
    fracs = []
    for d, col in ((16, GRAY), (64, BLUE), (256, RED)):
        alpha = mp.mpf(d) / 2 - 1
        l20 = mp.laguerre(2, alpha, 0)

        def g(r):
            q = 2 * mp.pi * r * r
            return (mp.laguerre(2, alpha, q) / l20 - 1) \
                * mp.e**(-mp.pi * r * r)

        def logdens(r):
            gv = g(r)
            if gv == 0:
                return mp.mpf(-9e9)
            return (d - 1) * mp.log(r) + mp.log(abs(gv))

        xs = [0.02 + 0.7 * i / 300 for i in range(301)]
        ld = [logdens(x * mp.sqrt(d)) + (d / 2) * mp.log(d) for x in xs]
        mx = max(ld)
        ys = [float(mp.e**(v - mx)) for v in ld]
        curves.append((d, col, list(zip(xs, ys))))
        # interior mass fraction inside c*sqrt(d), c=0.25
        c = mp.mpf("0.25")
        tot = mp.quad(lambda r: r**(d - 1) * abs(g(r)), [0, mp.inf])
        inn = mp.quad(lambda r: r**(d - 1) * abs(g(r)),
                      [0, c * mp.sqrt(d)])
        fracs.append((d, float(inn / tot)))
    ratio_ok = fracs[0][1] > fracs[1][1] > fracs[2][1]
    check("c11 interior mass fraction decays (c=0.25)", ratio_ok,
          " ".join(f"d={d}:{f:.3g}" for d, f in fracs))
    ax = Ax(0, 0.72, 0, 1.12, W=460, H=240)
    body = [ax.axes(xticks=(0, 0.2, 0.4, 0.6), yticks=(0, 0.5, 1),
                    xlab="r / √d", ylab="radial mass of |g| (norm.)")]
    body.append(ax.vline(1 / math.pi, color=GREEN, w=1.6, dash="6 4",
                         y0=0, y1=1.02))
    body.append(ax.text(1 / math.pi + 0.012, 1.09, "1/π", fs=12,
                        color=GREEN))
    for d, col, pts in curves:
        body.append(ax.polyline(pts, col, 2.0))
    # legend stacked in the top-left dead zone (all curves < 0.4 there)
    labels = [("d=16 (broadest)", GRAY, 0.045, 0.99),
              ("d=64", BLUE, 0.045, 0.87),
              ("d=256 (narrowest)", RED, 0.045, 0.75)]
    for name, col, lx, ly in labels:
        body.append(ax.text(lx, ly, name, fs=11, color=col))
    write_fig(OUT, "fig_mass", svg_wrap(
        "\n".join(body), 460, 240,
        "Mass concentration of Fourier eigenfunctions"))
    return fracs


# ------------------------------------------------------------ simple plots
def gauss_fig():
    f = lambda x: math.exp(-2 * math.pi * x * x)  # noqa: E731
    fh = lambda x: math.exp(-math.pi * x * x / 2) / math.sqrt(2)  # noqa
    ax = Ax(-2.2, 2.2, 0, 1.1, W=460, H=210)
    body = [ax.axes(xticks=(-2, -1, 0, 1, 2), yticks=(0.5, 1), xlab="x")]
    body.append(ax.polyline(sample(f, -2.2, 2.2), RED, 2.2))
    body.append(ax.polyline(sample(fh, -2.2, 2.2), BLUE, 2.2, dash="7 4"))
    body.append(ax.text(0.02, 1.045, "f(x) = e^", fs=0, color=RED))
    body = body[:-1]
    body.append(ax.text(0.28, 0.83, "narrow f", fs=12, color=RED))
    body.append(ax.text(1.15, 0.42, "wide  f̂", fs=12, color=BLUE))
    write_fig(OUT, "fig_gauss", svg_wrap(
        "\n".join(body), 460, 210, "Gaussian and its Fourier transform"))


def tent_figs():
    tent = lambda x: max(1 - abs(x), 0.0)  # noqa: E731

    def sinc2(x):
        if x == 0:
            return 1.0
        return (math.sin(math.pi * x) / (math.pi * x))**2

    ax = Ax(-2.6, 2.6, -0.28, 1.1, W=460, H=240)
    body = [ax.axes(xticks=(-2, -1, 0, 1, 2), yticks=(0.5, 1), xlab="x",
                    x_axis_at=0.0)]
    body.append(ax.polyline(sample(tent, -2.6, 2.6, 400), RED, 2.4))
    body.append(ax.polyline(sample(sinc2, -2.6, 2.6, 400), BLUE, 2.0,
                            dash="7 4"))
    body.append(ax.vline(1, color=GRAY, w=1, dash="3 3", y0=0.0, y1=1.0))
    body.append(ax.vline(-1, color=GRAY, w=1, dash="3 3", y0=0.0, y1=1.0))
    body.append(ax.text(-0.62, 0.94, "f = tent", fs=12, color=RED))
    body.append(ax.text(1.28, 0.5,
                        'f̂ = sinc² ≥ 0', fs=12, color=BLUE))
    body.append(ax.text(1.03, -0.18, "f ≤ 0 for |x| ≥ 1", fs=11,
                        color=MUT))
    write_fig(OUT, "fig_tent", svg_wrap(
        "\n".join(body), 460, 240, "The d=1 certificate pair"))
    # bridge figure: g = sinc^2 - tent
    g = lambda x: sinc2(x) - tent(x)  # noqa: E731
    ax = Ax(0, 3.0, -0.12, 0.12, W=460, H=230)
    pts = sample(g, 0.001, 3.0, 500)
    body = [ax.axes(xticks=(0, 2, 3), yticks=(-0.1, 0, 0.1),
                    xlab="r", x_axis_at=0.0)]
    neg = [(x, y) for x, y in pts if x <= 1.0]
    body.append(ax.path_area(neg, RED, opacity=0.25))
    pos = [(x, y) for x, y in pts if x >= 1.0]
    body.append(ax.path_area(pos, BLUE, opacity=0.25))
    body.append(ax.polyline(pts, DARK, 2.2))
    body.append(ax.vline(1.0, color=GRAY, w=1.2, dash="4 3"))
    body.append(ax.dot(0, 0, 4, GREEN))
    body.append(ax.text(0.10, 0.095, "g(0) = 0", fs=11, color=GREEN))
    body.append(ax.polyline([(1.38, -0.075), (0.64, -0.083)], RED, 1.0))
    body.append(ax.text(1.42, -0.075, "negative mass trapped in r < 1",
                        fs=11, color=RED))
    body.append(ax.text(1.7, 0.055, "g ≥ 0 outside", fs=11, color=BLUE))
    body.append(ax.text(1.06, 0.1, "r = 1/a", fs=11, color=MUT))
    write_fig(OUT, "fig_bridge", svg_wrap(
        "\n".join(body), 460, 230,
        "The anti-self-Fourier witness g = h-hat minus h"))
    # int_0^inf sinc^2 = 1/2 = int_0^inf tent exactly; numeric with tail
    intg = simpson(g, 0, 60, 40000)
    tail = 1 / (2 * math.pi**2 * 60)  # avg of sinc^2 tail beyond 60
    check("c20a int g = 0 (tent bridge)", abs(intg + tail) < 2e-4,
          f"{intg} (+tail {tail:.2e})")
    # ghat = -g numerically (cos transform)
    gh = lambda xi: 2 * simpson(  # noqa: E731
        lambda x: g(x) * math.cos(2 * math.pi * x * xi), 0, 40, 20000)
    err = max(abs(gh(xi) + g(xi)) for xi in (0.3, 0.8, 1.7))
    check("c20b ghat = -g (tent bridge)", err < 2e-3, f"{err}")
    negmass = -simpson(lambda x: min(g(x), 0.0), 0, 1, 2000)
    tot = simpson(lambda x: abs(g(x)), 0, 60, 40000) \
        + 1 / (2 * math.pi**2 * 60)
    check("c20c neg part = ||g||1/2 inside r<1",
          abs(2 * negmass / tot - 1) < 5e-3, f"ratio={2*negmass/tot:.5f}")
    return g


def td_fig(g):
    """T_d in d=1 (lambda=1/2) applied to g = sinc^2 - tent."""
    lam = 0.5

    def Tg(x):
        if x == 0:
            return 0.0
        # substitute s = t*x: (lam/2) x^{-lam} int_x^inf s^{lam-1} g(s) ds
        val = simpson(lambda s: s**(lam - 1) * g(s), x, 60, 8000)
        return lam / 2 * x**(-lam) * val

    n1g = simpson(lambda x: abs(g(x)), 0.0005, 60, 40000)
    n1T = simpson(lambda x: abs(Tg(x)), 0.0005, 60, 4000)
    check("c21a ||Tg||1 <= ||g||1/2", n1T <= n1g / 2 * 1.001,
          f"{n1T:.5f} vs {n1g/2:.5f}")
    # last sign change of Tg
    xs = [0.02 * i for i in range(1, 100)]
    vals = [Tg(x) for x in xs]
    rT = 0.0
    for i in range(len(xs) - 1):
        if vals[i] < 0 <= vals[i + 1]:
            rT = xs[i + 1]
    check("c21b r(Tg) < r(g)=1", 0 < rT < 1, f"r(Tg)~{rT:.2f}")
    # self-Fourier numerically
    Tgh = lambda xi: 2 * simpson(  # noqa: E731
        lambda x: Tg(x) * math.cos(2 * math.pi * x * xi), 0.0005, 40, 8000)
    err = max(abs(Tgh(xi) - Tg(xi)) for xi in (0.4, 1.1))
    check("c21c (Tg)hat = Tg", err < 5e-3, f"{err}")
    ax = Ax(0, 2.6, -0.12, 0.12, W=460, H=230)
    ptsg = sample(g, 0.001, 2.6, 400)
    ptsT = [(x, Tg(x)) for x, _ in sample(lambda x: 0, 0.02, 2.6, 200)]
    body = [ax.axes(xticks=(0, 2), yticks=(-0.1, 0, 0.1), xlab="r",
                    x_axis_at=0.0)]
    body.append(ax.polyline(ptsg, GRAY, 1.8, dash="5 4"))
    body.append(ax.polyline(ptsT, PUR, 2.4))
    body.append(ax.vline(1.0, color=GRAY, w=1, dash="3 3"))
    body.append(ax.vline(rT, color=PUR, w=1.4, dash="6 3"))
    body.append(ax.text(1.03, 0.1, "r(g) = 1", fs=11, color=MUT))
    body.append(ax.text(0.06, -0.109, f"r(Tg) ≈ {rT:.2f} →", fs=11,
                        color=PUR))
    body.append(ax.text(1.55, -0.05, "g (anti-self-Fourier)", fs=11,
                        color=GRAY))
    body.append(ax.polyline([(1.30, 0.083), (1.05, 0.014)], PUR, 1.0))
    body.append(ax.text(1.32, 0.086, "Tg (self-Fourier)", fs=11,
                        color=PUR))
    write_fig(OUT, "fig_Td", svg_wrap(
        "\n".join(body), 460, 230, "Tail integration shrinks the sign radius"))


def hex_fig():
    """Hexagonal packing patch, unit circles r=0.5, centers dist 1."""
    pts = []
    for i in range(0, 5):
        for j in range(0, 4):
            x = i + (j % 2) * 0.5
            y = j * math.sqrt(3) / 2
            pts.append((x, y))
    # isotropic scales so tangent unit disks are actually tangent
    ax = Ax(-0.62, 5.12, -0.62, 3.22, W=460, H=314, ml=10, mr=10,
            mt=10, mb=10)
    body = []
    for x, y in pts:
        body.append(f'<circle cx="{fmt(ax.X(x))}" cy="{fmt(ax.Y(y))}" '
                    f'r="{fmt(ax.X(0.5)-ax.X(0))}" fill="{BLUE}" '
                    f'opacity="0.25" stroke="{BLUE}" stroke-width="1.5"/>')
    for x, y in pts:
        body.append(f'<circle cx="{fmt(ax.X(x))}" cy="{fmt(ax.Y(y))}" '
                    f'r="2.4" fill="{DARK}"/>')
    # highlight one center pair at distance 1
    body.append(ax.polyline([(2, math.sqrt(3)), (2.5, math.sqrt(3) / 2 * 3)],
                            RED, 2.2))
    body.append(ax.text(2.62, 1.55 * math.sqrt(3) / 1.5 + 0.28, "1", fs=14,
                        color=RED, cls="v"))
    write_fig(OUT, "fig_hex", svg_wrap(
        "\n".join(body), 460, 314, "Hexagonal circle packing"))
    dens = mp.pi / mp.sqrt(12)
    check("c00 hex density pi/sqrt(12)", abs(dens - 0.9069) < 1e-3,
          f"{float(dens):.4f}")


def timeline_fig():
    data = [(1929, 0.5, "Blichfeldt 1929"),
            (1958, 0.5096, "Rogers 1958"),
            (1978, 0.5990, "Kabatianskii–Levenshtein 1978"),
            (2025, 0.6044, "this chapter")]
    ax = Ax(1900, 2036, 0.42, 1.1, W=460, H=250, ml=52)
    body = [ax.axes(xticks=(1920, 1960, 2000),
                    yticks=(0.5, 0.6, 1.0), xlab="year",
                    ylab='exponent α in Δd ≤ 2<tspan baseline-shift="super" font-size="8">−αd</tspan>',
                    xtickfmt=lambda v: str(int(v)))]
    # lower-bound line at alpha=1 (Minkowski: density >= 2^-d)
    body.append(ax.hline(1.0, color=GREEN, w=1.6, dash="6 4"))
    body.append(ax.text(1904, 1.035, "best packings known: α = 1",
                        fs=11, color=GREEN))
    # staircase of upper bounds
    steps = []
    for i, (yr, al, _) in enumerate(data):
        if i:
            steps.append((yr, data[i - 1][1]))
        steps.append((yr, al))
    steps.append((2036, data[-1][1]))
    body.append(ax.polyline(steps, RED, 2.4))
    for yr, al, name in data:
        body.append(ax.dot(yr, al, 3.5, RED))
    body.append(ax.text(1901, 0.458, "Blichfeldt ½", fs=10, color=MUT))
    body.append(ax.text(1959, 0.535, "Rogers .5096", fs=10, color=MUT))
    body.append(ax.text(1944, 0.625, "Kabatianskii–Levenshtein .5990",
                        fs=10, color=MUT))
    body.append(ax.text(1987, 0.674, "this chapter: .60440",
                        fs=11, color=RED))
    body.append(ax.text(1968, 0.9,
                        "the gap: nobody knows where Δd really lives",
                        fs=10, color=MUT))
    write_fig(OUT, "fig_timeline", svg_wrap(
        "\n".join(body), 460, 250, "Upper-bound exponent timeline"))


def expmap_fig():
    """Number line for Delta_d^{1/d} window."""
    lo, hi = 0.46, 0.72
    ax = Ax(lo, hi, 0, 1, W=460, H=150, ml=20, mr=20, mt=8, mb=30)
    y = 0.42
    body = []
    body.append(f'<line x1="{fmt(ax.X(lo))}" y1="{fmt(ax.Y(y))}" '
                f'x2="{fmt(ax.X(hi))}" y2="{fmt(ax.Y(y))}" '
                f'stroke="#333" stroke-width="1.6"/>')
    marks = [(0.5, GREEN, "½", "Minkowski packings"),
             (0.65774, RED, "0.6577", "new upper bound = LP limit"),
             (0.66017, ORA, "0.6602", "KL 1978"),
             (0.70711, GRAY, "0.7071", "Blichfeldt")]
    # shaded truth interval
    body.append(f'<rect x="{fmt(ax.X(0.5))}" y="{fmt(ax.Y(y)-9)}" '
                f'width="{fmt(ax.X(0.65774)-ax.X(0.5))}" height="18" '
                f'fill="{BLUE}" opacity="0.15"/>')
    lab_y = [0.72, 0.72, 0.16, 0.72]
    for (x, col, lab, note), ly in zip(marks, lab_y):
        body.append(f'<line x1="{fmt(ax.X(x))}" y1="{fmt(ax.Y(y)-10)}" '
                    f'x2="{fmt(ax.X(x))}" y2="{fmt(ax.Y(y)+10)}" '
                    f'stroke="{col}" stroke-width="2.4"/>')
        body.append(ax.text(x, ly, lab, fs=11, anchor="middle", color=col))
        body.append(ax.text(x, ly - 0.14, note, fs=9.5, anchor="middle",
                            color=MUT))
    body.append(ax.text(0.578, 0.30,
                        'Δd<tspan baseline-shift="super" font-size="7">1/d</tspan> lives somewhere in here', fs=10.5,
                        anchor="middle", color=BLUE))
    body.append(ax.text(hi, 0.06, 'Δd<tspan baseline-shift="super" font-size="7">1/d</tspan> →', fs=11,
                        anchor="end", color=DARK))
    write_fig(OUT, "fig_expmap", svg_wrap(
        "\n".join(body), 460, 150, "Where the packing rate lives"))


def vd_fig():
    ds = list(range(1, 65))
    ys = [float(mp.log10(mp.pi**(mp.mpf(d) / 2)
                / mp.gamma(mp.mpf(d) / 2 + 1) / 2**d)) for d in ds]
    ax = Ax(0, 64, -100, 5, W=460, H=230)
    body = [ax.axes(xticks=(8, 24, 48, 64), yticks=(0, -40, -80),
                    xlab="dimension d")]
    body.append(ax.text(2, -88, 'log₁₀ of ball volume vd/2<tspan '
                        'baseline-shift="super" font-size="8">d</tspan>',
                        fs=12, color=DARK))
    body.append(ax.polyline(list(zip(ds, ys)), BLUE, 2.2))
    body.append(ax.dot(8, ys[7], 3.5, RED))
    body.append(ax.dot(24, ys[23], 3.5, RED))
    body.append(ax.polyline([(19.2, -31.5), (9.5, -5.5)], GRAY, 1.0))
    body.append(ax.polyline([(33.4, -52.5), (26.2, -14.5)], GRAY, 1.0))
    body.append(ax.text(20, -35, "d=8: 0.0025", fs=10.5, color=RED))
    body.append(ax.text(34, -57, "d=24: 1.9×10⁻¹⁰",
                        fs=10.5, color=RED))
    write_fig(OUT, "fig_vd", svg_wrap(
        "\n".join(body), 460, 230, "A unit-separation ball's volume"))


def phi_fig():
    lam = 4.0
    f = lambda v: math.exp(lam * v - math.pi * math.exp(2 * v))  # noqa
    mx = max(f(v) for v, _ in sample(lambda v: 0, -4, 1.5, 300))
    g = lambda v: f(v) / mx  # noqa: E731
    ax = Ax(-4, 1.5, 0, 1.12, W=460, H=210)
    body = [ax.axes(xticks=(-4, -3, -2, -1, 0, 1), yticks=(0.5, 1),
                    xlab="v = log r")]
    body.append(ax.polyline(sample(g, -4, 1.5, 400), PUR, 2.4))
    body.append(ax.text(-3.6, 0.75,
                        'Φg(v) = e<tspan baseline-shift="super" font-size="8">λv</tspan> g(e<tspan baseline-shift="super" font-size="8">v</tspan>)', fs=12, color=PUR))
    body.append(ax.text(-3.6, 0.60, "(Gaussian g, λ = 4)", fs=10.5,
                        color=MUT))
    body.append(ax.text(-3.5, 0.17, 'e<tspan baseline-shift="super" font-size="8">λv</tspan> decay', fs=10.5, color=MUT))
    body.append(ax.text(0.4, 0.17, "super-exp decay", fs=10.5, color=MUT))
    write_fig(OUT, "fig_phi", svg_wrap(
        "\n".join(body), 460, 210, "Log-radius profile of a Gaussian"))


def psigma_fig():
    ax = Ax(-4, 4, 0, 0.65, W=460, H=220)
    body = [ax.axes(xticks=(-4, -2, 0, 2, 4), yticks=(0.2, 0.4, 0.6),
                    xlab="T")]
    for s, col in ((-0.5, GRAY), (0.0, BLUE), (0.7, PUR), (0.95, RED)):
        th = math.pi * (1 + s) / 2

        def P(T, th=th):
            return math.sin(th) / (4 * (math.cosh(math.pi * T / 2)
                                        - math.cos(th)))

        body.append(ax.polyline(sample(P, -4, 4, 400), col, 2.0))
    body.append(ax.text(0.42, 0.57, "σ = −0.5", fs=11, color=GRAY))
    body.append(ax.text(0.75, 0.27, "σ = 0", fs=11, color=BLUE))
    body.append(ax.text(1.35, 0.115, "σ = 0.7", fs=11, color=PUR))
    body.append(ax.text(2.5, 0.05, "σ = 0.95", fs=11, color=RED))
    write_fig(OUT, "fig_Psigma", svg_wrap(
        "\n".join(body), 460, 220, "Lower-edge harmonic measure kernels"))


def pconv_fig():
    ax = Ax(-3, 3, 0, 0.9, W=460, H=220)
    body = [ax.axes(xticks=(-3, -2, -1, 0, 1, 2, 3),
                    yticks=(0.4, 0.8), xlab="u")]
    for s, col in ((0.5, GRAY), (0.9, BLUE), (0.995, PUR)):
        th = math.pi * (1 + s) / 2

        def ps(u, th=th, s=s):
            return math.sin(th) / ((1 - s)
                                   * (math.cosh(math.pi * u) - math.cos(th)))

        body.append(ax.polyline(sample(ps, -3, 3, 400), col, 1.8))
    plim = lambda u: math.pi / 4 / math.cosh(math.pi * u / 2)**2  # noqa
    body.append(ax.polyline(sample(plim, -3, 3, 400), RED, 2.6, dash="2 3"))
    body.append(ax.text(0.0, 0.845, "limit p(u) = (π/4) sech²"
                        "(πu/2)", fs=11.5, color=RED, anchor="middle"))
    body.append(ax.text(-2.4, 0.42, "σ = 0.5", fs=11, color=GRAY))
    body.append(ax.text(1.62, 0.28, "σ = 0.9", fs=11, color=BLUE))
    body.append(ax.text(0.9, 0.62, "σ = 0.995", fs=11, color=PUR))
    write_fig(OUT, "fig_pconv", svg_wrap(
        "\n".join(body), 460, 220,
        "Normalized kernels converge to the sech-squared density"))


def thresh_fig():
    f = lambda c: math.log(math.pi**2 * c * c)  # noqa: E731
    ax = Ax(0.1, 0.55, -2.4, 1.2, W=460, H=230)
    body = [ax.axes(xticks=(0.2, 0.3, 0.4, 0.5), yticks=(-2, -1, 0, 1),
                    xlab="c")]
    pts = sample(f, 0.1, 0.55, 300)
    neg = [(x, y) for x, y in pts if x <= 1 / math.pi]
    body.append(ax.path_area(neg, GREEN, opacity=0.18))
    body.append(ax.polyline(pts, DARK, 2.4))
    body.append(ax.vline(1 / math.pi, color=GREEN, w=1.6, dash="6 4",
                         y0=-2.4, y1=0.85))
    body.append(ax.dot(1 / math.pi, 0, 4, GREEN))
    body.append(ax.text(1 / math.pi + 0.008, 1.0, "c = 1/π ≈ 0.318",
                        fs=11.5, color=GREEN))
    body.append(ax.text(0.125, -0.32, "rate < 0:", fs=11.5, color=GREEN))
    body.append(ax.text(0.125, -0.64, "interior mass dies", fs=11,
                        color=MUT))
    body.append(ax.text(0.24, 0.78, "log(π²c²)", fs=12.5, anchor="end",
                        color=DARK))
    write_fig(OUT, "fig_thresh", svg_wrap(
        "\n".join(body), 460, 230, "The exponential-rate threshold"))


def hlam_fig():
    lam, c = 16.0, 0.25
    R = c * math.sqrt(2 * lam)

    def h(y):
        t1 = lam * math.log(math.pi * R * R)
        t2 = float(mp.re(mp.loggamma(-1j * y / 2)))
        t3 = float(mp.re(mp.loggamma(lam + 1j * y / 2)))
        return t1 + t2 - t3

    ax = Ax(0, 40, -32, 10, W=460, H=230)
    body = [ax.axes(xticks=(10, 20, 30, 40), yticks=(-30, -20, -10, 0),
                    xlab="y")]
    pts = [(y, v) for y, v in sample(h, 0.12, 40, 600) if v <= 9.5]
    body.append(ax.polyline(pts, RED, 2.2))
    body.append(ax.hline(0, color=GRAY, w=1, dash="3 3", x1=38))
    body.append(ax.text(1.6, 7.6, "log spike to +∞ at y = 0", fs=10.5,
                        color=MUT))
    body.append(ax.text(14, -24, "eventually −λ log|y| decay",
                        fs=11, color=MUT))
    body.append(ax.text(24, 4.0, "hλ(y),  λ = 16, c = ¼",
                        fs=12, color=RED))
    write_fig(OUT, "fig_hlam", svg_wrap(
        "\n".join(body), 460, 230, "The lower-boundary majorant"))


def stair_fig():
    T = 0.6
    fT = lambda x: math.log(math.sqrt(x * x + T * T / 4))  # noqa: E731
    n = 8
    ax = Ax(-0.02, 1.05, -1.35, 0.25, W=460, H=220)
    body = [ax.axes(xticks=(0, 0.5, 1), yticks=(-1, -0.5, 0), xlab="x",
                    x_axis_at=0.0)]
    for k in range(n):
        x0, x1 = k / n, (k + 1) / n
        y = fT(x0)
        body.append(
            f'<rect x="{fmt(ax.X(x0))}" y="{fmt(min(ax.Y(y), ax.Y(0)))}" '
            f'width="{fmt(ax.X(x1)-ax.X(x0))}" '
            f'height="{fmt(abs(ax.Y(y)-ax.Y(0)))}" fill="{BLUE}" '
            f'opacity="0.25" stroke="{BLUE}" stroke-width="1"/>')
    body.append(ax.polyline(sample(fT, 0, 1.02, 300), RED, 2.4))
    body.append(ax.text(0.52, -0.92, "fT(x) = log √(x²+T²/4)",
                        fs=12, color=RED))
    body.append(ax.text(0.32, -1.13, "left Riemann sum (k/n, n = 8)",
                        fs=11, color=BLUE))
    write_fig(OUT, "fig_stair", svg_wrap(
        "\n".join(body), 460, 220, "Riemann sum against the log kernel"))


# --------------------------------------------------- saddle machinery (S4)
class Shells:
    def __init__(self, eps):
        self.eps = eps
        self.a0 = eps * eps
        self.A = math.log(1 / eps)
        self.B = eps**-3
        self.u0 = 1 + eps / 4
        self.U = 1 + eps / 2
        self.q = ((self.u0 - 1) + (self.U - 1)) / 2
        self.logQ = -self.q * self.B
        self.beta = self.u0 - 1

    def b(self, a):
        return 1 - 2 * self.eps * (1 + a)

    def ws(self, a):
        if not (self.a0 <= a <= self.A):
            return 0.0
        return -self.b(a) * math.exp(-2 * a) / (2 * a * a * math.cosh(a))

    def ws_int(self, f, n=1200):
        """integral of ws(a)*f(a) over [a0,A], log substitution a=e^s."""
        s0, s1 = math.log(self.a0), math.log(self.A)
        return simpson(lambda s: self.ws(math.exp(s)) * f(math.exp(s))
                       * math.exp(s), s0, s1, n)

    def wB_moment(self, kind, u=0.0, p=0, T=None, n=200):
        """log-space integral over [B,B+1] of wB(a)*factor:
        kind 'cosh':  cosh(u a) * a^p
        kind 'asinh': a * sinh(u a)
        kind 'osc':   cosh(u a) * (1 - cos(a T))
        wB(a) = Q / cosh(a). Values may be huge; float up to ~e^700."""
        def lcosh(x):
            x = abs(x)
            return x + math.log1p(math.exp(-2 * x)) - math.log(2)

        def lsinh(x):
            x = abs(x)
            return x + math.log1p(-math.exp(-2 * x)) - math.log(2)

        def g(a):
            lw = self.logQ - lcosh(a)
            if kind == "cosh":
                lv = lw + lcosh(u * a) + p * math.log(a)
            elif kind == "asinh":
                lv = lw + lsinh(u * a) + math.log(a)
            else:
                lv = lw + lcosh(u * a)
            if lv < -700:
                return 0.0
            if lv > 700:
                lv = 700.0  # clamp; only hit far past the plotted range
            val = math.exp(lv)
            if kind == "osc":
                val *= (1 - math.cos(a * T))
            return val
        return simpson(g, self.B, self.B + 1, n)

    def h_eps(self, zeta):
        """h_eps(zeta) = int w(a) (cos(a zeta) - 1) da, zeta complex.
        wB part negligible for |Im zeta| < 2 (checked separately)."""
        return self.ws_int(lambda a: cmath.cos(a * zeta) - 1)


def vsaddle(sh, lam, u):
    m = lam * (1 + u) / 2
    core = -0.5 * math.log(math.pi) + 0.5 * float(mp.digamma(m))
    shell = sh.ws_int(lambda a: a * math.sinh(u * a))
    shell += sh.wB_moment("asinh", u)
    return core + shell


def saddle_table():
    """R_{eps,d}/sqrt(d) for growing d, then the eps->0 limit trend.
    For fixed eps, digamma asymptotics give the limit
    sqrt((1+u0)/4pi) * exp(int w a sinh(u0 a) da)  -- eq. (84)."""
    rows = []
    sep = []
    for eps in (0.1, 0.05):
        sh = Shells(eps)
        disp = sh.ws_int(lambda a: a * math.sinh(sh.u0 * a)) \
            + sh.wB_moment("asinh", sh.u0)
        pred = math.sqrt((1 + sh.u0) / (4 * math.pi)) * math.exp(disp)
        for d in (64, 256, 1024, 4096):
            lam = d / 2
            R = math.exp(vsaddle(sh, lam, sh.u0))
            rows.append((eps, d, R / math.sqrt(d)))
        conv = abs(rows[-1][2] - pred) < 2e-3
        check(f"c17{'a' if eps == 0.1 else 'b'} R/sqrt(d) -> (84) limit "
              f"eps={eps}", conv,
              f"d=4096: {rows[-1][2]:.4f} vs limit {pred:.4f}")
        if eps == 0.1:
            check("c15 shell displacement -0.5 log(pi/2)+O(eps)",
                  abs(disp + 0.5 * math.log(math.pi / 2)) < 3 * eps,
                  f"disp={disp:.5f} target={-0.5*math.log(math.pi/2):.5f}")
        s1 = sh.B * math.exp(sh.logQ + (sh.u0 - 1) * sh.B)
        C0 = sh.A + 1 / sh.a0
        s2 = C0 * math.exp(-sh.logQ - (sh.U - 1) * (sh.B - sh.A))
        sep.append((s1, s2))
    check("c16 shell separations shrink with eps",
          sep[1][0] < sep[0][0] < 5e-3 and sep[1][1] < sep[0][1] < 5e-3,
          f"eps=0.1: {sep[0][0]:.1e},{sep[0][1]:.1e}  "
          f"eps=0.05: {sep[1][0]:.1e},{sep[1][1]:.1e}")
    # eps -> 0 trend of the limiting ratio (84): must approach 1/pi
    lims = []
    for eps in (0.1, 0.05, 0.02, 0.01):
        sh = Shells(eps)
        disp = sh.ws_int(lambda a: a * math.sinh(sh.u0 * a), n=4000) \
            + sh.wB_moment("asinh", sh.u0)
        lims.append((eps, math.sqrt((1 + sh.u0) / (4 * math.pi))
                     * math.exp(disp)))
    mono = all(lims[i][1] > lims[i + 1][1] for i in range(3))
    check("c17c limit(eps) decreasing -> 1/pi", mono
          and abs(lims[-1][1] - 1 / math.pi) < 0.006,
          " ".join(f"eps={e}: {v:.4f}" for e, v in lims)
          + f"  (1/pi={1/math.pi:.4f})")
    return rows


def vu_fig():
    eps, lam = 0.1, 32.0
    sh = Shells(eps)
    ustar = -1 + math.log(lam) / (4 * lam)
    # monotonicity check on the wide range (v huge but finite past U)
    uchk = [ustar + (1.25 - ustar) * i / 260 for i in range(261)]
    vchk = [vsaddle(sh, lam, u) for u in uchk]
    V_ok = all(vchk[i + 1] > vchk[i] for i in range(len(vchk) - 1))
    check("c18 v(u) strictly increasing on [u*,1.25]", V_ok)
    # plot only where e^v is plot-scale; past that the positive shell
    # sends the radius up a near-vertical wall
    umax = 1.033
    us = [ustar + (umax - ustar) * i / 260 for i in range(261)]
    vs = [vsaddle(sh, lam, u) for u in us]
    rs = [math.exp(min(v, 5.0)) / math.sqrt(2 * lam) for v in vs]
    # truncate at the wall: keep only points inside the plotted range
    keep = [i for i, r in enumerate(rs) if r <= 0.47]
    us = [us[i] for i in keep]
    rs = [rs[i] for i in keep]
    ax = Ax(-1.05, 1.32, 0, 0.5, W=460, H=240)
    body = [ax.axes(xticks=(-1, -0.5, 0, 0.5, 1),
                    yticks=(0.1, 0.2, 0.3, 0.4), xlab="u",
                    ylab='e<tspan baseline-shift="super" font-size="8">v(u)</tspan> / √d')]
    body.append(ax.hline(1 / math.pi, color=GREEN, w=1.4, dash="6 4",
                         x1=0.88))
    body.append(ax.polyline(list(zip(us, rs)), RED, 2.4))
    iu0 = min(range(len(us)), key=lambda i: abs(us[i] - sh.u0))
    body.append(ax.vline(sh.u0, color=BLUE, w=1.2, dash="4 3", y0=0,
                         y1=rs[iu0]))
    body.append(ax.dot(sh.u0, rs[iu0], 4, BLUE))
    body.append(ax.dot(ustar, rs[0], 4, PUR))
    body.append(ax.text(-1.02, 0.36, "1/π", fs=12, color=GREEN))
    body.append(ax.text(sh.u0 + 0.04, 0.06, "u₀", fs=13, color=BLUE))
    body.append(ax.text(-0.33, 0.038, "u∗: interior regime",
                        fs=11, color=PUR))
    body.append(ax.polyline([(-0.36, 0.038), (ustar + 0.02, 0.028)],
                            PUR, 1.0))
    body.append(ax.text(-0.6, 0.42,
                        "saddle radius sweeps out all radii", fs=11.5,
                        color=MUT))
    body.append(ax.text(-0.6, 0.385, "(d = 64, ε = 0.1)", fs=10.5,
                        color=MUT))
    body.append(ax.text(1.06, 0.335, "past U the", fs=10, color=ORA))
    body.append(ax.text(1.06, 0.30, "wB wall:", fs=10, color=ORA))
    body.append(ax.text(1.06, 0.265, "r → ∞", fs=10, color=ORA))
    # SMIL: dot sweeping the curve (relative path so the static
    # position — used by getBBox — is the first curve point)
    x0, y0 = ax.X(us[0]), ax.Y(rs[0])
    rel = " ".join(f"l {fmt(ax.X(u) - px)} {fmt(ax.Y(r) - py)}"
                   for (u, r, px, py) in
                   [(us[i], rs[i], ax.X(us[i - 10]), ax.Y(rs[i - 10]))
                    for i in range(10, len(us), 10)])
    body.append(
        f'<circle cx="{fmt(x0)}" cy="{fmt(y0)}" r="5" fill="{ORA}" '
        f'opacity="0.9"><animateMotion dur="7s" '
        f'repeatCount="indefinite" path="M 0 0 {rel}"/></circle>')
    write_fig(OUT, "fig_vu", svg_wrap(
        "\n".join(body), 460, 240, "The saddle radius as u sweeps"))
    return sh, lam, ustar


def damping_pieces(sh, lam, u):
    m = lam * (1 + u) / 2

    def Dg(T):
        return float(mp.re(mp.loggamma(m)) - mp.re(mp.loggamma(
            m - 1j * lam * T / 2)))

    def Ds(T):
        return -lam * sh.ws_int(lambda a: math.cosh(u * a)
                                * (1 - math.cos(a * T)))

    def DB(T):
        return lam * sh.wB_moment("osc", u, T=T)
    return Dg, Ds, DB


def damp_fig(sh, lam):
    u = sh.u0
    Dg, Ds, DB = damping_pieces(sh, lam, u)
    ax = Ax(0, 6, 0, 40, W=460, H=240)
    body = [ax.axes(xticks=(0, 2, 4, 6), yticks=(10, 20, 30, 40), xlab="T",
                    ylab="damping at u = u₀  (d = 64)")]
    ptsg = [(T, y) for T, y in sample(Dg, 0, 6, 200) if y <= 39.2]
    ptsu = [(T, Dg(T) - Ds(T) + DB(T)) for T, _ in sample(lambda T: 0,
                                                          0, 6, 120)]
    ptsu_pl = [(T, y) for T, y in ptsu if y <= 39.2]
    ptss = [(T, Ds(T)) for T, _ in sample(lambda T: 0, 0, 6, 120)]
    ptss_pl = [(T, y) for T, y in ptss if y <= 39.2]
    body.append(ax.polyline(ptsg, BLUE, 2.2))
    body.append(ax.polyline(ptsu_pl, RED, 2.4))
    body.append(ax.polyline(ptss_pl, GRAY, 1.8, dash="5 4"))
    body.append(ax.text(0.35, 34.5, "Dγ: gamma damping", fs=11.5,
                        color=BLUE))
    body.append(ax.text(4.4, 8.0, "Du = Dγ − Ds + DB", fs=11.5,
                        color=RED))
    body.append(ax.text(3.6, 5.2, "Ds: negative shell", fs=11, color=GRAY))
    Du_ok = all(y > 0 for T, y in ptsu if T > 0.001)
    check("c14b Du(T) > 0 at u0 on grid", Du_ok)
    write_fig(OUT, "fig_damp", svg_wrap(
        "\n".join(body), 460, 240, "Damping budget on the target contour"))
    # resonance panel: 1 - sinc(T/2) cos((B+1/2)T), illustrative B=6
    Bil = 6.0

    def res(T):
        s = 1.0 if T == 0 else math.sin(T / 2) / (T / 2)
        return 1 - s * math.cos((Bil + 0.5) * T)

    ax2 = Ax(0, 4, 0, 2.4, W=460, H=200)
    body2 = [ax2.axes(xticks=(0, 1, 2, 3, 4), yticks=(1, 2), xlab="T")]
    body2.append(ax2.polyline(sample(res, 0, 4, 700), PUR, 1.8))
    q = lambda T: min(T * T, 1.0) * 0.12  # noqa: E731
    body2.append(ax2.polyline(sample(q, 0, 4, 200), GREEN, 2.0, dash="6 4"))
    body2.append(ax2.text(1.15, 2.22, "interval shell: never returns to 0",
                          fs=11.5, color=PUR))
    body2.append(ax2.text(3.02, 0.29, "floor ≫ min(T²,1)", fs=11,
                          color=GREEN))
    write_fig(OUT, "fig_reson", svg_wrap(
        "\n".join(body2), 460, 200, "Interval support avoids resonances"))


def lemma42_check():
    eps = 0.1
    sh = Shells(eps)
    lam = 32.0
    worst = 1e9
    for iu in range(0, 22):
        u = -0.9 + (sh.U + 0.9) * iu / 21
        eta = 1 + u
        for ia in range(1, 40):
            a = sh.a0 + (sh.A - sh.a0) * ia / 40
            lhs = lam * abs(sh.ws(a)) * math.cosh(u * a)
            mu = math.exp(-eta * a) / (a * (1 - math.exp(-2 * a / lam)))
            worst = min(worst, 1 - lhs / mu)
    check("c14 (54) margin 1 - ratio >= c*eps", worst > 0.01,
          f"min margin={worst:.4f}")


def laplace_fig(sh_unused):
    # eps = 0.05: the positive shell is genuinely dormant at u0 here, so
    # the variance normalization sqrt(lam V) matches the active damping
    # (at eps = 0.1 and small d, V_B inflates V without damping -- the
    # lambda_eps finite-size effect discussed in S11).
    sh = Shells(0.05)
    u = sh.u0
    curves = []
    for lam, col in ((32.0, BLUE), (512.0, RED)):
        Dg, Ds, DB = damping_pieces(sh, lam, u)
        m = lam * (1 + u) / 2
        Vs = -sh.ws_int(lambda a: a * a * math.cosh(u * a))
        VB = sh.wB_moment("cosh", u, p=2)
        V = lam / 4 * float(mp.polygamma(1, m)) - Vs + VB
        sc = math.sqrt(lam * V)

        def env(x, Dg=Dg, Ds=Ds, DB=DB, sc=sc):
            T = x / sc
            return math.exp(-(Dg(T) - Ds(T) + DB(T)))

        curves.append((lam, col, sample(env, -3.5, 3.5, 160)))
    gsn = lambda x: math.exp(-x * x / 2)  # noqa: E731
    ax = Ax(-3.5, 3.5, 0, 1.12, W=460, H=220)
    body = [ax.axes(xticks=(-2, 0, 2), yticks=(0.5, 1),
                    xlab="x = T √(λV)")]
    body.append(ax.polyline(sample(gsn, -3.5, 3.5, 200), GREEN, 3.2,
                            opacity=0.45))
    for lam, col, pts in curves:
        body.append(ax.polyline(pts, col, 1.8))
    body.append(ax.text(-3.3, 1.02, 'Gaussian e<tspan baseline-shift="super" font-size="8">−x²/2</tspan>', fs=11,
                        color=GREEN))
    body.append(ax.text(1.0, 0.72, "d = 64", fs=11, color=BLUE))
    body.append(ax.text(1.62, 0.5, "d = 1024", fs=11, color=RED))
    write_fig(OUT, "fig_laplace", svg_wrap(
        "\n".join(body), 460, 220,
        "The centered integrand collapses onto a Gaussian"))


def poly_fig():
    beta = 0.025  # eps=0.1
    u0 = 1.025
    Pp = lambda u: beta + (1 - u)**2 * (1 + u)  # noqa: E731
    Pm = lambda u: beta + (1 - u) * (1 + u)**2  # noqa: E731
    P0 = lambda u: u * u - 1  # noqa: E731
    ax = Ax(-1.05, 1.45, -1.6, 2.2, W=460, H=250)
    body = [ax.axes(xticks=(-1, -0.5, 0, 0.5, 1), yticks=(-1, 0, 1, 2),
                    xlab="u", x_axis_at=0.0)]
    um = 1.45
    while Pm(um) < -1.58:
        um -= 0.005
    body.append(ax.polyline(sample(Pp, -1.05, 1.45, 300), BLUE, 2.4))
    body.append(ax.polyline(sample(Pm, -1.05, um, 300), RED, 2.4))
    body.append(ax.polyline(sample(P0, -1.05, 1.45, 300), GREEN, 2.4,
                            dash="6 4"))
    body.append(ax.vline(u0, color=GRAY, w=1.2, dash="4 3", y0=-1.6,
                         y1=1.75))
    body.append(ax.text(u0 - 0.06, 1.98, "u₀ = 1+ε/4", fs=11.5,
                        color=MUT, anchor="end"))
    body.append(ax.text(-0.88, 1.62, "P₊(iu) > 0 always", fs=12,
                        color=BLUE))
    body.append(ax.text(-0.40, -0.62, "P₋(iu)", fs=12, color=RED))
    body.append(ax.text(0.32, -1.44, "P₀(iu) = u²−1", fs=12,
                        color=GREEN))
    body.append(ax.dot(u0, Pm(u0), 3.5, RED))
    body.append(ax.dot(u0, P0(u0), 3.5, GREEN))
    check("c13 polynomial signs at u0",
          Pp(u0) > 0 and Pm(u0) < 0 and P0(u0) > 0,
          f"P+={Pp(u0):.4f} P-={Pm(u0):.5f} P0={P0(u0):.5f}")
    z = 0.3 + 0.7j
    Ppz = 1 + z * z + beta + 1j * z * (1 + z * z)
    Pmz = 1 + z * z + beta - 1j * z * (1 + z * z)
    Pmz_ref = 1 + z * z + beta - 1j * (-z) * (1 + z * z)
    _ = Pmz
    check("c13b reflection P-(-z) = P+(z)",
          abs((1 + z * z + beta - 1j * (-z) * (1 + z * z)) - Ppz) < 1e-15,
          f"{abs(Pmz_ref - Ppz):.1e}")
    write_fig(OUT, "fig_poly", svg_wrap(
        "\n".join(body), 460, 250, "Sign polynomials on the imaginary axis"))


def shells_fig():
    eps = 0.1
    sh = Shells(eps)
    wstar = lambda a: -math.exp(-2 * a) / (2 * a * a * math.cosh(a))  # noqa
    ymin = -13.0
    # start each curve where it enters the plotted range (avoid clipping:
    # both densities behave like -1/2a^2 near 0)
    a1 = 0.05
    while wstar(a1) < ymin:
        a1 += 0.002
    a2 = sh.a0
    while sh.ws(a2) < ymin:
        a2 += 0.002
    ax = Ax(0, 2.7, ymin, 1.5, W=460, H=250)
    body = [ax.axes(xticks=(0.5, 1, 1.5, 2, 2.5), yticks=(-10, -5, 0),
                    xlab="a")]
    body.append(ax.hline(0, color=GRAY, w=1, dash="2 3", x1=2.65))
    body.append(ax.polyline(sample(wstar, a1, 2.7, 500), GRAY, 1.8,
                            dash="5 4"))
    body.append(ax.polyline(sample(sh.ws, a2, sh.A, 600), RED, 2.4))
    body.append(ax.vline(sh.A, color=MUT, w=1, dash="3 3", y0=ymin,
                         y1=0.4))
    body.append(ax.text(0.30, -12.3, "both dive like −1/2a², cut at a₀ = ε²", fs=10.5, color=MUT))
    body.append(ax.text(sh.A + 0.03, 0.75, "A = log(1/ε)", fs=11,
                        color=MUT))
    body.append(ax.text(0.62, -8.5, "ideal w∗ (infinite mass at 0)",
                        fs=11.5, color=GRAY))
    body.append(ax.text(0.55, -3.2, "negative shell ws:", fs=11.5,
                        color=RED))
    body.append(ax.text(0.55, -4.6, "truncated + tapered", fs=11,
                        color=RED))
    body.append(ax.text(1.05, 1.05, "positive shell wB lives on [B, B+1],",
                        fs=10.5, color=BLUE))
    body.append(ax.text(1.05, 0.5, "B = ε⁻³ = 1000 →", fs=10.5,
                        color=BLUE))
    write_fig(OUT, "fig_shells", svg_wrap(
        "\n".join(body), 460, 250, "The two shells against the ideal"))


# --------------------------------------- inverse-Mellin profiles (S4 final)
# Demo parameters d = 4000 (lambda = 2000), eps = 0.15: large enough that
# the positive shell's residue corrections are genuinely negligible
# (log lam + logQ + 2NB/lam < 0), so the construction is in its asymptotic
# regime and the saddle formula / interior series are honest.

def lcosh(x):
    x = abs(x)
    return x + math.log1p(math.exp(-2 * x)) - math.log(2)


def lsinh(x):
    x = abs(x)
    if x < 1e-8:
        return math.log(x) if x > 0 else -746.0
    return x + math.log1p(-math.exp(-2 * x)) - math.log(2)


class Saddle:
    def __init__(self, eps=0.1, lam=2000.0):
        self.sh = Shells(eps)
        self.lam = lam

    def h_shell_c(self, zeta, n=1200):
        """h_eps(zeta) for complex zeta = T + iu, both shells.
        ws by Simpson (log grid); wB in log space (safe for u <= 1.06)."""
        sh = self.sh
        tot = sh.ws_int(lambda a: cmath.cos(a * zeta) - 1, n=n)
        T, u = zeta.real, zeta.imag
        B = sh.B

        def gB(a):
            lw = sh.logQ - lcosh(a)
            re = math.exp(lw + lcosh(u * a)) * math.cos(a * T) \
                - math.exp(lw)
            im = -math.exp(lw + lsinh(abs(u) * a)) * math.sin(a * T) \
                * (1 if u >= 0 else -1)
            return complex(re, im)
        tot += simpson(gB, B, B + 1, 100)
        return tot

    def V(self, u):
        sh, lam = self.sh, self.lam
        m = lam * (1 + u) / 2
        Vs = -sh.ws_int(lambda a: a * a * math.cosh(u * a))
        VB = sh.wB_moment("cosh", u, p=2)
        return lam / 4 * float(mp.polygamma(1, m)) - Vs + VB

    def Lu_grid(self, u, Ts, nws=300):
        """centered phase (46) on a T-grid, shared work hoisted per u."""
        sh, lam = self.sh, self.lam
        m = lam * (1 + u) / 2
        lgm = mp.loggamma(m)
        psim = float(mp.digamma(m))
        # precompute shell quadrature nodes/weights (log grid for ws)
        s0, s1 = math.log(sh.a0), math.log(sh.A)
        if nws % 2:
            nws += 1
        hs = (s1 - s0) / nws
        ws_nodes = []
        for i in range(nws + 1):
            s = s0 + i * hs
            a = math.exp(s)
            wgt = (1 if i in (0, nws) else (4 if i % 2 else 2)) * hs / 3
            ws_nodes.append((a, sh.ws(a) * a * wgt,
                             math.cosh(u * a), math.sinh(u * a)))
        B = sh.B
        nB = 60
        hB = 1.0 / nB
        wB_nodes = []
        for i in range(nB + 1):
            a = B + i * hB
            wgt = (1 if i in (0, nB) else (4 if i % 2 else 2)) * hB / 3
            lw = sh.logQ - lcosh(a)
            wB_nodes.append((a, math.exp(lw + lcosh(u * a)) * wgt,
                             math.exp(lw + lsinh(u * a)) * wgt))
        out = []
        for T in Ts:
            G = complex(mp.loggamma(m - 1j * lam * T / 2) - lgm) \
                + 1j * lam * T / 2 * psim
            cosr = sinr = 0.0
            for a, wsw, ch, sh_ in ws_nodes:
                cosr += wsw * ch * (math.cos(a * T) - 1)
                sinr += wsw * sh_ * (a * T - math.sin(a * T))
            for a, chB, shB in wB_nodes:
                cosr += chB * (math.cos(a * T) - 1)
                sinr += shB * (a * T - math.sin(a * T))
            out.append(G + lam * cosr + 1j * lam * sinr)
        return out

    def P(self, j, z):
        beta = self.sh.beta
        if j == 0:
            return -(1 + z * z)
        s = 1 if j > 0 else -1
        return 1 + z * z + beta + s * 1j * z * (1 + z * z)

    def I_norm_all(self, u, n=400):
        """sqrt(lam V / 2 pi) * I_{lam,P_j}(u) for j in (-1,0,1);
        Lemma 4.8 says these tend to P_j(iu)."""
        lam = self.lam
        V = self.V(u)
        Tr = 12 / math.sqrt(lam * V)
        h = 2 * Tr / n
        Ts = [-Tr + i * h for i in range(n + 1)]
        Ls = self.Lu_grid(u, Ts)
        tots = {-1: 0j, 0: 0j, 1: 0j}
        for i, (T, L) in enumerate(zip(Ts, Ls)):
            wgt = 1 if i in (0, n) else (4 if i % 2 else 2)
            eL = cmath.exp(L) * wgt
            z = complex(T, u)
            for j in (-1, 0, 1):
                tots[j] += eL * self.P(j, z)
        sc = math.sqrt(lam * V / (2 * math.pi)) * h / 3
        return {j: tots[j] * sc for j in (-1, 0, 1)}


def fprof_figs():
    # lam = 3e6: the positive-shell damping floor lam*Q*e^{(u0-1)B} must
    # exceed 1 for the saddle formula to bite near u0; for eps = 0.1 that
    # needs lam >> e^{12.5} ~ 3e5.  (At d = 4000 the formula is already
    # excellent for u <= 0.95 but wobbles near u0 -- a teachable point.)
    sd = Saddle(eps=0.1, lam=3.0e6)
    sh, lam = sd.sh, sd.lam
    # u grid: coarse below 0.9, fine 0.9..1.032 (positive-shell wall
    # sits at u = 1 + q_eps = 1.0375 for eps = 0.1)
    us = [-0.9 + (0.9 + 0.9) * i / 11 for i in range(12)] \
        + [0.9 + 0.132 * i / 28 for i in range(1, 29)]
    curves = {1: [], 0: [], -1: []}
    absdev = 0.0
    reldev = 0.0
    for u in us:
        vals = sd.I_norm_all(u)
        for j in (1, 0, -1):
            val = vals[j]
            curves[j].append((u, val.real))
            Pref = sd.P(j, complex(0, u)).real
            absdev = max(absdev, abs(val.real - Pref), abs(val.imag))
            if abs(Pref) > 0.5:
                reldev = max(reldev, abs(val.real - Pref) / abs(Pref))
    check("c19 saddle formula sqrt(lamV/2pi) I ~ P(iu)",
          absdev < 0.02 and reldev < 0.02,
          f"max abs dev = {absdev:.4f}, rel dev (|P|>0.5) = {reldev:.4f}"
          "  (d=6e6, eps=0.1)")
    # f- sign flip: happens at the root of P-(iu) = beta + (1-u)(1+u)^2,
    # which sits just below u0 (the paper evaluates at u0 for sign margin)
    uroot = float(mp.findroot(lambda u: sh.beta + (1 - u) * (1 + u)**2,
                              1.006))
    vm = curves[-1]
    ucross = None
    for i in range(len(vm) - 1):
        if vm[i][1] >= 0 > vm[i + 1][1]:
            ucross = 0.5 * (vm[i][0] + vm[i + 1][0])
    check("c19a f- sign flip at root of P-(iu) just below u0",
          ucross is not None and abs(ucross - uroot) < 0.006,
          f"cross at u={ucross:.4f}, root {uroot:.4f}, u0={sh.u0}")
    beta = sh.beta
    ax = Ax(-0.95, 1.12, -1.6, 2.3, W=460, H=260)
    body = [ax.axes(xticks=(-0.5, 0, 0.5, 1), yticks=(-1, 0, 1, 2),
                    x_axis_at=0.0)]
    body.append(ax.text(1.10, -1.45, "u", fs=13, anchor="end",
                        color=DARK, cls="v"))
    refs = {1: lambda u: beta + (1 - u)**2 * (1 + u),
            -1: lambda u: beta + (1 - u) * (1 + u)**2,
            0: lambda u: u * u - 1}
    cols = {1: BLUE, -1: RED, 0: GREEN}
    for j in (1, -1, 0):
        body.append(ax.polyline(sample(refs[j], -0.95, 1.12, 200),
                                cols[j], 1.4, dash="5 4", opacity=0.7))
        for u, y in curves[j]:
            body.append(ax.dot(u, y, 2.6, cols[j]))
    body.append(ax.vline(sh.u0, color=GRAY, w=1.2, dash="4 3"))
    body.append(ax.text(sh.u0 - 0.02, 2.12, "u₀", fs=12, color=MUT,
                        anchor="end"))
    body.append(ax.text(-0.85, 1.9, "dots: computed √(λV/2π)·I(u), "
                        "d = 6×10⁶", fs=11, color=MUT))
    body.append(ax.text(-0.85, 1.62, "dashes: predicted P(iu)", fs=11,
                        color=MUT))
    body.append(ax.text(-0.28, 1.32, "P₊", fs=13, color=BLUE))
    body.append(ax.text(-0.28, 0.44, "P₋", fs=13, color=RED))
    body.append(ax.text(-0.28, -0.70, "P₀", fs=13, color=GREEN))
    write_fig(OUT, "fig_fprofiles", svg_wrap(
        "\n".join(body), 460, 260,
        "The saddle integral lands on the sign polynomials"))
    return sd


def interior_fig(sd):
    """f+(r)/f+(0) vs e^{-y} on [0, r*], contour at u*, mpmath log-sum."""
    sh, lam = sd.sh, sd.lam
    ustar = -1 + math.log(lam) / (4 * lam)
    rstar = math.exp(vsaddle(sh, lam, ustar))
    h1p = sh.ws_int(lambda a: a * math.sinh(a)) + sh.wB_moment("asinh", 1.0)
    # log f+(0) = log(2 beta) + (lam/2) log pi + lam h_eps(i)
    h_i = sh.ws_int(lambda a: math.cosh(a) - 1) \
        + sh.wB_moment("cosh", 1.0) - sh.wB_moment("cosh", 0.0)
    logf0 = mp.log(2 * sh.beta) + lam / 2 * mp.log(mp.pi) + lam * h_i
    n = 500
    Tr = 12 / math.sqrt(lam * sd.V(ustar))
    hstep = 2 * Tr / n
    pre = []
    for i in range(n + 1):
        T = -Tr + i * hstep
        t = lam * complex(T, ustar)
        lg = mp.loggamma((lam - 1j * t) / 2)
        he = sd.h_shell_c(complex(T, ustar))
        Pp = sd.P(1, complex(T, ustar))
        pre.append((T, lg + 1j * t / 2 * mp.log(mp.pi)
                    + lam * mp.mpc(he), mp.mpc(Pp)))
    rows = []
    for k in range(0, 41):
        r = rstar * k / 40
        if r == 0:
            rows.append((0.0, 1.0, 1.0))
            continue
        lr = math.log(r)
        tot = mp.mpc(0)
        for i, (T, base, Pp) in enumerate(pre):
            wgt = 1 if i in (0, n) else (4 if i % 2 else 2)
            tot += wgt * Pp * mp.e**(base + 1j * lam * T * lr
                                     - lam * (1 + ustar) * lr - logf0)
        val = lam / (2 * mp.pi) * tot * hstep / 3
        y = math.pi * math.exp(2 * h1p) * r * r
        rows.append((r, float(mp.re(val)), math.exp(-y)))
    supd = max(abs(s - e) for _, s, e in rows)
    check("c19b f+/f+(0) ~ e^{-y} on [0,r*]", supd < 0.05,
          f"sup diff = {supd:.4f}, r* = {rstar:.3f}")
    check("c19c f+ > 0 on [0,r*]", all(s > 0 for _, s, _ in rows))
    ax = Ax(0, rstar, 0, 1.1, W=460, H=220)
    body = [ax.axes(xticks=(0, 0.3, 0.6), yticks=(0.5, 1), xlab="r")]
    body.append(ax.polyline([(r, e) for r, _, e in rows], GREEN, 3.4,
                            opacity=0.45))
    body.append(ax.polyline([(r, s) for r, s, _ in rows], PUR, 1.8))
    body.append(ax.text(rstar * 0.05, 0.5,
                        "f₊(r)/f₊(0) (purple, computed)", fs=11.5,
                        color=PUR))
    body.append(ax.text(rstar * 0.05, 0.36, "e⁻ʸ reference (green)",
                        fs=11.5, color=GREEN))
    body.append(ax.text(rstar * 0.05, 0.18,
                        f"agree to {supd:.1e} on [0, r∗];  d = 6×10⁶",
                        fs=10.5, color=MUT))
    write_fig(OUT, "fig_fplus_int", svg_wrap(
        "\n".join(body), 460, 220, "Interior positivity of f-plus"))
    # envelope symmetry checks at lam=2000 (all in mpmath: huge gammas)
    def mlam(t):
        return mp.pi**(1j * t) * mp.gamma((lam - 1j * t) / 2) \
            / mp.gamma((lam + 1j * t) / 2)

    def E(t):
        he = sh.ws_int(lambda a: math.cos(a * t / lam) - 1)
        return mp.pi**(1j * t / 2) * mp.gamma((lam - 1j * t) / 2) \
            * mp.e**(lam * he)

    err = max(float(abs(E(t) - mlam(t) * E(-t)) / abs(E(t)))
              for t in (0.7, 3.3))
    check("c12b envelope symmetry E(t)=m(t)E(-t)", err < 1e-10,
          f"{err:.1e}")
    err2 = max(float(abs(abs(mlam(t)) - 1)) for t in (0.5, 2.0, 11.0))
    check("c12 |m_lambda(t)| = 1", err2 < 1e-12, f"{err2:.1e}")
    # saddle radius for the demo parameters
    R = math.exp(vsaddle(sh, lam, sh.u0))
    print(f"INFO demo sign radius R/sqrt(d) = {R/math.sqrt(2*lam):.4f} "
          f"(1/pi = {1/math.pi:.4f}), r* = {rstar:.3f}")


# ---------------------------------------------------------------- SMIL figs
def poisson_anim():
    s = 0.55
    f = lambda x: math.exp(-math.pi * x * x / s)  # noqa: E731
    ax = Ax(-1.6, 1.6, 0, 2.1, W=460, H=230)
    body = [ax.axes(xticks=(-1, 0, 1), yticks=(1, 2), xlab="x")]
    total = lambda x: sum(f(x + n) for n in range(-4, 5))  # noqa: E731
    for k, n in enumerate((-2, -1, 0, 1, 2)):
        pts = sample(lambda x, n=n: f(x + n), -1.6, 1.6, 200)
        beg = 0.9 * k
        body.append(
            ax.polyline(pts, BLUE, 1.6).replace(
                "/>",
                f' opacity="0"><animate attributeName="opacity" '
                f'values="0;0;0.75;0.75" keyTimes="0;{beg/9:.3f};'
                f'{(beg+0.7)/9:.3f};1" dur="9s" '
                f'repeatCount="indefinite"/></polyline>'))
    ptsT = sample(total, -1.6, 1.6, 240)
    body.append(
        ax.polyline(ptsT, RED, 2.6).replace(
            "/>",
            ' opacity="0"><animate attributeName="opacity" '
            'values="0;0;1;1" keyTimes="0;0.55;0.72;1" dur="9s" '
            'repeatCount="indefinite"/></polyline>'))
    body.append(ax.hline(math.sqrt(s), color=GREEN, w=1.2, dash="5 4"))
    body.append(ax.text(-1.55, 2.0, "translates f(x+n) sum to a periodic "
                        "function", fs=11, color=MUT))
    body.append(ax.text(0.30, 1.22, "mean value = f̂(0)", fs=11,
                        color=GREEN))
    write_fig(OUT, "fig_poisson_anim", svg_wrap(
        "\n".join(body), 460, 230, "Periodization of a Gaussian"))


def dilate_anim():
    lam = 4.0
    g = lambda r: r * r * math.exp(-math.pi * r * r) * 8  # noqa: E731

    def avg(r, a):
        return 0.5 * (math.exp(a) * g(r * math.exp(a / lam))
                      + math.exp(-a) * g(r * math.exp(-a / lam)))

    ax = Ax(0, 2.4, 0, 1.88, W=460, H=230)
    body = [ax.axes(xticks=(0, 1, 2), yticks=(0.5, 1, 1.5), xlab="r")]
    body.append(ax.polyline(sample(g, 0, 2.4, 200), GRAY, 1.8, dash="5 4"))
    frames = []
    for a in (0.0, 0.3, 0.6, 0.9, 1.2, 0.9, 0.6, 0.3, 0.0):
        pts = sample(lambda r, a=a: avg(r, a), 0, 2.4, 120)
        d = "M" + " L".join(ax.P(x, y) for x, y in pts)
        frames.append(d)
    body.append(
        f'<path d="{frames[0]}" fill="none" stroke="{RED}" '
        f'stroke-width="2.4"><animate attributeName="d" '
        f'values="{";".join(frames)}" dur="8s" '
        f'repeatCount="indefinite"/></path>')
    body.append(ax.text(0.06, 1.76, '½(eᵃ g(re<tspan baseline-shift="super" font-size="8">a/λ</tspan>) + e⁻ᵃ g(re<tspan baseline-shift="super" font-size="8">−a/λ</tspan>))',
                        fs=12, color=RED))
    body.append(ax.text(0.06, 1.62, "a swings 0 → 1.2 → 0",
                        fs=10.5, color=MUT))
    body.append(ax.text(1.5, 0.32, "g (a = 0)", fs=11, color=GRAY))
    write_fig(OUT, "fig_dilate_anim", svg_wrap(
        "\n".join(body), 460, 230, "The dilation pair behind cos(at/lambda)"))


def main():
    os.makedirs(OUT, exist_ok=True)
    only = sys.argv[1] if len(sys.argv) > 1 else None

    def want(name):
        return only is None or only == name

    if want("const"):
        print("== constants ==")
        constants()
    if want("mellin"):
        print("== mellin + eigen (d=3) ==")
        mellin_checks()
    if want("toy"):
        print("== toy + bridge ==")
        g = tent_figs()
        td_fig(g)
    if want("simple"):
        print("== simple figures ==")
        hex_fig()
        timeline_fig()
        expmap_fig()
        vd_fig()
        gauss_fig()
        phi_fig()
        psigma_fig()
        pconv_fig()
        thresh_fig()
        hlam_fig()
        stair_fig()
        shells_fig()
        poly_fig()
        poisson_anim()
        dilate_anim()
    if want("mass"):
        print("== mass concentration ==")
        mass_fig()
    if want("saddle"):
        print("== saddle machinery ==")
        lemma42_check()
        saddle_table()
        sh, lam, ustar = vu_fig()
        damp_fig(sh, lam)
        laplace_fig(sh)
    if want("prof"):
        print("== inverse-Mellin profiles (slow) ==")
        sd = fprof_figs()
        interior_fig(sd)
    print()
    bad = [n for n, ok in CHECKS if not ok]
    print(f"TOTAL {len(CHECKS)} checks, {len(CHECKS)-len(bad)} OK, "
          f"{len(bad)} FAIL {bad if bad else ''}")


if __name__ == "__main__":
    main()
