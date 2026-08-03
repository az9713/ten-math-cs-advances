"""Figure-data generator + numeric validation for ch08 explainer.

Run: PYTHONIOENCODING=utf-8 python figgen/ch08_figs.py
Outputs: SVG coordinate strings and validation numbers quoted in the HTML.
"""
import math


# ---- 1. Triangle 3*Delta_2 - (1,1): interior lattice points ----
# K = {x > -1, y > -1, x + y < 1}, vertices (-1,-1), (2,-1), (-1,2)
def interior_pts(k):
    pts = []
    for x in range(-3 * k, 3 * k + 1):
        for y in range(-3 * k, 3 * k + 1):
            if x > -k and y > -k and x + y < k:
                pts.append((x, y))
    return pts


print("interior lattice pts of K (k=1):", interior_pts(1))
print("d_k = #int(kK) cap Z^2 vs 4.5 k^2:")
for k in (1, 2, 3, 5, 10, 20):
    dk = len(interior_pts(k))
    print(f"  k={k:2}  d_k={dk:5}  d_k/k^2={dk/k**2:.3f}  (vol=4.5)")

# ---- 2. n=1 explicit transport potential ----
# phi(x) = 2 log cosh(x/2) + log 2;  grad phi = tanh(x/2) maps R -> (-1,1)
# e^{-phi} = (1/2) sech^2(x/2); total mass = 2 = vol([-1,1]).  Check:
import numpy as np

xs = np.linspace(-40, 40, 400001)
dx = xs[1] - xs[0]
phi = 2 * np.log(np.cosh(xs / 2)) + math.log(2)
mass = np.trapezoid(np.exp(-phi), xs)
bary = np.trapezoid(xs * np.exp(-phi), xs)
print(f"n=1 check: mass={mass:.6f} (want 2), barycenter={bary:.2e}")
sup_dev = np.max(np.abs(phi - np.abs(xs)))
print(f"sup |phi - h_K| = {sup_dev:.4f} (bounded; h_K(x)=|x|)")


# ---- 3. I_k(u) and Bergman kernel curves for the n=1 toy ----
def I_k(k, u):
    f = np.exp(k * (u * xs - phi))
    return np.trapezoid(f, xs)


def log_Bk_over_k(k, xg):
    tot = np.zeros_like(xg)
    vals = []
    for m in range(-(k - 1), k):
        vals.append((m / k, I_k(k, m / k)))
    out = []
    for x in xg:
        terms = [k * u * x - math.log(Ik) for u, Ik in vals]
        mx = max(terms)
        s = sum(math.exp(t - mx) for t in terms)
        out.append((mx + math.log(s)) / k)
    return np.array(out)


xg = np.linspace(-4, 4, 33)
phig = 2 * np.log(np.cosh(xg / 2)) + math.log(2)
for k in (2, 8, 32):
    lb = log_Bk_over_k(k, xg)
    dev = np.max(np.abs(lb - phig))
    print(f"k={k:3}: max |(1/k)log B_k - phi| on [-4,4] = {dev:.4f}")

# SVG polylines for Fig (viewBox 0 0 480 260): x in [-4,4] -> [50,450],
# y values phi in [0, ~4.2] -> [230, 20]


def X(x):
    return 50 + (x + 4) / 8 * 400


def Y(v):
    return 230 - v / 4.2 * 210


print("phi_line:", " ".join(
    f"{X(x):.1f},{Y(v):.1f}" for x, v in zip(xg, phig)))
print("hK_line:", f"{X(-4):.1f},{Y(4):.1f} {X(0):.1f},{Y(0):.1f}"
      f" {X(4):.1f},{Y(4):.1f}")
for k in (2, 8, 32):
    lb = log_Bk_over_k(k, xg)
    print(f"B{k}_line:", " ".join(
        f"{X(x):.1f},{Y(max(v, -0.1)):.1f}" for x, v in zip(xg, lb)))

# ---- 4. sharp-bound rate chart: (1/n) log((n+1)^n/n!) ----
# viewBox 0 0 480 240, n in [1,40] -> [55,450], rate in [0.8, 1.5] -> [210,25]


def Xn(n):
    return 55 + (n - 1) / 39 * 395


def Yr(v):
    return 210 - (v - 0.8) / 0.7 * 185


pts = []
for n in range(1, 41):
    r = (n * math.log(n + 1) - math.lgamma(n + 1)) / n
    pts.append((n, r))
print("rate examples:", [(n, round(r, 3)) for n, r in pts if n in
                         (1, 2, 5, 10, 20, 40)])
print("rate_line:", " ".join(f"{Xn(n):.1f},{Yr(r):.1f}" for n, r in pts))
print("log4_y:", f"{Yr(math.log(4)):.1f}", " loge_y:", f"{Yr(1.0):.1f}")

# ---- 5. transport animation particle keyframes (n=1 toy) ----
targets = [i / 10 for i in range(-9, 10, 2)]           # -0.9 .. 0.9
starts = [2 * math.atanh(v) for v in targets]
print("particles (start x -> target v):")
print("  ", " ".join(f"{s:+.2f}->{v:+.1f}" for s, v in zip(starts, targets)))
# SVG coords: x in [-6,6] -> [40,460] for animation panel
print("  svg:", " ".join(
    f"({40 + (s + 6) / 12 * 420:.0f},{40 + (v + 6) / 12 * 420:.0f})"
    for s, v in zip(starts, targets)))

# ---- 6. lower-slope integral: area under vol - s^n/n! for n=2, vol=4.5 ----
# c_K = (n! vol)^{1/n} = 3;  integral / vol = (n/(n+1)) c_K = 2
n = 2
vol = 4.5
cK = (math.factorial(n) * vol) ** (1 / n)
val = cK - cK ** (n + 1) / (math.factorial(n + 1) * vol)
print(f"n=2 toy: c_K={cK}, integral={val} (want (n/(n+1))c_K="
      f"{n / (n + 1) * cK})")
