"""Figure-data generator for ch09 explainer (growth chart + pentagon).

Run: PYTHONIOENCODING=utf-8 python figgen/ch09_figs.py
Outputs SVG polyline coordinate strings to paste into the explainer.
"""
import math

# ---- Pentagon K5 (viewBox 0 0 220 210) ----
cx, cy, r = 110, 108, 85
pent = [(cx + r*math.cos(math.radians(-90 + 72*i)),
         cy + r*math.sin(math.radians(-90 + 72*i))) for i in range(5)]
print("pentagon:", [(round(x, 1), round(y, 1)) for x, y in pent])

# ---- Growth chart: y = (1/k) ln(bound), x = log10 k in [1,12] ----
# plot area x:60..440, y:20..230; y-value range [-0.5, 8]
YLO, YHI = -0.5, 8.0


def X(e10):
    return 60 + (e10 - 1) * (380 / 11)


def Y(v):
    return 230 - (v - YLO) / (YHI - YLO) * 210


yold = math.log(380) / 5          # exponential seed bound, constant
pts_new, pts_up = [], []
for i in range(2, 49):            # e10 = 0.5 .. 12 step .25
    e10 = i / 4
    k = 10.0 ** e10
    if e10 >= 1:
        pts_new.append((e10, math.log(k) / 3 - math.log(math.log(k))))
        pts_up.append((e10, (math.lgamma(k + 1)
                             + math.log(math.e - 1 / 6)) / k))

# truncate factorial curve at y = 8 (interpolate crossing)
up_t = []
for (e0, y0), (e1, y1) in zip(pts_up, pts_up[1:]):
    if y0 <= YHI:
        up_t.append((e0, y0))
        if y1 > YHI:
            fr = (YHI - y0) / (y1 - y0)
            up_t.append((e0 + fr * (e1 - e0), YHI))
            break
new_line = " ".join(f"{X(e):.1f},{Y(y):.1f}" for e, y in pts_new)
up_line = " ".join(f"{X(e):.1f},{Y(y):.1f}" for e, y in up_t)
print("newline:", new_line)
print("upline:", up_line)
print("oldline:", f"{X(1):.1f},{Y(yold):.1f} {X(12):.1f},{Y(yold):.1f}")

# crossover of new bound (c=1) over old constant
lo, hi = 10.0, 1e8
while hi / lo > 1.0001:
    mid = math.sqrt(lo * hi)
    if math.log(mid) / 3 - math.log(math.log(mid)) > yold:
        hi = mid
    else:
        lo = mid
kx = math.sqrt(lo * hi)
print("crossover:", f"k={kx:.3e}",
      f"X={X(math.log10(kx)):.1f}", f"Y={Y(yold):.1f}")
# y-axis ticks 0,2,4,6,8 and x ticks e10=1..12
print("yticks:", [(v, round(Y(v), 1)) for v in (0, 2, 4, 6, 8)])
print("xticks:", [(e, round(X(e), 1)) for e in (1, 3, 5, 7, 9, 11)])
