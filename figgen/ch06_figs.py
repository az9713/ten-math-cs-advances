"""Chapter 6 (quantum parallel repetition): figure data + checks.

Emits SVG figures to .ignore/ch06_figs/ and prints a validation battery.
Run:  PYTHONIOENCODING=utf-8 python figgen/ch06_figs.py
"""
import itertools
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch06_figlib import Ax, svg_wrap, sample, write_fig, fmt  # noqa: E402

import numpy as np  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), ".ignore", "ch06_figs")

RED = "#c0392b"
BLUE = "#2e6da4"
GREEN = "#2a7d2a"
GRAY = "#888"
DARK = "#333"
MUT = "#555"
PUR = "#7d3c98"
ORA = "#b06a00"

CHECKS = []
rng = np.random.default_rng(20260802)


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(f"CHECK {name}: {'OK' if ok else 'FAIL'}  {detail}")


# ===================================================================== CHSH
def chsh_win_classical():
    """Max over deterministic strategies of P(a xor b = x and y)."""
    best = 0.0
    for fa in itertools.product((0, 1), repeat=2):
        for fb in itertools.product((0, 1), repeat=2):
            w = sum(1 for x in (0, 1) for y in (0, 1)
                    if fa[x] ^ fb[y] == (x & y)) / 4.0
            best = max(best, w)
    return best


def meas_vec(theta, out):
    if out == 0:
        return np.array([math.cos(theta), math.sin(theta)])
    return np.array([-math.sin(theta), math.cos(theta)])


def chsh_win_quantum2(tA, tB):
    """Win prob of the rotation strategy on |Phi+> with angle maps."""
    w = 0.0
    for x in (0, 1):
        for y in (0, 1):
            for a in (0, 1):
                for b in (0, 1):
                    if a ^ b != (x & y):
                        continue
                    amp = meas_vec(tA[x], a) @ meas_vec(tB[y], b)
                    w += 0.25 * (amp * amp) / 2.0
    return w


def obs(theta):
    v0 = meas_vec(theta, 0)
    v1 = meas_vec(theta, 1)
    return np.outer(v0, v0) - np.outer(v1, v1)


def bell_norm(a0, a1, b0, b1):
    B = (np.kron(a0, b0) + np.kron(a0, b1)
         + np.kron(a1, b0) - np.kron(a1, b1))
    return np.linalg.norm(B, 2)


def checks_chsh():
    best = chsh_win_classical()
    check("c01 CHSH classical value = 3/4", abs(best - 0.75) < 1e-12,
          f"max over 16 deterministic pairs = {best}")
    tA = (0.0, math.pi / 4)
    tB = (math.pi / 8, -math.pi / 8)
    w = chsh_win_quantum2(tA, tB)
    tgt = math.cos(math.pi / 8) ** 2
    check("c02 CHSH rotation strategy = cos^2(pi/8)",
          abs(w - tgt) < 1e-12, f"w={w:.10f} target={tgt:.10f}")
    nrm = bell_norm(obs(0), obs(math.pi / 4),
                    obs(math.pi / 8), obs(-math.pi / 8))
    check("c03a Bell operator norm = 2 sqrt 2",
          abs(nrm - 2 * math.sqrt(2)) < 1e-12, f"|B|={nrm:.10f}")
    worst = 0.0
    for _ in range(300):
        th = rng.uniform(0, 2 * math.pi, 4)
        worst = max(worst, bell_norm(obs(th[0]), obs(th[1]),
                                     obs(th[2]), obs(th[3])))
    ok4 = True
    for _ in range(60):
        mats = []
        for _k in range(4):
            H = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
            H = H + H.conj().T
            w_, V = np.linalg.eigh(H)
            s = np.sign(rng.normal(size=4))
            s[s == 0] = 1
            mats.append((V * s) @ V.conj().T)
        B = (np.kron(mats[0], mats[2]) + np.kron(mats[0], mats[3])
             + np.kron(mats[1], mats[2]) - np.kron(mats[1], mats[3]))
        if np.linalg.norm(B, 2) > 2 * math.sqrt(2) + 1e-9:
            ok4 = False
    check("c03b Tsirelson: random observables never beat 2 sqrt 2",
          worst <= 2 * math.sqrt(2) + 1e-9 and ok4,
          f"max over 300 qubit + 60 dim-4 samples = {worst:.6f}")


# ==================================================================== Feige
# Answers 0..3 encode (c,d): c = idx//2 + 1 (player pointer), d = idx%2.
def feige_win(x, y, a, b):
    if a != b:
        return 0
    c, d = a // 2 + 1, a % 2
    return 1 if (c == 1 and d == x) or (c == 2 and d == y) else 0


def checks_feige():
    best = 0.0
    for fa in itertools.product(range(4), repeat=2):
        for fb in itertools.product(range(4), repeat=2):
            w = sum(feige_win(x, y, fa[x], fb[y])
                    for x in (0, 1) for y in (0, 1)) / 4.0
            best = max(best, w)
    check("c04 Feige game value = 1/2", abs(best - 0.5) < 1e-12,
          f"max over 256 deterministic pairs = {best}")

    # cross strategy on G x G: a=( (1,x1), (2,x1) ), b=( (1,y2), (2,y2) )
    tot = 0.0
    for x1 in (0, 1):
        for x2 in (0, 1):
            for y1 in (0, 1):
                for y2 in (0, 1):
                    a1, a2 = 0 + x1, 2 + x1
                    b1, b2 = 0 + y2, 2 + y2
                    tot += (feige_win(x1, y1, a1, b1)
                            * feige_win(x2, y2, a2, b2)) / 16.0
    check("c05 Feige cross strategy wins both = 1/2",
          abs(tot - 0.5) < 1e-12, f"P(win both)={tot}")

    # exact value of the two-fold repetition by best response
    T = np.zeros((16, 16, 4, 4), dtype=np.int8)
    for ap in range(16):
        a1, a2 = ap // 4, ap % 4
        for bp in range(16):
            b1, b2 = bp // 4, bp % 4
            for xp in range(4):
                x1, x2 = xp // 2, xp % 2
                for yp in range(4):
                    y1, y2 = yp // 2, yp % 2
                    T[ap, bp, xp, yp] = (feige_win(x1, y1, a1, b1)
                                         * feige_win(x2, y2, a2, b2))
    idx = np.arange(16 ** 4)
    sA = np.stack([(idx // (16 ** k)) % 16 for k in range(4)],
                  axis=1)  # (65536, 4): answer pair for each x-pair
    # value(s) = 1/16 sum_yp max_bp sum_xp T[sA[s,xp], bp, xp, yp]
    acc = np.zeros((16 ** 4, 16, 4))
    for xp in range(4):
        acc += T[sA[:, xp], :, xp, :]
    val = acc.max(axis=1).sum(axis=1) / 16.0
    v2 = float(val.max())
    check("c06 Feige two-fold value = 1/2 (not 1/4)",
          abs(v2 - 0.5) < 1e-12,
          f"exact max over all 65536 Alice strategies x best "
          f"Bob response = {v2}")


# ============================================== classical information tools
def dtv(P, Q):
    return 0.5 * float(np.abs(P - Q).sum())


def dkl(P, Q):
    m = P > 0
    return float((P[m] * np.log(P[m] / Q[m])).sum())


def rand_dist(k):
    v = rng.random(k) + 1e-3
    return v / v.sum()


def checks_info():
    ok = True
    worst = 0.0
    for _ in range(500):
        k = int(rng.integers(2, 7))
        P, Q = rand_dist(k), rand_dist(k)
        gap = math.sqrt(0.5 * dkl(P, Q)) - dtv(P, Q)
        worst = min(worst, gap) if _ else gap
        if gap < -1e-12:
            ok = False
    check("c07 Pinsker holds on 500 random pairs", ok)
    t = 1e-3
    P = np.array([0.5 + t, 0.5 - t])
    Q = np.array([0.5, 0.5])
    ratio = dtv(P, Q) / math.sqrt(0.5 * dkl(P, Q))
    check("c07b Pinsker near-tight for close binary pair",
          abs(ratio - 1) < 1e-4, f"ratio={ratio:.8f}")

    # correlated sampling: exact disagreement + Monte Carlo
    P = np.array([0.5, 0.3, 0.2])
    Q = np.array([0.3, 0.3, 0.4])
    d = dtv(P, Q)
    exact = float(np.abs(P - Q).sum()) / float(np.maximum(P, Q).sum())
    check("c08a rejection-coupling disagreement = 2dTV/(1+dTV)",
          abs(exact - 2 * d / (1 + d)) < 1e-12,
          f"exact={exact:.6f} 2dTV={2*d:.6f}")
    ra = random.Random(7)
    bad = runs = 0
    for _ in range(200000):
        while True:
            z = ra.randrange(3)
            u = ra.random()
            accP, accQ = u <= P[z], u <= Q[z]
            if accP or accQ:
                break
        runs += 1
        if not (accP and accQ):
            bad += 1
    mc = bad / runs
    check("c08b Monte Carlo matches and is <= 2dTV",
          abs(mc - exact) < 0.01 and mc <= 2 * d,
          f"mc={mc:.4f} exact={exact:.4f}")

    ok = True
    for _ in range(200):
        a = rng.random(4) + 1e-3
        b = rng.random(4) + 1e-3
        lhs = a.sum() * math.log(a.sum() / b.sum())
        rhs = float((a * np.log(a / b)).sum())
        if lhs > rhs + 1e-12:
            ok = False
    check("c13 log-sum inequality on 200 random tuples", ok)


# ===================================== resolvent purification (exact Gram)
def gram_s(s, t):
    if abs(s - t) < 1e-12 * max(s, t, 1e-30):
        return s
    return s * t * math.log(t / s) / (t - s)


def eig_pos(F, tol=1e-12):
    w, V = np.linalg.eigh(F)
    return [(w[k], V[:, k]) for k in range(len(w)) if w[k] > tol]


def gram_M(A, B):
    """M(A,B) with <Gamma(A)u, Gamma(B)w> = u^dag M(A,B) w."""
    MA = np.zeros(A.shape, dtype=complex)
    for s, a in eig_pos(A):
        Pa = np.outer(a, a.conj())
        for t, b in eig_pos(B):
            Pb = np.outer(b, b.conj())
            MA += gram_s(s, t) * (Pa @ Pb)
    return MA


def H1op(F):
    w, V = np.linalg.eigh(F)
    h = np.where(w > 1e-15, -w * np.clip(w, 1e-300, None)
                 * 0, 0.0)
    h = np.array([-x * math.log(x) if x > 1e-15 else 0.0 for x in w])
    return (V * h) @ V.conj().T


def rand_contraction(d, singular=False):
    H = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    H = H @ H.conj().T
    w, V = np.linalg.eigh(H)
    w = w / w.max()
    if singular:
        w[0] = 0.0
    return (V * w) @ V.conj().T


def checks_resolvent():
    import mpmath as mp
    for sg in (0.3, 1.0, 2.7):
        val = mp.quad(lambda u: (sg / (sg + u)) ** 2, [0, mp.inf])
        check(f"c09 int (s/(s+u))^2 du = s  (s={sg})",
              abs(val - sg) < 1e-12)
    s, t = 0.4, 1.3
    num = mp.quad(lambda u: (s / (s + u)) * (t / (t + u)), [0, mp.inf])
    check("c09b Gram formula st log(t/s)/(t-s) matches quadrature",
          abs(num - gram_s(s, t)) < 1e-10,
          f"quad={float(num):.8f} formula={gram_s(s, t):.8f}")
    ok = True
    for d in (2, 3):
        for sing in (False, True):
            F = rand_contraction(d, sing)
            if np.linalg.norm(gram_M(F, F) - F) > 1e-8:
                ok = False
    check("c10 Born property Gamma(F)^dag Gamma(F) = F "
          "(2x2 and 3x3, incl. singular)", ok)

    okA = True
    wmin = 1.0
    for trial in range(40):
        d = 2 if trial % 2 == 0 else 3
        kk = int(rng.integers(2, 5))
        Fs = [rand_contraction(d, singular=(trial % 5 == 0 and j == 0))
              for j in range(kk)]
        ps = rand_dist(kk)
        Fbar = sum(p * F for p, F in zip(ps, Fs))
        lhs = np.zeros((d, d), dtype=complex)
        for p, F in zip(ps, Fs):
            lhs += p * (gram_M(F, F) - gram_M(F, Fbar)
                        - gram_M(Fbar, F) + gram_M(Fbar, Fbar))
        rhs = H1op(Fbar) - sum(p * H1op(F) for p, F in zip(ps, Fs))
        gap = np.linalg.eigvalsh(rhs - lhs).min()
        wmin = min(wmin, gap)
        if gap < -1e-8:
            okA = False
    check("c11 Lemma 4.3 entropy control on 40 random ensembles",
          okA, f"min eig(RHS-LHS) over trials = {wmin:.2e}")


def checks_martingale44():
    # toy filtration: omega = 3 bits; F_j = E[F_omega | first j bits]
    dA = dB = 2
    Fo = [rand_contraction(dA) for _ in range(8)]
    K = rand_contraction(dB)
    psi = rng.normal(size=dA * dB) + 1j * rng.normal(size=dA * dB)
    psi = psi / np.linalg.norm(psi)
    N = 3

    def cond_F(bits):
        j = len(bits)
        sel = [o for o in range(8)
               if all((o >> (2 - k)) & 1 == bits[k] for k in range(j))]
        return sum(Fo[o] for o in sel) / len(sel)

    def sq_incr(A, B):
        DG = (gram_M(A, A) - gram_M(A, B)
              - gram_M(B, A) + gram_M(B, B))
        op = np.kron(DG, K)
        return float(np.real(psi.conj() @ (op @ psi)))

    tot = 0.0
    for j in range(N):
        for bits in itertools.product((0, 1), repeat=j + 1):
            prev = cond_F(bits[:-1])
            cur = cond_F(bits)
            tot += 2.0 ** (-(j + 1)) * sq_incr(cur, prev)
    F0 = cond_F(())
    p0 = float(np.real(psi.conj()
                       @ (np.kron(F0, K) @ psi)))
    Ksq = np.zeros((dB, dB), dtype=complex)
    wK, VK = np.linalg.eigh(K)
    Ksq = (VK * np.sqrt(np.clip(wK, 0, None))) @ VK.conj().T
    Psi = psi.reshape(dA, dB)
    rhoK = (Psi @ Ksq).conj().T
    rhoK = (Psi @ Ksq) @ (Psi @ Ksq).conj().T
    trH = float(np.real(np.trace(rhoK @ H1op(F0))))
    h1p0 = -p0 * math.log(p0)
    check("c12 Lemma 4.4 toy martingale: sum <= Tr(rho H1(F0)) <= H1(p0)",
          tot <= trH + 1e-9 and trH <= h1p0 + 1e-9,
          f"sum={tot:.6f} trH={trH:.6f} H1(p0)={h1p0:.6f}")
    check("c12b uniform-step bound sum/N <= H1(p0)/N",
          tot / N <= h1p0 / N + 1e-12)


# =================================== Lemma 4.1 reverse distribution (m=4)
def checks_reverse():
    m = 4
    M = list(range(m))
    fwd = {}
    for i in M:
        rest = [c for c in M if c != i]
        for mask in range(2 ** (m - 1)):
            LX = [rest[k] for k in range(m - 1) if (mask >> k) & 1]
            LY = [rest[k] for k in range(m - 1) if not (mask >> k) & 1]
            for pX in itertools.permutations(LX):
                for pY in itertools.permutations(LY):
                    for kX in range(len(LX) + 1):
                        for kY in range(len(LY) + 1):
                            pr = (1.0 / m) * 2.0 ** (-(m - 1))
                            pr /= math.factorial(len(LX))
                            pr /= math.factorial(len(LY))
                            pr /= (len(LX) + 1) * (len(LY) + 1)
                            key = (i, tuple(sorted(LX)), pX, pY, kX, kY)
                            fwd[key] = fwd.get(key, 0.0) + pr
    rev = {}
    for mask in range(2 ** m):
        LX = [M[k] for k in range(m) if (mask >> k) & 1]
        Lp = [M[k] for k in range(m) if not (mask >> k) & 1]
        N = len(Lp)
        if N == 0:
            continue
        ppart = 2.0 ** (-m) * 2 * N / m
        for pYfull in itertools.permutations(Lp):
            for kY in range(N):
                i = pYfull[kY]
                pY = tuple(c for c in pYfull[:kY]) + tuple(
                    c for c in pYfull[kY + 1:])
                for pX in itertools.permutations(LX):
                    for kX in range(len(LX) + 1):
                        pr = ppart / math.factorial(N) / N
                        pr /= math.factorial(len(LX))
                        pr /= (len(LX) + 1)
                        key = (i, tuple(sorted(LX)), pX, pY, kX, kY)
                        rev[key] = rev.get(key, 0.0) + pr
    keys = set(fwd) | set(rev)
    worst = max(abs(fwd.get(k, 0) - rev.get(k, 0)) for k in keys)
    check("c14 Lemma 4.1 reverse = forward distribution (m=4, exact)",
          worst < 1e-14,
          f"{len(keys)} outcomes, max prob deviation {worst:.1e}")
    # wait: reverse pY should keep the order structure; kY cut uniform
    tot_f = sum(fwd.values())
    tot_r = sum(rev.values())
    check("c14b both experiments have total mass 1",
          abs(tot_f - 1) < 1e-12 and abs(tot_r - 1) < 1e-12,
          f"fwd={tot_f:.12f} rev={tot_r:.12f}")
    # size-bias cancellation sum
    ssum = sum(2.0 ** (-m) for mask in range(2 ** m)
               if mask != 2 ** m - 1)
    check("c16 sum of fair-partition probs (N>0) = 1 - 2^-m <= 1",
          abs(ssum - (1 - 2.0 ** (-m))) < 1e-14, f"sum={ssum}")


# =============================================== greedy conditioning (toy)
def checks_greedy():
    n = 12
    ra = np.random.default_rng(5)
    comp = []
    for _ in range(3):
        comp.append(ra.uniform(0.6, 0.98, size=n))
    mixw = np.array([0.5, 0.3, 0.2])
    outs = np.array(list(itertools.product((0, 1), repeat=n)))
    probs = np.zeros(len(outs))
    for w, pv in zip(mixw, comp):
        probs += w * np.prod(np.where(outs == 1, pv, 1 - pv), axis=1)
    theta = probs[np.all(outs == 1, axis=1)].sum()
    delta = 0.2
    D = []
    while True:
        inD = np.all(outs[:, D] == 1, axis=1) if D else np.ones(
            len(outs), bool)
        pD = probs[inD].sum()
        rest = [i for i in range(n) if i not in D]
        cond = [(probs[inD & (outs[:, i] == 1)].sum() / pD, i)
                for i in rest]
        avg = sum(c for c, _ in cond) / len(cond)
        if 1 - avg <= delta:
            break
        cbad = min(cond)[1]
        D.append(cbad)
    bound = math.log(1 / theta) / delta
    check("c22 greedy conditioning toy: |D|, mass, avg all as in L3.1",
          len(D) <= bound and pD >= theta and avg >= 1 - delta,
          f"|D|={len(D)} bound={bound:.2f} pD={pD:.4f} "
          f"theta={theta:.4f} avg={avg:.4f}")


# ===================================================== embezzlement decay
def embezzle_error(b, lam=0.3):
    j = np.arange(1, b + 1, dtype=float)
    c2 = (1 / j) / (1 / j).sum()
    tgt = np.sort(np.concatenate([lam * c2, (1 - lam) * c2]))[::-1]
    src = np.concatenate([c2, np.zeros(b)])
    src = np.sort(src)[::-1]
    F = float(np.sqrt(src * tgt).sum())
    return math.sqrt(max(0.0, 2 - 2 * F))


def checks_embezzle():
    errs = {b: embezzle_error(2 ** k)
            for k, b in [(4, 2 ** 4), (8, 2 ** 8), (12, 2 ** 12),
                         (16, 2 ** 16)]}
    es = list(errs.values())
    dec = all(es[k + 1] < es[k] for k in range(len(es) - 1))
    prod = [e * math.log(b) for b, e in errs.items()]
    check("c17 van Dam-Hayden embezzlement error falls like 1/log b",
          dec and prod[-1] < 2.5 * prod[0] and prod[-1] > 0.2 * prod[0],
          "err(b): " + ", ".join(f"2^{k}:{errs[2**k]:.4f}"
                                 for k in (4, 8, 12, 16)))


# ================================================= constants and Sec. 3.5
def constants_chain():
    Kqs = 1.0
    Ustar = Kqs * (1 + 2 ** (1 / 6)) + 2
    Bqs = ((5 + 2 * Ustar) * math.sqrt(1.5)
           + 2 * Kqs * 32 ** (1 / 12) + 2 * math.sqrt(8))
    cqs = 1.0 / (8 * (4 * Bqs) ** 12)
    return Kqs, Ustar, Bqs, cqs


def checks_constants():
    Kqs, Ustar, Bqs, cqs = constants_chain()
    check("c18 constants (K=1): U*, B, c computed",
          Bqs > 1 and cqs > 0,
          f"U*={Ustar:.4f} B={Bqs:.4f} c={cqs:.3e}")
    ok = True
    for eps in (0.1, 0.5, 0.9):
        for ell in (math.log(4.0), 5.0):
            gam = cqs * eps ** 13 / (eps + ell)
            delta = eps / 4
            if not gam <= delta / 2 + 1e-18:
                ok = False
            etab = 2 * gam * (1 + 4 * ell / eps)
            if not etab <= (eps / (4 * Bqs)) ** 12 + 1e-18:
                ok = False
            alpha = (eps / (16 * Kqs)) ** 12
            q = 1 - eps / 4
            win = (q - Bqs * etab ** (1 / 12)
                   - 2 * Kqs * alpha ** (1 / 12))
            if not win > 1 - eps:
                ok = False
    check("c19 Sec 3.5 chain: gamma<=delta/2, eta<=(eps/4B)^12, "
          "win > 1-eps on grid", ok)
    ok = True
    for eta in np.linspace(1e-9, 1.0, 200):
        kap = math.sqrt(1.5 * eta)
        lhs = ((5 + 2 * Ustar) * kap
               + 2 * Kqs * (32 * eta) ** (1 / 12)
               + 2 * math.sqrt(8 * eta))
        if lhs > Bqs * eta ** (1 / 12) + 1e-12:
            ok = False
    check("c20 A.8 assembly: kappa/sqrt-eta terms fold into B eta^(1/12)",
          ok)


# ============================= Lemma A.1 refinement operators (sqrt Gamma)
def checks_refine():
    d = 3
    F = rand_contraction(d, singular=True)
    A1 = rand_contraction(d)
    A2 = np.eye(d) - A1
    w, V = np.linalg.eigh(F)
    sq = (V * np.sqrt(np.clip(w, 0, None))) @ V.conj().T
    Fa = [sq @ A1 @ sq, sq @ A2 @ sq]
    isq = (V * np.array([1 / math.sqrt(x) if x > 1e-12 else 0.0
                         for x in w])) @ V.conj().T
    # Gamma = sqrt purification on same space; U_F = sq . isq = supp proj
    U = sq @ isq
    P0 = U @ U.conj().T
    Ms = []
    for k, Fk in enumerate(Fa):
        Mk = U @ (isq @ Fk @ isq) @ U.conj().T
        if k == 0:
            Mk = Mk + (np.eye(d) - P0)
        Ms.append(Mk)
    okpos = all(np.linalg.eigvalsh(M).min() > -1e-10 for M in Ms)
    oksum = np.linalg.norm(sum(Ms) - np.eye(d)) < 1e-9
    okrep = all(np.linalg.norm(sq @ M @ sq - Fk) < 1e-9
                for M, Fk in zip(Ms, Fa))
    check("c21 Lemma A.1 refinement POVM (singular F): PSD, sums to I, "
          "Gamma^dag M Gamma = F^a", okpos and oksum and okrep)


# ================== end-to-end toy: n=2 entangled strategy, D={1}, i=2
def kron_meas(Us, base):
    return Us @ base


def checks_end_to_end():
    # Repeated game: CHSH x CHSH. Alice/Bob spaces C^2 x C^2.
    # Correlated strategy: Alice entangles her two qubits with a unitary
    # before product measurement; likewise Bob.  State: Phi+ x Phi+.
    tA = (0.0, math.pi / 4)
    tB = (math.pi / 8, -math.pi / 8)

    ra = np.random.default_rng(11)
    QA = np.linalg.qr(ra.normal(size=(4, 4)))[0]
    QB = np.linalg.qr(ra.normal(size=(4, 4)))[0]

    def povmA(x1, x2):
        E = {}
        for a1 in (0, 1):
            for a2 in (0, 1):
                v = QA.T @ np.kron(
                    meas_vec(tA[x1], a1), meas_vec(tA[x2], a2))
                E[(a1, a2)] = np.outer(v, v)
        return E

    def povmB(y1, y2):
        E = {}
        for b1 in (0, 1):
            for b2 in (0, 1):
                v = QB.T @ np.kron(
                    meas_vec(tB[y1], b1), meas_vec(tB[y2], b2))
                E[(b1, b2)] = np.outer(v, v)
        return E

    phip = np.array([1, 0, 0, 1]) / math.sqrt(2)
    # full state on A1 A2 (x) B1 B2, ordered (A1A2),(B1B2):
    # psi = sum_ij sum_kl (Phi+)_{ik} (Phi+)_{jl} |ij>_A |kl>_B
    psi = np.zeros(16)
    for i1 in (0, 1):
        for k1 in (0, 1):
            for j2 in (0, 1):
                for l2 in (0, 1):
                    aA = i1 * 2 + j2
                    aB = k1 * 2 + l2
                    amp = (phip[i1 * 2 + k1] * phip[j2 * 2 + l2]) * 2
                    amp = ((1 / math.sqrt(2)) if i1 == k1 else 0) * (
                        (1 / math.sqrt(2)) if j2 == l2 else 0)
                    psi[aA * 4 + aB] = amp
    check("c30a toy state is normalized",
          abs(np.linalg.norm(psi) - 1) < 1e-12)

    mu = {(x, y): 0.25 for x in (0, 1) for y in (0, 1)}
    # joint distribution P(x1 y1 x2 y2 a1 a2 b1 b2)
    P = {}
    for x1, y1 in mu:
        for x2, y2 in mu:
            EA = povmA(x1, x2)
            EB = povmB(y1, y2)
            for (a1, a2), Ea in EA.items():
                for (b1, b2), Fb in EB.items():
                    pr = float(psi @ (np.kron(Ea, Fb) @ psi))
                    P[(x1, y1, x2, y2, a1, a2, b1, b2)] = (
                        mu[(x1, y1)] * mu[(x2, y2)] * pr)
    tot = sum(P.values())
    check("c30b toy joint distribution sums to 1", abs(tot - 1) < 1e-10)

    def win(x, y, a, b):
        return 1 if (a ^ b) == (x & y) else 0

    # condition on winning coordinate 1 (D={1}, i=2, CX=CY={1})
    p = sum(v for k, v in P.items()
            if win(k[0], k[1], k[4], k[6]))
    # histories r = (x1, y1, a1, b1); live questions (x2, y2)
    oks4 = True
    oks32 = True
    for x1 in (0, 1):
        for y1 in (0, 1):
            for a1 in (0, 1):
                for b1 in (0, 1):
                    for x2 in (0, 1):
                        for y2 in (0, 1):
                            EA = povmA(x1, x2)
                            EB = povmB(y1, y2)
                            H = sum(EA[(a1, t)] for t in (0, 1))
                            K = sum(EB[(b1, t)] for t in (0, 1))
                            HK = np.kron(H, K)
                            prx = float(psi @ (HK @ psi))
                            # direct: P(Z=z | T, X2, Y2)
                            num = sum(
                                P[(x1, y1, x2, y2, a1, t, b1, u)]
                                for t in (0, 1) for u in (0, 1))
                            den = mu[(x1, y1)] * mu[(x2, y2)]
                            if abs(prx - num / den) > 1e-10:
                                oks4 = False
                            if prx < 1e-14:
                                continue
                            # Lemma 3.2 via sqrt purification + A.1
                            wH, VH = np.linalg.eigh(H)
                            sqH = (VH * np.sqrt(np.clip(wH, 0, None))
                                   ) @ VH.conj().T
                            wK, VK = np.linalg.eigh(K)
                            sqK = (VK * np.sqrt(np.clip(wK, 0, None))
                                   ) @ VK.conj().T
                            phiv = np.kron(sqH, sqK) @ psi
                            Psi = phiv / np.linalg.norm(phiv)
                            for a in (0, 1):
                                for b in (0, 1):
                                    Ha = EA[(a1, a)]
                                    Kb = EB[(b1, b)]
                                    lhs = float(
                                        psi @ (np.kron(Ha, Kb) @ psi)
                                    ) / prx
                                    num2 = P[(x1, y1, x2, y2,
                                              a1, a, b1, b)]
                                    rhs = num2 / den / (num / den)
                                    if abs(lhs - rhs) > 1e-9:
                                        oks32 = False
    check("c31 eq.(4): <psi|H x K|psi> = P(Z=z | T, X_i, Y_i) "
          "on all 64 branches", oks4)
    check("c32 Lemma 3.2 ratio identity = Q(A_i,B_i | r,x,y) "
          "on all branches", oks32)

    # ideal success probability = q  (average over posterior)
    num_q = sum(v for k, v in P.items()
                if win(k[0], k[1], k[4], k[6])
                and win(k[2], k[3], k[5], k[7]))
    q = num_q / p
    s = 0.0
    for k, v in P.items():
        if not win(k[0], k[1], k[4], k[6]):
            continue
        s += v * win(k[2], k[3], k[5], k[7])
    check("c33 E_Q[ideal win] = P(W2|W1) = q on the toy",
          abs(s / p - q) < 1e-12, f"q={q:.6f} p={p:.6f}")
    return q, p


# ========================================================== figure helpers
def arrow_def(pid, color=DARK):
    return (f'<marker id="{pid}" viewBox="0 0 10 10" refX="9" refY="5" '
            f'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>')


def box(x, y, w, h, fill="#fff", stroke=DARK, rx=7, sw=1.4, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')


def txt(x, y, s, fs=13, anchor="middle", color=DARK, cls=None, w=None):
    c = f' class="{cls}"' if cls else ""
    wt = f' font-weight="{w}"' if w else ""
    return (f'<text x="{x}" y="{y}" font-size="{fs}" '
            f'text-anchor="{anchor}" fill="{color}"{c}{wt}>{s}</text>')


def line(x1, y1, x2, y2, color=DARK, w=1.4, dash=None, marker=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = f' marker-end="url(#{marker})"' if marker else ""
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
            f'stroke="{color}" stroke-width="{w}"{d}{m}/>')


def sb(s):
    return (f'<tspan baseline-shift="sub" font-size="70%">{s}</tspan>')


def sp(s):
    return (f'<tspan baseline-shift="super" font-size="70%">{s}'
            '</tspan>')


# ------------------------------------------------------------ fig_game
def fig_game():
    b = []
    b.append(f'<defs>{arrow_def("c6gA")}'
             f'{arrow_def("c6gR", RED)}</defs>')
    b.append(box(200, 12, 120, 44, fill="#f3f1ec"))
    b.append(txt(260, 30, "referee", 13, w="bold"))
    b.append(txt(260, 47, "samples (x, y) ~ μ", 11, color=MUT))
    b.append(box(30, 120, 130, 56, fill="#eef4fa", stroke=BLUE))
    b.append(txt(95, 141, "Alice", 13, w="bold", color=BLUE))
    b.append(txt(95, 158, "measures with x", 11, color=MUT))
    b.append(box(360, 120, 130, 56, fill="#eef4fa", stroke=BLUE))
    b.append(txt(425, 141, "Bob", 13, w="bold", color=BLUE))
    b.append(txt(425, 158, "measures with y", 11, color=MUT))
    b.append(line(215, 56, 75, 118, marker="c6gA"))
    b.append(txt(112, 76, "question x", 11.5, color=DARK))
    b.append(line(305, 56, 445, 118, marker="c6gA"))
    b.append(txt(408, 76, "question y", 11.5))
    b.append(line(125, 118, 235, 60, color=RED, marker="c6gR"))
    b.append(txt(196, 106, "answer a", 11.5, color=RED,
                 anchor="start"))
    b.append(line(395, 118, 285, 60, color=RED, marker="c6gR"))
    b.append(txt(324, 106, "answer b", 11.5, color=RED,
                 anchor="end"))
    # entangled state
    b.append('<path d="M162,148 C 220,132 300,164 358,148" fill="none" '
             f'stroke="{PUR}" stroke-width="2.2" stroke-dasharray="1 5" '
             'stroke-linecap="round"/>')
    b.append(txt(260, 136, "shared entangled state ψ", 11.5,
                 color=PUR))
    b.append(txt(260, 197, "win iff V(x, y, a, b) = 1 "
                 "— no communication after the questions", 12))
    return svg_wrap("\n".join(b), 520, 210,
                    "two-player one-round entangled game")


# ------------------------------------------------------------ fig_born
def fig_born():
    cx, cy, R = 120, 108, 78
    th = math.pi / 5
    b = []
    b.append(f'<defs>{arrow_def("c6bA", BLUE)}'
             f'{arrow_def("c6bR", RED)}</defs>')
    b.append(f'<circle cx="{cx}" cy="{cy}" r="{R}" fill="none" '
             f'stroke="{GRAY}" stroke-width="1.2"/>')
    b.append(line(cx - R - 8, cy, cx + R + 8, cy, color="#bbb", w=1))
    b.append(line(cx, cy + R + 8, cx, cy - R - 8, color="#bbb", w=1))
    # state vector along e1
    b.append(line(cx, cy, cx + R, cy, color=RED, w=2.4, marker="c6bR"))
    b.append(txt(cx + R - 6, cy + 20, "state ψ", 12, color=RED,
                 anchor="end"))
    vx = cx + R * math.cos(th)
    vy = cy - R * math.sin(th)
    wx = cx - R * math.sin(th) * 0.72
    wy = cy - R * math.cos(th) * 0.72
    b.append(line(cx, cy, vx, vy, color=BLUE, w=2.2, marker="c6bA"))
    b.append(line(cx, cy, wx, wy, color=BLUE, w=2.2, marker="c6bA",
                  dash="5 3"))
    b.append(txt(vx + 8, vy + 2, "v₀", 12.5, color=BLUE,
                 anchor="start", cls="v"))
    b.append(txt(wx - 2, wy - 8, "v₁", 12.5, color=BLUE, cls="v"))
    # angle arc
    b.append(f'<path d="M {cx+34},{cy} A 34 34 0 0 0 '
             f'{cx+34*math.cos(th)},{cy-34*math.sin(th)}" fill="none" '
             f'stroke="{DARK}" stroke-width="1.2"/>')
    b.append(txt(cx + 46, cy - 12, "θ", 13))
    # projection dashed
    px = cx + R * math.cos(th) ** 2
    py = cy - R * math.cos(th) * math.sin(th)
    b.append(line(cx + R, cy, px, py, color=GRAY, w=1.2, dash="3 3"))
    x0 = 268
    b.append(txt(x0, 52, "outcome 0:  probability "
                 "|⟨v₀, ψ⟩|² = cos²θ",
                 12.5, anchor="start"))
    b.append(txt(x0, 78, "outcome 1:  probability "
                 "|⟨v₁, ψ⟩|² = sin²θ",
                 12.5, anchor="start"))
    b.append(txt(x0, 116, "a measurement is a choice of", 12,
                 anchor="start", color=MUT))
    b.append(txt(x0, 134, "orthonormal directions; the state", 12,
                 anchor="start", color=MUT))
    b.append(txt(x0, 152, "picks each with the squared", 12,
                 anchor="start", color=MUT))
    b.append(txt(x0, 170, "overlap (Born rule).", 12,
                 anchor="start", color=MUT))
    return svg_wrap("\n".join(b), 500, 210,
                    "Born rule for a rotated basis measurement")


# -------------------------------------------------------- fig_entangle
def fig_entangle():
    b = []

    def grid(x0, y0, cells, title, sub):
        out = [txt(x0 + 80, y0 - 26, title, 12.5, w="bold"),
               txt(x0 + 80, y0 - 10, sub, 11, color=MUT)]
        lab = ["00", "01", "10", "11"]
        for i in range(2):
            for j in range(2):
                v = cells[i][j]
                shade = 0.15 + 0.85 * v * 2
                col = (f'rgb({int(255-120*min(1,shade))},'
                       f'{int(255-60*min(1,shade))},255)')
                out.append(box(x0 + j * 80, y0 + i * 44, 78, 42,
                               fill=col, stroke="#999", rx=3, sw=1))
                out.append(txt(x0 + j * 80 + 39, y0 + i * 44 + 26,
                               fmt(v), 13))
        out.append(txt(x0 - 8, y0 + 26, "a=0", 11, anchor="end",
                       color=MUT))
        out.append(txt(x0 - 8, y0 + 70, "a=1", 11, anchor="end",
                       color=MUT))
        out.append(txt(x0 + 39, y0 + 100, "b=0", 11, color=MUT))
        out.append(txt(x0 + 119, y0 + 100, "b=1", 11, color=MUT))
        return out

    b += grid(56, 60, [[0.25, 0.25], [0.25, 0.25]],
              "independent coins",
              "P(a,b) = P(a)P(b): no correlation")
    b += grid(320, 60, [[0.5, 0.0], [0.0, 0.5]],
              "Φ⁺, both measure angle θ",
              "outcomes always equal — for every θ")
    b.append(txt(260, 196, "entanglement is stronger-than-shared-coin "
                 "correlation: same-basis outcomes agree with "
                 "probability 1", 11.5))
    return svg_wrap("\n".join(b), 520, 210,
                    "product versus entangled correlations")


# ------------------------------------------------------------ fig_chsh
def fig_chsh():
    b = []
    # left: rules table
    b.append(txt(120, 24, "CHSH: win iff a ⊕ b = x ∧ y", 13,
                 w="bold"))
    rows = [("x", "y", "need"), ("0", "0", "a = b"), ("0", "1", "a = b"),
            ("1", "0", "a = b"), ("1", "1", "a ≠ b")]
    for r, (c1, c2, c3) in enumerate(rows):
        y = 44 + r * 26
        w0 = "bold" if r == 0 else None
        b.append(txt(60, y, c1, 12.5, w=w0))
        b.append(txt(110, y, c2, 12.5, w=w0))
        b.append(txt(175, y, c3, 12.5, w=w0))
        if r == 0:
            b.append(line(35, y + 8, 225, y + 8, color="#999", w=1))
    b.append(box(30, 162, 200, 30, fill="#fdf7e3",
                 stroke="#c9b45a", rx=5, sw=1))
    b.append(txt(130, 182, "classical max: 3 of 4 ⇒ 3/4", 12))
    # right: measurement angles on circle
    cx, cy, R = 380, 108, 74
    b.append(f'<circle cx="{cx}" cy="{cy}" r="{R}" fill="none" '
             f'stroke="{GRAY}" stroke-width="1.2"/>')
    angs = [(0.0, "A₀ = 0°", BLUE, 16),
            (math.pi / 4, "A₁ = 45°", BLUE, 16),
            (math.pi / 8, "B₀ = 22.5°", RED, 30),
            (-math.pi / 8, "B₁ = −22.5°", RED, 30)]
    for th, lab, col, ext in angs:
        x2 = cx + (R + 2) * math.cos(th)
        y2 = cy - (R + 2) * math.sin(th)
        b.append(line(cx, cy, x2, y2, color=col, w=2.2))
        lx = cx + (R + ext) * math.cos(th)
        ly = cy - (R + ext) * math.sin(th) + 4
        anch = "start" if abs(th) < 0.7 else "middle"
        b.append(txt(lx + (4 if anch == "start" else 0), ly, lab, 11.5,
                     color=col, anchor=anch))
    b.append(txt(cx, cy + R + 26, "every needed pair is 22.5° "
                 "apart (67.5° for the ≠ case)", 11.5))
    b.append(txt(cx, cy + R + 44,
                 "each case wins with cos²(22.5°) "
                 "≈ 0.854", 12, w="bold", color=GREEN))
    return svg_wrap("\n".join(b), 540, 260,
                    "CHSH rules and the optimal measurement angles")


# ------------------------------------------------------- fig_chshcurve
def fig_chshcurve():
    def winb(beta):
        return chsh_win_quantum2((0.0, math.pi / 4), (beta, -beta))

    ax = Ax(0, 45, 0.45, 1.0, W=480, H=250, ml=52, mb=40)
    pts = sample(lambda d: winb(math.radians(d)), 0, 45, 80)
    b = ['<defs></defs>']
    b.append(ax.axes(xticks=(0, 11.25, 22.5, 33.75, 45),
                     yticks=(0.5, 0.75, 0.854, 1.0),
                     xlab="Bob offset β (degrees)",
                     ylab="",
                     xtickfmt=lambda t: fmt(t),
                     ytickfmt=lambda t: f"{t:g}"))
    b.append(ax.text(1.5, 0.955, "win probability", 11.5))
    b.append(ax.hline(0.75, color=ORA, dash="5 3"))
    b.append(ax.text(33.5, 0.705, "best classical 3/4", 11,
                     anchor="middle", color=ORA))
    b.append(ax.hline(math.cos(math.pi / 8) ** 2, color=GREEN,
                      dash="5 3"))
    b.append(ax.text(2, 0.868, "Tsirelson ceiling cos²(π/8)",
                     11, color=GREEN))
    b.append(ax.polyline(pts, BLUE, w=2.2))
    b.append(ax.vline(22.5, color=GRAY, dash="3 3"))
    b.append(ax.dot(22.5, winb(math.pi / 8), color=RED, r=4))
    b.append(ax.text(23.5, 0.9, "β = 22.5°", 11.5, color=RED))
    return svg_wrap("\n".join(b), 480, 250,
                    "entangled win probability of the rotation family")


# ---------------------------------------------------------- fig_repeat
def fig_repeat():
    b = []
    b.append(f'<defs>{arrow_def("c6rA", BLUE)}'
             f'{arrow_def("c6rB", RED)}</defs>')
    n = 5
    for i in range(n):
        x = 60 + i * 88
        b.append(box(x, 78, 72, 50, fill="#f7f7f7", stroke="#aaa",
                     rx=5, sw=1))
        b.append(txt(x + 36, 98, f"copy {i+1}", 11.5, color=MUT))
        b.append(txt(x + 36, 116, f"(x{chr(0x2080+i+1)}, "
                     f"y{chr(0x2080+i+1)})", 11.5))
    b.append(box(40, 20, 460, 34, fill="#eef4fa", stroke=BLUE))
    b.append(txt(270, 42, "Alice: ONE joint measurement on all of "
                 "x₁ … xₙ  →  answers a₁ … "
                 "aₙ", 12.5, color=BLUE))
    b.append(box(40, 152, 460, 34, fill="#fdeeee", stroke=RED))
    b.append(txt(270, 174, "Bob: ONE joint measurement on all of "
                 "y₁ … yₙ  →  answers b₁ … "
                 "bₙ", 12.5, color=RED))
    for i in range(n):
        x = 60 + i * 88 + 36
        b.append(line(x, 56, x, 76, color=BLUE, w=1.2, marker="c6rA"))
        b.append(line(x, 130, x, 150, color=RED, w=1.2, marker="c6rB"))
    b.append(txt(270, 208, "accept iff every copy is won — but "
                 "answer i may depend on all n questions of that "
                 "player", 12))
    return svg_wrap("\n".join(b), 540, 222,
                    "the repeated game and its cross-coordinate freedom")


# ----------------------------------------------------------- fig_feige
def fig_feige():
    b = []
    b.append(txt(270, 22, "Feige's counterexample game G"
                 + sb("F"), 13.5, w="bold"))
    b.append(txt(270, 40, "x, y independent fair bits; answers are "
                 "claims (player 1 got d) or (player 2 got d)",
                 11.5))
    b.append(txt(270, 55, "win iff both players make the same "
                 "true claim", 11.5))
    # two coordinate panels
    for k, (x0, ttl, aa, bb, cond) in enumerate([
            (50, "coordinate 1", "a₁ = (1, x₁)",
             "b₁ = (1, y₂)", "agree ⇔ y₂ = x₁"),
            (300, "coordinate 2", "a₂ = (2, x₁)",
             "b₂ = (2, y₂)", "agree ⇔ x₁ = y₂")]):
        b.append(box(x0, 62, 190, 96, fill="#f7f7f7", stroke="#aaa",
                     rx=6, sw=1))
        b.append(txt(x0 + 95, 82, ttl, 12, w="bold"))
        b.append(txt(x0 + 95, 106, aa, 12.5, color=BLUE))
        b.append(txt(x0 + 95, 128, bb, 12.5, color=RED))
        b.append(txt(x0 + 95, 150, cond, 11.5, color=MUT))
    b.append('<path d="M 148,166 C 210,190 330,190 392,166" '
             f'fill="none" stroke="{PUR}" stroke-width="1.6" '
             'stroke-dasharray="4 3"/>')
    b.append(txt(270, 196, "one shared event x₁ = y₂ makes "
                 "BOTH coordinates win at once", 12, color=PUR))
    b.append(box(90, 210, 360, 30, fill="#fdf7e3", stroke="#c9b45a",
                 rx=5, sw=1))
    gf = "G" + sb("F")
    b.append(txt(270, 230, "ω(" + gf + " ⊗ " + gf + ") = 1/2 "
                 "= ω(" + gf + ")  —  not 1/4", 12.5, w="bold"))
    return svg_wrap("\n".join(b), 540, 252,
                    "the cross strategy that defeats naive repetition")


# ------------------------------------------------------- fig_landscape
def fig_landscape():
    rows = [
        ("XOR games", "perfect: ω* (G⊗ⁿ) = ω*(G)"
         "ⁿ", "CSUU 2008", GREEN, 0),
        ("unique games", "exponential", "KRT 2010", GREEN, 0),
        ("projection games", "exponential", "DSV 2015", GREEN, 0),
        ("free games (product μ)", "exponential", "JPY14, CS15",
         GREEN, 0),
        ("anchored games (modified G)", "exponential", "BVY 2017",
         ORA, 1),
        ("ALL games", "polynomial n" + sp("−1/4"),
         "Yuen 2016", ORA, 2),
        ("ALL games", "exponential", "this chapter", RED, 3),
    ]
    b = []
    b.append(txt(60, 24, "class of entangled games", 12, w="bold",
                 anchor="start"))
    b.append(txt(310, 24, "decay proved", 12, w="bold", anchor="start"))
    b.append(txt(455, 24, "who", 12, w="bold", anchor="start"))
    for r, (cls, decay, who, col, mark) in enumerate(rows):
        y = 50 + r * 30
        if mark == 3:
            b.append(box(48, y - 19, 478, 27, fill="#fdf7e3",
                         stroke="#c9b45a", rx=4, sw=1.2))
        b.append(f'<circle cx="{56}" cy="{y-6}" r="4.5" '
                 f'fill="{col}"/>')
        b.append(txt(68, y - 1, cls, 12, anchor="start"))
        b.append(txt(310, y - 1, decay, 12, anchor="start",
                     color=(RED if mark == 3 else DARK)))
        b.append(txt(455, y - 1, who, 11.5, anchor="start", color=MUT))
    b.append(txt(270, 268, "green: special structure · orange: "
                 "general but weaker · the open cell was "
                 "“all games, exponential”", 11.5))
    return svg_wrap("\n".join(b), 540, 280,
                    "the landscape of parallel-repetition results")


# ----------------------------------------------------------- fig_decay
def fig_decay():
    ax = Ax(1, 60, -6.2, 0.15, W=500, H=250, ml=56, mb=40)
    b = []
    poly = sample(lambda n: math.log(min(1.0, 1.2 * math.log(max(
        n, 2.0)) / n ** 0.25)), 1, 60, 80)
    expo = sample(lambda n: -0.1 * n, 1, 60, 2)
    b.append(ax.axes(xticks=(1, 20, 40, 60),
                     yticks=(0, -2, -4, -6),
                     xlab="number of copies n",
                     ylab="",
                     ytickfmt=lambda t: f"{t:g}"))
    b.append(ax.polyline(poly, ORA, w=2.2))
    b.append(ax.polyline(expo, RED, w=2.2))
    b.append(ax.text(2.5, -5.4, "vertical axis: log of the bound "
                     "on ω*(G" + sp("⊗n") + ")", 11))
    b.append(ax.text(34, -1.35, "polynomial (Yuen 2016):", 11.5,
                     color=ORA))
    b.append(ax.text(34, -1.85, "log-bound flattens", 11.5, color=ORA))
    b.append(ax.text(8, -3.45, "exponential (Thm 1.1):", 11.5,
                     color=RED))
    b.append(ax.text(8, -3.95, "log-bound is a straight line", 11.5,
                     color=RED))
    return svg_wrap("\n".join(b), 500, 250,
                    "shape of polynomial versus exponential decay")


# --------------------------------------------------------- fig_roadmap
def fig_roadmap():
    b = []
    b.append(f'<defs>{arrow_def("c6mA")}</defs>')

    def node(x, y, w, h, lines, fill="#f7f7f7", stroke="#888",
             fs=11.5):
        out = [box(x, y, w, h, fill=fill, stroke=stroke, rx=6, sw=1.3)]
        for k, s in enumerate(lines):
            out.append(txt(x + w / 2, y + 16 + k * 15, s, fs))
        return out

    b += node(20, 14, 226, 46, ["assume the bound fails at n:",
                                "ϑ = P(win all n) ≫ bound"],
              fill="#fdeeee", stroke=RED)
    b += node(20, 82, 226, 46, ["Lemma 3.1 greedy conditioning:",
                                "core D, q ≥ 1−δ, "
                                "info cost η small"])
    b += node(20, 150, 226, 46, ["Lemma 3.2 ideal experiment:",
                                 "state Ψ wins live copy w.p. q"])
    b += node(20, 218, 226, 60, ["Lemma 3.3 (NEW): local states",
                                 "8η-close, histories "
                                 "κ-sampleable,", "NO 1/p loss"],
              fill="#fdf7e3", stroke="#c9b45a")
    b += node(294, 82, 226, 46, ["Lemma 3.4 classical correlated",
                                 "sampling (proved)"])
    b += node(294, 150, 226, 46, ["Lemma 3.5 = Import 1:",
                                  "embezzlement (DSV15)"])
    b += node(294, 218, 226, 60, ["Lemma 3.6 rounding: real",
                                  "strategy for G wins at least",
                                  "q − Bη" + sp("1/12")
                                  + " − 2Kα" + sp("1/12")])
    b += node(148, 302, 246, 44, ["contradiction: that exceeds",
                                  "ω*(G) = 1 − ε"],
              fill="#fdeeee", stroke=RED)
    b.append(line(133, 60, 133, 80, marker="c6mA"))
    b.append(line(133, 128, 133, 148, marker="c6mA"))
    b.append(line(133, 196, 133, 216, marker="c6mA"))
    b.append(line(246, 248, 292, 248, marker="c6mA"))
    b.append(line(407, 128, 407, 148, marker="c6mA"))
    b.append(line(407, 196, 407, 216, marker="c6mA"))
    b.append(line(350, 278, 300, 300, marker="c6mA"))
    b.append(txt(270, 360, "proof of Lemma 3.3 = Section 4: reveal "
                 "martingales, resolvent purification,", 11))
    b.append(txt(270, 374, "entropy budget, classical sampleability "
                 "(this document: §11–§14)", 11))
    return svg_wrap("\n".join(b), 540, 380,
                    "the proof roadmap")


# --------------------------------------------------------- fig_analogy
def fig_analogy():
    b = []
    b.append(txt(135, 24, "classical player", 12.5, w="bold"))
    b.append(txt(405, 24, "entangled player", 12.5, w="bold"))
    b.append(box(28, 36, 214, 118, fill="#f2f7f2", stroke=GREEN, rx=6))
    b.append(box(298, 36, 214, 118, fill="#fdeeee", stroke=RED, rx=6))
    for k, s in enumerate([
            "strategy = a function; its", "randomness can be resampled",
            "consistently with any", "conditioning — Holenstein",
            "embeds a fresh question", "into a conditioned transcript"]):
        b.append(txt(135, 56 + 16 * k, s, 11))
    for k, s in enumerate([
            "strategy = measurements on", "a state; conditioning on",
            "a rare win POSTSELECTS the", "state — measuring "
            "destroyed", "it, and a naive redo pays a",
            "factor 1/p ≈ eᵞⁿ"]):
        b.append(txt(405, 56 + 16 * k, s, 11))
    b.append(txt(270, 180, "the new lemma makes quantum conditioning "
                 "survivable: its estimates never divide by p", 12,
                 w="bold"))
    return svg_wrap("\n".join(b), 540, 196,
                    "why quantum conditioning is the hard part")


# ---------------------------------------------------- fig_greedy (SMIL)
def fig_greedy():
    # toy conditional-success bars over 3 greedy rounds (computed here)
    n = 10
    ra = np.random.default_rng(9)
    comp = [ra.uniform(0.55, 0.99, size=n) for _ in range(3)]
    mixw = np.array([0.45, 0.35, 0.2])
    outs = np.array(list(itertools.product((0, 1), repeat=n)))
    probs = np.zeros(len(outs))
    for w, pv in zip(mixw, comp):
        probs += w * np.prod(np.where(outs == 1, pv, 1 - pv), axis=1)
    stages = []
    D = []
    for _ in range(3):
        inD = (np.all(outs[:, D] == 1, axis=1) if D
               else np.ones(len(outs), bool))
        pD = probs[inD].sum()
        cond = []
        for i in range(n):
            if i in D:
                cond.append(None)
            else:
                cond.append(probs[inD & (outs[:, i] == 1)].sum() / pD)
        stages.append(list(cond))
        worst = min((c, i) for i, c in enumerate(cond)
                    if c is not None)[1]
        D.append(worst)
    picks = D[:]
    W, H = 520, 260
    x0, y0, bw, bh = 56, 210, 40, 150
    b = []
    b.append(line(x0 - 10, y0, x0 + n * (bw + 4), y0, color=DARK,
                  w=1.2))
    b.append(txt(x0 - 16, y0 - bh + 4, "1.0", 10.5, anchor="end"))
    b.append(txt(x0 - 16, y0 + 4, "0.0", 10.5, anchor="end"))
    dur = 9
    for i in range(n):
        x = x0 + i * (bw + 4)
        b.append(txt(x + bw / 2, y0 + 16, str(i + 1), 10.5))
        vals = []
        for s in range(3):
            v = stages[s][i]
            vals.append(v if v is not None else 1.0)
        hseq = [bh * v for v in vals]
        fills = []
        for s in range(3):
            if stages[s][i] is None:
                fills.append(RED)
            elif i == picks[s]:
                fills.append(ORA)
            else:
                fills.append(BLUE)
        hs = ";".join(f"{h:.1f}" for h in
                      [hseq[0], hseq[0], hseq[1], hseq[1],
                       hseq[2], hseq[2]])
        ys = ";".join(f"{y0-h:.1f}" for h in
                      [hseq[0], hseq[0], hseq[1], hseq[1],
                       hseq[2], hseq[2]])
        fs = ";".join([fills[0], fills[0], fills[1], fills[1],
                       fills[2], fills[2]])
        kt = "0;0.32;0.34;0.66;0.68;1"
        b.append(
            f'<rect x="{x}" y="{y0-hseq[0]:.1f}" width="{bw}" '
            f'height="{hseq[0]:.1f}" fill="{fills[0]}" opacity="0.85">'
            f'<animate attributeName="height" values="{hs}" '
            f'keyTimes="{kt}" dur="{dur}s" repeatCount="indefinite"/>'
            f'<animate attributeName="y" values="{ys}" keyTimes="{kt}" '
            f'dur="{dur}s" repeatCount="indefinite"/>'
            f'<animate attributeName="fill" values="{fs}" '
            f'keyTimes="{kt}" dur="{dur}s" repeatCount="indefinite"/>'
            f'</rect>')
    b.append(txt(258, 246, "coordinate i — bars: P(Wᵢ | "
                 "win the red core so far); orange = next greedy pick, "
                 "red = already in D", 11))
    b.append(txt(258, 30, "greedy conditioning: absorb the worst "
                 "coordinate into the core D; the rest rise", 12,
                 w="bold"))
    return svg_wrap("\n".join(b), W, H,
                    "greedy conditioning animation")


# ---------------------------------------------------------- fig_reveal
def fig_reveal():
    b = []
    n = 12
    Dset = {2, 7}
    live = 5
    cX = {0, 1, 3, 8, 10}   # extra alice-revealed
    x0, bw = 40, 38
    b.append(txt(20, 58, "Alice", 12, anchor="start", color=BLUE,
                 w="bold"))
    b.append(txt(20, 108, "Bob", 12, anchor="start", color=RED,
                 w="bold"))
    for i in range(n):
        x = x0 + 32 + i * bw
        for row, y in ((0, 40), (1, 90)):
            inC = (i in Dset or (row == 0 and i in cX)
                   or (row == 1 and i == 0)
                   or (row == 1 and i not in cX and i != live
                       and i not in Dset))
            if i == live:
                fill = "#fdf7e3"
                stroke = "#c9b45a"
            elif i in Dset:
                fill = "#fdeeee"
                stroke = RED
            elif inC:
                fill = "#e8f0e8"
                stroke = GREEN
            else:
                fill = "#fff"
                stroke = "#bbb"
            b.append(box(x, y, bw - 5, 32, fill=fill, stroke=stroke,
                         rx=4, sw=1.2))
            lab = ("?" if not inC and i != live else
                   ("★" if i == live else "✓"))
            if i == live:
                lab = "live"
            b.append(txt(x + (bw - 5) / 2, y + 21,
                         lab, 10.5 if i == live else 12,
                         color=(MUT if lab == "?" else DARK)))
        b.append(txt(x + (bw - 5) / 2, 30, str(i + 1), 10, color=MUT))
    b.append(txt(270, 143, "✓ revealed by the history T · "
                 "? unrevealed", 11))
    b.append(txt(270, 157, "red: core D (both revealed) · gold: "
                 "live coordinate i (neither revealed)", 11))
    b.append(txt(270, 174, "every non-live column has at least one "
                 "✓ — that is C" + sb("X") + " ∪ C"
                 + sb("Y") + " = [n] ∖ {i}", 11.5, w="bold"))
    return svg_wrap("\n".join(b), 540, 182,
                    "the reveal pattern of the history")


# ---------------------------------------------------------- fig_factor
def fig_factor():
    b = []
    b.append(box(30, 34, 220, 40, fill="#eef4fa", stroke=BLUE, rx=6))
    b.append(txt(140, 52, "still-random Alice questions", 11.5,
                 color=BLUE))
    b.append(txt(140, 66, "live only on [n] ∖ (C" + sb("X")
                 + " ∪ {i})", 11, color=BLUE))
    b.append(box(290, 34, 220, 40, fill="#fdeeee", stroke=RED, rx=6))
    b.append(txt(400, 52, "still-random Bob questions", 11.5,
                 color=RED))
    b.append(txt(400, 66, "live only on [n] ∖ (C" + sb("Y")
                 + " ∪ {i})", 11, color=RED))
    b.append(txt(270, 100, "C" + sb("X") + " ∪ C" + sb("Y")
                 + " = [n] ∖ {i} ⇒ these two index sets "
                 "are DISJOINT", 12.5, w="bold"))
    b.append(txt(270, 126, "so, given (T, Xᵢ, Yᵢ), no jointly "
                 "sampled pair (Xⱼ, Yⱼ) remains:", 12))
    b.append(box(60, 138, 420, 32, fill="#fdf7e3", stroke="#c9b45a",
                 rx=5))
    b.append(txt(270, 159, "P(xⁿ, yⁿ | t, x, y) = "
                 "P(xⁿ | t, x) · P(yⁿ | t, y)   "
                 "under the prior P", 12.5))
    return svg_wrap("\n".join(b), 540, 184,
                    "why the unrevealed questions factorize")


# ---------------------------------------------------------- fig_purify
def fig_purify():
    b = []
    b.append(f'<defs>{arrow_def("c6pA")}</defs>')
    b.append(box(30, 20, 130, 46, fill="#f7f7f7", stroke="#888", rx=6))
    b.append(txt(95, 40, "effect F", 12.5))
    b.append(txt(95, 56, "0 ⪯ F ⪯ I on ℋ", 11,
                 color=MUT))
    b.append(line(160, 43, 220, 43, marker="c6pA"))
    b.append(txt(190, 34, "Γ", 13, cls="v"))
    b.append(box(222, 20, 176, 46, fill="#eef4fa", stroke=BLUE, rx=6))
    b.append(txt(310, 40, "cross operator Γ(F)", 12.5,
                 color=BLUE))
    b.append(txt(310, 56, "ℋ → ℋ ⊗ 𝒜",
                 11.5, color=BLUE))
    b.append(box(70, 90, 400, 34, fill="#fdf7e3", stroke="#c9b45a",
                 rx=5))
    b.append(txt(270, 112, "Γ(F)†Γ(F) = F  ⇒  "
                 "‖Γ(F)ψ‖² = ⟨ψ| F "
                 "|ψ⟩:  branch keeps the exact Born "
                 "probability", 12))
    b.append(txt(270, 150, "a branch records “this outcome "
                 "happened” without renormalizing yet; Γ "
                 "chooses HOW the state updates", 11.5, color=MUT))
    return svg_wrap("\n".join(b), 540, 164,
                    "purification of a POVM effect")


# -------------------------------------------------------- fig_triangle
def fig_triangle():
    b = []
    tx, ty = 230, 44
    lx, ly = 80, 180
    rx_, ry = 380, 180
    b.append(line(tx, ty + 10, lx + 14, ly - 12, color=GRAY, w=1.6))
    b.append(line(tx, ty + 10, rx_ - 14, ry - 12, color=GRAY, w=1.6))
    b.append(line(lx + 22, ly, rx_ - 22, ry, color=GRAY, w=1.6,
                  dash="5 3"))
    b.append(f'<circle cx="{tx}" cy="{ty}" r="7" fill="{PUR}"/>')
    b.append(f'<circle cx="{lx}" cy="{ly}" r="7" fill="{BLUE}"/>')
    b.append(f'<circle cx="{rx_}" cy="{ry}" r="7" fill="{RED}"/>')
    b.append(txt(tx, ty - 16, "ideal Ψ (needs r, x AND y)", 12.5,
                 color=PUR, w="bold"))
    b.append(txt(lx, ly + 24, "Ψᴬ: Alice can", 12,
                 color=BLUE))
    b.append(txt(lx, ly + 40, "describe from (r, x)", 12, color=BLUE))
    b.append(txt(rx_, ry + 24, "Ψᴮ: Bob can", 12, color=RED))
    b.append(txt(rx_, ry + 40, "describe from (r, y)", 12, color=RED))
    b.append(txt(136, 102, "E‖Ψ − Ψᴬ‖"
                 "² ≤ 8η", 12.5, color=DARK))
    b.append(txt(330, 102, "E‖Ψ − Ψᴮ‖"
                 "² ≤ 8η", 12.5))
    b.append(txt(230, 168, "E‖Ψᴬ − Ψᴮ"
                 "‖² ≤ 32η (triangle)", 11.5,
                 color=MUT))
    return svg_wrap("\n".join(b), 460, 232,
                    "the ideal state and its two local descriptions")


# --------------------------------------------------- fig_distributions
def fig_distributions():
    b = []
    b.append(f'<defs>{arrow_def("c6dA")}</defs>')
    b.append(f'<ellipse cx="250" cy="66" rx="104" ry="34" '
             f'fill="#f3ecf7" stroke="{PUR}" stroke-width="1.5"/>')
    b.append(txt(250, 62, "posterior Q = P( · | W" + sb("D")
                 + " )", 12.5, color=PUR))
    b.append(txt(250, 80, "nobody can sample it alone", 10.5,
                 color=MUT))
    b.append(f'<ellipse cx="96" cy="180" rx="86" ry="32" '
             f'fill="#eef4fa" stroke="{BLUE}" stroke-width="1.5"/>')
    b.append(txt(96, 176, "J" + sb("A") + ": Alice samples", 12,
                 color=BLUE))
    b.append(txt(96, 193, "r from Q(R | i, Xᵢ = x)", 11,
                 color=BLUE))
    b.append(f'<ellipse cx="404" cy="180" rx="86" ry="32" '
             f'fill="#fdeeee" stroke="{RED}" stroke-width="1.5"/>')
    b.append(txt(404, 176, "J" + sb("B") + ": Bob samples", 12,
                 color=RED))
    b.append(txt(404, 193, "r from Q(R | i, Yᵢ = y)", 11,
                 color=RED))
    b.append(line(180, 108, 128, 148, color=DARK, w=1.3,
                  marker="c6dA"))
    b.append(line(320, 108, 372, 148, color=DARK, w=1.3,
                  marker="c6dA"))
    b.append(txt(122, 122, "d" + sb("TV") + " ≤ κ", 12))
    b.append(txt(382, 122, "d" + sb("TV") + " ≤ κ", 12))
    b.append(txt(250, 235, "correlated sampling couples them: "
                 "P(r" + sb("A") + " ≠ r" + sb("B")
                 + ") ≤ 4κ,   κ = √(3η/2)", 12.5,
                 w="bold"))
    return svg_wrap("\n".join(b), 500, 250,
                    "one ideal posterior, two locally samplable copies")


# -------------------------------------------------------- fig_pipeline
def fig_pipeline():
    b = []
    b.append(f'<defs>{arrow_def("c6qA")}</defs>')

    def stage(x, w, lines, err, fill="#f7f7f7", stroke="#888"):
        out = [box(x, 44, w, 74, fill=fill, stroke=stroke, rx=6)]
        for k, s in enumerate(lines):
            out.append(txt(x + w / 2, 62 + k * 15, s, 10.8))
        if err:
            out.append(box(x, 132, w, 30, fill="#fdeeee", stroke=RED,
                           rx=5, sw=1))
            out.append(txt(x + w / 2, 151, err, 10.5, color=RED))
        return out

    rA = "r" + sb("A")
    rB = "r" + sb("B")
    b += stage(16, 122, ["shared randomness:", "uniform live i,",
                         "correlated-sampling", "table (L3.4)"],
               "history mismatch 4κ")
    b += stage(154, 122, ["histories " + rA + ", " + rB,
                          "from own questions",
                          rA + " = " + rB + " = r w.h.p.",
                          "distribution J" + sb("A")],
               "bias vs Q: κ each")
    b += stage(292, 122, ["embezzlement (L3.5):", "U(Ψᴬ) "
                          "⊗ V(Ψᴮ)", "on |E⟩ "
                          "prepares ≈ Ψᴬ",
                          "no communication"],
               "K(α" + sp("1/12") + " + ‖Δ‖"
               + sp("1/6") + ")")
    b += stage(430, 96, ["measure with", "Mᵃ ⊗ Nᵇ",
                         "of Lemma 3.2", "→ answers"],
               "state error 2e")
    b.append(line(138, 81, 152, 81, marker="c6qA"))
    b.append(line(276, 81, 290, 81, marker="c6qA"))
    b.append(line(414, 81, 428, 81, marker="c6qA"))
    b.append(txt(270, 24, "the rounded single-game strategy S"
                 + sb("α") + " on fresh questions (x, y) ~ μ",
                 12.5, w="bold"))
    b.append(txt(270, 190, "total: win(S" + sb("α") + ") ≥ q − "
                 "B" + sb("qs") + "η" + sp("1/12") + " − 2K"
                 + sb("qs") + "α" + sp("1/12") + " — every loss "
                 "is a power of η or α, never 1/p", 12))
    return svg_wrap("\n".join(b), 540, 204,
                    "the rounding pipeline and its error budget")


# ------------------------------------------------------- fig_embezzle
def fig_embezzle():
    bvals = [2 ** k for k in range(2, 17)]
    errs = [embezzle_error(v) for v in bvals]
    # left: coefficients of E_b for b=16
    bpx = []
    bN = 16
    j = np.arange(1, bN + 1, dtype=float)
    c2 = (1 / j) / (1 / j).sum()
    c = np.sqrt(c2)
    x0, y0, bw = 46, 180, 14
    bpx.append(txt(150, 26, "Schmidt weights of |E₁₆⟩",
                   12, w="bold"))
    for k in range(bN):
        h = 150 * c[k] / c[0]
        bpx.append(f'<rect x="{x0 + k * bw}" y="{y0 - h:.1f}" '
                   f'width="{bw - 3}" height="{h:.1f}" fill="{BLUE}" '
                   f'opacity="0.8"/>')
    bpx.append(line(x0 - 6, y0, x0 + bN * bw + 4, y0, w=1.2))
    bpx.append(txt(150, y0 + 18, "j = 1 … 16:  coefficient "
                   "∝ 1/√j", 11.5))
    ax = Ax(2, 16, 0, 0.62, W=220, H=210, ml=40, mr=10, mt=30, mb=40)
    right = []
    right.append(txt(420, 26, "embezzlement error", 12, w="bold"))
    right.append(ax.axes(xticks=(4, 8, 12, 16),
                         yticks=(0, 0.2, 0.4, 0.6),
                         xlab="log₂ b", ylab=""))
    pts = list(zip([math.log2(v) for v in bvals], errs))
    right.append(ax.polyline(pts, RED, w=2.2))
    for lb, e in zip([math.log2(v) for v in bvals], errs):
        right.append(ax.dot(lb, e, r=2.6, color=RED))
    right.append(ax.text(6.4, 0.14, "error ≈ c / log b → 0",
                         11.5, color=RED))
    g = ('<g transform="translate(300,0)">'
         + "\n".join(right[1:]) + "</g>")
    return svg_wrap("\n".join(bpx) + "\n" + right[0] + "\n" + g,
                    540, 240,
                    "the van Dam-Hayden embezzlement state")


# --------------------------------------------------- fig_sampling (SMIL)
def fig_sampling():
    P = [0.5, 0.3, 0.2]
    Q = [0.3, 0.3, 0.4]
    x0, y0, cw, sc = 70, 200, 110, 150
    b = []
    for k in range(3):
        x = x0 + k * cw
        hP = sc * P[k]
        hQ = sc * Q[k]
        b.append(box(x, y0 - sc, 84, sc, fill="#fafafa", stroke="#ccc",
                     rx=0, sw=1))
        b.append(f'<rect x="{x}" y="{y0-hP:.0f}" width="40" '
                 f'height="{hP:.0f}" fill="{BLUE}" opacity="0.55"/>')
        b.append(f'<rect x="{x+44}" y="{y0-hQ:.0f}" width="40" '
                 f'height="{hQ:.0f}" fill="{RED}" opacity="0.55"/>')
        b.append(txt(x + 42, y0 + 16, f"z = {k+1}", 11))
        b.append(txt(x + 20, y0 - hP - 6, f"P={P[k]:g}", 10.5,
                     color=BLUE))
        b.append(txt(x + 64, y0 - hQ - 6, f"Q={Q[k]:g}", 10.5,
                     color=RED))
    # darts: (z,u) proposals; animate a falling dot cycling positions
    seq = [(0, 0.42, "onlyP"), (2, 0.34, "onlyQ"), (0, 0.25, "both"),
           (1, 0.18, "both"), (2, 0.36, "onlyQ")]
    dur = 10
    xs = []
    ys = []
    cols = []
    for z, u, kind in seq:
        xs.append(x0 + z * cw + 42)
        ys.append(y0 - sc * u)
        cols.append(GREEN if kind == "both" else ORA)
    n = len(seq)
    kt = ";".join(f"{k/n:.3f}" for k in range(n)) + ";1"
    xv = ";".join(f"{x}" for x in xs + [xs[0]])
    yv = ";".join(f"{y:.0f}" for y in ys + [ys[0]])
    cv = ";".join(cols + [cols[0]])
    b.append(f'<circle cx="{xs[0]}" cy="{ys[0]:.0f}" r="6" '
             f'fill="{GREEN}">'
             f'<animate attributeName="cx" values="{xv}" '
             f'keyTimes="{kt}" dur="{dur}s" calcMode="discrete" '
             f'repeatCount="indefinite"/>'
             f'<animate attributeName="cy" values="{yv}" '
             f'keyTimes="{kt}" dur="{dur}s" calcMode="discrete" '
             f'repeatCount="indefinite"/>'
             f'<animate attributeName="fill" values="{cv}" '
             f'keyTimes="{kt}" dur="{dur}s" calcMode="discrete" '
             f'repeatCount="indefinite"/></circle>')
    b.append(txt(270, 26, "shared darts (z, u): accept when u falls "
                 "under your own curve", 12.5, w="bold"))
    b.append(txt(270, 232, "green dart: accepted by BOTH → same "
                 "output · orange: accepted by one only → "
                 "possible disagreement", 11))
    dtv_ = "d" + sb("TV")
    b.append(txt(270, 252, "P(outputs differ) = 2" + dtv_ + "/(1 + "
                 + dtv_ + ") ≤ 2" + dtv_ + "(P, Q)", 12,
                 w="bold"))
    return svg_wrap("\n".join(b), 540, 266,
                    "correlated sampling by shared darts")


# --------------------------------------------------------- fig_pinsker
def fig_pinsker():
    def Dkl(p):
        out = 0.0
        for a, bq in ((p, 0.5), (1 - p, 0.5)):
            if a > 0:
                out += a * math.log(a / bq)
        return out

    ax = Ax(0, 1, 0, 0.8, W=480, H=240, ml=52, mb=40)
    b = []
    b.append(ax.axes(xticks=(0, 0.25, 0.5, 0.75, 1),
                     yticks=(0, 0.2, 0.4, 0.6, 0.8),
                     xlab="p (binary distribution vs fair coin)",
                     ylab=""))
    b.append(ax.polyline(sample(lambda p: abs(p - 0.5), 0, 1, 80),
                         BLUE, w=2.2))
    b.append(ax.polyline(
        sample(lambda p: math.sqrt(max(0.0, 0.5 * Dkl(p))),
               1e-6, 1 - 1e-6, 100),
        RED, w=2.2))
    b.append(ax.text(0.03, 0.62, "√(D/2): the Pinsker "
                     "majorant", 11.5, color=RED))
    b.append(ax.text(0.72, 0.115, "d" + sb("TV")
                     + " = |p − 1/2|", 11.5, color=BLUE))
    return svg_wrap("\n".join(b), 480, 240,
                    "Pinsker's inequality on the binary family")


# -------------------------------------------------------------- fig_h1
def fig_h1():
    ax = Ax(0, 1, 0, 0.4, W=480, H=240, ml=52, mb=40)
    b = []
    b.append(ax.axes(xticks=(0, 0.25, 0.5, 0.75, 1),
                     yticks=(0, 0.1, 0.2, 0.3, 0.37),
                     xlab="v", ylab=""))
    b.append(ax.polyline(
        sample(lambda v: -v * math.log(v) if v > 0 else 0.0,
               0, 1, 300), BLUE, w=2.2))
    va, vb_ = 0.15, 0.75
    ha = -va * math.log(va)
    hb = -vb_ * math.log(vb_)
    vm = 0.5 * (va + vb_)
    hm = -vm * math.log(vm)
    b.append(ax.polyline([(va, ha), (vb_, hb)], ORA, w=1.8,
                         dash="5 3"))
    b.append(ax.dot(va, ha, color=ORA))
    b.append(ax.dot(vb_, hb, color=ORA))
    b.append(ax.dot(vm, hm, color=RED, r=4))
    b.append(ax.dot(vm, 0.5 * (ha + hb), color=ORA, r=4))
    b.append(ax.vline(vm, color=GRAY, dash="2 3", y0=0.5 * (ha + hb),
                      y1=hm))
    b.append(ax.text(vm + 0.02, hm + 0.02,
                     "H₁(E v) ≥ E H₁(v)", 11.5,
                     color=RED))
    b.append(ax.text(0.62, 0.115, "chord: E H₁(v)", 11,
                     color=ORA))
    b.append(ax.text(0.05, 0.36, "H₁(v) = −v log v", 12,
                     color=BLUE))
    return svg_wrap("\n".join(b), 480, 240,
                    "the entropy function and its concavity gap")


# ------------------------------------------------- fig_martingale (SMIL)
def fig_martingale():
    # genuine Doob martingale: f(omega) on 8 bits revealed one at a time
    ra = np.random.default_rng(3)
    f = ra.uniform(0.05, 0.95, size=256)
    omega = 179  # fixed realized path
    path = []
    for j in range(9):
        mask = 256
        sel = [o for o in range(256)
               if (o >> (8 - j)) == (omega >> (8 - j))]
        path.append(float(np.mean(f[sel])))
    ax = Ax(0, 8, 0, 1, W=520, H=250, ml=50, mb=42)
    b = []
    b.append(ax.axes(xticks=list(range(9)), yticks=(0, 0.5, 1),
                     xlab="questions revealed j",
                     ylab=""))
    pts = [(j, v) for j, v in enumerate(path)]
    b.append(ax.polyline(pts, BLUE, w=2.2))
    for j, v in pts:
        b.append(ax.dot(j, v, r=3, color=BLUE))
    kcut = 4
    b.append(ax.vline(kcut, color=GRAY, dash="4 3"))
    b.append(ax.vline(kcut + 1, color=GRAY, dash="4 3"))
    x1, y1 = ax.X(kcut), ax.Y(path[kcut])
    x2, y2 = ax.X(kcut + 1), ax.Y(path[kcut + 1])
    b.append(f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x1)}" '
             f'y2="{fmt(y2)}" stroke="{RED}" stroke-width="3"/>')
    b.append(ax.text(kcut + 1.25, 0.16, "live step k: the", 11.5,
                     color=RED))
    b.append(ax.text(kcut + 1.25, 0.09, "increment we must bound",
                     11.5, color=RED))
    b.append(ax.text(0.15, 0.93, "the prediction E[F | first j "
                     "reveals]", 12, color=BLUE))
    b.append(ax.text(0.15, 0.86, "is a Doob martingale", 12,
                     color=BLUE))
    # sweeping dot
    xv = ";".join(fmt(ax.X(j)) for j, _ in pts)
    yv = ";".join(fmt(ax.Y(v)) for _, v in pts)
    b.append(f'<circle cx="{fmt(ax.X(0))}" cy="{fmt(ax.Y(path[0]))}" '
             f'r="6" fill="{ORA}" opacity="0.9">'
             f'<animate attributeName="cx" values="{xv}" dur="8s" '
             f'repeatCount="indefinite"/>'
             f'<animate attributeName="cy" values="{yv}" dur="8s" '
             f'repeatCount="indefinite"/></circle>')
    return svg_wrap("\n".join(b), 520, 250,
                    "a reveal martingale and its live increment")


# ---------------------------------------------------------- fig_gsigma
def fig_gsigma():
    ax = Ax(0, 6, 0, 1.05, W=500, H=240, ml=52, mb=42)
    b = []
    b.append(ax.axes(xticks=(0, 2, 4, 6), yticks=(0, 0.5, 1),
                     xlab="u", ylab=""))
    cols = [(0.5, BLUE), (1.5, RED), (4.0, GREEN)]
    for sg, col in cols:
        b.append(ax.polyline(
            sample(lambda u, s=sg: s / (s + u), 0, 6, 80), col,
            w=2.2))
    b.append(ax.path_area(
        sample(lambda u: 0.5 / (0.5 + u), 0, 6, 80), BLUE,
        opacity=0.14))
    b.append(ax.text(1.05, 0.18, "σ = 1/2", 11.5, color=BLUE))
    b.append(ax.text(1.9, 0.52, "σ = 3/2", 11.5, color=RED))
    b.append(ax.text(4.1, 0.62, "σ = 4", 11.5, color=GREEN))
    b.append(ax.text(0.85, 0.99, "g" + sb("σ") + "(u) = "
                     "σ/(σ + u):  ∫ g" + sb("σ")
                     + "² du = σ exactly", 12))
    b.append(ax.text(1.35, 0.88, "the squared curve encloses "
                     "area σ — Born survives", 11, color=MUT))
    return svg_wrap("\n".join(b), 500, 240,
                    "the resolvent profile functions")


# -------------------------------------------------------- fig_telescope
def fig_telescope():
    vals = [0.36, 0.30, 0.22, 0.11]
    ax = Ax(-0.5, 3.5, 0, 0.47, W=500, H=230, ml=56, mb=40)
    b = []
    b.append(ax.axes(xticks=(0, 1, 2, 3), yticks=(0, 0.2, 0.37),
                     xlab="martingale step j",
                     ylab="", xtickfmt=lambda t: str(int(t))))
    for j, v in enumerate(vals):
        x = ax.X(j) - 26
        y = ax.Y(v)
        h = ax.Y(0) - y
        b.append(f'<rect x="{fmt(x)}" y="{fmt(y)}" width="52" '
                 f'height="{fmt(h)}" fill="{BLUE}" opacity="0.7"/>')
        if j:
            b.append(f'<line x1="{fmt(x-8)}" y1="{fmt(ax.Y(vals[j-1]))}"'
                     f' x2="{fmt(x+60)}" y2="{fmt(ax.Y(vals[j-1]))}" '
                     f'stroke="{ORA}" stroke-width="1.4" '
                     f'stroke-dasharray="4 3"/>')
            ym = 0.5 * (ax.Y(vals[j - 1]) + ax.Y(v))
            b.append(txt(ax.X(j) + 34, ym + 4,
                         f"cost {j}", 10.5, color=RED, anchor="start"))
    b.append(ax.text(-0.4, 0.455, "entropy budget Tr(ρ H₁"
                     "(Fⱼ)): each reveal step SPENDS some", 11.5))
    b.append(ax.text(-0.4, 0.415, "budget; total spent ≤ H₁"
                     "(p₀), so a random step costs ≤ "
                     "H₁(p₀)/N", 11.5))
    return svg_wrap("\n".join(b), 500, 230,
                    "the telescoping entropy budget")


# --------------------------------------------------------- fig_sizebias
def fig_sizebias():
    b = []
    b.append(txt(260, 26, "the size-bias cancellation in Lemma 4.6",
                 13, w="bold"))
    b.append(box(30, 44, 210, 58, fill="#eef4fa", stroke=BLUE, rx=6))
    b.append(txt(135, 66, "reverse partition picks the", 11.5))
    b.append(txt(135, 82, "live block: weight 2N/m", 11.5))
    b.append(box(280, 44, 210, 58, fill="#fdeeee", stroke=RED, rx=6))
    b.append(txt(385, 66, "uniform cut inside the block:", 11.5))
    b.append(txt(385, 82, "each step weight 1/N", 11.5))
    b.append(txt(260, 130, "(2N/m) × (1/N) = 2/m  —  the "
                 "block size N cancels exactly", 13, w="bold",
                 color=GREEN))
    b.append(txt(260, 158, "summing the fair-partition weights "
                 "2⁻ᵐ over all partitions gives ≤ 1, so "
                 "I" + sb("A") + ", I" + sb("B")
                 + " ≤ 2p(τ + s)/m", 12))
    return svg_wrap("\n".join(b), 520, 176,
                    "why no block-size factor survives")


# -------------------------------------------------------- fig_stability
def fig_stability():
    b = []
    b.append(f'<defs>{arrow_def("c6sA")}</defs>')
    b.append(txt(135, 24, "naive route (lossy)", 12.5, w="bold",
                 color=ORA))
    b.append(box(28, 36, 214, 96, fill="#fff8ee", stroke=ORA, rx=6))
    for k, s in enumerate([
            "prove a bound on the PRIOR", "experiment, then condition:",
            "posterior error ≤ prior/p.", "with p ≈ e⁻"
            "ᵞⁿ this explodes"]):
        b.append(txt(135, 58 + 17 * k, s, 11))
    b.append(txt(405, 24, "postselection-stable route", 12.5, w="bold",
                 color=GREEN))
    b.append(box(298, 36, 214, 96, fill="#f2f7f2", stroke=GREEN, rx=6))
    for k, s in enumerate([
            "keep branches UNNORMALIZED:", "budget itself carries a "
            "factor p", "(Lemma 4.5).  Posterior weight",
            "divides by p — exact cancellation"]):
        b.append(txt(405, 58 + 17 * k, s, 11))
    b.append(txt(135, 158, "cost ∝ (τ+s)/(m·p)", 12.5,
                 color=ORA))
    b.append(txt(405, 158, "cost ≤ 8(τ+s)/m = 8η", 12.5,
                 color=GREEN, w="bold"))
    b.append(txt(270, 186, "this single cancellation is what upgrades "
                 "polynomial decay (Yuen 2016) to exponential decay",
                 12))
    return svg_wrap("\n".join(b), 540, 200,
                    "the postselection-stability idea")


# -------------------------------------------------------- fig_waterfall
def fig_waterfall():
    Kqs, Ustar, Bqs, cqs = constants_chain()
    eps, ell = 0.3, math.log(4.0)
    gam = cqs * eps ** 13 / (eps + ell)
    delta = eps / 4
    etab = 2 * gam * (1 + 4 * ell / eps)
    alpha = (eps / (16 * Kqs)) ** 12
    vals = [("ε", eps), ("δ = ε/4", delta),
            ("α", alpha), ("η bound", etab),
            ("γ", gam)]
    ax = Ax(-0.6, 4.6, -34, 1.5, W=520, H=290, ml=58, mb=52, mt=48)
    b = []
    b.append(ax.axes(xticks=(), yticks=(0, -10, -20, -30),
                     xlab="", ylab="",
                     ytickfmt=lambda t: str(int(t))))
    for k, (lab, v) in enumerate(vals):
        lg = math.log10(v)
        y = ax.Y(lg)
        h = ax.Y(-34) - y
        b.append(f'<rect x="{fmt(ax.X(k)-30)}" y="{fmt(y)}" '
                 f'width="60" height="{fmt(h)}" fill="{BLUE}" '
                 f'opacity="0.75"/>')
        b.append(txt(ax.X(k), y - 7, f"{v:.1e}", 10.5))
        b.append(txt(ax.X(k), ax.Y(-34) + 18, lab, 11.5))
    b.append(txt(260, 20, "ε = 0.3, ℓ = log 4, K"
                 + sb("qs") + " = 1 (illustrative floor): "
                 "log₁₀ of each parameter", 12))
    b.append(ax.text(1.62, -6, "the cascade spans 30 orders", 11.5,
                     color=RED))
    b.append(ax.text(1.62, -9.6, "of magnitude — the theorem is",
                     11.5, color=RED))
    b.append(ax.text(1.62, -13.2, "qualitative, not practical", 11.5,
                     color=RED))
    return svg_wrap("\n".join(b), 520, 290,
                    "the parameter cascade at work")


# ----------------------------------------------------------- fig_eps13
def fig_eps13():
    b = []
    b.append(f'<defs>{arrow_def("c6eA")}</defs>')
    b.append(txt(260, 24, "where ε¹³ comes from", 13,
                 w="bold"))
    steps = [
        ("rounding must beat ε/4:", "Bη" + sp("1/12")
         + " ≤ ε/4 needs η ≲ ε¹²", 46),
        ("the 12th power is the embezzlement", "exponent 1/12 of "
         "Import 1 (DSV15)", 88),
        ("η ≈ γ(1 + ℓ/ε)·const: one more "
         "ε", "and the alphabet term ε + ℓ", 130),
        ("γ ∝ ε¹² · ε/(ε + "
         "ℓ) = ε¹³/(ε + ℓ)", "", 172),
    ]
    for t1, t2, y in steps:
        b.append(box(60, y - 16, 400, 40 if t2 else 28, fill="#f7f7f7",
                     stroke="#999", rx=5, sw=1.1))
        b.append(txt(260, y, t1, 11.8))
        if t2:
            b.append(txt(260, y + 16, t2, 11.8))
        if y < 172:
            b.append(line(260, y + (24 if t2 else 12), 260,
                          y + (40 if t2 else 28), marker="c6eA"))
    b.append(txt(260, 216, "the exponent 13 is an artifact of the "
                 "tools, not conjectured optimal", 11.5, color=MUT))
    return svg_wrap("\n".join(b), 520, 230,
                    "anatomy of the epsilon exponent")


FIGS = [
    ("fig_game", fig_game), ("fig_born", fig_born),
    ("fig_entangle", fig_entangle), ("fig_chsh", fig_chsh),
    ("fig_chshcurve", fig_chshcurve), ("fig_repeat", fig_repeat),
    ("fig_feige", fig_feige), ("fig_landscape", fig_landscape),
    ("fig_decay", fig_decay), ("fig_roadmap", fig_roadmap),
    ("fig_analogy", fig_analogy), ("fig_greedy", fig_greedy),
    ("fig_reveal", fig_reveal), ("fig_factor", fig_factor),
    ("fig_purify", fig_purify), ("fig_triangle", fig_triangle),
    ("fig_distributions", fig_distributions),
    ("fig_pipeline", fig_pipeline), ("fig_embezzle", fig_embezzle),
    ("fig_sampling", fig_sampling), ("fig_pinsker", fig_pinsker),
    ("fig_h1", fig_h1), ("fig_martingale", fig_martingale),
    ("fig_gsigma", fig_gsigma),
    ("fig_telescope", fig_telescope), ("fig_sizebias", fig_sizebias),
    ("fig_stability", fig_stability), ("fig_waterfall", fig_waterfall),
    ("fig_eps13", fig_eps13),
]


def main():
    checks_chsh()
    checks_feige()
    checks_info()
    checks_resolvent()
    checks_martingale44()
    checks_reverse()
    checks_greedy()
    checks_embezzle()
    checks_constants()
    checks_refine()
    checks_end_to_end()
    for name, fn in FIGS:
        write_fig(OUT, name, fn())
        print(f"FIG {name} written")
    nfail = sum(1 for _, ok in CHECKS if not ok)
    print(f"\n{len(CHECKS)} checks, {nfail} failures, "
          f"{len(FIGS)} figures")
    if nfail:
        sys.exit(1)


if __name__ == "__main__":
    main()
