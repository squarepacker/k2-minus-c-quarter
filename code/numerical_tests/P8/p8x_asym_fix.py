# -*- coding: utf-8 -*-
"""
p8x_asym_fix.py - follow-up to p8x_asym.py (4a): the 25-digit optimiser there reported sup rho = 0.1696603 > c_inf
at y0=1e30 with eps = 1.28e-22, in apparent contradiction with Proposition 11.6(b).  Suspected artefact of
cancellation in 1-(1-eps)^{3/4} at 25 digits in p8x_asym.py itself (not in the manuscript).
This script
  (i)  re-evaluates the reported point at 25, 50, 80 digits and with 50-digit interval arithmetic,
  (ii) redoes the (4a) maximisation at y0 = 1e12, 1e20, 1e30 and 1e40 with a cancellation-free invF
       (1-(1-eps)^{3/4} = -expm1(0.75*log1p(-eps))) at 60 digits,
  (iii) checks with interval arithmetic (120 digits) that rho(pi) < c_inf at the optimiser's best points.
usage:  python p8x_asym_fix.py        (about 20 s)
Output: out/p8x_asym_fix_out.txt, out/p8x_asym_fix_out.json
"""
import os, json, math
from mpmath import mp, iv, mpf, nstr, sqrt, expm1, log1p, ceil
from p8x_core import constants, params, derived, lo, hi, pw

LINES = []


def out(*a):
    s = ' '.join(str(x) for x in a)
    LINES.append(s)
    print(s, flush=True)


Kr = constants(mp, 'round9')
mp.dps = 50
lam = 2 * Kr['delta'] + mpf('2e-12')
CINF = mpf(16) * mpf(2) ** mpf(-0.25) / (5 * sqrt(Kr['A'] * Kr['Ap'])) * (1 - lam) ** mpf(1.25) * mpf('5.0001') ** mpf(-0.25)
out('c_inf(round9) =', nstr(CINF, 20))

pt = ('0.447204', '0.334366', '1' + '0' * 30, '1.28e-22', '9.11e-16')
for d in (25, 50, 80):
    mp.dps = d
    K = constants(mp, 'round9')
    P = params(mp, *pt)
    D = derived(mp, P, K)
    out('dps=%d: invF=%s rho=%s  rho<c_inf: %s' % (d, nstr(D['invF'], 20), nstr(D['rho'], 15), D['rho'] < CINF))
mp.dps = 50
Ki = constants(iv, 'round9')
Pi = params(iv, *pt)
Di = derived(iv, Pi, Ki)
out('interval(50 digits, outward): invF in [%s, %s], rho in [%s, %s]' % (
    nstr(lo(Di['invF']), 12), nstr(hi(Di['invF']), 12), nstr(lo(Di['rho']), 15), nstr(hi(Di['rho']), 15)))


def rho_stable(c0, c1, y0, eps, om0, K):
    """rho with cancellation-free invF; all args mpf"""
    delta, bW, A, Ap, CL = K['delta'], K['bW'], K['A'], K['Ap'], K['CL']
    alpha0 = c0 / sqrt(y0)
    w0 = mpf('0.50001') * c0 ** 2 * (1 + mpf('1.0001') / y0) + delta + mpf('1.0001') * c1 * ((1 - eps) * y0) ** mpf(-0.75) + mpf('1e-12')
    h0 = 1 - 2 * w0
    one_minus_e34 = -expm1(mpf(0.75) * log1p(-eps))
    e34 = 1 - one_minus_e34
    invF = 4 * e34 * one_minus_e34 / (3 * eps * (1 + bW / (eps * y0)))
    Gamma = mpf('1.00002') * (2 * c1 * A / sqrt(y0) + mpf(4) / 5 * c1 ** 2 / delta * y0 ** mpf(-1.25))
    B = 1 / (om0 * y0 ** mpf(0.75)) + A / (2 * c1) + 2 * c1 / (c0 * invF) * (Ap + alpha0 * CL / delta) + Gamma
    return 8 * mpf(2) ** mpf(-0.25) * h0 * (1 - om0) * (1 - eps) ** mpf(0.25) / B, invF, h0


def nelder_mead(f, x0, step, iters=4000, tol=1e-16):
    n = len(x0)
    pts = [list(x0)] + [[x0[j] + (step[j] if j == i else 0) for j in range(n)] for i in range(n)]
    vals = [f(x) for x in pts]
    for it in range(iters):
        order = sorted(range(n + 1), key=lambda i: vals[i])
        pts = [pts[i] for i in order]
        vals = [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) < tol and it > 50:
            break
        cen = [sum(pts[i][j] for i in range(n)) / n for j in range(n)]
        xr = [2 * cen[j] - pts[-1][j] for j in range(n)]
        fr = f(xr)
        if fr < vals[0]:
            xe = [3 * cen[j] - 2 * pts[-1][j] for j in range(n)]
            fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = [cen[j] + 0.5 * (pts[-1][j] - cen[j]) for j in range(n)]
            fc = f(xc)
            if fc < vals[-1]:
                pts[-1], vals[-1] = xc, fc
            else:
                for i in range(1, n + 1):
                    pts[i] = [pts[0][j] + 0.5 * (pts[i][j] - pts[0][j]) for j in range(n)]
                    vals[i] = f(pts[i])
    i = min(range(n + 1), key=lambda i: vals[i])
    return pts[i], vals[i]


RES = {}
for y0e in (12, 20, 30, 40):
    mp.dps = 60
    K = constants(mp, 'round9')
    y0 = mpf(10) ** y0e
    c0max = mpf('1.5e-6') * sqrt(y0) * (1 - mpf('1e-12'))

    def f(x):
        c0 = min(mpf(abs(x[0])), c0max)
        c1 = mpf(abs(x[1]))
        eps = mpf(10) ** mpf(x[2])
        om0 = mpf(10) ** mpf(x[3])
        if not (0 < eps < 0.5 and 0 < om0 < 0.5 and c0 > 0 and c1 > 0):
            return 1e9
        r, invF, h0 = rho_stable(c0, c1, y0, eps, om0, K)
        if h0 <= 0:
            return 1e9
        # return the GAP to c_inf in units of 1e-12 to keep double-precision resolution of the optimiser
        return float((CINF - r) * mpf(10) ** 12)

    best = None
    for st in ([0.4472, 0.3344, -y0e / 2, -y0e / 4], [0.44, 0.33, -10.0, -8.0], [0.447204, 0.334366, -22.0, -15.0]):
        x, v = nelder_mead(f, st, [0.003, 0.003, 1.0, 1.0], iters=5000)
        x, v = nelder_mead(f, x, [0.0003, 0.0003, 0.3, 0.3], iters=5000)
        if best is None or v < best[1]:
            best = (x, v)
    x, v = best
    c0 = min(mpf(abs(x[0])), c0max)
    pt2 = (nstr(c0, 20), nstr(mpf(abs(x[1])), 20), str(int(10 ** y0e)), nstr(mpf(10) ** mpf(x[2]), 20), nstr(mpf(10) ** mpf(x[3]), 20))
    iv.dps = 120
    Pi = params(iv, *pt2)
    Di = derived(iv, Pi, constants(iv, 'round9'))
    iv.dps = 50
    out('y0=1e%d: best rho = c_inf - %s  (c0=%s c1=%s eps=%s om0=%s) ; interval rho_hi=%s < c_inf: %s' % (
        y0e, nstr(mpf(v) / mpf(10) ** 12, 6), nstr(c0, 8), nstr(abs(x[1]), 8), nstr(mpf(10) ** mpf(x[2]), 4),
        nstr(mpf(10) ** mpf(x[3]), 4), nstr(hi(Di['rho']), 18), hi(Di['rho']) < CINF))
    RES['1e%d' % y0e] = dict(gap=nstr(mpf(v) / mpf(10) ** 12, 8), rho_hi=nstr(hi(Di['rho']), 20), below=bool(hi(Di['rho']) < CINF))
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
with open(os.path.join(OUTDIR, 'p8x_asym_fix_out.txt'), 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(LINES) + '\n')
with open(os.path.join(OUTDIR, 'p8x_asym_fix_out.json'), 'w', encoding='utf-8') as fh:
    json.dump(RES, fh, indent=1)
