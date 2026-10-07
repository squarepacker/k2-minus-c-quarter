# -*- coding: utf-8 -*-
"""
p8x_asym.py - the asymptotic constant (Section 11.4, Table 1):
 (1) closed form c_inf of (11.1) (exact / rounded constants, K = 9 and K = 13)
 (2) independent 2-D numerical maximisation of the limit functional over (c0, c1) (Lemma 11.4)
 (3) the family pi_n of the proof of Proposition 11.5 (y0 -> infinity): rho(pi) = Psi_infty(pi) values
     ('stated' values 0.169124, 0.169630, 0.169649, 0.16965029 are those of the research draft; the
     manuscript does not list them)
 (4) Proposition 11.6: (4a) sup of rho(pi) over admissible (c0,c1,eps,om0) at fixed y0 (own Nelder-Mead,
     objective evaluated at 25 digits), (4b) random admissible (pi,k): Psi_pi(k) < c_inf k^{1/4}
 (5) c0 = 0.4472 vs constraint (C1) alpha(y0) <= 1.5e-6 (rows k = 1e30 of Table 4)
 (6) exploratory probe of the thresholds k0 (Section 12 does not claim that they are the least ones):
     own Nelder-Mead, verified at 50 digits with intervals
Own code only (p8x_core.py).
CAUTION: in full mode, line (4a) at y0 = 10^30 reports 'False'.  This is an artefact of this script
(cancellation in 1-(1-eps)^{3/4} at 25 digits for eps ~ 1e-22, giving invF > 1); see p8x_asym_fix.py and
README.md.

usage:  python p8x_asym.py [--quick]
  --quick : (4a) only y0 = 10^12 and 10^20 with 600 instead of 3000 Nelder-Mead iterations; (4b) 2000
            instead of 20000 random samples per variant; (6) only k = 4.62e12 (round9), one start,
            1000 instead of 6000 iterations per stage.  Sections (1), (2), (3), (5) unchanged.
Output: out/p8x_asym_out.txt, out/p8x_asym_out.json
"""
import os, json, time, sys, math, random
from mpmath import mp, iv, mpf, nstr, sqrt, log, ceil, floor
from p8x_core import num, pw, lo, hi, constants, params, derived, at_k

mp.dps = 50
iv.dps = 50
T0 = time.time()
QUICK = '--quick' in sys.argv[1:]
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
LINES = []
JS = {}


def out(*a):
    s = ' '.join(str(x) for x in a)
    LINES.append(s)
    print(s)
    sys.stdout.flush()


def golden_max(f, a, b, tol):
    g = (sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    while b - a > tol:
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
    return (a + b) / 2


out('=' * 100)
out('(1)+(2) closed form vs independent numerical maximisation of the limit functional')
for var in ('exact9', 'round9', 'exact13', 'round13'):
    K = constants(mp, var)
    A, Ap, delta = K['A'], K['Ap'], K['delta']
    lam = 2 * delta + mpf('2e-12')
    kap = mpf('1.00002')
    # limit functional: y0->inf, eps->0, eps*y0->inf, om0->0  =>  w0 -> 0.50001 c0^2 + delta + 1e-12, invF->1, Gamma->0,
    #   1/(om0 y0^{3/4}) -> 0 (with om0 = y0^{-1/4}), alpha(y0)->0
    lim = lambda c0, c1: 8 * mpf(2) ** mpf(-0.25) * (1 - lam - kap * c0 ** 2) / (A / (2 * c1) + 2 * c1 * Ap / c0)
    mp.dps = 40
    best_c1 = lambda c0: golden_max(lambda c1: lim(c0, c1), mpf('0.01'), mpf('2'), mpf('1e-30'))
    c0n = golden_max(lambda c0: lim(c0, best_c1(c0)), mpf('0.05'), mpf('0.99'), mpf('1e-30'))
    c1n = best_c1(c0n)
    vn = lim(c0n, c1n)
    mp.dps = 50
    c0o = sqrt((1 - lam) / (5 * kap))
    c1o = sqrt(A * c0o / (4 * Ap))
    closed = mpf(16) * mpf(2) ** mpf(-0.25) / (5 * sqrt(A * Ap)) * (1 - lam) ** mpf(1.25) * mpf('5.0001') ** mpf(-0.25)
    out('%-8s numeric: c0=%s c1=%s value=%s | closed: c0=%s c1=%s c_inf=%s | |diff|=%s' % (
        var, nstr(c0n, 12), nstr(c1n, 12), nstr(vn, 15), nstr(c0o, 12), nstr(c1o, 12), nstr(closed, 15),
        nstr(abs(vn - closed), 3)))
    JS['cinf_' + var] = dict(numeric=nstr(vn, 20), closed=nstr(closed, 20), c0=nstr(c0o, 15), c1=nstr(c1o, 15))
    if var == 'exact9':
        dd = 4 * mpf(2) ** mpf(-0.25) * (1 - mpf(1) / 5) * sqrt(1 / sqrt(mpf(5))) / sqrt(A * Ap)
        out('   draft formula (c0=1/sqrt5, h0->1-c0^2): %s ; relative difference to c_inf: %s' % (
            nstr(dd, 12), nstr(dd / closed - 1, 4)))
        out('   value of the true limit functional at c0=1/sqrt5 (with lambda, kappa): %s' % nstr(lim(1 / sqrt(mpf(5)), sqrt(A / sqrt(mpf(5)) / (4 * Ap))), 12))

out('=' * 100)
out('(3) family pi_n of Prop 11.5: c0=c0_opt, c1=sqrt(A c0/(4A\')), eps=y0^{-1/2}, om0=y0^{-1/4} (round9 constants)')
K = constants(mp, 'round9')
lam = 2 * K['delta'] + mpf('2e-12')
c0o = sqrt((1 - lam) / mpf('5.0001'))
c1o = sqrt(K['A'] * c0o / (4 * K['Ap']))
STATED = {'1e12': '0.169124', '1e16': '0.169630', '1e20': '0.169649', '1e30': '0.16965029', '1e40': '0.16965029'}
fam = {}
prev = None
for e in (12, 14, 16, 18, 20, 25, 30, 40, 60):
    y0 = mpf(10) ** e
    P = dict(c0=c0o, c1=c1o, y0=y0, eps=y0 ** mpf(-0.5), om0=y0 ** mpf(-0.25))
    D = derived(mp, P, K)
    ok_c1 = D['alpha0'] <= mpf('1.5e-6') and D['beta0'] <= mpf('1.5e-6')
    s = STATED.get('1e%d' % e)
    out('y0=1e%-3d rho=%s  stated=%s  alpha(y0)=%s  C1,C2 ok=%s  increasing=%s' % (
        e, nstr(D['rho'], 12), s, nstr(D['alpha0'], 4), ok_c1, (prev is None or D['rho'] > prev)))
    fam['1e%d' % e] = nstr(D['rho'], 15)
    prev = D['rho']
JS['family'] = fam
cinf_r = mpf(16) * mpf(2) ** mpf(-0.25) / (5 * sqrt(K['A'] * K['Ap'])) * (1 - lam) ** mpf(1.25) * mpf('5.0001') ** mpf(-0.25)
out('c_inf(round9) = %s ; last family value below it: %s' % (nstr(cinf_r, 15), prev < cinf_r))

out('=' * 100)
out('(4a) sup over admissible (c0,c1,eps,om0) of rho(pi) at fixed y0 (round9), own Nelder-Mead [Prop 11.6(b)]')
Y0E_LIST = (12, 20) if QUICK else (10.6, 11, 12, 14, 16, 20, 30)
NM_ITERS_4A = 600 if QUICK else 3000
N_RANDOM_4B = 2000 if QUICK else 20000


def nelder_mead(f, x0, step, iters=4000, tol=1e-15):
    n = len(x0)
    pts = [list(x0)]
    for i in range(n):
        x = list(x0)
        x[i] += step[i]
        pts.append(x)
    vals = [f(x) for x in pts]
    for it in range(iters):
        order = sorted(range(n + 1), key=lambda i: vals[i])
        pts = [pts[i] for i in order]
        vals = [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) < tol * (1 + abs(vals[0])) and it > 50:
            break
        cen = [sum(pts[i][j] for i in range(n)) / n for j in range(n)]
        xr = [cen[j] + (cen[j] - pts[-1][j]) for j in range(n)]
        fr = f(xr)
        if fr < vals[0]:
            xe = [cen[j] + 2 * (cen[j] - pts[-1][j]) for j in range(n)]
            fe = f(xe)
            if fe < fr:
                pts[-1], vals[-1] = xe, fe
            else:
                pts[-1], vals[-1] = xr, fr
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


BB = 1.5e-6
supr = {}
for y0e in Y0E_LIST:
    y0 = mpf(10) ** mpf(y0e)
    y0 = ceil(y0)
    c0max = float(mpf('1.5e-6') * sqrt(y0)) * (1 - 1e-12)

    def negrho(x):
        c0 = min(abs(x[0]), c0max)
        c1 = abs(x[1])
        eps = 10 ** x[2]
        om0 = 10 ** x[3]
        if not (0 < eps < 0.5 and 0 < om0 < 0.5 and c0 > 1e-6 and c1 > 1e-6):
            return 1e9
        P = dict(c0=mpf(c0), c1=mpf(c1), y0=y0, eps=mpf(eps), om0=mpf(om0))
        mp.dps = 25
        D = derived(mp, P, K)
        mp.dps = 50
        if D['h0'] <= 0 or D['beta0'] > mpf('1.5e-6') or not (P['c1'] * (1 - P['eps']) ** mpf(-0.75) < P['c0'] * y0 ** mpf(0.25)):
            return 1e9
        return -float(D['rho'])

    best = None
    for start in ([min(0.447, c0max), 0.33, -5, -4], [min(0.3, c0max), 0.27, -4, -3], [min(0.2, c0max), 0.25, -6, -5]):
        x, v = nelder_mead(negrho, start, [0.02, 0.02, 0.5, 0.5], iters=NM_ITERS_4A)
        if best is None or v < best[1]:
            best = (x, v)
    out('y0=10^%-5s sup rho ~ %.10f  (c0=%.6f c1=%.6f eps=%.3g om0=%.3g)  < c_inf=%.10f : %s' % (
        y0e, -best[1], min(abs(best[0][0]), c0max), abs(best[0][1]), 10 ** best[0][2], 10 ** best[0][3], float(cinf_r),
        -best[1] < float(cinf_r)))
    supr[str(y0e)] = -best[1]
JS['sup_rho_at_y0'] = supr

out('(4b) random admissible (pi,k): check Psi_pi(k) < c_inf k^{1/4} (Prop 11.6(a)); also K_w=13')
random.seed(12345)
for var in ('round9', 'round13'):
    Kv = constants(mp, var)
    lamv = 2 * Kv['delta'] + mpf('2e-12')
    cv = mpf(16) * mpf(2) ** mpf(-0.25) / (5 * sqrt(Kv['A'] * Kv['Ap'])) * (1 - lamv) ** mpf(1.25) * mpf('5.0001') ** mpf(-0.25)
    mp.dps = 30
    nt, nviol, maxr = 0, 0, mpf(0)
    for _ in range(N_RANDOM_4B):
        y0 = mpf(int(10 ** random.uniform(4, 30)))
        c0 = mpf(random.uniform(0.001, 0.99)) * min(1, mpf('1.5e-6') * sqrt(y0))
        c1 = mpf(random.uniform(0.001, 2)) * min(1, mpf('1.5e-6') * y0 ** mpf(0.75))
        eps = mpf(10) ** random.uniform(-12, -0.31)
        om0 = mpf(10) ** random.uniform(-12, -0.31)
        P = dict(c0=c0, c1=c1, y0=y0, eps=eps, om0=om0)
        D = derived(mp, P, Kv)
        if D['h0'] <= 0 or not (c1 * (1 - eps) ** mpf(-0.75) < c0 * y0 ** mpf(0.25)) or D['alpha0'] > mpf('1.5e-6'):
            continue
        kmin = 2 * ((y0 + 1) / (1 - eps) + 3)
        k = kmin * mpf(10) ** random.uniform(0, 30)
        R = at_k(mp, P, Kv, D, k)
        nt += 1
        r = R['Psi'] / k ** mpf(0.25)
        if r > maxr:
            maxr = r
        if not r < cv:
            nviol += 1
    mp.dps = 50
    out('   %s: admissible samples=%d  violations of Psi < c_inf k^{1/4}: %d  max Psi/k^{1/4} = %s (c_inf=%s)' % (
        var, nt, nviol, nstr(maxr, 10), nstr(cv, 10)))
    JS['random_' + var] = dict(samples=nt, violations=nviol, max_ratio=nstr(maxr, 12))

out('=' * 100)
out('(5) c0 = 0.4472 vs alpha(y0) <= 1.5e-6  [constraint (C1); rows k=1e30 of Table 4]')
ymin = (mpf('0.4472') / mpf('1.5e-6')) ** 2
out('   y0 must be >= (0.4472/1.5e-6)^2 = %s ; c0_opt=%s needs y0 >= %s' % (
    nstr(ymin, 12), nstr(c0o, 10), nstr((c0o / mpf('1.5e-6')) ** 2, 12)))
for k, y0s in (('1e30 (K_w*=9 table)', '547300000000000'), ('1e30 (K_w=13 table)', '428700000000000')):
    a = mpf('0.4472') / sqrt(mpf(y0s))
    out('   table row k=%s: y0=%s alpha(y0)=%s <= 1.5e-6: %s' % (k, y0s, nstr(a, 6), a <= mpf('1.5e-6')))
out('   at k0=4.62e12: (y0+1)^{1/4} for y0=(c0_opt/1.5e-6)^2 is %s vs ((k0/2-3))^{1/4}=%s' % (
    nstr(((c0o / mpf('1.5e-6')) ** 2 + 1) ** mpf(0.25), 6), nstr((mpf('4.62e12') / 2 - 3) ** mpf(0.25), 6)))

out('=' * 100)
out('(6) probe of the thresholds k0 (not part of the proof): maximise Psi_pi(k)/k^{1/4} at given k')
NM_ITERS_6 = 1000 if QUICK else 6000


def best_ratio(var, kstr, starts):
    Kv = constants(mp, var)
    kk = mpf(kstr)

    def f(x):
        c0, c1 = x[0], x[1]
        if not (0.01 < c0 < 0.9 and 0.01 < c1 < 2):
            return 1e9
        y0 = max(10 ** x[2], (c0 / BB) ** 2 * (1 + 1e-9))
        y0 = mpf(math.ceil(y0))
        eps, om0 = 10 ** x[3], 10 ** x[4]
        if not (0 < eps < 0.5 and 0 < om0 < 0.5):
            return 1e9
        mp.dps = 30
        P = dict(c0=mpf(c0), c1=mpf(c1), y0=y0, eps=mpf(eps), om0=mpf(om0))
        D = derived(mp, P, Kv)
        if D['alpha0'] > mpf('1.5e-6') or (1 - P['eps']) * (kk / 2 - 3) < y0 + 1:
            mp.dps = 50
            return 1e9
        R = at_k(mp, P, Kv, D, kk)
        mp.dps = 50
        return -float(R['ratio'])

    best = None
    for s in starts:
        x, v = nelder_mead(f, s, [0.01, 0.01, 0.05, 0.1, 0.1], iters=NM_ITERS_6, tol=1e-14)
        x, v = nelder_mead(f, x, [0.002, 0.002, 0.01, 0.02, 0.02], iters=NM_ITERS_6, tol=1e-15)
        if best is None or v < best[1]:
            best = (x, v)
    return best


PROBES = [('round9', '4.62e12', '0.1'), ('round9', '4.61e12', '0.1'), ('round9', '4.60e12', '0.1'),
          ('round13', '2.17e13', '0.1'), ('round13', '2.16e13', '0.1')]
if QUICK:
    PROBES = PROBES[:1]
probe_res = {}
for var, kstr, cs in PROBES:
    starts = ([0.3078, 0.2754, math.log10(4.2108e10), math.log10(5.8e-6), math.log10(1.7e-5)],
              [0.3674, 0.3012, math.log10(6e10), math.log10(5e-6), math.log10(1.4e-5)],
              [0.30, 0.28, 10.6, -5.0, -4.6])
    if QUICK:
        starts = starts[:1]
    x, v = best_ratio(var, kstr, starts)
    c0, c1 = x[0], x[1]
    y0 = math.ceil(max(10 ** x[2], (c0 / BB) ** 2 * (1 + 1e-9)))
    eps, om0 = 10 ** x[3], 10 ** x[4]
    # round to decimal strings conservatively and re-verify rigorously (interval arithmetic)
    c0s = '%.6f' % (math.floor(c0 * 1e6) / 1e6)
    y0i = int(math.ceil(max(y0, (float(c0s) / BB) ** 2 * (1 + 1e-9)) / 1000.0) * 1000)
    ps = (c0s, '%.6f' % c1, str(y0i), '%.4g' % eps, '%.4g' % om0)
    Ki = constants(iv, var)
    Pi = params(iv, *ps)
    Di = derived(iv, Pi, Ki)
    Ri = at_k(iv, Pi, Ki, Di, kstr, cs)
    from p8x_core import constraints as _cons
    okC1 = all(lo(mg) > 0 if (nm.startswith('C3') or nm.startswith('C4') or nm.startswith('C5') or nm.startswith('C0'))
               else lo(mg) >= 0 for nm, mg in _cons(iv, Pi, Ki, Di, kstr))
    rig = lo(Ri['Psi']) - hi(Ri['target'])
    out('   %s k=%s: optimiser best Psi/k^1/4 = %.9f ; rounded pi=%s ; interval Psi - c*k^1/4 >= %s ; C1 ok=%s' % (
        var, kstr, -v, ps, nstr(rig, 8), okC1))
    probe_res[var + '_' + kstr] = dict(best=-v, pi=ps, rigorous_margin=nstr(rig, 10), C1=okC1)
JS['probe'] = probe_res

out('time %.1fs  mode=%s' % (time.time() - T0, 'quick' if QUICK else 'full'))
with open(os.path.join(OUTDIR, 'p8x_asym_out.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(LINES) + '\n')
with open(os.path.join(OUTDIR, 'p8x_asym_out.json'), 'w', encoding='utf-8') as f:
    json.dump(JS, f, indent=1)
