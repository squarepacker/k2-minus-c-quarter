# -*- coding: utf-8 -*-
"""
p8x_scan.py - attack on 'Psi_pi(k) >= c* k^{1/4} for all k >= k0' (Proposition 11.2; Sections 11.3, 11.5,
Table 3):
  (a) fine log grid k0..1e40 at 50 digits (Psi, Psi/k^{1/4}, phi', the lower bound for phi' of the proof of
      Proposition 11.2)
  (b) integer resolution k0-50 .. k0+3000
  (c) crossing point k_c (phi(k_c)=0) for the fixed parameters
  (d) rigorous adaptive interval cover of [k0, 1e100] with monotone bounds (mpmath.iv, 50 digits)
  (e) analytic derivative vs central finite difference
for Theorem A (pi*, k0=4.62e12), Theorem A13 (pi*_13, k0=2.17e13), and the other (c*,k0) pairs.
Own code only (p8x_core.py).

usage:  python p8x_scan.py [--quick]
  --quick : only Theorem A and Theorem A13; 2001 instead of 40001 grid points; integers k0-50..k0+300
            instead of k0-50..k0+3000; the rigorous cover (d) is unchanged.
Output: out/p8x_scan_out.txt, out/p8x_scan_out.json
"""
import os, json, time, sys
from mpmath import mp, iv, mpf, nstr, log, floor, ceil
from p8x_core import num, pw, lo, hi, constants, params, derived, at_k, psi_deriv

mp.dps = 50
iv.dps = 50
T0 = time.time()
QUICK = '--quick' in sys.argv[1:]
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
LINES = []


def out(*a):
    s = ' '.join(str(x) for x in a)
    LINES.append(s)
    print(s)
    sys.stdout.flush()


CASES = [
    ('Thm A pi* c*=0.1', 'round9', ('0.3078', '0.2754', '42108000000', '5.8e-6', '1.7e-5'), '4.62e12', '0.1', True),
    ('Thm A13 c*=0.1', 'round13', ('0.3674', '0.3012', '60000000000', '5e-6', '1.4e-5'), '2.17e13', '0.1', True),
    ('pair9 0.05', 'round9', ('0.0989', '0.1561', '4349000000', '1.7e-5', '2.9e-5'), '1.51e11', '0.05', False),
    ('pair9 0.08', 'round9', ('0.2279', '0.237', '23084000000', '7.6e-6', '1.9e-5'), '1.25e12', '0.08', False),
    ('pair9 0.12', 'round9', ('0.3674', '0.3009', '60000000000', '5e-6', '1.5e-5'), '2.17e13', '0.12', False),
    ('pair9 0.14', 'round9', ('0.4086', '0.3173', '74220000000', '4.7e-6', '1.5e-5'), '2.01e14', '0.14', False),
    ('pair9 0.15', 'round9', ('0.4243', '0.3234', '80030000000', '4.6e-6', '1.4e-5'), '1.18e15', '0.15', False),
    ('pair9 0.16', 'round9', ('0.4376', '0.3284', '85109000000', '4.5e-6', '1.4e-5'), '2.77e16', '0.16', False),
    ('pair13 0.05', 'round13', ('0.1395', '0.1857', '8652000000', '1.2e-5', '2.3e-5'), '3.28e11', '0.05', False),
    ('pair13 0.08', 'round13', ('0.2934', '0.2692', '38260000000', '6e-6', '1.6e-5'), '3.53e12', '0.08', False),
    ('pair13 0.12', 'round13', ('0.4153', '0.3203', '76670000000', '4.6e-6', '1.3e-5'), '3.78e14', '0.12', False),
]

SUMMARY = {}
NINT = 301 if QUICK else 3001          # integers k0-50 .. k0+NINT-1
NG_FULL = 2001 if QUICK else 40001     # log-grid points for Theorems A and A13
out('mode=%s' % ('quick' if QUICK else 'full'))
for (label, var, p, k0s, cs, full) in CASES:
    if QUICK and not full:
        continue
    out('=' * 100)
    out('CASE', label, var, p, 'k0=', k0s, 'c*=', cs)
    Km, Ki = constants(mp, var), constants(iv, var)
    Pm, Pi = params(mp, *p), params(iv, *p)
    Dm, Di = derived(mp, Pm, Km), derived(iv, Pi, Ki)
    k0 = mpf(k0s)
    c = mpf(cs)
    phi = lambda k: at_k(mp, Pm, Km, Dm, k)['Psi'] - c * k ** mpf(0.25)
    S = {}
    # (c) crossing point for these fixed parameters (phi(k0) >= 0 expected, phi(k_c) = 0 for some k_c < k0)
    a, b = k0 / 4, k0
    if phi(a) < 0 <= phi(b):
        for _ in range(200):
            m = (a + b) / 2
            if phi(m) < 0:
                a = m
            else:
                b = m
        S['k_cross'] = nstr(b, 15)
        out('  crossing point of Psi_pi(k) = c* k^{1/4} for these fixed parameters: k_c = %s  (k0/k_c = %s)' % (
            nstr(b, 12), nstr(k0 / b, 8)))
    else:
        out('  no crossing found in [k0/4, k0]: phi(k0/4)=%s phi(k0)=%s' % (nstr(phi(a), 8), nstr(phi(b), 8)))
    # (b) integer resolution near k0
    if full:
        mn, kmn, prev, mono = None, None, None, True
        neg_below = []
        for i in range(-50, NINT):
            kk = k0 + i
            v = phi(kk)
            if i >= 0:
                if mn is None or v < mn:
                    mn, kmn = v, kk
                if prev is not None and v <= prev:
                    mono = False
                prev = v
            elif v < 0:
                neg_below.append(i)
        out('  integers k0..k0+%d: min phi = %s at k=k0+%s ; strictly increasing: %s ; integers below k0 with phi<0: %d of 50' % (
            NINT - 1, nstr(mn, 12), nstr(kmn - k0, 6), mono, len(neg_below)))
        S['int_min_phi'] = nstr(mn, 15)
        S['int_monotone'] = mono
    # (a) log grid
    NG = NG_FULL if full else 4001
    L0, L1 = log(k0, 10), mpf(40)
    mn_phi, k_mn, mn_ratio_margin, mono_ratio, min_dphi, viol_bound, negs = None, None, None, True, None, 0, 0
    prev_ratio = None
    D1 = Dm['P'] - c / (4 * mpf(2) ** mpf(0.75))
    KW2B = Dm['KW'] / (2 * Dm['B'])
    for i in range(NG):
        kk = mpf(10) ** (L0 + (L1 - L0) * i / (NG - 1))
        if i == 0:
            kk = k0
        R = at_k(mp, Pm, Km, Dm, kk)
        v = R['Psi'] - c * R['k14']
        if v < 0:
            negs += 1
        if mn_phi is None or v < mn_phi:
            mn_phi, k_mn = v, kk
        rt = R['ratio']
        if prev_ratio is not None and rt <= prev_ratio:
            mono_ratio = False
        prev_ratio = rt
        dphi = psi_deriv(mp, Pm, Dm, kk, cs)
        u = kk / 2 - 3
        lb = (u ** mpf(0.25) * D1 - KW2B) / u
        if not dphi > lb:
            viol_bound += 1
        if min_dphi is None or dphi < min_dphi:
            min_dphi = dphi
    out('  log grid (%d pts, k0..1e40): negatives=%d ; min phi=%s at k=%s ; Psi/k^1/4 strictly increasing: %s ; '
        'min phi\'=%s ; Prop-11.2 derivative lower bound violated at %d pts' % (
            NG, negs, nstr(mn_phi, 10), nstr(k_mn, 8), mono_ratio, nstr(min_dphi, 6), viol_bound))
    out('  Psi/k^1/4 at k0 = %s ; at 1e40 = %s ; rho(pi) = %s' % (
        nstr(at_k(mp, Pm, Km, Dm, k0)['ratio'], 10), nstr(prev_ratio, 10), nstr(Dm['rho'], 10)))
    S.update(grid_negatives=negs, grid_min_phi=nstr(mn_phi, 15), grid_min_at=nstr(k_mn, 10),
             ratio_increasing=mono_ratio, min_dphi=nstr(min_dphi, 8), deriv_bound_violations=viol_bound)
    # (e) analytic derivative vs central difference
    for kk in (k0, 10 * k0, mpf('1e20'), mpf('1e35')):
        h = kk * mpf('1e-12')
        fd = (phi(kk + h) - phi(kk - h)) / (2 * h)
        an = psi_deriv(mp, Pm, Dm, kk, cs)
        out('  derivative check k=%s: analytic=%s  central diff=%s  rel.err=%s' % (
            nstr(kk, 6), nstr(an, 12), nstr(fd, 12), nstr(abs(fd / an - 1), 3)))
    # (d) rigorous interval cover of [k0, 1e100]
    cI = iv.mpf(cs)
    A8 = (1 - Pi['om0']) * 8 * Di['h0']
    y0p = pw(iv, Pi['y0'] + 1, 1, 4)
    KO = Di['Kcoef'] * num(iv, '1.000001') * Pi['c0'] ** 2
    one_m_eps = 1 - Pi['eps']

    def lower_phi(ka, kb):
        ua = iv.mpf(ka) / 2 - 3
        ub = iv.mpf(kb) / 2 - 3
        val = (A8 * (pw(iv, one_m_eps * ua, 1, 4) - y0p) - KO * (iv.log(ub / Pi['y0']) + 2 / Pi['y0'])) / Di['B'] \
            - cI * pw(iv, iv.mpf(kb), 1, 4)
        return lo(val)

    KMAX = mpf('1e100')
    stack = [(k0, KMAX, 0)]
    acc = []
    fail = None
    nev = 0
    while stack:
        ka, kb, dep = stack.pop()
        nev += 1
        if lower_phi(ka, kb) > 0:
            acc.append((ka, kb))
            continue
        if dep > 80:
            fail = (ka, kb)
            break
        mid = mp.sqrt(ka * kb) if kb / ka > 1.001 else (ka + kb) / 2
        stack.append((mid, kb, dep + 1))
        stack.append((ka, mid, dep + 1))
    acc.sort()
    contiguous = acc[0][0] == k0 and acc[-1][1] == KMAX and all(acc[i][1] == acc[i + 1][0] for i in range(len(acc) - 1))
    out('  RIGOROUS interval cover [k0, 1e100]: %s ; pieces=%d ; evaluations=%d ; contiguous=%s ; smallest piece ratio kb/ka=%s' % (
        'PASS' if fail is None else 'FAIL at %s' % str(fail), len(acc), nev, contiguous,
        nstr(min(b / a for a, b in acc), 10)))
    S.update(cover_pass=(fail is None and contiguous), cover_pieces=len(acc))
    SUMMARY[label] = S

out('=' * 100)
out('time %.1fs' % (time.time() - T0))
with open(os.path.join(OUTDIR, 'p8x_scan_out.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(LINES) + '\n')
with open(os.path.join(OUTDIR, 'p8x_scan_out.json'), 'w', encoding='utf-8') as f:
    json.dump(SUMMARY, f, indent=1)
