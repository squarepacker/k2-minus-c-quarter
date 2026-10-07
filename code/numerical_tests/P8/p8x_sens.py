# -*- coding: utf-8 -*-
"""
p8x_sens.py - slack of Theorem A (pi*, k0=4.62e12; Section 11.3) and Theorem A13 (pi*_13, k0=2.17e13; Section 11.5)
against changes of the assumed constants (compare Section 12.6): C_Lambda (Corollary 7.2), A' (Corollary 8.13),
A ([R, Lemmas 4.13 and 4.15]), b_W (Theorem 9.1), or an extra loss C' W/delta in Theorem 9.1, at FIXED parameters.
Also: smallest k0' (fixed parameters) under a few perturbed scenarios.  Own code; 50 digits.
usage:  python p8x_sens.py        (a few seconds)
Output: out/p8x_sens_out.txt, out/p8x_sens_out.json
"""
import os, json
from mpmath import mp, mpf, nstr, findroot
from p8x_core import constants, params, derived, at_k

mp.dps = 50
LINES = []


def out(*a):
    s = ' '.join(str(x) for x in a)
    LINES.append(s)
    print(s, flush=True)


RES = {}
for label, var, p, k0, cs in (('Thm A', 'round9', ('0.3078', '0.2754', '42108000000', '5.8e-6', '1.7e-5'), '4.62e12', '0.1'),
                              ('Thm A13', 'round13', ('0.3674', '0.3012', '60000000000', '5e-6', '1.4e-5'), '2.17e13', '0.1')):
    K = constants(mp, var)
    P = params(mp, *p)
    D = derived(mp, P, K)
    R = at_k(mp, P, K, D, k0, cs)
    Nnum = R['Psi'] * D['B']           # numerator of (10.3)
    Bmax = Nnum / R['target']
    dB = Bmax - D['B']
    coefT3 = 2 * P['c1'] / (P['c0'] * D['invF'])
    alpha0, delta = D['alpha0'], K['delta']
    out('=' * 90)
    out('%s: B=%s, max B keeping Psi(k0)>=c*k0^1/4: %s  -> slack dB=%s (relative %s)' % (
        label, nstr(D['B'], 12), nstr(Bmax, 12), nstr(dB, 6), nstr(dB / D['B'], 4)))
    out('   max C_Lambda (Cor 7.2)         : %s' % nstr(K['CL'] + dB / (coefT3 * alpha0 / delta), 8))
    out("   max A' (Cor 8.13)              : %s  (now %s)" % (nstr(K['Ap'] + dB / coefT3, 8), nstr(K['Ap'], 8)))
    out('   max A  ([R, L4.13/L4.15])      : %s  (now %s)' % (
        nstr(K['A'] + dB / (1 / (2 * P['c1']) + mpf('1.00002') * 2 * P['c1'] / mp.sqrt(P['y0'])), 8), nstr(K['A'], 8)))
    # b_W: T3 = coefT3*(...) with 1/invF proportional to (1+bW/(eps y0)); also K_W changes (wall) -> solve exactly
    def psi_bw(bw):
        K2 = dict(K)
        K2['bW'] = bw
        D2 = derived(mp, P, K2)
        R2 = at_k(mp, P, K2, D2, k0, cs)
        return R2['Psi'] - R2['target']
    bwmax = findroot(psi_bw, mpf(40))
    out('   max b_W (Thm 9.1)              : %s  (now 2.0002)' % nstr(bwmax, 8))
    # extra loss Lambda_1(s) <= C' W/delta in Theorem 9.1 acts like C_Lambda + C'
    out("   max extra C' in Lambda_1<=C'W/delta : %s" % nstr(dB / (coefT3 * alpha0 / delta), 8))
    # smallest k with Psi >= c* k^{1/4} for fixed parameters under perturbed scenarios
    for scen, mod in (('as stated', {}), ('C_Lambda=2', {'CL': mpf(2)}), ('C_Lambda=3', {'CL': mpf(3)}),
                      ("A'=A'+0.01", {'Ap': K['Ap'] + mpf('0.01')}), ('b_W=4', {'bW': mpf(4)})):
        K2 = dict(K)
        K2.update(mod)
        D2 = derived(mp, P, K2)
        f = lambda kk: at_k(mp, P, K2, D2, kk)['Psi'] - mpf(cs) * kk ** mpf(0.25)
        a, b = mpf(k0) / 2, mpf(k0) * 4
        if f(b) < 0:
            out('   scenario %-12s: fixed parameters fail even at 4*k0' % scen)
            continue
        for _ in range(200):
            m = (a + b) / 2
            if f(m) < 0:
                a = m
            else:
                b = m
        out('   scenario %-12s: crossing k = %s  (holds at k0=%s: %s)' % (scen, nstr(b, 6), k0, f(mpf(k0)) >= 0))
    RES[label] = dict(B=nstr(D['B'], 15), dB=nstr(dB, 8), rel=nstr(dB / D['B'], 6))
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
with open(os.path.join(OUTDIR, 'p8x_sens_out.txt'), 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(LINES) + '\n')
with open(os.path.join(OUTDIR, 'p8x_sens_out.json'), 'w', encoding='utf-8') as fh:
    json.dump(RES, fh, indent=1)
