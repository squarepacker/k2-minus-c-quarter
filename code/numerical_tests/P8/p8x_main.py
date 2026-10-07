# -*- coding: utf-8 -*-
"""
p8x_main.py - recomputation of the constants of Sections 10-12 of the manuscript (Table 1, Table 2,
Table 3, Table 4), of the constraints (C0)-(C7), (Q1), (Q2) and of conditions (i)-(iii) of
Proposition 11.2 for every parameter set, plus independent quadrature checks of (10.1), (10.2), Gamma
and invF, the asymptotic constants of Section 11.4, the integer claims of Corollary A, and a comparison
with the bounds of [R] quoted in Section 1.
Own code (p8x_core.py). 50-digit mpmath + 50-digit interval arithmetic (mpmath.iv).

'stated' values: the decimal strings below are the values as written in the research draft of
Sections 10-12 against which this test was written (2026-10-06).  Most of them coincide with the
manuscript; eleven of them were found to be rounded in the unsafe direction and appear in the
manuscript with safe rounding (see README.md).  A line 'UNSF' refers to the research-draft value.

usage:  python p8x_main.py [--quick]
  --quick : 30 instead of 300 random small exact examples of (10.1); everything else unchanged.
Output: out/p8x_main_out.txt, out/p8x_main_out.json
"""
import os, json, sys, time
from decimal import Decimal
from mpmath import mp, iv, mpf, nstr, quad, log, tan, sqrt, zeta, floor, ceil
from p8x_core import (num, pw, lo, hi, constants, params, derived, at_k,
                      psi_deriv, constraints)

mp.dps = 50
iv.dps = 50
T0 = time.time()
QUICK = '--quick' in sys.argv[1:]
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
LINES = []
JS = {'discrepancies': [], 'unsafe_rounding': [], 'checks': 0, 'fails': 0}


def out(*a):
    s = ' '.join(str(x) for x in a)
    LINES.append(s)
    print(s)


def chk(name, ok, info=''):
    JS['checks'] += 1
    if not ok:
        JS['fails'] += 1
    out(('OK   ' if ok else 'FAIL ') + name + ('   ' + info if info else ''))
    return ok


def ulp(s):
    e = Decimal(s).as_tuple().exponent
    return mpf(10) ** e


def stated(name, true, s, direction):
    """compare a 50-digit value with a stated decimal string.
    direction: 'lower' -> stated must be <= true (stated used as a lower bound)
               'upper' -> stated must be >= true
               'approx' -> only closeness"""
    st = mpf(s)
    u = ulp(s)
    diff = true - st
    close = abs(diff) <= u          # truncation or rounding to the shown digits
    near = abs(diff) <= u / 2 * (1 + mpf('1e-30'))
    if direction == 'lower':
        safe = st <= true
    elif direction == 'upper':
        safe = st >= true
    else:
        safe = True
    tag = 'OK  '
    if not close:
        tag = 'DISC'
        JS['discrepancies'].append(dict(name=name, true=nstr(true, 20), stated=s))
    elif not safe:
        tag = 'UNSF'
        JS['unsafe_rounding'].append(dict(name=name, true=nstr(true, 20), stated=s, direction=direction))
    JS['checks'] += 1
    if tag == 'DISC':
        JS['fails'] += 1
    out('%s %-46s true=%-26s stated=%-16s dir=%-6s %s' % (
        tag, name, nstr(true, 18), s, direction, 'nearest' if near else ('within 1 ulp' if close else 'OFF')))
    return tag


VARIANTS = {}
for v in ('exact9', 'round9', 'exact13', 'round13'):
    VARIANTS[v] = (constants(mp, v), constants(iv, v))

# ---------------------------------------------------------------- constants
out('NOTE: "stated" = value as written in the research draft of Sections 10-12 (see README.md);'
    ' mode=%s' % ('quick' if QUICK else 'full'))
out('=' * 100)
out('SECTION C. constants (Section 12.2, Table 1)')
Km9, _ = VARIANTS['exact9']
Km13, _ = VARIANTS['exact13']
stated('A (K_w*=9)', Km9['A'], '10.60661762772', 'approx')
stated("A' (K_w*=9)", Km9['Ap'], '10.60661763773', 'approx')
stated('C_Lambda', Km9['CL'], '1.00000011856', 'approx')
stated('A_13', Km13['A'], '12.74756790533', 'approx')
stated("A'_13", Km13['Ap'], '12.74756791534', 'approx')
stated('C_Lambda^(13)', Km13['CL'], '1.00000017126', 'approx')
Ki9 = VARIANTS['exact9'][1]
Ki13 = VARIANTS['exact13'][1]
chk('interval: A < 10.6067', hi(Ki9['A']) < mpf('10.6067'), nstr(hi(Ki9['A']), 20))
chk("interval: A' < 10.6068", hi(Ki9['Ap']) < mpf('10.6068'), nstr(hi(Ki9['Ap']), 20))
chk('interval: C_L <= 1.0000002', hi(Ki9['CL']) <= mpf('1.0000002'))
chk('interval: A_13 < 12.7476', hi(Ki13['A']) < mpf('12.7476'), nstr(hi(Ki13['A']), 20))
chk("interval: A'_13 < 12.7476", hi(Ki13['Ap']) < mpf('12.7476'), nstr(hi(Ki13['Ap']), 20))
chk('interval: C_L13 <= 1.0000002', hi(Ki13['CL']) <= mpf('1.0000002'))
# expression of Theorem 8.3(c) bounded by A'_13 (last paragraph of Section 12.2)
alt = 4 * sqrt(mpf(13) / mpf('0.32')) / ((2 - mpf('1.5e-6')) * mp.cos(mpf('1.5e-6'))) + mpf('1e-8')
chk("Sec 12.2: 4 sqrt(13/Q*)/((2-bb)cos bb)+1e-8 <= A'_13 (=A13/cos bb+1e-8)", alt <= Km13['Ap'],
    'alt=%s' % nstr(alt, 16))

# elementary inequalities (spot checks): tan x <= 1.000001 x (Step 11 of the proof of Theorem 10.3),
# sec x - 1 <= 0.50001 x^2 for 0 <= x <= 1e-3 (Section 2, used in Section 4), and
# sin a + 1 - cos a <= 1.0001 a for small a (Lemma 9.3(e), Remark 9.4)
x = mpf('1.5e-6')
chk('tan x <= 1.000001 x at x=1.5e-6 (tan x/x increasing)', tan(x) <= mpf('1.000001') * x, nstr(tan(x) / x - 1, 5))
chk('sec x - 1 <= 0.50001 x^2 at x=1e-3', 1 / mp.cos(mpf('1e-3')) - 1 <= mpf('0.50001') * mpf('1e-3') ** 2)
for a in (mpf('1e-30'), mpf('1e-12'), mpf('1e-9'), mpf('1.5e-6')):
    chk('sin a + 1 - cos a <= 1.0001 a at a=%s' % nstr(a, 3), mp.sin(a) + 1 - mp.cos(a) <= mpf('1.0001') * a)
# polynomial identity of the proof of Lemma 10.2 (invF < 1)
import random
random.seed(1)
okp = True
for _ in range(200):
    v = mpf(random.random())
    lhs = 3 * (1 - v ** 4) - 4 * v ** 3 * (1 - v ** 3)
    rhs = (v - 1) ** 2 * (4 * v ** 4 + 8 * v ** 3 + 9 * v ** 2 + 6 * v + 3)
    okp = okp and abs(lhs - rhs) < mpf('1e-45')
chk('identity 3(1-v^4)-4v^3(1-v^3) = (v-1)^2(4v^4+8v^3+9v^2+6v+3) (200 random v)', okp)

# ---------------------------------------------------------------- parameter sets
# 'Thm A'  = explicit part of Theorem A (Section 11.3, Table 2);  'Thm A13' = Theorem A13 (Section 11.5,
# Table 2);  'pair9'/'pair13' = pairs (c*, k0) (Table 3 contains c* = 0.12, 0.14, 0.15, 0.16 for K = 9 and
# c* = 0.12 for K = 13; the pairs with c* = 0.05, 0.08 of the research draft are not in the manuscript);
# 'table9'/'table13' = the rows of Table 4 (lower bounds for W at given k).
PI_STAR = ('0.3078', '0.2754', '42108000000', '5.8e-6', '1.7e-5')
SETS = []
# (label, variant, params, k, cstar, stated dict)
SETS.append(('Thm A pi* k0=4.62e12', 'round9', PI_STAR, '4.62e12', '0.1', {
    'alpha0': ('1.4999829e-6', 'approx'), 'beta0': ('2.96e-9', 'approx'),
    'C3rhs': ('139.43', 'approx'), 'w0': ('0.0473814', 'upper'), 'h0': ('0.9052372', 'lower'),
    'epsy0': ('244226', 'lower'), 'invF': ('0.99998819', 'lower'), 'Gamma': ('2.85e-5', 'upper'),
    'T1': ('0.000633', 'upper'), 'T2': ('19.25690', 'upper'), 'T3': ('19.24924', 'upper'),
    'Gamma2': ('0.0000285', 'upper'), 'B': ('38.50680', 'upper'), 'Kcoef': ('3.57899', 'upper'),
    'OmUB': ('0.379415', 'upper'), 'wall': ('1.35793', 'upper'), 'wall_short': ('1.358', 'upper'),
    'lH': ('5647.4855', 'lower'), 'Psi': ('146.62429', 'lower'), 'target': ('146.60895', 'upper'),
    'D1': ('0.0086430', 'lower'), 'u14D1': ('10.655', 'lower'), 'KW2B': ('0.0044029', 'upper'),
    'rho': ('0.158143', 'approx'), 'ratio': ('0.1000104', 'lower')}))
OTHER9 = [('0.05', '1.51e11', ('0.0989', '0.1561', '4349000000', '1.7e-5', '2.9e-5'), '0.0500115', '0.09805'),
          ('0.08', '1.25e12', ('0.2279', '0.237', '23084000000', '7.6e-6', '1.9e-5'), '0.0800208', '0.14251'),
          ('0.10', '4.62e12', PI_STAR, '0.1000104', '0.15814'),
          ('0.12', '2.17e13', ('0.3674', '0.3009', '60000000000', '5e-6', '1.5e-5'), '0.1200448', '0.16509'),
          ('0.14', '2.01e14', ('0.4086', '0.3173', '74220000000', '4.7e-6', '1.5e-5'), '0.1400019', '0.16767'),
          ('0.15', '1.18e15', ('0.4243', '0.3234', '80030000000', '4.6e-6', '1.4e-5'), '0.1500058', '0.16818'),
          ('0.16', '2.77e16', ('0.4376', '0.3284', '85109000000', '4.5e-6', '1.4e-5'), '0.1600075', '0.16840')]
for cs, k0, p, r, rh in OTHER9:
    SETS.append(('pair9 c*=%s k0=%s' % (cs, k0), 'round9', p, k0, cs, {'ratio': (r, 'lower'), 'rho': (rh, 'approx')}))
TAB9 = [('1e13', ('0.3421', '0.2903', '52015000000', '5.3e-6', '1.6e-5'), '36.5254', '196.77', 197, '0.11065'),
        ('1e14', ('0.3993', '0.3136', '70863000000', '4.8e-6', '1.5e-5'), '33.8081', '426.16', 427, '0.13476'),
        ('1e16', ('0.4345', '0.3273', '83930000000', '4.5e-6', '1.4e-5'), '32.4097', '1575.76', 1576, '0.15757'),
        ('1e20', ('0.4457', '0.3324', '245800000000', '2.7e-6', '9.5e-6'), '31.9102', '16751.28', 16752, '0.16751'),
        ('1e30', ('0.4472', '0.3343', '547300000000000', '5.7e-8', '5.3e-7'), '31.7249', '5363347.52', 5363348, '0.16960')]
for k, p, Bs, Ps, Wn, rs in TAB9:
    SETS.append(('table9 k=%s' % k, 'round9', p, k, None, {'B': (Bs, 'upper'), 'Psi': (Ps, 'lower'),
                                                         'ratio': (rs, 'lower'), 'W': (Wn, 'int')}))
PI13 = ('0.3674', '0.3012', '60000000000', '5e-6', '1.4e-5')
SETS.append(('Thm A13 k0=2.17e13', 'round13', PI13, '2.17e13', '0.1', {
    'B': ('42.30946', 'upper'), 'Psi': ('215.83530', 'lower'), 'target': ('215.83155', 'upper'),
    'D1': ('0.0055791', 'lower'), 'u14D1': ('10.125', 'lower'), 'KW2B': ('0.00524', 'upper'),
    'rho': ('0.13753', 'approx'), 'ratio': ('0.1000017', 'lower')}))
OTHER13 = [('0.05', '3.28e11', ('0.1395', '0.1857', '8652000000', '1.2e-5', '2.3e-5'), '0.0500175'),
           ('0.08', '3.53e12', ('0.2934', '0.2692', '38260000000', '6e-6', '1.6e-5'), '0.0800070'),
           ('0.12', '3.78e14', ('0.4153', '0.3203', '76670000000', '4.6e-6', '1.3e-5'), '0.1200086')]
for cs, k0, p, r in OTHER13:
    SETS.append(('pair13 c*=%s k0=%s' % (cs, k0), 'round13', p, k0, cs, {'ratio': (r, 'lower')}))
TAB13 = [('1e13', ('0.3421', '0.2907', '52015000000', '5.3e-6', '1.4e-5'), '43.8461', '163.91', 164, '0.09217'),
         ('1e14', ('0.3993', '0.314', '70863000000', '4.8e-6', '1.3e-5'), '40.5842', '355.01', 356, '0.11226'),
         ('1e16', ('0.4345', '0.3276', '83930000000', '4.5e-6', '1.3e-5'), '38.9055', '1312.66', 1313, '0.13126'),
         ('1e20', ('0.4458', '0.3325', '192800000000', '3e-6', '9.5e-6'), '38.3369', '13948.51', 13949, '0.13948'),
         ('1e30', ('0.4472', '0.3343', '428700000000000', '6.4e-8', '5.3e-7'), '38.1280', '4462691.45', 4462692, '0.14112')]
for k, p, Bs, Ps, Wn, rs in TAB13:
    SETS.append(('table13 k=%s' % k, 'round13', p, k, None, {'B': (Bs, 'upper'), 'Psi': (Ps, 'lower'),
                                                           'ratio': (rs, 'lower'), 'W': (Wn, 'int')}))

RESULTS = {}
for label, var, p, k, cs, st in SETS:
    out('=' * 100)
    out('SET', label, ' variant', var, ' pi=', p, ' k=', k, ' c*=', cs)
    Km, Ki = VARIANTS[var]
    Pm, Pi = params(mp, *p), params(iv, *p)
    Dm, Di = derived(mp, Pm, Km), derived(iv, Pi, Ki)
    Rm, Ri = at_k(mp, Pm, Km, Dm, k, cs), at_k(iv, Pi, Ki, Di, k, cs)
    # y0 integer
    chk('C0 y0 is a positive integer', Pm['y0'] == floor(Pm['y0']) and Pm['y0'] > 0)
    for (nm, mg) in constraints(iv, Pi, Ki, Di, k):
        strict = nm.startswith('C3') or nm.startswith('C4') or nm.startswith('C5') or nm.startswith('C0')
        okc = lo(mg) > 0 if strict else lo(mg) >= 0
        chk('[iv] ' + nm, okc, 'margin_lo=%s' % nstr(lo(mg), 8))
    vals = dict(alpha0=Dm['alpha0'], beta0=Dm['beta0'], C3rhs=Pm['c0'] * pw(mp, Pm['y0'], 1, 4),
                w0=Dm['w0'], h0=Dm['h0'], epsy0=Dm['epsy0'], invF=Dm['invF'], Gamma=Dm['Gamma'],
                Gamma2=Dm['Gamma'], T1=Dm['T1'], T2=Dm['T2'], T3=Dm['T3'], B=Dm['B'], Kcoef=Dm['Kcoef'],
                OmUB=Rm['OmUB'], wall=Rm['wall'], wall_short=Rm['wall'], lH=Rm['lH'], Psi=Rm['Psi'],
                rho=Dm['rho'], ratio=Rm['ratio'])
    if cs is not None:
        vals.update(target=Rm['target'], D1=Rm['D1'], u14D1=Rm['u14D1'], KW2B=Rm['KW2B'])
    out('  50-digit values: B=%s  Psi=%s  Psi/k^1/4=%s  rho=%s  h0=%s  invF=%s' % (
        nstr(Dm['B'], 15), nstr(Rm['Psi'], 15), nstr(Rm['ratio'], 12), nstr(Dm['rho'], 12),
        nstr(Dm['h0'], 12), nstr(Dm['invF'], 12)))
    out('  interval  Psi in [%s, %s]' % (nstr(lo(Ri['Psi']), 20), nstr(hi(Ri['Psi']), 20)))
    for key, (s, d) in st.items():
        if d == 'int':
            chk('W >= %d  (interval Psi_lo > %d)' % (s, s - 1), lo(Ri['Psi']) > s - 1,
                'Psi_lo=%s ceil=%s' % (nstr(lo(Ri['Psi']), 12), nstr(ceil(Rm['Psi']), 12)))
            chk('stated W equals ceil(Psi) (not understated)', ceil(Rm['Psi']) == s)
        else:
            stated(key, vals[key], s, d)
    if cs is not None:
        chk('[iv] Prop 11.2 (i) Psi(k0) >= c* k0^{1/4}', lo(Ri['Psi']) >= hi(Ri['target']),
            'margin=%s' % nstr(lo(Ri['Psi']) - hi(Ri['target']), 8))
        chk('[iv] Prop 11.2 (ii) D1 > 0', lo(Ri['D1']) > 0, nstr(lo(Ri['D1']), 10))
        chk('[iv] Prop 11.2 (iii) u0^{1/4} D1 >= K_W/(2B)', lo(Ri['u14D1']) >= hi(Ri['KW2B']),
            '%s >= %s' % (nstr(lo(Ri['u14D1']), 8), nstr(hi(Ri['KW2B']), 8)))
        chk('rho(pi) > c* (needed for asymptotics)', Dm['rho'] > mpf(cs))
    RESULTS[label] = dict(B=nstr(Dm['B'], 20), Psi=nstr(Rm['Psi'], 20), Psi_lo_iv=nstr(lo(Ri['Psi']), 20),
                          ratio=nstr(Rm['ratio'], 15), rho=nstr(Dm['rho'], 15), h0=nstr(Dm['h0'], 15),
                          invF=nstr(Dm['invF'], 15), alpha0=nstr(Dm['alpha0'], 15))

# ---------------------------------------------------------------- section Q: quadrature / structural checks at pi*
out('=' * 100)
out('SECTION Q. independent quadrature checks of (10.1), (10.2), Gamma, invF at pi*, k0 (round9)')
Km, Ki = VARIANTS['round9']
Pm = params(mp, *PI_STAR)
Dm = derived(mp, Pm, Km)
c0, c1, y0, eps, om0 = Pm['c0'], Pm['c1'], Pm['y0'], Pm['eps'], Pm['om0']
k0 = mpf('4.62e12')
y1 = k0 / 2 - 3
mp.dps = 30
f_om = lambda t: (lambda s: Pm['c0'] * (s + 2) * tan(Pm['c0'] / sqrt(s)) * s ** mpf(-1.5) * s)(mp.exp(t))
pts = [log(y0) + (log(y1) - log(y0)) * i / 16 for i in range(17)]
OmQ = quad(f_om, pts)
mp.dps = 50
Rm = at_k(mp, Pm, Km, Dm, '4.62e12', '0.1')
chk('Omega_W (quadrature) <= (10.2) bound', OmQ <= Rm['OmUB'], 'quad=%s UB=%s' % (nstr(OmQ, 12), nstr(Rm['OmUB'], 12)))
stated('Omega_W quadrature', OmQ, '0.3794148', 'approx')
stated('Omega_W bound', Rm['OmUB'], '0.3794152', 'upper')
# Gamma by quadrature of 2*0.50001*int_{y0}^inf (A c1 y^{-3/2} + c1^2/delta y^{-9/4}) dy (proof of Lemma 10.5)
A, delta = Km['A'], Km['delta']
mp.dps = 30
gq = 2 * mpf('0.50001') * quad(lambda t: (A * c1 * mp.exp(t) ** mpf(-1.5) + c1 ** 2 / delta * mp.exp(t) ** mpf(-2.25)) * mp.exp(t),
                               [log(y0), log(y0) + 5, log(y0) + 20, log(y0) + 60, log(y0) + 200])
mp.dps = 50
chk('Gamma closed form vs quadrature (rel 1e-20)', abs(gq / Dm['Gamma'] - 1) < mpf('1e-20'), nstr(gq, 15))
# invF inner integral (Step 10 of the proof of Theorem 10.3):
#   int_y^{y/(1-eps)} s^{-3/4}/(eps s + bW) ds >= 4(1-(1-eps)^{3/4})/(3 eps (1+bW/(eps y0))) y^{-3/4}
bW = Km['bW']
e34 = (1 - eps) ** mpf(0.75)
for yy in (y0, 10 * y0, 1000 * y0, (1 - eps) * y1):
    mp.dps = 40
    I = quad(lambda s: s ** mpf(-0.75) / (eps * s + bW), [yy, yy / (1 - eps)])
    R = 4 * (1 - e34) / (3 * eps * (1 + bW / (eps * y0))) * yy ** mpf(-0.75)
    mp.dps = 50
    chk('invF step: inner integral >= bound at y=%s' % nstr(yy, 6), I >= R, 'ratio=%s' % nstr(I / R, 12))

# (10.1) (Step 2 of the proof of Theorem 10.3): integer-part version and exact l(H) via Hurwitz zeta, at pi*, k0
w0, h0 = Dm['w0'], Dm['h0']
j0 = y0
j1 = floor((1 - eps) * y1)
LH_int = 8 * h0 * ((j1 + 1) ** mpf(0.25) - (j0 + 1) ** mpf(0.25))
LH_cont = 8 * h0 * (((1 - eps) * y1) ** mpf(0.25) - (y0 + 1) ** mpf(0.25))
chk('(10.1) integer-part >= continuous', LH_int >= LH_cont, 'int=%s cont=%s' % (nstr(LH_int, 18), nstr(LH_cont, 18)))
stated('l(H) integer-part version', LH_int, '5647.48557862663', 'approx')
stated('l(H) continuous version', LH_cont, '5647.48557862566', 'approx')
# exact l(H_b) for the full cells m=j0..j1-1 (lower bound of l(H_b) itself):
#   sum_m 4((m+1-w0)^{1/4} - (m+w0)^{1/4}) = 4[ S(1-w0) - S(w0) ],  S(c)=sum_{m=j0}^{j1-1}(m+c)^{1/4}
mp.dps = 60
S = lambda c: zeta(mpf(-0.25), j0 + c) - zeta(mpf(-0.25), j1 + c)
cells = 4 * (S(1 - w0) - S(w0))
mp.dps = 50
chk('exact cell sum 2*sum_m int_{m+w0}^{m+1-w0} y^{-3/4} >= (10.1) bound', 2 * cells >= LH_int,
    'exact=%s bound=%s rel.slack=%s' % (nstr(2 * cells, 15), nstr(LH_int, 15), nstr(2 * cells / LH_int - 1, 5)))
# cell inequality at extreme cells
for m in (j0, j1 - 1):
    lhs = 4 * ((m + 1 - w0) ** mpf(0.25) - (m + w0) ** mpf(0.25))
    rhs = h0 * 4 * ((m + 2) ** mpf(0.25) - (m + 1) ** mpf(0.25))
    chk('cell inequality at m=%s' % nstr(m, 15), lhs >= rhs, nstr(lhs / rhs - 1, 5))

# small exact examples of (10.1), random, including odd k and non-integer y0 (j0=ceil y0)
random.seed(7)
bad = 0
NTRIAL = 30 if QUICK else 300
for trial in range(NTRIAL):
    kk = random.randint(60, 4000)
    yy0 = mpf(random.randint(1, 20)) + (mpf(random.random()) if trial % 3 == 0 else 0)
    ee = mpf(random.choice(['0.001', '0.01', '0.05', '0.2']))
    ww = mpf(random.choice(['0.01', '0.1', '0.3', '0.45']))
    yy1 = mpf(kk) / 2 - 3
    top = (1 - ee) * yy1
    if top < yy0 + 1:
        continue
    # exact l(H) = 2 l(H_b), H_b = {y in [yy0, top]: frac(y) in [ww, 1-ww]}
    tot = mpf(0)
    m = int(floor(yy0)) - 1
    while m <= int(top) + 1:
        a, b = max(m + ww, yy0), min(m + 1 - ww, top)
        if b > a:
            tot += 4 * (b ** mpf(0.25) - a ** mpf(0.25))
        m += 1
    lHx = 2 * tot
    jj0, jj1 = ceil(yy0), floor(top)
    bnd_int = 8 * (1 - 2 * ww) * ((jj1 + 1) ** mpf(0.25) - (jj0 + 1) ** mpf(0.25))
    bnd = 8 * (1 - 2 * ww) * (top ** mpf(0.25) - (jj0 + 1) ** mpf(0.25))
    # also reflection check: l(H_t) computed directly
    tot_t = mpf(0)
    lo_t, hi_t = kk - top, kk - yy0
    m = int(floor(lo_t)) - 1
    while m <= int(hi_t) + 1:
        a, b = max(m + ww, lo_t), min(m + 1 - ww, hi_t)
        if b > a:
            tot_t += 4 * ((kk - a) ** mpf(0.25) - (kk - b) ** mpf(0.25))
        m += 1
    if not (lHx >= bnd_int - mpf('1e-40') and bnd_int >= bnd and abs(tot_t - tot) < mpf('1e-35')):
        bad += 1
chk('(10.1) and l(H_t)=l(H_b) on %d random small exact examples (odd/even k, non-integer y0)' % NTRIAL, bad == 0,
    'bad=%d' % bad)

# Step 10 of the proof of Theorem 10.3 (second use of Tonelli's theorem; the constant invF) on a toy instance
# (moderate numbers so quadrature is exact enough)
mp.dps = 30
tc0, tc1, ty0, ty1, te, tbW = mpf('0.4'), mpf('0.3'), mpf(1000), mpf(50000), mpf('0.05'), mpf('2.0002')
Yb = [(mpf(1000), mpf(1500)), (mpf(3000.5), mpf(3001.25)), (mpf(9000), mpf(20000)), (mpf(40000), mpf(47500))]
assert Yb[-1][1] <= (1 - te) * ty1


def meas_inter(a, b):
    return sum(max(mpf(0), min(b, q) - max(a, p)) for (p, q) in Yb)


def lhs_int(s):
    G = ((1 - te) * s) ** mpf(0.75) / tc1 * meas_inter((1 - te) * s, s)
    return G / (te * s + tbW) * tc0 / 2 * s ** mpf(-1.5)


brk = sorted(set([ty0, ty1] + [p for (p, q) in Yb] + [q for (p, q) in Yb] + [p / (1 - te) for (p, q) in Yb]
                 + [q / (1 - te) for (p, q) in Yb]))
brk = [b for b in brk if ty0 <= b <= ty1]
LHS = quad(lhs_int, brk)
te34 = (1 - te) ** mpf(0.75)
tinvF = 4 * te34 * (1 - te34) / (3 * te * (1 + tbW / (te * ty0)))
RHS = tc0 * tinvF / (2 * tc1) * sum(4 * (q ** mpf(0.25) - p ** mpf(0.25)) for (p, q) in Yb)
mp.dps = 50
chk('toy Tonelli II: M >= (c0 invF/(2c1)) int f dl', LHS >= RHS, 'LHS=%s RHS=%s ratio=%s' % (nstr(LHS, 12), nstr(RHS, 12), nstr(LHS / RHS, 10)))

# ---------------------------------------------------------------- section A: asymptotic claims (short; full study in p8x_asym.py)
# c_infty of (11.1) / Theorem A, Lemma 11.4 (two closed forms), Table 1
out('=' * 100)
out('SECTION A. asymptotic constants (Section 11.4, Table 1)')
for var, cstate, lab in (('exact9', '0.16965241', 'c_inf exact A'), ('round9', '0.16965029', 'c_inf rounded A'),
                         ('exact13', '0.14115933', 'c13_inf exact'), ('round13', '0.14115898', 'c13_inf rounded')):
    Km = VARIANTS[var][0]
    lam = 2 * Km['delta'] + mpf('2e-12')
    c0o = sqrt((1 - lam) / mpf('5.0001'))
    cinf = 4 * mpf(2) ** mpf(-0.25) * (1 - lam - mpf('1.00002') * c0o ** 2) * sqrt(c0o) / sqrt(Km['A'] * Km['Ap'])
    closed = mpf(16) * mpf(2) ** mpf(-0.25) / (5 * sqrt(Km['A'] * Km['Ap'])) * (1 - lam) ** mpf(1.25) * mpf('5.0001') ** mpf(-0.25)
    chk('%s: two closed forms agree' % lab, abs(cinf - closed) < mpf('1e-45'))
    stated(lab, cinf, cstate, 'lower')
    if var == 'exact9':
        stated('c0_opt', c0o, '0.44720465', 'approx')
        stated('c_inf (0.169652)', cinf, '0.169652', 'lower')
        # an earlier formula of the research draft (c0 = 1/sqrt5, h0 -> 1 - c0^2); not used in the manuscript
        cdraft = 4 * mpf(2) ** mpf(-0.25) * (1 - mpf(1) / 5) * sqrt(1 / sqrt(mpf(5))) / sqrt(Km['A'] * Km['Ap'])
        stated('draft formula c (c0=1/sqrt5)', cdraft, '0.16965750', 'approx')
    if var == 'exact13':
        stated('c13_inf (0.141159)', cinf, '0.141159', 'lower')

# ---------------------------------------------------------------- section D: integer claims of Corollary A
out('=' * 100)
out('SECTION D. integer claims of Corollary A (Section 11.6, Remark 11.7) and Section 12.5')
chk('10^4*146^4+1 <= 4.62e12 < 10^4*147^4+1', 10 ** 4 * 146 ** 4 + 1 <= 4620000000000 < 10 ** 4 * 147 ** 4 + 1)
chk('10^4*215^4+1 <= 2.17e13 < 10^4*216^4+1', 10 ** 4 * 215 ** 4 + 1 <= 21700000000000 < 10 ** 4 * 216 ** 4 + 1)
t = mpf('0.1') * mpf('4.62e12') ** mpf(0.25)
chk('ceil(0.1 k0^{1/4}) = 147 at k0=4.62e12', ceil(t) == 147, nstr(t, 15))
t = mpf('0.1') * mpf('2.17e13') ** mpf(0.25)
chk('ceil(0.1 k0^{1/4}) = 216 at k0=2.17e13', ceil(t) == 216, nstr(t, 15))
chk('c>=2065 => (c/0.16)^4 >= 2.77e16 (and 2064 fails)', (mpf(2065) / mpf('0.16')) ** 4 >= mpf('2.77e16') > (mpf(2064) / mpf('0.16')) ** 4)
# first c where pair (0.16, 2.77e16) beats pair (0.1, 4.62e12): max(2.77e16,(c/.16)^4) < max(4.62e12, 1e4 c^4 + 1)
cfirst = None
for c in range(1, 5000):
    if max(mpf('2.77e16'), (mpf(c) / mpf('0.16')) ** 4) < max(mpf('4.62e12'), mpf(10) ** 4 * c ** 4 + 1):
        cfirst = c
        break
out('INFO first integer c for which the (0.16, 2.77e16) pair gives a smaller k-threshold than (0.1, 4.62e12):', cfirst)

# ---------------------------------------------------------------- section V: comparison with the bounds of [R]
# quoted in Section 1 of the manuscript ('[R]' columns); '0.033 ln k' is the constant of an earlier version of [R]
out('=' * 100)
out('SECTION V. comparison with [R] (Section 1): W >= max(1, F(k)), F(k)=1.99954(log k-13.06675)/30.418; kappa=0.0353 (W>=0.0353 log k)')
Km, Ki = VARIANTS['round9']
COMP = {}
for k, p in [('4.62e12', PI_STAR)] + [(k, p) for (k, p, *_r) in TAB9]:
    Pm = params(mp, *p)
    Dm = derived(mp, Pm, Km)
    Rm = at_k(mp, Pm, Km, Dm, k)
    kk = mpf(k)
    F = mpf('1.99954') * (log(kk) - mpf('13.06675')) / mpf('30.418')
    kap = mpf('0.0353') * log(kk)
    pub = mpf('0.033') * log(kk)
    out('k=%-8s  [R] F(k)=%-10s 0.0353 ln k=%-10s (W>=%d)   0.033 ln k=%-8s   new Psi=%-14s (W>=%d)  ratio new/old=%s' % (
        k, nstr(F, 6), nstr(kap, 6), int(ceil(max(1, F, kap))), nstr(pub, 5), nstr(Rm['Psi'], 10), int(ceil(Rm['Psi'])),
        nstr(Rm['Psi'] / max(F, kap), 5)))
    COMP[k] = dict(F=nstr(F, 10), kappa=nstr(kap, 10), new=nstr(Rm['Psi'], 12))
stated('[R] F(1e13) quoted 1.1087', mpf('1.99954') * (log(mpf('1e13')) - mpf('13.06675')) / mpf('30.418'), '1.1087', 'approx')
th4 = mpf('13.06675') + mpf('30.418') * 4 / mpf('1.99954')
chk('[R] bound gives s(k^2-4)=k only for log k > 73.92 (exact threshold %s)' % nstr(th4, 8), mpf('73.91') < th4 < mpf('73.92'))
chk('log(4.62e12) >= 29.16 (quoted)', log(mpf('4.62e12')) >= mpf('29.16'), nstr(log(mpf('4.62e12')), 10))
stated('log10(4.62e12) quoted 12.665', log(mpf('4.62e12'), 10), '12.665', 'approx')

out('=' * 100)
out('TOTAL checks=%d  fails=%d  discrepancies=%d  unsafe_roundings=%d   time=%.1fs' % (
    JS['checks'], JS['fails'], len(JS['discrepancies']), len(JS['unsafe_rounding']), time.time() - T0))
for d in JS['discrepancies']:
    out('  DISC', d)
for d in JS['unsafe_rounding']:
    out('  UNSAFE', d)
JS['results'] = RESULTS
JS['comparison'] = COMP
with open(os.path.join(OUTDIR, 'p8x_main_out.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(LINES) + '\n')
with open(os.path.join(OUTDIR, 'p8x_main_out.json'), 'w', encoding='utf-8') as f:
    json.dump(JS, f, indent=1, ensure_ascii=False)
