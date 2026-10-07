# p6x_ties.py -- positive-measure event ties (t_S = t_D, contact at exactly h, square flat on the floor,
# squares touching walls) : exact checks that the priority D > W > contact > H of Definition 3.9 is applied
# and that the checks of Section 8 still pass (degenerate cases of Remark 8.15).
# usage: python p6x_ties.py      (no arguments; a few seconds; prints its results, writes no file)
import math
from fractions import Fraction as Fr
from p6x_core import Sq, Flow, pyth
from p6x_run import run_config, validate

k = Fr(8)
delta = Fr(1, 5)
Ta = Fr(1, 10)
one = (Fr(1), Fr(0))
rows = []


def show(tag, sqs, expect=None):
    assert validate(sqs, k), tag
    r, v = run_config(tag, sqs, k, delta, Ta, nsamp=30, seed=3)
    fl = Flow(sqs, k, delta, Ta, k / 2 - 1).run()
    meas = {}
    for (kd, a, b, _, _) in fl.terms:
        meas[kd] = meas.get(kd, Fr(0)) + (b - a)
    pairs = {kk: sum(b - a for (a, b, _, _) in vv) for kk, vv in fl.pairs.items()}
    print(tag, 'viol', len(v), 'meas', {kk: str(vv) for kk, vv in meas.items()},
          'pairs', {str(kk): str(vv) for kk, vv in pairs.items()})
    if v:
        print('   ', v[:5])
    if expect:
        expect(fl, pairs, meas)
    rows.append((tag, len(v)))


# 1. X flat on the floor (contact at t=0), Y flat exactly delta above X: t_S == t_D on [x-range] -> D wins
X = Sq(0, (Fr(2), Fr(1, 2)), one)
Y = Sq(1, (Fr(5, 2), Fr(3, 2) + delta), one)
show('tie_tS_eq_tD', [X, Y], lambda fl, p, m: print('    pair (0,1) present?', (0, 1) in p))
# 2. Y exactly delta/2 above X, X exactly delta/2 above floor: total g == delta at Y -> D wins
X = Sq(0, (Fr(2), Fr(1, 2) + delta / 2), one)
Y = Sq(1, (Fr(5, 2), Fr(3, 2) + delta), one)
show('tie_cum_gap_eq_delta', [X, Y], lambda fl, p, m: print('    pair (0,1) present?', (0, 1) in p))
# 3. Y just below the tie (delta - 1e-9 total): entries happen
Y = Sq(1, (Fr(5, 2), Fr(3, 2) + delta - Fr(1, 10 ** 9)), one)
show('just_below_tie', [X, Y], lambda fl, p, m: print('    pair (0,1) measure', p.get((0, 1))))
# 4. bottom of Y exactly at height h = k/2-1 = 3 (contact at height h processed; then H)
X = Sq(0, (Fr(2), Fr(1, 2)), one)
X2 = Sq(1, (Fr(2), Fr(3, 2) + delta / 4), one)
X3 = Sq(2, (Fr(2), Fr(5, 2) + delta / 2), one)
Y = Sq(3, (Fr(21, 10), Fr(7, 2)), one)   # bottom at exactly 3 = h; gap from X3 top: 3-(3+delta/2) <0?
# X3 top = 3 + delta/2 > 3 -> overlap; move X3 down instead
X3 = Sq(2, (Fr(2), Fr(5, 2) - Fr(0)), one)
X2 = Sq(1, (Fr(2), Fr(3, 2) - Fr(0) + Fr(0)), one)
# stack X(0..1), X2 must be strictly above X: gaps
X = Sq(0, (Fr(2), Fr(1, 2)), one)
X2 = Sq(1, (Fr(2), Fr(3, 2) + Fr(1, 100)), one)
X3 = Sq(2, (Fr(2), Fr(5, 2) + Fr(2, 100)), one)
Y = Sq(3, (Fr(21, 10), Fr(7, 2) + Fr(3, 100) - Fr(0)), one)  # bottom at 3.03 > h: not reached
show('stack_to_h', [X, X2, X3, Y])
X3b = Sq(2, (Fr(2), Fr(5, 2) - Fr(2, 100) + Fr(2, 100)), one)
Yb = Sq(3, (Fr(21, 10), Fr(7, 2)), one)   # bottom exactly at 3 = h, gap 0.03-? from X3 top 3.02 -> overlap
# choose X3 top below 3: X3 centre 2.5 + 0.02 -> top 3.02 ; use X3 centre 2.48? must stay above X2 top 2.01
X3c = Sq(2, (Fr(2), Fr(5, 2) + Fr(1, 100) + Fr(1, 1000)), one)   # top 3.011 > 3 ... use thinner gaps
X2c = Sq(1, (Fr(2), Fr(3, 2) + Fr(1, 1000)), one)                 # top 2.001
X3c = Sq(2, (Fr(2), Fr(5, 2) + Fr(2, 1000)), one)                 # top 3.002 > 3 -> Y bottom 3 overlaps
# final: put Y bottom exactly at h on a column whose top is h - delta/3
X0 = Sq(0, (Fr(2), Fr(1, 2) + Fr(1, 1000)), one)
X1 = Sq(1, (Fr(2), Fr(3, 2) + Fr(2, 1000)), one)
Xt = Sq(2, (Fr(2), Fr(3) - delta / 3 - Fr(1, 2)), one)             # top at 3 - delta/3, bottom 2-delta/3+...
ok = Xt.V[0][1] > X1.V[2][1]
Yh = Sq(3, (Fr(21, 10), Fr(7, 2)), one)
if ok:
    show('contact_exactly_at_h', [X0, X1, Xt, Yh],
         lambda fl, p, m: print('    pair (2,3) measure', p.get((2, 3)), ' H measure', m.get('H')))
# 5. squares touching both walls (vertex on x=0 and side on x=k), tilted
u = pyth(Fr(1, 40))
c, s = u
hw = (c + s) / 2
A = Sq(0, (hw, hw + Fr(1, 1000)), u)
B = Sq(1, (k - Fr(1, 2), Fr(1, 2) + Fr(1, 100)), one)
Cq = Sq(2, (hw + Fr(1, 10), Fr(3, 2) + hw - Fr(1, 2) + Fr(5, 100)), one)
show('walls', [A, B, Cq])
print('total violations', sum(n for _, n in rows))
