"""p7x_targeted2.py -- fixed T1 (death exactly at contact on a positive-measure set) and T4 (floor/wall contacts).
usage: python p7x_targeted2.py        output: out/targeted2.jsonl (appended)"""
import os, json, random
from fractions import Fraction as Fr
from p7x_core import Sq, validate, trace_point
from p7x_meas import Params
from p7x_eval import evaluate
from p7x_configs import tq, cs

rng = random.Random(778)
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
out = open(os.path.join(OUTDIR, 'targeted2.jsonl'), 'a', encoding='utf-8')


def emit(tag, r):
    r['tag'] = tag
    out.write(json.dumps(r, default=str) + '\n'); out.flush()
    ps = r.get('per_s', [])
    brief = [(c['s'], c['nZ'], round(c['a_ratio'], 4), len(c['b_live_on_Z']), c['d_Tstar_symdiff'],
              c['pointwise_fail'][:2], c['bad_pieces'][:1]) for c in ps]
    print(tag, r.get('invalid', ''), r.get('F0_tot'), 'xbad', r.get('xcheck', {}).get('nbad'), r.get('anom'), brief,
          'nested', r.get('d_Ts_nested'))


for eq in (True, False):
    k = 12
    delta = Fr(1, 10)
    gap = delta if eq else delta - Fr(1, 10 ** 9)
    sq = []; sid = 0
    for i in range(k - 1):
        sq.append(Sq(sid, Fr(1, 2) + i * Fr(10001, 10000), Fr(1, 2), 0)); sid += 1
    for j in range(1, 5):
        for i in range(k - 1):
            # first gap exactly `gap`; later gaps 1e-12 each
            sq.append(Sq(sid, Fr(1, 2) + i * Fr(10001, 10000), Fr(3, 2) + gap + (j - 1) * (1 + Fr(1, 10 ** 12)), 0)); sid += 1
    # decoration: tilted squares above with fractional heights (Z candidates)
    t = tq(0.004, 40000)
    c, s = cs(t)
    for i in range(k - 2):
        sq.append(Sq(sid, (c + s) / 2 + i * Fr(10101, 10000), Fr(11, 2) + Fr(3, 10) + (c + s) / 2, t)); sid += 1
    P = Params(k=k, delta=delta, c0=Fr(1, 10), c1=Fr(3, 100), y0=2, eps=Fr(1, 2), omega0=k + 1)
    print('T1 validate', validate(sq, k)[:3])
    r = evaluate('T1_gap_eq_delta=%s' % eq, sq, P, [P.y0, P.y1, Fr(5, 2)], rng, nys=10, xcheck_n=40)
    emit('T1', r)
    for x in (Fr(3, 10), Fr(57, 10)):
        tp = trace_point(x, sq, k, delta, P.alpha_max, P.hmax)
        print('   point x=%s -> %s g=%s end=%s' % (x, tp[0], float(tp[3]), [float(v) for v in tp[2]]))

# T4: floor / wall / ceiling contacts
k = 10
P = Params(k=k, delta=Fr(1, 10), c0=Fr(3, 10), c1=Fr(3, 100), y0=Fr(3, 2), eps=Fr(1, 2), omega0=k + 1)
sq = []; sid = 0
sq.append(Sq(sid, Fr(1, 2), Fr(1, 2), 0)); sid += 1
sq.append(Sq(sid, k - Fr(1, 2), Fr(1, 2), 0)); sid += 1
sq.append(Sq(sid, Fr(1, 2), k - Fr(1, 2), 0)); sid += 1
t = tq(0.15, 20000); c, s = cs(t)
sq.append(Sq(sid, Fr(3), (c + s) / 2, t)); sid += 1                  # lowest vertex exactly on the floor
sq.append(Sq(sid, (c + s) / 2, Fr(3), -t)); sid += 1                 # touching the left wall
sq.append(Sq(sid, k - (c + s) / 2, Fr(3), t)); sid += 1              # touching the right wall
for i in range(2, 7):
    sq.append(Sq(sid, Fr(1, 2) + i * Fr(11, 10) + Fr(1, 3), Fr(5, 2) + Fr(1, 100), 0)); sid += 1
print('T4 validate', validate(sq, k))
r = evaluate('T4_floor_wall', sq, P, [P.y0, P.y1, Fr(7, 4)], rng, nys=10, xcheck_n=30)
emit('T4', r)
out.close()
