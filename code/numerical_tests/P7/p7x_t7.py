"""p7x_t7.py -- P2 tightness from the gap side: axis-parallel stack with gaps summing to just below delta, then a tilted
square touched near its lowest vertex.  Reports max maxdist(R(Z))/w(s,Z,q) over live F_s contacts (must be <= 1;
the width test of Theorem 4.11 with w(s, Z) of Lemma 9.3(e)).
usage: python p7x_t7.py        output: out/t7.jsonl (appended)"""
import os, json, random
from fractions import Fraction as Fr
from p7x_core import Sq, validate
from p7x_meas import Params
from p7x_eval import evaluate
from p7x_configs import tq, cs

rng = random.Random(9)
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
out = open(os.path.join(OUTDIR, 't7.jsonl'), 'a', encoding='utf-8')
best = None; nviol = 0; n = 0
for c0 in (Fr(1, 50), Fr(1, 20), Fr(1, 10)):
    for dl in (Fr(1, 10), Fr(1, 5), Fr(3, 10)):
        for m in (1, 2, 3, 4):
            for aZ in (0.005, 0.02, 0.05):
                k = 12
                P = Params(k=k, delta=dl, c0=c0, c1=Fr(1, 20), y0=2, eps=Fr(1, 2), omega0=k + 1)
                sq = []; sid = 0
                g = (dl - Fr(2, 10 ** 6)) / m                     # floor gap + (m-1) internal gaps = delta - 2e-6
                y = g
                for i in range(m):
                    sq.append(Sq(sid, Fr(11, 2), y + Fr(1, 2), 0)); sid += 1
                    y += 1 + g
                y -= g
                # tilted square whose lowest vertex is 1e-7 above the stack top, over the stack
                t = tq(aZ, 100000); c, s = cs(t)
                top = y
                BLx = Fr(11, 2) - Fr(2, 5); BLy = top + Fr(1, 10 ** 7)
                sq.append(Sq(sid, BLx + c / 2 - s / 2, BLy + s / 2 + c / 2, t)); sid += 1
                if validate(sq, k):
                    continue
                r = evaluate('T7', sq, P, [P.y0], rng, nys=4, xcheck_n=8, ovg_grid=4)
                n += 1
                c_ = r['per_s'][0]
                w = c_.get('p2_worst')
                if w and (best is None or w[0] > best[0]):
                    best = (w[0], float(c0), float(dl), m, aZ, w)
                nviol += len(c_['p2_viol']) + len(c_['pointwise_fail'])
                out.write(json.dumps(r, default=str) + '\n')
print('T7 configs', n, 'violations', nviol, 'best P2 ratio', best)
out.write(json.dumps({'tag': 'T7_best', 'n': n, 'viol': nviol, 'best': best}) + '\n')
out.close()
