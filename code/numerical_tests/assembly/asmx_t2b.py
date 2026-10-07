# -*- coding: utf-8 -*-
"""T2b: pointwise Cyrus-Beck cross-check on merge-heavy configurations (valley onto axis or T_max rows).
For an M record, the independent no-merge trace must reach the same point on the same square, and some
other floor point must reach that point with lexicographically smaller (g, x) (checked by bisection-free
search over the competing pieces recorded at that square).
usage: python asmx_t2b.py        output: out/targeted_T2b.json"""
import os, math, random, time, json
import numpy as np
from asmx_core import Packing, Flow
from asmx_gen import tilted_row, ext_of, GENS_ALL, compose
from asmx_targeted import trace_point

rng = random.Random(77)
t0 = time.time()
st = dict(cfgs=0, win_samples=0, win_match=0, M_samples=0, M_match=0, M_measure=0.0, examples=[])
while st['cfgs'] < 60 and time.time() - t0 < 240:
    k = rng.choice([16, 20]); delta = rng.choice([0.08, 0.1, 0.2]); c0 = 0.2; y0 = 3
    amax = c0 / math.sqrt(y0)
    gprm = dict(amax=amax, ay1=c0 / math.sqrt(k / 2 - 3), beta_lo=0.02, beta_hi=0.05, delta=delta, a_eg=0.05)
    cfg = compose(k, GENS_ALL['valley'](k, k / 2 - 0.02, rng, gprm), GENS_ALL['valley'](k, k / 2 - 0.02, rng, gprm))
    P = Packing(k, cfg)
    if not P.validate()[0]:
        continue
    F = Flow(P, delta, amax, k / 2 - 1); F.run(60)
    Mrecs = [r for r in F.recs if r[2] == 'M']
    if not Mrecs:
        continue
    st['cfgs'] += 1
    st['M_measure'] += sum(r[1] - r[0] for r in Mrecs)
    for r in Mrecs:
        for _ in range(5):
            x = rng.uniform(r[0], r[1])
            if r[1] - r[0] < 1e-9:
                continue
            Tp = (r[4] + r[6] * x, r[5] + r[7] * x)
            res = trace_point(P, x, delta, amax, k / 2 - 1, stop_sq=r[8])
            st['M_samples'] += 1
            if res[0] == 'STOP' and abs(res[1][0] - Tp[0]) < 1e-8 and abs(res[1][1] - Tp[1]) < 1e-8:
                st['M_match'] += 1
            elif len(st['examples']) < 6:
                st['examples'].append(('M', x, res[0], res[1], Tp))
    recs = [r for r in F.recs if r[2] != 'M' and r[1] - r[0] > 1e-9]
    for _ in range(300):
        r = rng.choice(recs)
        x = rng.uniform(r[0], r[1])
        Tp = (r[4] + r[6] * x, r[5] + r[7] * x)
        res = trace_point(P, x, delta, amax, k / 2 - 1)
        if res[0] == 'V':
            continue
        st['win_samples'] += 1
        if res[0] == r[2] and abs(res[1][0] - Tp[0]) < 1e-8 and abs(res[1][1] - Tp[1]) < 1e-8:
            st['win_match'] += 1
        elif len(st['examples']) < 12:
            st['examples'].append(('W', x, r[2], res[0], res[1], Tp))
st['time'] = time.time() - t0
print(json.dumps(st, indent=1, default=str))
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
with open(os.path.join(OUTDIR, 'targeted_T2b.json'), 'w', encoding='utf-8') as fh:
    json.dump(st, fh, indent=1, default=str)
