# p6x_negctl.py -- negative controls / mutation tests: the checkers must FIRE when the rules
# (R1 of Definition 3.14, the order "M before T", R2 of Definition 3.10) are removed or the constructions
# of Lemmas 8.6/8.7 are deliberately broken ('mut_noDtrim': forget the exceptional end interval of
# Remark 8.8; 'mut_widen': widen the hull [w1, w2] by 0.15 on each side).
# usage: python p6x_negctl.py [time_limit_s] [n_configs]   (defaults 300 and 40, as in the original run)
# output: out/negctl_out.json
import os, sys, json, time, math, random
from collections import Counter
from fractions import Fraction as Fr
from p6x_core import Flow, reflect_squares
import p6x_configs as C
from p6x_run import gen, validate
from p6x_checks import check_pairs, pointwise_audit

T0 = time.time()
tlimit = float(sys.argv[1]) if len(sys.argv) > 1 else 300
nconf = int(sys.argv[2]) if len(sys.argv) > 2 else 40
out = {}
cases = Counter()
for i in range(nconf):
    if time.time() - T0 > tlimit:
        break
    s = 555000 + i
    rng = random.Random(s)
    kind = ['merge', 'bridge', 'label45', 'valley', 'merge'][i % 5]
    delta, amax = [(0.2, 0.1), (0.1, 0.05)][i % 2]
    k = rng.choice((8, 10))
    sqs = gen(kind, rng, k, delta, amax)
    dl = Fr(delta).limit_denominator(1000)
    Ta = Fr(math.tan(amax)).limit_denominator(10 ** 6)
    kk = Fr(k)
    if not validate(sqs, kk):
        continue
    for var in ('none', 'no_R1', 'T_before_M', 'no_R2', 'mut_noDtrim', 'mut_widen'):
        fl = Flow(sqs, kk, dl, Ta, kk / 2 - 1, 'floor')
        if var in ('no_R1', 'T_before_M', 'no_R2'):
            fl.variant.add(var)
        fl.run()
        mut = {'mut_noDtrim': 'noDtrim', 'mut_widen': 'widen'}.get(var)
        st, v = check_pairs(fl, mut=mut)
        c = Counter(x[0] for x in v)
        if var in ('none', 'no_R1'):
            d, nM, v3 = pointwise_audit(fl, 20, random.Random(s))
            c.update(x[0] for x in v3)
        cases[var] += 1
        for key, n in c.items():
            out.setdefault(var, Counter())[key] += n
        out.setdefault(var + '_configs_with_violation', Counter())['n'] += (1 if c else 0)
print('configs per variant', dict(cases))
for var, c in out.items():
    print(var, dict(c))
outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(outdir, exist_ok=True)
with open(os.path.join(outdir, 'negctl_out.json'), 'w') as fh:
    json.dump({k: dict(v) for k, v in out.items()} | {'cases': dict(cases)}, fh, indent=1)
print('time', time.time() - T0)
