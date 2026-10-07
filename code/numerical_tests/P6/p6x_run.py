# p6x_run.py -- driver: generate configs, trace floor+ceiling flows exactly, run checks.
#
# usage:  python p6x_run.py <seed> <n_configs> <time_limit_s> <out_file.jsonl> [scaled|tiny]
#   seed        campaign seed (configuration i uses random.Random(seed*100000 + i));
#               seeds < 4 use k in {6,8,10,12}, seeds >= 4 use k in {10,12,16,20}
#   scaled      delta, alpha_max from PARAMS below (delta <= 0.2)            [default]
#   tiny        delta = 1e-5, alpha_max = 1.5e-6 (the values of the manuscript), exact rationals
# The output file is written to the subfolder out/ of this script's folder (appended, one JSON line per
# configuration).  Quick run (about 1-2 min):  python p6x_run.py 3 25 120 quick_scaled_3.jsonl scaled
# Original campaign, e.g.:                     python p6x_run.py 3 100000 2400 camp_scaled_3.jsonl scaled
import os
import sys
import json
import time
import math
import random
from fractions import Fraction as Fr
import mpmath as mp
from p6x_core import Flow, reflect_squares, closed_disjoint
import p6x_configs as C
from p6x_checks import check_pairs, check_L7, check_chain, pointwise_audit


def validate(sqs, k):
    for S in sqs:
        for v in S.V:
            if not (0 <= v[0] <= k and 0 <= v[1] <= k):
                return False
    for i in range(len(sqs)):
        for j in range(i + 1, len(sqs)):
            A, B = sqs[i], sqs[j]
            if A.bb[2] < B.bb[0] - 1e-6 or B.bb[2] < A.bb[0] - 1e-6 or A.bb[3] < B.bb[1] - 1e-6 or B.bb[3] < A.bb[1] - 1e-6:
                continue
            if not closed_disjoint(A.V, B.V):
                return False
    return True


def run_config(tag, sqs, k, delta, Ta, nsamp=60, seed=0, tiny=False):
    t0 = time.time()
    k = Fr(k)
    h = k / 2 - 1
    delta = Fr(delta)
    Ta = Fr(Ta)
    res = dict(tag=tag, k=float(k), N=len(sqs), delta=float(delta), amax=float(math.atan(Ta)))
    if not validate(sqs, k):
        res['invalid'] = True
        return res, []
    W = float(k * k - len(sqs))
    res['W'] = W
    viol = []
    polys = []
    sqsets = {}
    agg = {}
    rng = random.Random(seed)
    for fam, S in (('floor', sqs), ('ceil', reflect_squares(sqs, k))):
        fl = Flow(S, k, delta, Ta, h, fam).run()
        if not fl.conservation:
            viol.append(('conservation', fam))
        for a in fl.anom:
            viol.append(('anomaly', fam) + tuple(a))
        st, v = check_pairs(fl, tiny=tiny, polys_out=polys, k=k)
        viol += [(fam,) + tuple(x) for x in v]
        cnt, mn, v2, Sm = check_chain(fl)
        viol += [(fam,) + tuple(x) for x in v2]
        done, nM, v3 = pointwise_audit(fl, nsamp, rng)
        viol += [(fam,) + tuple(x) for x in v3]
        used = set()
        for (X, Y) in fl.pairs:
            used.add(Y)
            if X >= 0:
                used.add(X)
        sqsets[fam] = used
        meas = {}
        for (kd, a, b, _, _) in fl.terms:
            meas[kd] = meas.get(kd, 0) + float(b - a)
        st['chain_entries'] = cnt
        st['chain_min'] = None if mn is None else float(mn)
        st['Sm'] = float(Sm)
        st['pw_done'] = done
        st['pw_M'] = nM
        st['meas'] = meas
        st['degen_mid'] = fl.degen_mid
        st['merge_meas'] = getattr(fl, 'merge_meas', None)
        st['merge_at_T'] = getattr(fl, 'merge_at_T', None)
        # squares of this family's pairs must reach height <= h (in the family's own frame)
        for sid in used:
            if min(v[1] for v in fl.byid[sid].V) > h:
                viol.append(('family_square_above_h', fam, sid))
        for kk in ('L6_sharp_min', 'L6_weak_min'):
            if st[kk] is not None:
                st[kk] = float(st[kk])
        agg[fam] = st
    tests, v4 = check_L7(polys, sqsets)
    viol += v4
    res['L7_tests'] = tests
    res['L7_polys'] = len(polys)
    res['fams'] = agg
    res['nviol'] = len(viol)
    res['time'] = time.time() - t0
    return res, viol


def gen_half(kind, rng, k, delta, amax):
    lim = k / 2 - 0.2
    if kind == 'pile':
        return C.pile(rng, k, 5 * k, delta, amax, ymax=lim).sqs
    if kind == 'lattice':
        return C.pile(rng, k, 4 * k, delta, amax, mix=(0.2, 0.7, 0.07, 0.03), ymax=lim, lattice=True).sqs
    if kind == 'jam':
        return C.pile(rng, k, 8 * k, delta * 0.3, amax, mix=(0.1, 0.75, 0.1, 0.05), ymax=lim).sqs
    if kind == 'colA':
        return C.columns(rng, k, delta, amax, rng.choice((0.5, 0.9, 1.3)) * amax, alternate=rng.random() < 0.5)
    if kind == 'valley':
        return C.valleys(rng, k, delta, rng.uniform(0.5, 1.5) * amax, rows=rng.randint(2, 6))
    if kind == 'merge':
        return C.valleys2(rng, k, delta, amax, layers=rng.randint(2, 4))
    if kind == 'bridge':
        return C.bridge(rng, k, delta, amax)
    if kind == 'fan':
        return C.fan(rng, k, delta, amax, rows=rng.randint(2, 5))
    if kind == 'lstrip':
        return C.lstrip(rng, k, delta, amax)
    if kind == 'label45':
        return C.pile(rng, k, 5 * k, delta, amax, mix=(0.1, 0.45, 0.1, 0.35), ymax=lim).sqs
    raise ValueError(kind)


def gen(kind, rng, k, delta, amax):
    from p6x_core import Sq
    bot = gen_half(kind, rng, k, delta, amax)
    kind2 = kind if rng.random() < 0.6 else rng.choice(KINDS)
    top = gen_half(kind2, rng, k, delta, amax)
    sqs = list(bot)
    for S in top:
        sqs.append(Sq(len(sqs), (S.c[0], k - S.c[1]), (S.u[0], -S.u[1])))
    return sqs


PARAMS = [(0.2, 0.1), (0.2, 0.03), (0.1, 0.05), (0.05, 0.03), (0.02, 0.01)]
KINDS = ['pile', 'lattice', 'jam', 'colA', 'valley', 'merge', 'fan', 'lstrip', 'label45', 'merge', 'bridge', 'bridge']

if __name__ == '__main__':
    seed0 = int(sys.argv[1])
    nconf = int(sys.argv[2])
    tlimit = float(sys.argv[3])
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, os.path.basename(sys.argv[4]))
    mode = sys.argv[5] if len(sys.argv) > 5 else 'scaled'
    tiny = (mode == 'tiny')
    if tiny:
        C.CFG['rotden'] = 10 ** 9
        C.CFG['posden'] = 10 ** 12
        params = [(Fr(1, 100000), 1.5e-6)]
    else:
        params = [(Fr(d).limit_denominator(1000), a) for (d, a) in PARAMS]
    T0 = time.time()
    tot = 0
    with open(out, 'a', encoding='utf-8') as f:
        for i in range(nconf):
            if time.time() - T0 > tlimit:
                break
            s = seed0 * 100000 + i
            rng = random.Random(s)
            kind = KINDS[i % len(KINDS)]
            dl, amax = params[(i // len(KINDS)) % len(params)]
            delta = float(dl)
            k = rng.choice((6, 8, 10, 12) if seed0 < 4 else (10, 12, 16, 20))
            sqs = gen(kind, rng, k, delta, amax)
            Ta = Fr(math.tan(amax)).limit_denominator(10 ** 6 if not tiny else 10 ** 12)
            try:
                res, viol = run_config('%s_s%d' % (kind, s), sqs, k, dl, Ta, nsamp=40, seed=s, tiny=tiny)
            except Exception as e:
                res, viol = dict(tag='%s_s%d' % (kind, s), error=repr(e)), [('exception', repr(e))]
            res['viol'] = [repr(v)[:300] for v in viol[:20]]
            f.write(json.dumps(res) + '\n')
            f.flush()
            tot += 1
            print(res['tag'], res.get('N'), 'nviol', len(viol), 'pairs',
                  sum(res['fams'][x]['pairs'] for x in res.get('fams', {})), '%.1fs' % res.get('time', 0), flush=True)
    print('done', tot, time.time() - T0)
