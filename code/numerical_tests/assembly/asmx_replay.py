# -*- coding: utf-8 -*-
"""asmx_replay.py -- replay campaign configurations (same RNG consumption as asmx_run.run_one) and recompute the
r_shadow diagnostic with a fine quadrature, plus exact pointwise checks, for flagged cases.
usage: python asmx_replay.py <campaign> <seed> <families> <k_values> <flagged.json>
  e.g. python asmx_replay.py campC 13 valley,hill,columns,floortouch,rows,lshape 20,24,28 results/original/flagged.json
  (campaigns campA and campB are replayed with gen_hill_v1)          output: out/replay_<campaign>.json"""
import os, sys, json, math, random, time
import numpy as np
from asmx_core import Packing
from asmx_chain import Params, LineStruct, FamilyFlows, family_links
from asmx_gen import GENS, GENS_ALL, compose
from asmx_run import pick_flow_params, line_variants


def replay_prefix(fam, k, rng):
    eg = (fam == 'nestedEG')
    if eg:
        delta = rng.choice([0.02, 0.03, 0.05]); c0 = rng.choice([0.08, 0.1, 0.15]); y0 = rng.choice([3, 4])
        c1d = rng.choice([0.25, 0.3, 0.35]); epsd = rng.choice([0.3, 0.35])
    else:
        delta, c0, y0 = pick_flow_params(rng, k)
        c1d = rng.choice([0.04, 0.06, 0.1, 0.15, 0.25, 0.35]); epsd = 0.35
    y1 = k / 2 - 3
    top = (1 - epsd) * y1
    gprm = dict(amax=c0 / math.sqrt(y0), ay1=c0 / math.sqrt(y1), beta_lo=c1d * top ** -0.75,
                beta_hi=c1d * y0 ** -0.75, delta=delta, a_eg=0.97 * c1d * top ** -0.75)
    Hh = k / 2.0 - 0.02
    if eg:
        fam2 = 'nestedEG'
    else:
        fam2 = rng.choice(list(GENS.keys())) if rng.random() < 0.5 else fam
        if fam2 == 'random' and k > 20:
            fam2 = 'rows'
    bot = GENS_ALL[fam](k, Hh, rng, gprm)
    topd = GENS_ALL[fam2](k, Hh, rng, gprm)
    P = Packing(k, compose(k, bot, topd))
    ok, why = P.validate()
    if not ok:
        return None
    s_grid = sorted(set([float(y0), float(y1)] + list(np.linspace(y0, y1, 15)[1:-1] + 1e-7 * math.pi) +
                        [rng.uniform(y0, y1) for _ in range(5)]))
    if eg:
        variants = []
        for om in (0.9, 0.95, 0.97, 0.99, float(k)):
            p = Params(k, delta, c0, y0, c1d, epsd, om)
            if not p.ok():
                variants.append(p)
    else:
        variants = line_variants(rng, k, delta, c0, y0)
        for om in (0.97, float(k)):
            p = Params(k, delta, c0, y0, c1d, 0.35, om)
            if not p.ok():
                variants.append(p)
    return P, delta, c0, y0, s_grid, variants


def use_hill_v1():
    """campaigns A and B were started before gen_hill was rewritten: replay them with the old one."""
    import asmx_gen
    asmx_gen.GENS['hill'] = asmx_gen.gen_hill_v1
    asmx_gen.GENS_ALL['hill'] = asmx_gen.gen_hill_v1
    GENS['hill'] = asmx_gen.gen_hill_v1
    GENS_ALL['hill'] = asmx_gen.gen_hill_v1


def main():
    camp = sys.argv[1]; seed = int(sys.argv[2]); fams = sys.argv[3].split(','); ks = [int(x) for x in sys.argv[4].split(',')]
    if camp in ('campA', 'campB'):
        use_hill_v1()
    with open(sys.argv[5], encoding='utf-8') as fh:   # flagged json file {camp: {"i": [[vi, side], ...]}}
        targets = json.load(fh)[camp]
    targets = {int(a): b for a, b in targets.items()}
    rng = random.Random(seed)
    imax = max(targets)
    out = []
    for i in range(imax + 1):
        fam = fams[i % len(fams)]
        k = rng.choice(ks)
        if fam == 'random' and k > 20:
            k = 20
        res = replay_prefix(fam, k, rng)
        if i not in targets or res is None:
            continue
        P, delta, c0, y0, s_grid, variants = res
        for (vi, side) in targets[i]:
            prm = variants[vi]
            PP = P if side == 'floor' else P.reflect()
            FF = FamilyFlows(PP, k, delta, c0, y0, s_grid, nov=100)
            LS = LineStruct(PP, prm)
            fl = family_links(FF, LS, prm, nypts=24)
            fl_f = family_links(FF, LS, prm, nypts=3000)
            for d24, dfi in zip(fl['per_s'], fl_f['per_s']):
                if d24.get('r_shadow', 0) > 1 + 1e-9 or dfi.get('r_shadow', 0) > 1 + 1e-9:
                    out.append(dict(cid='%s-%d' % (camp, i), fam=fam, k=k, vi=vi, side=side, s=d24['s'], nZ=d24['nZ'],
                                    r_shadow_24=d24['r_shadow'], r_shadow_3000=dfi['r_shadow'],
                                    step3min_3000=dfi['step3min'], r_sharp=dfi['r_sharp'], r_meas=dfi['r_meas']))
    print(json.dumps(out, indent=1))
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'replay_%s.json' % camp), 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1)


if __name__ == '__main__':
    main()
