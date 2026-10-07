# -*- coding: utf-8 -*-
"""asmx_replay2.py -- two-stage replay (each stage < 1 min, < 100 MB):
  stage 'gen'  : replay the RNG of a campaign and pickle the flagged configurations + parameters
  stage 'eval' : for each pickled target, find the scales with 24-point r_shadow > 1 and recompute
                 int_{V_s} Sh dy with 3000 midpoint samples (pointwise Sh >= sum c_Z also checked).
usage: python asmx_replay2.py gen <campaign> <seed> <families> <k_values> <flagged.json> <out.pkl>
       python asmx_replay2.py genpart <campaign> <seed> <families> <k_values> <flagged.json> <out.pkl> <i_from> <i_to> <state.pkl>
       python asmx_replay2.py eval <campaign> <in.pkl> <lo> <hi>          (output: out/replay2_<campaign>_<lo>_<hi>.json)
(campaigns campA and campB are replayed with gen_hill_v1)"""
import os, sys, json, math, random, time, pickle
import numpy as np
from asmx_chain import Params, LineStruct, FamilyFlows, family_links
from asmx_replay import replay_prefix, use_hill_v1
from asmx_core import Packing


def stage_gen(camp, seed, fams, ks, flagged_file, out_file):
    if camp in ('campA', 'campB'):
        use_hill_v1()
    with open(flagged_file, encoding='utf-8') as fh:
        targets = {int(a): b for a, b in json.load(fh)[camp].items()}
    rng = random.Random(seed)
    saved = []
    for i in range(max(targets) + 1):
        fam = fams[i % len(fams)]
        k = rng.choice(ks)
        if fam == 'random' and k > 20:
            k = 20
        res = replay_prefix(fam, k, rng)
        if i in targets and res is not None:
            P, delta, c0, y0, s_grid, variants = res
            cfg = [(s.cx, s.cy, s.phi) for s in P.S]
            saved.append(dict(i=i, fam=fam, k=k, cfg=cfg, delta=delta, c0=c0, y0=y0, s_grid=s_grid,
                              variants=[(p.c1, p.eps, p.omega0) for p in variants], targets=targets[i]))
    with open(out_file, 'wb') as fh:
        pickle.dump(saved, fh)
    print('saved', len(saved))


def stage_eval(camp, in_file, lo, hi):
    with open(in_file, 'rb') as fh:
        saved = pickle.load(fh)
    out = []
    for item in saved[lo:hi]:
        k = item['k']
        P = Packing(k, item['cfg'])
        for (vi, side) in item['targets']:
            c1, eps, om = item['variants'][vi]
            prm = Params(k, item['delta'], item['c0'], item['y0'], c1, eps, om)
            PP = P if side == 'floor' else P.reflect()
            FF = FamilyFlows(PP, k, item['delta'], item['c0'], item['y0'], item['s_grid'], nov=60)
            LS = LineStruct(PP, prm)
            fl = family_links(FF, LS, prm, nypts=24)
            for d, rec in zip(fl['per_s'], FF.S):
                if d.get('r_shadow', 0) <= 1 + 1e-9:
                    continue
                s = d['s']; Z = np.array(LS.Zset(s), dtype=int)
                ya = max((1 - eps) * s - prm.h, 1e-9); yb = min(s + prm.h, s + 2 - 1e-9)
                n = 3000
                ys = ya + (np.arange(n) + 0.5) / n * (yb - ya)
                shint = 0.0; zint = 0.0; mn = 1e9
                for y in ys:
                    r = rec['LE'].at(float(y), want_ov=False)
                    zc = float(r['c'][Z].sum())
                    shint += r['Sh'] * (yb - ya) / n; zint += zc * (yb - ya) / n
                    mn = min(mn, r['Sh'] - zc)
                out.append(dict(cid='%s-%d' % (camp, item['i']), fam=item['fam'], k=k, vi=vi, side=side, s=s,
                                nZ=int(len(Z)), r_shadow_24=d['r_shadow'], r_shadow_3000=len(Z) / shint,
                                sumcZ_3000=zint, step3min_3000=mn, r_sharp=d.get('r_sharp'), r_meas=d['r_meas']))
    print(json.dumps(out))
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'replay2_%s_%d_%d.json' % (camp, lo, hi)), 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1)


def stage_genpart(camp, seed, fams, ks, flagged_file, out_file, i_from, i_to, state_file):
    """same as stage_gen but for i in [i_from, i_to], resuming the RNG state from state_file (if i_from>0)."""
    import os
    if camp in ('campA', 'campB'):
        use_hill_v1()
    with open(flagged_file, encoding='utf-8') as fh:
        targets = {int(a): b for a, b in json.load(fh)[camp].items()}
    rng = random.Random(seed)
    if i_from > 0:
        with open(state_file, 'rb') as fh:
            rng.setstate(pickle.load(fh))
    saved = []
    if os.path.exists(out_file) and i_from > 0:
        with open(out_file, 'rb') as fh:
            saved = pickle.load(fh)
    for i in range(i_from, i_to + 1):
        fam = fams[i % len(fams)]
        k = rng.choice(ks)
        if fam == 'random' and k > 20:
            k = 20
        res = replay_prefix(fam, k, rng)
        if i in targets and res is not None:
            P, delta, c0, y0, s_grid, variants = res
            saved.append(dict(i=i, fam=fam, k=k, cfg=[(s.cx, s.cy, s.phi) for s in P.S], delta=delta, c0=c0, y0=y0,
                              s_grid=s_grid, variants=[(p.c1, p.eps, p.omega0) for p in variants], targets=targets[i]))
    with open(out_file, 'wb') as fh:
        pickle.dump(saved, fh)
    with open(state_file, 'wb') as fh:
        pickle.dump(rng.getstate(), fh)
    print('saved', len(saved), 'up to', i_to)


if __name__ == '__main__':
    if sys.argv[1] == 'genpart':
        stage_genpart(sys.argv[2], int(sys.argv[3]), sys.argv[4].split(','), [int(x) for x in sys.argv[5].split(',')],
                      sys.argv[6], sys.argv[7], int(sys.argv[8]), int(sys.argv[9]), sys.argv[10])
    elif sys.argv[1] == 'gen':
        stage_gen(sys.argv[2], int(sys.argv[3]), sys.argv[4].split(','), [int(x) for x in sys.argv[5].split(',')],
                  sys.argv[6], sys.argv[7])
    else:
        stage_eval(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
