"""At exceptional heights (Lemma 6.7(iv); the inequality of Theorem 6.9(f) may fail there): which type of
exceptional height actually breaks the inequality?
Types: 'D' (constant death height only), 'side' (axis-parallel bottom/top only), 'mixed'.

Usage:  python exceptional_types.py [budget_s]   (default 420 s; writes out/res_exclab.json)"""
import json
import random
import sys
import time
from collections import defaultdict
from core import Tracer, analyse, exc_set, validate, rss_mb, out_path
import generators as GEN

BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 420


def main():
    T = time.time()
    rng = random.Random(9090)
    gens = [lambda r: GEN.fam_deaths(r, False), lambda r: GEN.fam_floor(r, False), lambda r: GEN.fam_rows(r, False),
            lambda r: GEN.fam_valley(r, False), lambda r: GEN.fam_zigzag(r, False), lambda r: GEN.fam_deaths(r, True)]
    tab = defaultdict(lambda: [0, 0, None])     # type -> [n, n_fail, min margin]
    ex = {}
    ncfg = 0
    while time.time() - T < BUDGET:
        cfg = gens[ncfg % len(gens)](rng)
        ncfg += 1
        ok, _ = validate(cfg.k, cfg.sq)
        if not ok:
            continue
        P = Tracer(cfg).run().final
        exc = exc_set(P, cfg.hF)
        for e, labs in exc.items():
            kinds = {('D' if lab[0] == 'D' else 'side') for lab in labs}
            ty = 'mixed' if len(kinds) > 1 else kinds.pop()
            r = analyse(cfg, P, e)
            t = tab[ty]
            t[0] += 1
            if r['margin'] < 0:
                t[1] += 1
                if ty not in ex:
                    ex[ty] = {'cfg': cfg.name, 'y': str(e), 'margin': str(r['margin']), 'B': str(r['B']),
                              'U': str(r['U']), 'labels': sorted(str(l) for l in labs)}
            if t[2] is None or r['margin'] < t[2]:
                t[2] = r['margin']
    out = {ty: {'n': v[0], 'n_fail': v[1], 'min_margin': float(v[2])} for ty, v in tab.items()}
    print('configs', ncfg, json.dumps(out), 'time', round(time.time() - T, 1), 'rss', round(rss_mb(), 1))
    for ty, e in ex.items():
        print('example', ty, e)
    with open(out_path('res_exclab.json'), 'w', encoding='utf-8') as fo:
        json.dump({'configs': ncfg, 'table': out, 'examples': ex}, fo, indent=1)


if __name__ == '__main__':
    main()
