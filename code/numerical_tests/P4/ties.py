"""handcrafted exact-tie configurations (priority rules D > W > contact > H on whole intervals;
degenerate configurations, cf. Remark 6.16), each run through all checks of run_tests.run_config.

Usage:  python ties.py        (writes out/res_ties.jsonl; a few seconds)"""
import json
import random
import time
from fractions import Fraction as Fr
from core import Cfg, Sq, validate, rss_mb, out_path
from run_tests import run_config
import generators as GEN

F = Fr
H = F(1, 2)


def configs():
    out = []
    d = F(1, 5)
    # T-a: floor gap exactly delta under an axis square (D/contact tie on an interval, death ON the side)
    sq = [(F(1, 5) + H, d + H, 0), (F(3, 2) + H, F(1, 10) + H, 0), (F(3) + H, H, 0)]
    out.append(Cfg(4, sq, d, F(1, 4), F(7, 2), 'tie_a_floorgap=delta'))
    # T-b: floor-lying square then a gap exactly delta to the next one
    sq = [(F(13, 10) + H, H, 0), (F(5, 4) + H, 1 + d + H, 0), (F(1, 10) + H, F(1, 20) + H, 0),
          (F(29, 10) + H, H, 0), (F(29, 10) + H, 1 + d / 2 + H, 0)]
    out.append(Cfg(4, sq, d, F(1, 4), F(7, 2), 'tie_b_gap=delta_after_pass'))
    # T-c: contact exactly at height h_F (contact beats H), squares stacked to reach it
    hF = F(3, 2)
    sq = [(F(1, 2), H, 0), (F(1, 2) + F(1, 1000), hF + H, 0), (F(5, 2), H + F(1, 20), 0),
          (F(5, 2), hF + H, 0)]
    out.append(Cfg(4, sq, F(3, 5), F(1, 4), hF, 'tie_c_contact_at_hF'))
    # T-d: delta == h_F (D/H tie for unobstructed vertical paths)
    sq = [(F(2), H + F(1, 50), 0), (F(2) + F(1, 3), F(3, 2) + H + F(1, 25), GEN.t_of(0.1))]
    out.append(Cfg(4, sq, F(3, 2), F(1, 4), F(3, 2), 'tie_d_delta=hF'))
    # T-e: squares flush with walls / corners, tilted squares with vertex on floor and touching wall
    t = GEN.t_of(0.2, 10 ** 9)
    hw = GEN.halfwidth(t)
    sq = [(H, H, 0), (F(4) - hw, -Sq(F(4) - hw, 0, -t).ymin, -t), (H, F(1) + F(1, 10) + H, 0),
          (F(2), F(1, 10) + H, t / 3)]
    out.append(Cfg(4, sq, F(1, 2), F(1, 4), F(7, 2), 'tie_e_walls_corners'))
    # T-f: two axis squares sharing bottom height, one above at exactly the death height of the gap
    sq = [(F(1, 2) + F(1, 100), H, 0), (F(3, 2) + F(2, 100), H, 0), (F(1) + F(1, 100), 1 + F(3, 10) + H, 0),
          (F(3), F(3, 10) + H, 0), (F(3), F(13, 10) + F(3, 10) + H, 0)]
    out.append(Cfg(4, sq, F(3, 10), F(1, 4), F(7, 2), 'tie_f_death_height_equals_bottom'))
    # T-g: same as T-a but with tilt so the D/contact tie is only at isolated points
    sq = [(F(1, 5) + H, d + H + F(1, 100), GEN.t_of(0.01, 10 ** 9)), (F(2) + H, d + H, 0)]
    out.append(Cfg(4, sq, d, F(1, 4), F(7, 2), 'tie_g_mixed'))
    return out


def main():
    T = time.time()
    rng = random.Random(77)
    res = []
    for cfg in configs():
        ok, msg = validate(cfg.k, cfg.sq)
        if not ok:
            res.append({'name': cfg.name, 'skip': msg})
            print('SKIP', cfg.name, msg)
            continue
        # F_s variants at the tie heights as well
        extra = [(cfg.tF, h) for h in sorted({S.ymin for S in cfg.sq} | {S.ymax for S in cfg.sq} | {cfg.delta})
                 if 0 < h <= cfg.hF]
        o = run_config(cfg, rng, {'tl': 120, 'nrand': 60, 'nexc': 60, 'ntgt': 120, 'ns': 60, 'npt': 120,
                                  'nrare': 60, 'extra_variants': extra, 'many_var': True}, 'ties')
        res.append(o)
        print(cfg.name, 'nviol', o['nviol'], 'exc', o['n_exc'], o['exc_labels'], 'F0', o['F0'], 'point', o['point'],
              'Fs flows', len(o['Fs']), 'Fs exc fail', sum(s['n_exc_fail'] for s in o['Fs']),
              'meas', {k: round(v, 6) for k, v in o['meas'].items()})
        for v in o['viol'][:5]:
            print('   V', v)
    with open(out_path('res_ties.jsonl'), 'w', encoding='utf-8') as fo:
        for o in res:
            fo.write(json.dumps(o, default=str) + '\n')
    print('time', round(time.time() - T, 1), 'rss', round(rss_mb(), 1))


if __name__ == '__main__':
    main()
