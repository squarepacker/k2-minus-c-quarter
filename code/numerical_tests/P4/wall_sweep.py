"""wall term sharpness (Proposition 6.10): L_W^{F_s}(y) / (2 y tan alpha_s) on the explicit drift configuration
(one square per wall, tilted toward it, vertex on the floor, touching the wall; delta = 4k).
Predicted (hand computation): per side L_W(y) = (y - cos a) sin a cos a  for cos a < y < cos a + 1/sin a.

Usage:  python wall_sweep.py        (writes out/res_wallsweep.json; a few seconds)"""
import json
import math
import time
from fractions import Fraction as Fr
from core import Sq, Tracer, analyse, identity_checks, exc_set, truncate, validate, rss_mb, out_path
import generators as GEN


def build(a, k):
    den = 10 ** 12
    t = GEN.t_of(a, den)
    tF = t * Fr(1000001, 1000000)
    b = GEN.Builder(k)
    hw = GEN.halfwidth(t)
    b.add(hw, -Sq(hw, 0, t).ymin, t)
    b.add(Fr(k) - hw, -Sq(Fr(k) - hw, 0, -t).ymin, -t)
    return b.cfg(Fr(4 * k), tF, Fr(2 * k - 1, 2), 'wallsweep a=%g k=%d' % (a, k)), t, tF


def main():
    T = time.time()
    rows = []
    worst = 0.0
    nviol = 0
    for k in (6, 11, 20):
        for a in (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7):
            cfg, t, tF = build(a, k)
            ok, msg = validate(cfg.k, cfg.sq)
            if not ok or len(cfg.sq) != 2:
                rows.append({'a': a, 'k': k, 'skip': msg})
                continue
            P0 = Tracer(cfg).run().final
            exc0 = exc_set(P0, cfg.hF)
            for (ts, hs) in ((tF, cfg.hF), (tF, cfg.hF / 2)):
                Ps = truncate(cfg, P0, ts, hs)
                tanas = 2 * ts / (1 - ts * ts)
                best = (0.0, None)
                for i in list(range(1, 40)) + [39.999]:
                    y = hs * Fr(int(i * 1000), 40000)
                    r = analyse(cfg, Ps, y)
                    bad = identity_checks(r, cfg.k)
                    if bad or (y not in exc0 and (r['B'] != 0 or r['margin'] < 0 or r['margin'] != r['U'])):
                        nviol += 1
                    wb = 2 * y * tanas
                    if r['L']['W'] > wb:
                        nviol += 1
                    ratio = float(r['L']['W'] / wb)
                    if ratio > best[0]:
                        best = (ratio, float(y))
                aa = 2 * math.atan(float(t))
                yb = best[1] if best[1] else float(hs)
                pred = 2 * max(0.0, min(yb - math.cos(aa), 1 / math.sin(aa))) * math.sin(aa) * math.cos(aa) \
                    / (2 * yb * float(tanas))
                rows.append({'k': k, 'a': a, 'hs': float(hs), 'max_ratio': best[0], 'at_y': best[1],
                             'pred_at_y': pred})
                worst = max(worst, best[0])
    with open(out_path('res_wallsweep.json'), 'w', encoding='utf-8') as fo:
        json.dump(rows, fo, indent=1)
    for r in rows:
        print(r)
    print('max ratio', worst, 'violations', nviol, 'time', round(time.time() - T, 1), 'rss', round(rss_mb(), 1))


if __name__ == '__main__':
    main()
