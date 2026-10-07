"""Negative control for Theorem 6.9(d),(f): rule R1 switched off.
Valley X1(-a), X2(+a) on the floor; row Z_L, Z1, Z_R (axis-parallel) above it covering the line l_y
except gaps of 1e-6, so the unused waste U on l_y is tiny.  Without R1 the doubled beams in Z1 give
Ov > 0 and, since margin = U - Ov - B, the final inequality itself should fail.  With R1 it must hold.

Usage:  python negative_control.py        (writes out/res_negctl.json; about 10 s)"""
import json
import sys
import time
from fractions import Fraction as Fr
from core import Tracer, analyse, identity_checks, exc_set, validate, rss_mb, HALF, out_path
import generators as GEN


def build(a, w, gam, delta, k=4):
    b = GEN.Builder(k)
    t = GEN.t_of(a, 10 ** 9)
    hw = GEN.halfwidth(-t)
    b.drop(Fr(k, 2) - hw - w / 2, -t, 0, 0)
    b.drop(Fr(k, 2) + hw + w / 2, t, 0, 0)
    Z1 = b.drop(Fr(k, 2), 0, gam)
    wl = Fr(1, 10 ** 6)
    b.drop(Z1.xmin - wl - HALF, 0, gam)
    b.drop(Z1.xmax + wl + HALF, 0, gam)
    return b.cfg(delta, Fr(1, 3), Fr(k) - HALF, 'negctl a=%g w=%s gam=%s delta=%s' % (a, w, gam, delta))


def main():
    out = []
    T = time.time()
    for a in (0.1, 0.2, 0.3, 0.4, 0.5):
        for w in (Fr(1, 1000), Fr(1, 100)):
            for gam in (Fr(1, 100), Fr(1, 20)):
                for delta in (Fr(3), Fr(1, 2)):
                    cfg = build(a, w, gam, delta)
                    ok, msg = validate(cfg.k, cfg.sq)
                    if not ok or len(cfg.sq) != 5:
                        out.append({'name': cfg.name, 'skip': msg, 'nsq': len(cfg.sq)})
                        continue
                    lo = max(S.ymin for S in cfg.sq[2:])
                    hi = min(S.ymax for S in cfg.sq[2:])
                    res = {'name': cfg.name, 'band': [float(lo), float(hi)]}
                    for r1 in (False, True):
                        P = Tracer(cfg, r1=r1).run().final
                        exc = exc_set(P, cfg.hF)
                        st = {'n': 0, 'n_exc': 0, 'n_fail_nonexc': 0, 'min_margin_nonexc': None,
                              'max_Ov': 0.0, 'id_bad': 0, 'M': float(sum(p.hi - p.lo for p in P if p.kind == 'M'))}
                        for i in range(1, 40):
                            y = lo + (hi - lo) * Fr(i, 40) + Fr(1, 10 ** 11)
                            r = analyse(cfg, P, y)
                            st['n'] += 1
                            if identity_checks(r, cfg.k):
                                st['id_bad'] += 1
                            st['max_Ov'] = max(st['max_Ov'], float(r['Ov']))
                            if y in exc:
                                st['n_exc'] += 1
                                continue
                            if r['margin'] < 0:
                                st['n_fail_nonexc'] += 1
                                if 'example' not in st:
                                    st['example'] = {'y': str(y), 'margin': str(r['margin']), 'Ov': str(r['Ov']),
                                                     'U': str(r['U']), 'B': str(r['B']), 'Sh': str(r['Sh']),
                                                     'LT+Lam': str(r['LT'] + r['Lam'])}
                            if st['min_margin_nonexc'] is None or r['margin'] < st['min_margin_nonexc']:
                                st['min_margin_nonexc'] = r['margin']
                        st['min_margin_nonexc'] = None if st['min_margin_nonexc'] is None else float(st['min_margin_nonexc'])
                        res['R1' if r1 else 'noR1'] = st
                    out.append(res)
    with open(out_path('res_negctl.json'), 'w', encoding='utf-8') as fo:
        json.dump(out, fo, indent=1)
    nfail = sum(1 for o in out if 'noR1' in o and o['noR1']['n_fail_nonexc'] > 0)
    nok = sum(1 for o in out if 'R1' in o and o['R1']['n_fail_nonexc'] == 0 and o['R1']['id_bad'] == 0)
    nbuilt = sum(1 for o in out if 'R1' in o)
    print('built', nbuilt, '| noR1 configs with final-inequality failure at non-exc y:', nfail,
          '| R1 configs with no failure:', nok, '| time', round(time.time() - T, 1), 'rss', round(rss_mb(), 1))
    for o in out:
        if 'R1' in o:
            print(o['name'], 'noR1:', {kk: o['noR1'][kk] for kk in ('n_fail_nonexc', 'min_margin_nonexc', 'max_Ov', 'id_bad')},
                  'R1:', {kk: o['R1'][kk] for kk in ('n_fail_nonexc', 'min_margin_nonexc', 'max_Ov', 'M', 'id_bad')})
        else:
            print(o)
    ex = [o for o in out if 'noR1' in o and 'example' in o['noR1']]
    if ex:
        print('example:', ex[0]['name'], ex[0]['noR1']['example'])


if __name__ == '__main__':
    main()
