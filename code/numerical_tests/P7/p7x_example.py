"""p7x_example.py -- re-trace one stored drift case where F^0 reaches a square of Z(s) but F_s does not; print exact data.
Also: Remark 9.4 (sin a + 1 - cos a <= 1.0001 a, used in Lemma 9.3(e), holds only for small a).
usage: python p7x_example.py [records.jsonl]
  default: results/original/example_case.jsonl (the record of campaign B from which results/original/example_out.txt
  was printed); any output file of p7x_run.py can be given instead (the smallest such case in it is used)."""
import os, json, sys, random
from fractions import Fraction as Fr
import mpmath
from p7x_core import Sq, Flow, mpf_fr
from p7x_meas import Params, trunc_map, fs_pieces, Lines, p2_test

mpmath.mp.dps = 40
f = mpmath.findroot(lambda a: mpmath.sin(a) + 1 - mpmath.cos(a) - mpmath.mpf('1.0001') * a, mpmath.mpf('0.0002'))
print('Remark 9.4: sin a + 1 - cos a <= 1.0001 a holds exactly for a <= %s; at a=pi/4: lhs=%s rhs=%s' % (
    mpmath.nstr(f, 12), mpmath.nstr(mpmath.sin(mpmath.pi / 4) + 1 - mpmath.cos(mpmath.pi / 4), 8),
    mpmath.nstr(mpmath.mpf('1.0001') * mpmath.pi / 4, 8)))

fn = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results',
                                                          'original', 'example_case.jsonl')
best = None
for line in open(fn, encoding='utf-8'):
    r = json.loads(line)
    if 'per_s' not in r or not r.get('sq'):
        continue
    for c in r['per_s']:
        if c.get('F0_contact_on_Z') and (best is None or len(r['sq']) < len(best[0]['sq'])):
            best = (r, c)
if best is None:
    print('no case'); sys.exit()
r, c = best
pp = r['params']
sq = [Sq(i, Fr(a), Fr(b), Fr(t)) for i, (a, b, t) in enumerate(r['sq'])]
P = Params(k=Fr(pp['k']).limit_denominator(1), delta=Fr(pp['delta']).limit_denominator(1000), c0=Fr(pp['c0']).limit_denominator(100),
           c1=Fr(pp['c1']).limit_denominator(100), y0=Fr(pp['y0']).limit_denominator(10), eps=Fr(pp['eps']).limit_denominator(10),
           omega0=Fr(pp['omega0']).limit_denominator(10))
s = Fr(c['s']).limit_denominator(10 ** 9)
fl = Flow(sq, P.k, P.delta, P.alpha_max, P.hmax).run()
tm, _ = trunc_map(fl, P, s)
L = Lines(sq, P)
ivs, pts, bp = L.Yb_W(s)
Z = L.Zset(s, ivs, pts, bp)
print('case', r['name'], 'k', P.k, 'delta', P.delta, 'c0', P.c0, 'c1', P.c1, 'y0', P.y0, 'eps', P.eps, 's~', float(s))
print('alpha_max', mpmath.nstr(P.alpha_max, 10), 'alpha(s)', mpmath.nstr(P.alpha(s), 10), 'w0', float(P.w0))
for zid, yz in Z.items():
    S = sq[zid]
    print('Z id', zid, 'centre', (str(S.cx), str(S.cy)), 't', str(S.t), 'a', mpmath.nstr(S.a, 8), 'ymin', float(S.ymin),
          'ymax', float(S.ymax), 'witness y', float(yz), 'frac', float(yz - int(yz)))
    for cd in fl.cands:
        if cd.sid == zid:
            prior = [(sid, mpmath.nstr(sq[sid].a, 6), bool(tm[sid])) for (sid, idx) in cd.ent]
            print('   F0 contact on Z: x in (%s, %s) ~ (%.6f, %.6f), contact height ~ %.6f, prior R1-won entries (id, tilt, truncates F_s):'
                  % (cd.xlo, cd.xhi, float(cd.xlo), float(cd.xhi), float(cd.q[2] + cd.q[3] * (cd.xlo + cd.xhi) / 2)), prior[:12])
col = [S for S in sq if S.a > 0.05]
if col:
    print('column tilt', mpmath.nstr(col[0].a, 10), 't', str(col[0].t), 'number', len(col))
p2 = p2_test(fl, P, s, tm)
print('live F_s contacts on Z(s):', [z for z in Z if z in p2['live']])
