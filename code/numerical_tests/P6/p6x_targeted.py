# p6x_targeted.py -- hand-built exact configurations aimed at the bound of the exceptional end interval
# in Lemma 8.7 (P6.V) / Remark 8.8 (288 'sharp' cases with delta = 1e-5, alpha_max = 1.5e-6), and merge
# configurations whose cap square Y is a tilt terminator, run with three orders of rules (M before T as in
# Definition 3.14; T before M; R1 replaced by 'leftmost wins').
# usage: python p6x_targeted.py [time_limit_s]   (default 900; the original run finished in about 30 s)
# output: out/targeted_out.json
import os, sys, json, time, math, random
from collections import Counter
from fractions import Fraction as Fr
from p6x_core import Sq, pyth, rot_for_angle, Flow, HALF
import p6x_configs as C
from p6x_run import run_config, validate, gen
from p6x_checks import check_pairs, pointwise_audit

T0 = time.time()
TL = float(sys.argv[1]) if len(sys.argv) > 1 else 900
res_all = []


def sharp_case(k, delta, Ta, tX, phY, wv, epsv, g0, mirror=False):
    """X (tilt 2 atan tX > 0) almost on the floor (lowest vertex at height g0), Y steep (phY) with its
    lowest vertex at o_X + wv u_X + epsv n_X  ->  hull end where t ~ delta - g ~ delta."""
    uX = pyth(tX)
    c, s = uX
    cX = (Fr(3, 2), (c + s) / 2 + g0)
    X = Sq(0, cX, uX)
    oX = (X.c[0] + HALF * X.n[0], X.c[1] + HALF * X.n[1])
    uY = rot_for_angle(phY, 10 ** 9)
    Y0 = Sq(1, (Fr(0), Fr(0)), uY)
    # lowest vertex of Y0 relative to its centre
    lv = min(Y0.V, key=lambda v: v[1])
    P = (oX[0] + wv * X.u[0] + epsv * X.n[0], oX[1] + wv * X.u[1] + epsv * X.n[1])
    Y = Sq(1, (P[0] - lv[0], P[1] - lv[1]), uY)
    sqs = [X, Y]
    if mirror:
        sqs = [Sq(S.id, (k - S.c[0], S.c[1]), (-S.u[0], S.u[1])) for S in sqs]
    return sqs


def run_sharp():
    k = Fr(6)
    delta = Fr(1, 100000)
    Ta = Fr(math.tan(1.5e-6)).limit_denominator(10 ** 12)
    out = []
    for tX in (Fr(7, 10 ** 7), Fr(74, 10 ** 8), Fr(3, 10 ** 7)):
        for phY in (0.775, 0.6, 0.3, math.pi / 4 - 1e-4, math.pi / 4 + 1e-4, -0.5):
            for wv in (Fr(-4999, 10000), Fr(-49999, 100000), Fr(-2, 5), Fr(0)):
                for epsv in (delta / 100, delta / 3):
                    for mirror in (False, True):
                        if time.time() - T0 > TL / 3:
                            return out
                        sqs = sharp_case(k, delta, Ta, tX, phY, wv, epsv, Fr(1, 10 ** 9), mirror)
                        if not validate(sqs, k):
                            out.append(dict(case=(str(tX), phY, str(wv), str(epsv), mirror), invalid=True))
                            continue
                        r, v = run_config('sharp', sqs, k, delta, Ta, nsamp=10, seed=1, tiny=True)
                        f = r['fams']['floor']
                        out.append(dict(case=(str(tX), phY, str(wv), str(epsv), mirror), nviol=len(v),
                                        viol=[repr(x)[:200] for x in v[:5]], pairs=f['pairs'],
                                        exc=f['L6_exc_ratio_max'], sharp=f['L6_sharp_min'], weak=f['L6_weak_min'],
                                        gv=f['L6_gv_max'], L5gap=f['L5_gap_min']))
    return out


def run_merge_T():
    out = []
    i = 0
    for (delta, amax) in ((0.2, 0.1), (0.1, 0.05), (1e-5, 1.5e-6)):
        tiny = delta < 1e-3
        C.CFG['rotden'] = 10 ** 9 if tiny else 10 ** 4
        C.CFG['posden'] = 10 ** 12 if tiny else 10 ** 7
        dl = Fr(1, 100000) if tiny else Fr(delta).limit_denominator(1000)
        Ta = Fr(math.tan(amax)).limit_denominator(10 ** 12 if tiny else 10 ** 6)
        for yt in (0.3, -0.3, 0.775, -0.775, 0.5 * amax):
            for rep in range(4):
                if time.time() - T0 > TL:
                    return out
                i += 1
                rng = random.Random(777000 + i)
                k = rng.choice((8, 10))
                sqs = C.valleys2(rng, k, delta, amax, layers=rng.randint(2, 4), ytilt=yt)
                kk = Fr(k)
                if not validate(sqs, kk):
                    continue
                row = dict(delta=delta, ytilt=yt, k=k, N=len(sqs))
                for var in ('none', 'T_before_M', 'R1_leftmost'):
                    fl = Flow(sqs, kk, dl, Ta, kk / 2 - 1, 'floor')
                    if var != 'none':
                        fl.variant.add(var)
                    fl.run()
                    st, v = check_pairs(fl, tiny=tiny)
                    v = list(v) + [('anom',) + tuple(a) for a in fl.anom]
                    if var == 'none':
                        d, nM, v3 = pointwise_audit(fl, 20, random.Random(i))
                        v += v3
                        row['pw_M'] = nM
                    row[var] = dict(nviol=len(v), kinds=dict(Counter(x[0] for x in v)),
                                    merge=fl.merge_meas, merge_at_T=fl.merge_at_T, pairs=st['pairs'],
                                    exc=st['L6_exc_ratio_max'])
                out.append(row)
    return out


if __name__ == '__main__':
    A = run_sharp()
    B = run_merge_T()
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'targeted_out.json'), 'w') as fh:
        json.dump(dict(sharp=A, mergeT=B), fh, indent=1, default=str)
    nv = sum(1 for a in A if a.get('nviol'))
    print('sharp cases', len(A), 'invalid', sum(1 for a in A if a.get('invalid')), 'with viol', nv)
    ex = [a['exc'] for a in A if a.get('exc') is not None]
    print('max exc ratio', max(ex) if ex else None)
    sh = [a['sharp'] for a in A if a.get('sharp') is not None]
    print('min sharp margin/(delta bbar)', min(sh) if sh else None)
    df = [a['weak'] for a in A if a.get('weak') is not None]
    print('min weak-form (2.1 delta bbar) margin/(delta bbar)', min(df) if df else None)
    gv = [a['gv'] for a in A if a.get('gv') is not None]
    print('max gv/delta', max(gv) if gv else None)
    print('mergeT rows', len(B))
    for var in ('none', 'T_before_M', 'R1_leftmost'):
        print(var, 'viol', sum(r[var]['nviol'] for r in B if var in r),
              'merge', sum(r[var]['merge'] for r in B if var in r),
              'merge_at_T', sum(r[var]['merge_at_T'] for r in B if var in r),
              Counter(k for r in B if var in r for k in r[var]['kinds']))
    print('pw_M', sum(r.get('pw_M', 0) for r in B))
    print('time', time.time() - T0)
