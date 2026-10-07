"""p7x_targeted.py -- degenerate / tie configurations (measure-zero-type situations random fuzzing misses).
T1 death exactly at contact on a positive-measure set (gap == delta - g) vs. gap slightly smaller
T2 exits exactly at height s+2 (H at exit) and contacts exactly at height s+2
T3 ramp exactly on the H boundary (frac == w0) placed where drift columns end
T4 squares resting on the floor (t_S = 0) and touching the walls; squares touching the ceiling
T5 P2 tightness: drift column with tilt just below alpha(s) at s = y0, tiny delta (width test of Theorem 4.11
   with w(s, Z) of Lemma 9.3(e))
T6 ceiling version: tilde Z(s) of the mirrored configuration equals the mirror image of Z(s)
(In the run of 2026-10-06 the T1 and T4 configurations built here were rejected as invalid (overlapping squares);
p7x_targeted2.py repeats T1 and T4 with corrected configurations.)
usage: python p7x_targeted.py        output: out/targeted.jsonl (appended)
"""
import os, sys, json, random, time
from fractions import Fraction as Fr
import mpmath
from p7x_core import Sq, Flow, validate, trace_point, mpf_fr
from p7x_meas import Params
from p7x_eval import evaluate
from p7x_configs import cfg_drift, cfg_blocked, tq, bbox, cs, drop

rng = random.Random(777)
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
out = open(os.path.join(OUTDIR, 'targeted.jsonl'), 'a', encoding='utf-8')


def emit(tag, r):
    r['tag'] = tag
    out.write(json.dumps(r, default=str) + '\n'); out.flush()
    ps = r.get('per_s', [])
    brief = [(c['s'], c['nZ'], round(c['a_ratio'], 4), c.get('p2_worst', [None])[0] if c.get('p2_worst') else None,
              len(c['b_live_on_Z']), c['d_Tstar_symdiff'], c['pointwise_fail'][:2], c['bad_pieces'][:1]) for c in ps]
    print(tag, 'inv' if 'invalid' in r else '', r.get('F0_tot'), r.get('xcheck', {}).get('nbad'), r.get('anom'), brief)


# ---------------- T1: grid with a horizontal gap exactly delta between row 0 and row 1
for eq in (True, False):
    k = 12
    delta = Fr(1, 10)
    gap = delta if eq else delta - Fr(1, 10 ** 9)
    sq = []; sid = 0
    for i in range(k - 1):
        sq.append(Sq(sid, Fr(1, 2) + i * Fr(10001, 10000), Fr(1, 2), 0)); sid += 1
    for j in range(1, 6):
        for i in range(k - 1):
            sq.append(Sq(sid, Fr(1, 2) + i * Fr(10001, 10000), Fr(1, 2) + j + gap, 0)); sid += 1
    # tilted squares (Z candidates) at fractional heights near the top (unreachable decoration)
    P = Params(k=k, delta=delta, c0=Fr(1, 10), c1=Fr(3, 100), y0=2, eps=Fr(1, 2), omega0=k + 1)
    r = evaluate('T1_gap_eq_delta=%s' % eq, sq, P, [P.y0, P.y1], rng, nys=10, xcheck_n=40)
    emit('T1', r)
    # point tracer at a few x: dies exactly at contact when eq
    for x in (Fr(3, 10), Fr(57, 10)):
        t = trace_point(x, sq, k, delta, P.alpha_max, P.hmax)
        print('   point x=%s -> %s g=%s end=%s' % (x, t[0], float(t[3]), [float(v) for v in t[2]]))

# ---------------- T2: axis-parallel columns with tops exactly at s+2; tilted T_s squares with contacts near s+2
k = 14
P = Params(k=k, delta=Fr(1, 20), c0=Fr(1, 5), c1=Fr(3, 100), y0=2, eps=Fr(1, 2), omega0=k + 1)
sq = []; sid = 0
for i in range(k - 1):
    for j in range(4):
        sq.append(Sq(sid, Fr(1, 2) + i * Fr(10001, 10000), Fr(1, 2) + j * Fr(1000001, 1000000) if j else Fr(1, 2), 0)); sid += 1
s_vals = []
tops = sorted({S.ymax for S in sq})
for t in tops:
    s = t - 2
    if P.y0 <= s <= P.y1:
        s_vals.append(s)
s_vals += [P.y0, P.y1]
r = evaluate('T2_tops_at_s+2', sq, P, s_vals, rng, nys=8, xcheck_n=30)
emit('T2', r)

# tilted terminators sitting on an integer column so that contacts happen near height s+2 (both sides)
sq2 = list(sq)
sid = len(sq2)
tT = tq(0.2, 20000)
for i in range(0, k - 2, 2):
    T = drop(sq2, sid, Fr(1, 2) + i * Fr(10001, 10000) + Fr(1, 2), tT, Fr(1, 10 ** 4))
    if T.ymax <= k:
        sq2.append(T); sid += 1
ent_heights = sorted({S.ymin for S in sq2[len(sq):]})
s_vals = [h - 2 for h in ent_heights if P.y0 <= h - 2 <= P.y1] + [h - 2 + Fr(1, 10 ** 6) for h in ent_heights if P.y0 <= h - 2 <= P.y1]
r = evaluate('T2b_Ts_contacts_near_s+2', sq2, P, s_vals[:6] + [P.y0], rng, nys=8, xcheck_n=30)
emit('T2b', r)

# ---------------- T3: square with ramp exactly at the H boundary on top of drift columns
for c0, y0, dl in ((Fr(1, 2), Fr(2), Fr(1, 20)), (Fr(3, 5), Fr(3, 2), Fr(1, 50))):
    for ncol in (2, 3, 4, 6, 8):
        P0 = Params(k=40, delta=dl, c0=c0, c1=Fr(1, 100), y0=y0, eps=Fr(3, 10), omega0=41)
        am = float(P0.alpha_max)
        tcol = tq(am * 0.995, 40000)
        sq = cfg_drift(40, tcol, ncol, Fr(1, 10 ** 5), Fr(1, 200), Fr(1, 4000), Fr(20), g0=Fr(0), nZ=1, zdx=Fr(-45, 100))
        Z = sq[-1]
        # shift Z vertically so that its lowest vertex has fractional part exactly w0 (or 1-w0), if this keeps validity
        for target in (P0.w0, 1 - P0.w0):
            import math
            base = Fr(math.floor(Z.ymin))
            for b in (base, base + 1):
                v = b + target
                if v < Z.ymin:
                    continue
                Z2 = Sq(Z.id, Z.cx, Z.cy + (v - Z.ymin), Z.t)
                sqx = sq[:-1] + [Z2]
                if validate(sqx, 40):
                    continue
                ss = [z for z in (Z2.ymin, Z2.ymin + Fr(1, 10 ** 5), Z2.ymax, Z2.ymax - Fr(1, 10 ** 5), P0.y0)
                      if P0.y0 <= z <= P0.y1]
                r = evaluate('T3_rampEdge_ncol=%d_target=%.4f' % (ncol, float(target)), sqx, P0, ss, rng, nys=6, xcheck_n=20)
                emit('T3', r)
                break

# ---------------- T4: floor/wall/ceiling contacts
k = 10
P = Params(k=k, delta=Fr(1, 10), c0=Fr(3, 10), c1=Fr(3, 100), y0=Fr(3, 2), eps=Fr(1, 2), omega0=k + 1)
sq = []; sid = 0
sq.append(Sq(sid, Fr(1, 2), Fr(1, 2), 0)); sid += 1                  # corner (0,0)
sq.append(Sq(sid, k - Fr(1, 2), Fr(1, 2), 0)); sid += 1              # corner (k,0)
sq.append(Sq(sid, Fr(1, 2), k - Fr(1, 2), 0)); sid += 1              # ceiling corner
t = tq(0.15, 20000); c, s = cs(t)
sq.append(Sq(sid, Fr(3), (c + s) / 2, t)); sid += 1                  # tilted, lowest vertex on the floor
sq.append(Sq(sid, (c + s) / 2, Fr(3), -t)); sid += 1                 # tilted touching the left wall
sq.append(Sq(sid, k - (c + s) / 2, Fr(3), t)); sid += 1              # tilted touching the right wall
for i in range(2, 8):
    sq.append(Sq(sid, Fr(1, 2) + i * Fr(11, 10) + Fr(1, 3), Fr(3, 2) + Fr(1, 100), 0)); sid += 1
probs = validate(sq, k)
print('T4 validate', probs)
r = evaluate('T4_floor_wall', sq, P, [P.y0, P.y1], rng, nys=8, xcheck_n=30)
emit('T4', r)

# ---------------- T5: P2 tightness (drift column, tilt just below alpha(s), s = y0)
best = None
for c0 in (Fr(1, 2), Fr(3, 5), Fr(7, 10)):
    for y0 in (Fr(3, 2), Fr(2), Fr(3)):
        for dl in (Fr(1, 100), Fr(1, 50), Fr(1, 20)):
            for ncol in range(1, 7):
                k = 2 * (int(y0) + 3) + 10
                P0 = Params(k=k, delta=dl, c0=c0, c1=Fr(1, 100), y0=y0, eps=Fr(3, 10), omega0=k + 1)
                am = float(P0.alpha_max)
                for frac in (0.98, 0.995, 0.999):
                    tcol = tq(am * frac, 100000)
                    sq = cfg_drift(k, tcol, ncol, Fr(1, 10 ** 6), Fr(1, 10 ** 4), Fr(1, 8000), Fr(k, 2), g0=Fr(0),
                                   nZ=1, zdx=Fr(-499, 1000))
                    if validate(sq, k):
                        continue
                    ss = [P0.y0]
                    r = evaluate('T5_p2tight', sq, P0, ss, rng, nys=4, xcheck_n=10, ovg_grid=6)
                    c = r['per_s'][0]
                    w = c.get('p2_worst')
                    if w and (best is None or w[0] > best[0]):
                        best = (w[0], float(c0), float(y0), float(dl), ncol, frac, w)
                    if c['p2_viol'] or c['b_live_on_Z'] or c['pointwise_fail']:
                        emit('T5_FLAG', r)
print('T5 best P2 ratio', best)
out.write(json.dumps({'tag': 'T5_best', 'best': best}, default=str) + '\n')

# ---------------- T6: ceiling version -- tilde Z(s) on the mirrored configuration == rho Z(s)
from p7x_meas import Lines
from p7x_configs import cfg_corner
nT6 = 0; badT6 = []
cfgs = []
k = 24
sq, rows, bZ = cfg_blocked(k, tq(0.075, 20000), tq(0.008, 40000), Fr(5, 4) + Fr(1, 100), 5, Fr(1, 1000), Fr(1, 1000), Fr(1, 1000))
cfgs.append(('blocked', sq, Params(k=k, delta=Fr(1, 10), c0=Fr(1, 10), c1=Fr(3, 100), y0=2, eps=Fr(1, 2), omega0=k + 1)))
sqc = cfg_corner(20, 3, Fr(1450, 1000), tq(0.006, 20000), Fr(7, 10), nslots=2, slot_shift=1)
cfgs.append(('corner', sqc, Params(k=20, delta=Fr(2, 100), c0=Fr(1, 10), c1=Fr(5, 100), y0=Fr(3, 2), eps=Fr(3, 10), omega0=Fr(1, 2))))
for name, sqa, Pa in cfgs:
    mir = [Sq(S.id, S.cx, Pa.k - S.cy, -S.t) for S in sqa]
    La = Lines(sqa, Pa); Lm = Lines(mir, Pa)
    for s in [Pa.y0 + (Pa.y1 - Pa.y0) * Fr(j, 7) for j in range(8)]:
        ivs, pts, bp = La.Yb_W(s)
        Zb = set(La.Zset(s, ivs, pts, bp))
        Zt = Lm.Zset_ceiling(s)
        nT6 += 1
        if Zb != Zt:
            badT6.append((name, float(s), len(Zb), len(Zt)))
print('T6 ceiling/reflection checks', nT6, 'mismatches', badT6[:5])
out.write(json.dumps({'tag': 'T6', 'n': nT6, 'bad': badT6}, default=str) + '\n')
out.close()
