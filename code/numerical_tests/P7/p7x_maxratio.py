"""p7x_maxratio.py -- deterministic attempt to push the (a)-ratio |Z(s)| / (|V_s|(N + Lambda0 + wall)) towards 1.
blocked Z-columns with tiny c1 (small b = beta^((1-eps)y0)), tiny c0 (small wall term and E-loss), many rows m,
stack raised so that the H upper limit (1-eps) y1 >= s needs a smaller k.
usage: python p7x_maxratio.py [time_limit_s]   (default 1200; the run of 2026-10-06 took about 3.5 min)
output: out/maxratio.jsonl (appended)"""
import os, sys, json, random, time, math
from fractions import Fraction as Fr
from p7x_core import Sq
from p7x_meas import Params
from p7x_eval import evaluate
from p7x_configs import cfg_blocked, tq, bbox, cs

rng = random.Random(4242)
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
out = open(os.path.join(OUTDIR, 'maxratio.jsonl'), 'a', encoding='utf-8')
t_start = time.time()
TL = float(sys.argv[1]) if len(sys.argv) > 1 else 1200
cases = []
for m in (6, 8, 10, 12):
    for c1 in (Fr(4, 1000), Fr(1, 100)):
        for c0 in (Fr(2, 100), Fr(4, 100)):
            for live in (False, True):
                for raise_stack in (True, False):
                    cases.append((m, c1, c0, live, raise_stack))
rng.shuffle(cases)
cases.sort(key=lambda z: z[0])
best = None
for (m, c1, c0, live, raise_stack) in cases:
    if time.time() - t_start > TL:
        break
    y0 = Fr(2)
    am = float(c0) / math.sqrt(float(y0))
    tT = tq(am * 1.004, 100000)
    while 2 * math.atan(float(tT)) < am * 1.0005:
        tT += Fr(1, 100000)
    sinT = float(cs(tT)[1])
    g0 = Fr(1, 2000)
    delta = Fr(int(1000 * (sinT + 0.0005 + 0.004)) + 1, 1000)
    grow = Fr(1, 10 ** 4); gZ = Fr(1, 10 ** 4); gT = Fr(1, 10 ** 4)
    D_est = (m - 2) * 1.0 + 0.01
    A_target = max(2.3, D_est) if raise_stack else 2.3
    s_est = A_target + D_est
    aZ = 0.6 * float(c1) * s_est ** -0.75
    tZ = tq(aZ, 10 ** 6)
    bZ = bbox(tZ)
    c, s_ = cs(tZ); sca = abs(s_) * c
    # first Z row: its top ramp just above A_target, fractional part ~0.5
    v1 = Fr(int(math.floor(A_target)) - 1) + Fr(1, 2)
    if v1 < g0 + bbox(tT) + Fr(1, 100):
        v1 += 1
    sq, rows, bZ = cfg_blocked(200, tT, tZ, v1, m + 1, gT, gZ, grow, g0=g0)
    A = rows[0] + bZ - sca / 20
    sv = rows[m - 1] + sca / 20
    eps = 1 - A / sv
    k = int(math.ceil(2 * (float(sv) ** 2 / float(A) + 3)))
    if k > 130:
        continue
    sq, rows, bZ = cfg_blocked(k, tT, tZ, v1, m + 1, gT, gZ, grow, g0=g0, live_col=live)
    P = Params(k=k, delta=delta, c0=c0, c1=c1, y0=y0, eps=eps, omega0=k + 1)
    t0 = time.time()
    r = evaluate('maxratio(m=%d,c1=%s,c0=%s,live=%s,raise=%s)' % (m, c1, c0, live, raise_stack), sq, P, [sv], rng,
                 nys=4, xcheck_n=15, ovg_grid=8, tlimit=600)
    r['wall_time'] = time.time() - t0
    out.write(json.dumps(r, default=str) + '\n'); out.flush()
    c_ = r['per_s'][0] if r.get('per_s') else None
    if c_:
        print('m=%d c1=%s c0=%s live=%s raise=%s k=%d nsq=%d nZ=%d N=%.3f L0m=%.3f wall=%.3f Vlen=%.4f ratio=%.5f '
              'ratioZp=%.5f p2=%s live=%s fails=%s xbad=%s t=%.1f' % (
                  m, c1, c0, live, raise_stack, k, len(sq), c_['nZ'], c_['N'], c_['L0m'], c_['wall'], c_['Vlen'],
                  c_['a_ratio'], c_['a_ratio_Zprime'], c_['p2_worst'][0] if c_['p2_worst'] else None,
                  c_['b_live_on_Z'], c_['pointwise_fail'][:2], r['xcheck']['nbad'], r['wall_time']), flush=True)
        if best is None or c_['a_ratio'] > best[0]:
            best = (c_['a_ratio'], m, str(c1), str(c0), live, raise_stack, k)
    else:
        print('no per_s', r.get('invalid'), flush=True)
print('BEST', best)
out.write(json.dumps({'best': best}) + '\n')
out.close()
