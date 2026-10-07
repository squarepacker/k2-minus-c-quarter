"""Cross-check of the beam tracer by an independent per-point float tracer (no beams, no R1) on configurations
where the exact tracer reports no merge: event measures must agree up to sampling error.

Usage:  python cross_check.py        (writes out/cross_check_results.json; under a minute)"""
import math, random, json, time
from fractions import Fraction as Fr
from tracer import trace, check_config, out_path
import configs as C

def ftrace(sqs, k, delta, h, sF, x):
    P = [(float(S.c[0]), float(S.c[1]), float(S.u[0]), float(S.u[1]), float(S.n[0]), float(S.n[1]), abs(S.s) >= sF) for S in sqs]
    p = (x, 0.0); d = (0.0, 1.0); src = -1; g = 0.0
    for _ in range(1000):
        if p[1] >= h: return 'H'
        best = None
        for i, (cx, cy, ux, uy, nx, ny, tilt) in enumerate(P):
            if i == src: continue
            # ray p + t d vs square: slab intersection in local frame
            rx, ry = p[0] - cx, p[1] - cy
            a0, a1 = rx*ux + ry*uy, d[0]*ux + d[1]*uy
            b0, b1 = rx*nx + ry*ny, d[0]*nx + d[1]*ny
            tlo, thi = -1e18, 1e18; side = None
            for (o, s_, nm) in ((a0, a1, 'u'), (b0, b1, 'n')):
                if abs(s_) < 1e-18:
                    if abs(o) > 0.5: tlo, thi = 1, 0
                    continue
                t1, t2 = (-0.5 - o)/s_, (0.5 - o)/s_
                if t1 > t2: t1, t2 = t2, t1
                if t1 > tlo: tlo = t1; side = nm
                thi = min(thi, t2)
            if tlo <= thi and thi >= 0 and tlo >= -1e-12:
                if best is None or tlo < best[0]: best = (max(tlo, 0.0), i, side)
        tD = delta - g
        tW = (-p[0]/d[0]) if d[0] < 0 else ((k - p[0])/d[0] if d[0] > 0 else 1e18)
        tH = (h - p[1]) / d[1]
        tS = best[0] if best else 1e18
        m = min(tD, tW, tS, tH)
        if tD == m: return 'D'
        if tW == m: return 'W'
        if tS == m:
            t, i, side = best
            q = (p[0] + t*d[0], p[1] + t*d[1])
            cx, cy, ux, uy, nx, ny, tilt = P[i]
            bl = (q[0]-cx)*nx + (q[1]-cy)*ny
            if side == 'n' and bl < 0:     # bottom edge
                if tilt: return 'T'
                g += t
                p = (q[0] + nx, q[1] + ny); d = (nx, ny); src = i
                continue
            return 'E'
        return 'H'
    return '?'

out = []
rnd = random.Random(5)
t0 = time.time()
cands = []
for delta in (Fr(1, 100), Fr(1, 20), Fr(1, 10)):
    for seed in range(6):
        mode = seed % 3
        samp = [lambda r: r.uniform(-0.05, 0.05), lambda r: r.uniform(-0.3, 0.3), lambda r: r.choice([0.0, r.uniform(-0.6, 0.6)])][mode]
        sq = C.rowjam(6, delta, 3, samp, 300 + seed)
        cands.append((sq, delta, Fr(10) if mode < 2 else C.s_of(0.4)))
for sq, delta, sF in cands:
    if time.time() - t0 > 45: break
    if check_config(sq, 6): continue
    rec = trace(sq, 6, delta, 2, sF)
    if rec.M: continue
    n = 4000
    cnt = {}
    for j in range(n):
        x = (j + rnd.random()) * 6 / n
        e = ftrace(sq, 6, float(delta), 2.0, sF, x)
        cnt[e] = cnt.get(e, 0) + 1
    ex = dict(D=float(rec.D), W=float(rec.W), H=float(rec.H), T=float(rec.T), E=float(sum(rec.E.values(), Fr(0))))
    est = {k_: 6 * cnt.get(k_, 0) / n for k_ in ex}
    out.append(dict(delta=float(delta), N=len(sq), exact=ex, sampled=est, maxdiff=max(abs(ex[k_] - est[k_]) for k_ in ex)))
    print(out[-1], flush=True)
json.dump(out, open(out_path('cross_check_results.json'), 'w'), indent=1)
print('max |exact - sampled| =', max(o['maxdiff'] for o in out), 'over', len(out), 'configs (stratified 4000 samples, step 0.0015)')
