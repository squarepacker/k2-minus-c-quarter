# -*- coding: utf-8 -*-
"""asmx_targeted.py -- targeted checks of the assembly test (written from scratch, 2026-10-06).
T1 negative control "T before M" (order B): with merges at T_max squares the equality {T* <= s} = T_s of
   Lemma 9.6(c) breaks under order B, while it holds with M before T (Definition 3.14),
T2 independent pointwise (Cyrus-Beck) tracer cross-check,
T3 P7 stress (Theorem 9.1) with a fine s-grid, T4 wall-loss stress (Proposition 6.10) with 'hill' configurations.
usage: python asmx_targeted.py T1|T2|T3|T4|all        output: out/targeted_<T..>.json"""
import os, sys, math, time, random, json
import numpy as np
from asmx_core import Packing, Flow, walk, TYP, Sq
from asmx_chain import Params, LineStruct, FamilyFlows, family_links
from asmx_gen import GENS_ALL, compose, tilted_row, ext_of

OUT = {}


# ------------------------------------------------------------------ T2: pointwise Cyrus-Beck tracer
def ray_hit(p, d, Y):
    """smallest t>=0 with p+t d in Y (closed), and the index of the entering edge (or -1)."""
    tin, tout, ein = -1e300, 1e300, -1
    for e in range(4):
        nx, ny = Y.nrm[e]
        vx, vy = Y.V[e]
        c = nx * vx + ny * vy
        nd = nx * d[0] + ny * d[1]
        np_ = nx * p[0] + ny * p[1]
        if abs(nd) < 1e-15:
            if np_ > c + 1e-12:
                return None
            continue
        t = (c - np_) / nd
        if nd < 0:
            if t > tin:
                tin, ein = t, e
        else:
            tout = min(tout, t)
    if tin > tout + 1e-12 or tout < -1e-12:
        return None
    if tin < 0:
        return (0.0, ein)
    return (tin, ein)


def trace_point(P, x, delta, alphaF, hF, stop_sq=None):
    k = P.k
    p = (x, 0.0); d = (0.0, 1.0); src = -1; g = 0.0
    for step in range(100000):
        reach = delta - g
        x0 = min(p[0], p[0] + reach * d[0]); x1 = max(p[0], p[0] + reach * d[0])
        cand = P.query(x0, x1, p[1], p[1] + reach * d[1])
        bt = 1e300; best = None
        for j in cand:
            if j == src:
                continue
            h = ray_hit(p, d, P.S[j])
            if h is not None and h[0] < bt:
                bt, best = h[0], (j, h[1])
        tD = delta - g
        if d[0] > 1e-15:
            tW = (k - p[0]) / d[0]
        elif d[0] < -1e-15:
            tW = -p[0] / d[0]
        else:
            tW = 1e300
        tH = (hF - p[1]) / d[1]
        tm = min(tD, tW, bt, tH)
        tol = 1e-12 * max(1.0, abs(tm))
        if tD <= tm + tol:
            return ('D', (p[0] + tD * d[0], p[1] + tD * d[1]), -1)
        if tW <= tm + tol:
            return ('W', (p[0] + tW * d[0], p[1] + tW * d[1]), -1)
        if bt <= tm + tol:
            j, e = best
            Y = P.S[j]
            q = (p[0] + bt * d[0], p[1] + bt * d[1])
            g += bt
            if e != 0:
                return ('E', q, j)
            xi = (q[0] - Y.cx) * Y.u[0] + (q[1] - Y.cy) * Y.u[1]
            if abs(abs(xi) - 0.5) < 1e-9:
                return ('V', q, j)       # vertex: measure zero, skip
            if stop_sq is not None and j == stop_sq:
                return ('STOP', q, j)
            if Y.a >= alphaF:
                return ('T', q, j)
            p = (q[0] + Y.n[0], q[1] + Y.n[1]); d = Y.n; src = j
            if p[1] >= hF:
                return ('H', p, -1)
            continue
        return ('H', (p[0] + tH * d[0], p[1] + tH * d[1]), -1)
    return ('LOOP', p, -1)


def t2_crosscheck(ncfg, seed, tl):
    rng = random.Random(seed)
    t0 = time.time()
    stats = dict(cfgs=0, samples=0, match=0, mismatch=0, skipped=0, Mchecked=0, Mbad=0, examples=[])
    fams = ['rows', 'columns', 'valley', 'floortouch', 'random', 'hill', 'nested', 'lshape']
    while stats['cfgs'] < ncfg and time.time() - t0 < tl:
        fam = fams[stats['cfgs'] % len(fams)]
        k = 16 if fam == 'random' else rng.choice([16, 20])
        delta = rng.choice([0.05, 0.1, 0.2]); c0 = rng.choice([0.15, 0.2, 0.25]); y0 = 3
        amax = c0 / math.sqrt(y0)
        gprm = dict(amax=amax, ay1=c0 / math.sqrt(k / 2 - 3), beta_lo=0.02, beta_hi=0.05, delta=delta, a_eg=0.05)
        cfg = compose(k, GENS_ALL[fam](k, k / 2 - 0.02, rng, gprm), GENS_ALL['rows'](k, k / 2 - 0.02, rng, gprm))
        P = Packing(k, cfg)
        if not P.validate()[0]:
            continue
        stats['cfgs'] += 1
        F = Flow(P, delta, amax, k / 2 - 1); F.run(60)
        recs = sorted(F.recs, key=lambda r: r[0])
        xs = [rng.uniform(0, k) for _ in range(400)]
        for x in xs:
            r = None
            for rr in recs:
                if rr[0] + 1e-9 < x < rr[1] - 1e-9:
                    r = rr; break
            if r is None:
                stats['skipped'] += 1; continue
            typ = r[2]
            Tp = (r[4] + r[6] * x, r[5] + r[7] * x)
            if typ == 'M':
                res = trace_point(P, x, delta, amax, k / 2 - 1, stop_sq=r[8])
                stats['Mchecked'] += 1
                if res[0] != 'STOP' or abs(res[1][0] - Tp[0]) > 1e-8 or abs(res[1][1] - Tp[1]) > 1e-8:
                    stats['Mbad'] += 1
                    if len(stats['examples']) < 5:
                        stats['examples'].append(('M', fam, x, res[0], res[1], Tp))
                continue
            res = trace_point(P, x, delta, amax, k / 2 - 1)
            if res[0] == 'V':
                stats['skipped'] += 1; continue
            stats['samples'] += 1
            ok = (res[0] == typ) and abs(res[1][0] - Tp[0]) < 1e-8 and abs(res[1][1] - Tp[1]) < 1e-8
            if ok:
                stats['match'] += 1
            else:
                # a winner path whose pointwise (no-merge) trace differs only because it would have
                # merged: never possible for winners; record
                stats['mismatch'] += 1
                if len(stats['examples']) < 8:
                    stats['examples'].append((fam, x, typ, res[0], res[1], Tp, r[8], res[2]))
    stats['time'] = time.time() - t0
    return stats


# ------------------------------------------------------------------ T1: order-B negative control
def t1_orderB(seed, tl):
    rng = random.Random(seed)
    t0 = time.time()
    out = []
    tries = 0
    while time.time() - t0 < tl and len(out) < 12 and tries < 200:
        tries += 1
        k = 16; delta = 0.1; c0 = 0.2; y0 = 3
        amax = c0 / math.sqrt(y0)
        a = rng.uniform(0.2, 0.6) * min(amax, delta / 2.5)
        cfg = []
        cfg += tilted_row(k, 0.5, 0.0, k - 1, 1e-6 + rng.uniform(0, 0.9), 1e-5)
        base = 1.0 + 1e-6
        w = ext_of(a)
        nl = k // 2 - 1
        span = (nl - 1) * (1 + 1e-5) / math.cos(a) + w
        mid = k / 2.0 + rng.uniform(-0.3, 0.3)
        cfg += tilted_row(k, base + w / 2, -a, nl, mid - 5e-7 - span, 1e-5)
        cfg += tilted_row(k, base + w / 2, a, nl, mid + 5e-7, 1e-5)
        base += w + rng.uniform(0.05, 0.4) * delta
        b = rng.uniform(1.05, 2.5) * amax * rng.choice([-1, 1])
        wb = ext_of(b)
        nb = int((k - wb) / ((1 + 1e-4) / math.cos(abs(b)))) + 1
        spb = (nb - 1) * (1 + 1e-4) / math.cos(abs(b)) + wb
        xb0 = min(max(1e-6, mid - spb / 2 + rng.uniform(-0.5, 0.5)), k - spb - 1e-6)
        cfg += tilted_row(k, base + wb / 2, b, nb, xb0, 1e-4)
        P = Packing(k, cfg)
        if not P.validate()[0]:
            continue
        sg = list(np.linspace(y0, k / 2 - 3, 12))
        FA = FamilyFlows(P, k, delta, c0, y0, sg, nov=50)
        FB = FamilyFlows(P, k, delta, c0, y0, sg, nov=50, orderB=True)
        out.append(dict(M_A=float(FA.meas0[3]), T_A=float(FA.meas0[4]), M_B=float(FB.meas0[3]), T_B=float(FB.meas0[4]),
                        symd_A=max(d['symd'] for d in FA.S), symd_B=max(d['symd'] for d in FB.S),
                        N_B_minus_n=max(d['N'] - d['n'] for d in FB.S)))
    return dict(runs=out, time=time.time() - t0)


# ------------------------------------------------------------------ T3: P7 stress (Theorem 9.1) with a fine s-grid
def t3_p7stress(seed, tl, nconf):
    rng = random.Random(seed)
    t0 = time.time()
    best = []
    cnt = 0; viol = 0; worst = (0, None)
    allr = []
    while time.time() - t0 < tl and cnt < nconf:
        k = rng.choice([20, 24, 28])
        delta = rng.choice([0.02, 0.05]); c0 = rng.choice([0.08, 0.1, 0.15]); y0 = rng.choice([3, 4])
        c1 = rng.choice([0.25, 0.3, 0.35]); eps = rng.choice([0.2, 0.3, 0.35, 0.5])
        y1 = k / 2 - 3; top = (1 - eps) * y1
        prm = Params(k, delta, c0, y0, c1, eps, float(k))
        if prm.ok():
            continue
        gprm = dict(amax=c0 / math.sqrt(y0), ay1=c0 / math.sqrt(y1), beta_lo=c1 * top ** -0.75,
                    beta_hi=c1 * y0 ** -0.75, delta=delta, a_eg=0.97 * c1 * top ** -0.75)
        fam = rng.choice(['nestedEG', 'nested'])
        cfg = compose(k, GENS_ALL[fam](k, k / 2 - 0.02, rng, gprm), GENS_ALL[fam](k, k / 2 - 0.02, rng, gprm))
        P = Packing(k, cfg)
        if not P.validate()[0]:
            continue
        cnt += 1
        sg = sorted(set(list(np.linspace(y0, y1, 70)) + [rng.uniform(y0, y1) for _ in range(10)]))
        FF = FamilyFlows(P, k, delta, c0, y0, sg, nov=100)
        LS = LineStruct(P, prm)
        fl = family_links(FF, LS, prm, nypts=32)
        for d in fl['per_s']:
            if d['nZ'] > 0:
                rr = (d.get('r_sharp', 0), d['r_meas'], d['r_paper'], d['nZ'], d['s'], k, eps, fam)
                allr.append(rr)
                for kk in ('r_paper', 'r_meas', 'r_sharp', 'r_shadow'):
                    if d.get(kk, 0) > 1 + 1e-9:
                        viol += 1
                # diagnostic: half b_W
                Vh = eps * d['s'] + prm.h
                den = Vh * (d['N'] + FF.Lam0 + d['LW'])
                if den > 0 and d['nZ'] / den > worst[0]:
                    worst = (d['nZ'] / den, (d['s'], d['nZ'], k, eps, fam))
    allr.sort(reverse=True)
    return dict(configs=cnt, viol=viol, top=allr[:8], halfbW_worst=worst, time=time.time() - t0,
                n_cases=len(allr))


# ------------------------------------------------------------------ T4: wall-loss stress (hill)
def t4_wall(seed, tl, nconf):
    rng = random.Random(seed)
    t0 = time.time(); cnt = 0
    best = (0, None); viol = 0; halfviol = 0
    while time.time() - t0 < tl and cnt < nconf:
        k = rng.choice([16, 20, 24])
        delta = rng.choice([0.05, 0.1, 0.2]); c0 = rng.choice([0.1, 0.15, 0.2, 0.25]); y0 = rng.choice([2, 3])
        y1 = k / 2 - 3
        gprm = dict(amax=c0 / math.sqrt(y0), ay1=c0 / math.sqrt(y1), beta_lo=0.02, beta_hi=0.05, delta=delta, a_eg=0.05)
        cfg = compose(k, GENS_ALL['hill'](k, k / 2 - 0.02, rng, gprm), GENS_ALL['hill'](k, k / 2 - 0.02, rng, gprm))
        P = Packing(k, cfg)
        if not P.validate()[0]:
            continue
        cnt += 1
        sg = list(np.linspace(y0, y1, 25))
        FF = FamilyFlows(P, k, delta, c0, y0, sg, nov=50)
        for d in FF.S:
            s = d['s']; wb = 2 * (s + 2) * math.tan(c0 / math.sqrt(s))
            r = d['LW'] / wb
            if r > 1 + 1e-9:
                viol += 1
            if r > 0.5 + 1e-9:
                halfviol += 1
            if r > best[0]:
                best = (r, (s, d['LW'], wb, k, c0, y0))
    return dict(configs=cnt, viol=viol, gt_half=halfviol, best=best, time=time.time() - t0)


if __name__ == '__main__':
    which = sys.argv[1]
    res = {}
    if which in ('T1', 'all'):
        res['T1'] = t1_orderB(21, 150)
    if which in ('T2', 'all'):
        res['T2'] = t2_crosscheck(40, 22, 400)
    if which in ('T3', 'all'):
        res['T3'] = t3_p7stress(23, 420, 60)
    if which in ('T4', 'all'):
        res['T4'] = t4_wall(24, 200, 40)
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'targeted_%s.json' % which), 'w', encoding='utf-8') as fh:
        json.dump(res, fh, indent=1, default=str)
    print(json.dumps(res, indent=1, default=str)[:6000])
