"""p7x_merge.py <seed> <time_limit> -- merge-heavy configurations; check EVERY M piece and every winner entry:
  (1) point tracer (no merging) reproduces the loser's path up to its contact point,
  (2) a recorded winner reaches the same point with lexicographically smaller (g, x) (point-traced g),
  (3) at sampled winner entry points, no logged live arrival has smaller (g, x)  (R1 completeness),
plus the full P7 evaluation.  Original runs (2026-10-06): seeds 71 and 72 with time limit 480 s.
output: out/merge_<seed>.jsonl (appended)"""
import os, sys, json, random, time, math
from fractions import Fraction as Fr
from p7x_core import Sq, Flow, trace_point, pev
from p7x_meas import Params
from p7x_eval import evaluate
from p7x_configs import tq, bbox, drop

seed, TL = int(sys.argv[1]), float(sys.argv[2])
rng = random.Random(seed)
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
os.makedirs(OUTDIR, exist_ok=True)
out = open(os.path.join(OUTDIR, 'merge_%d.jsonl' % seed), 'a', encoding='utf-8')
t0 = time.time()
tot = dict(cfg=0, M=0, M_ok=0, M_bad=0, win=0, win_bad=0, Mmeas=0.0, eval_viol=0)
while time.time() - t0 < TL:
    k = rng.choice([12, 14, 16])
    y0 = rng.choice([Fr(3, 2), Fr(2)])
    c0 = Fr(rng.choice([30, 40, 50, 60]), 100)
    delta = Fr(rng.choice([10, 20, 30]), 100)
    P = Params(k=k, delta=delta, c0=c0, c1=Fr(2, 100), y0=y0, eps=Fr(1, 2), omega0=k + 1)
    am = float(P.alpha_max)
    th = am * rng.uniform(0.2, 0.97)
    g = Fr(rng.choice([1, 5, 20, 100, 500]), 10 ** 5)
    sq = []; sid = 0
    for r in range(rng.randint(2, 6)):
        x = Fr(rng.randrange(0, 700), 1000); i = 0
        while True:
            t = tq((th if (i + r) % 2 else -th) * rng.uniform(0.7, 1.0), 20000); b = bbox(t)
            if x + b > k:
                break
            S = drop(sq, sid, x + b / 2, t, g)
            if S.ymax > k / 2 - 1:
                break
            sq.append(S); sid += 1
            x += b + Fr(rng.randrange(1, 60), 10 ** 5); i += 1
    fl = Flow(sq, P.k, P.delta, P.alpha_max, P.hmax).run()
    tot['cfg'] += 1
    Ms = [r for r in fl.recs if r.typ == 'M']
    tot['Mmeas'] += float(sum((r.xhi - r.xlo) for r in Ms))
    for r in Ms[:60]:
        tot['M'] += 1
        x = (r.xlo + r.xhi) / 2
        typ, passed, endp, gx, contacts = trace_point(x, sq, P.k, P.delta, P.alpha_max, P.hmax)
        rp = [s_ for s_ in r.seg if s_ != -1]
        q = pev(r.way[-1], x)
        j = len(rp)
        if not (passed[:j] == rp and len(contacts) > j and contacts[j][0] == r.tsq and contacts[j][1] == q):
            tot['M_bad'] += 1
            out.write(json.dumps({'bad': 'M_path', 'x': str(x)}) + '\n'); continue
        gx = contacts[j][2]
        ok = False
        for w in fl.recs:
            for (s_, idx) in w.ent:
                if s_ != r.tsq:
                    continue
                Pw = w.way[idx]
                xw = (q[0] - Pw[0]) / Pw[1] if Pw[1] != 0 else ((q[1] - Pw[2]) / Pw[3] if Pw[3] != 0 else None)
                if xw is None or not (w.xlo < xw < w.xhi) or pev(Pw, xw) != q:
                    continue
                t2 = trace_point(xw, sq, P.k, P.delta, P.alpha_max, P.hmax)
                gw = [cg for (cs_, cq, cg) in t2[4] if cs_ == r.tsq and cq == q]
                ok = bool(gw) and (gw[0], xw) < (gx, x)
                break
            if ok:
                break
        if ok:
            tot['M_ok'] += 1
        else:
            tot['M_bad'] += 1
            out.write(json.dumps({'bad': 'R1_no_winner', 'x': str(x), 'sq': r.tsq}) + '\n')
    # (3) winner completeness at sampled winner entry points
    wins = [(w, s_, idx) for w in fl.recs for (s_, idx) in w.ent]
    rng.shuffle(wins)
    for (w, s_, idx) in wins[:60]:
        xw = (w.xlo + w.xhi) / 2
        q = pev(w.way[idx], xw)
        t2 = trace_point(xw, sq, P.k, P.delta, P.alpha_max, P.hmax)
        gw = [cg for (cs_, cq, cg) in t2[4] if cs_ == s_ and cq == q]
        if not gw:
            tot['win_bad'] += 1; continue
        tot['win'] += 1
        for cd in fl.cands:
            if cd.sid != s_:
                continue
            Q = cd.q
            xc = (q[0] - Q[0]) / Q[1] if Q[1] != 0 else ((q[1] - Q[2]) / Q[3] if Q[3] != 0 else None)
            if xc is None or not (cd.xlo < xc < cd.xhi) or pev(Q, xc) != q:
                continue
            gc = cd.g[0] + cd.g[1] * xc
            if (gc, xc) < (gw[0], xw):
                tot['win_bad'] += 1
                out.write(json.dumps({'bad': 'R1_not_lexmin', 'xw': str(xw), 'xc': str(xc), 'sq': s_}) + '\n')
    if Ms:
        ss = [P.y0, P.y1, (P.y0 + P.y1) / 2]
        r = evaluate('merge', sq, P, ss, rng, nys=6, xcheck_n=10, ovg_grid=10, tlimit=60)
        for c in r.get('per_s', []):
            tot['eval_viol'] += int(c['a_ratio'] > 1) + len(c['p2_viol']) + len(c['b_live_on_Z']) + len(c['pointwise_fail']) \
                + int(not c['c_ok']) + int(c['d_Tstar_symdiff'] != 0)
        r['merge_measure'] = tot['Mmeas']
        out.write(json.dumps(r, default=str) + '\n'); out.flush()
print('MERGE', seed, json.dumps(tot))
out.write(json.dumps({'tot': tot}) + '\n')
out.close()
