# -*- coding: utf-8 -*-
"""Non-vacuity check for the scaled analogue of Theorem 4.13 (see consequence.py).
Forbidden squares Z (Z(s)-analogue and/or P2-forbidden) are placed directly beside the top of a
column whose squares are tilted TOWARD Z (max tilt), with tiny horizontal gaps, so that live F_s
paths do reach Z.  Statement checked (analogue of Theorem 4.13(b),(c)): they touch Z only on a
side (E) or die/wall there, never on closed bot(Z), never pass through Z.  Also: Z shifted so
that it is NOT forbidden is reached on bot (control).

Usage:  python side_contacts.py [budget_s]      (default 600; writes out/side_contacts_out.json)
"""
import sys, time, json, random, math
from fractions import Fraction as Fr
import mpmath as mp
from tracer import Sq, Cfg, trace, rot, HALF, mpf, atan2t, Kcoef, maxdist_interval
from common import out_path, build_column, ok_with, rss_mb
from consequence import classify_Z, Back, f_of

T0 = time.time()
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 600
out = dict(cases=[], n_cases=0, n_forb=0, side_hits=0, bot_hits_forb=0, pass_forb=0,
           ctrl_cases=0, ctrl_bot_hits=0, viol=[], kinds={}, anomalies=[])
r = random.Random(12)
for c0 in (Fr(3, 10), Fr(9, 20), Fr(3, 5)):
    for s in (4, 9, 16):
        for delta in (Fr(1, 100), Fr(3, 100)):
            for sgn in (1, -1):
                for gapx in (Fr(1, 10 ** 6), Fr(1, 10 ** 4)):
                    if time.time() - T0 > BUDGET:
                        break
                    alpha_true = mp.mpf(c0.numerator) / c0.denominator / mp.sqrt(s)
                    tau = Fr(mp.nstr(mp.tan(alpha_true / 2), 15)).limit_denominator(10 ** 9)
                    alpha = atan2t(tau)
                    b = Fr(1, 100)
                    eps = Fr(3, 10)
                    w0 = Kcoef(alpha) * (s + 1 + mpf(b)) + mpf(delta) + f_of(mpf(b)) + mp.mpf("1e-12")
                    if w0 >= 0.5:
                        continue
                    k = Fr(2 * s + 8)
                    m = s - 1
                    # column tilted toward +x if sgn=+1: phi<0 gives n_x>0 (moves right)
                    tt = -sgn * tau * Fr(999, 1000)
                    x0 = Fr(s + 4)
                    specs = [(tt, Fr(0))] * m
                    gaps = [Fr(1, 10 ** 9)] * m
                    c1, s1 = rot(tt)
                    gaps[0] = abs(s1) / 2 + Fr(1, 10 ** 9)
                    col = build_column(x0, specs, gaps, k)
                    if col is None:
                        out["anomalies"].append("column build failed")
                        continue
                    sq, p, d, used = col
                    L = sq[-1]
                    # Z beside L on the side the paths drift to; Z spans the exit height of L's top
                    cornerT = L.pt(sgn * HALF, HALF)
                    ytop = cornerT[1]
                    ftop = ytop - math.floor(ytop)
                    w0f = Fr(mp.nstr(w0, 12))
                    opts = [("beside", w0f + Fr(1, 10 ** 6)), ("beside", (w0f + ftop) / 2),
                            ("beside", ftop - Fr(1, 100)), ("above_ctrl", None)]
                    for trial, (where, frc) in enumerate(opts):
                        aZ = Fr(1, 200)
                        tZ = Fr(math.tan(float(aZ) / 2)).limit_denominator(10 ** 9) * (-sgn)
                        cz, sz = rot(tZ)
                        best = None
                        if where == "beside":
                            if frc is None or frc <= 0 or frc >= ftop:
                                continue
                            v = Fr(math.floor(ytop)) + frc
                            cy = v + (cz + abs(sz)) / 2
                            cx = cornerT[0] + sgn * (HALF + gapx + abs(sz))
                            for it in range(80):
                                S = Sq(cx, cy, tZ, tag="Zbeside")
                                if ok_with(S, sq, k):
                                    best = S
                                    break
                                cx += sgn * Fr(1, 10 ** 4)
                        else:
                            # control: a translate of L along n_L by 1 + delta/10 (reachable, near-integer)
                            S = Sq(L.cx + (1 + delta / 10) * L.n()[0], L.cy + (1 + delta / 10) * L.n()[1],
                                   L.t, tag="Zabove")
                            if ok_with(S, sq, k):
                                best = S
                        if best is None:
                            continue
                        Z = best
                        cfg = Cfg(k, sq + [Z], "side_c%s_s%d_d%s_g%s_%d_%d" % (c0, s, delta, gapx, sgn, trial))
                        ok, msg = cfg.validate()
                        if not ok:
                            out["anomalies"].append("invalid " + msg)
                            continue
                        zs, wit, p2f, slack = classify_Z(Z, s, eps, w0, b, alpha, delta)
                        forb = zs or p2f
                        h = Fr(s + 2)
                        # forward: dense x over the column width, plus a fine band near the drift edge
                        xs = [x0 + Fr(j, 200) - HALF for j in range(201)]
                        xs += [x0 + sgn * (HALF - Fr(j, 20000)) for j in range(200)]
                        hits = dict(side=0, bot=0, vbot=0, vtop=0, passed=0)
                        for x in xs:
                            res = trace(cfg, x, delta, tau, h)
                            for i, rec in enumerate(res["contacts"]):
                                if rec["Z"] != Z.idx:
                                    continue
                                kk = rec["kind"] + ("/" + rec["flag"] if rec["flag"] else "")
                                out["kinds"][kk] = out["kinds"].get(kk, 0) + 1
                                onbot = rec["kind"] in ("bot", "vbot") and rec["eta"] == -HALF
                                if onbot:
                                    hits["bot"] += 1
                                else:
                                    hits["side"] += 1
                                last = (i == len(res["contacts"]) - 1)
                                if forb and (onbot or not last or res["term"] not in ("E", "D", "W")):
                                    out["viol"].append(dict(cfg=cfg.name, x=str(x), kind=rec["kind"], term=res["term"]))
                            if Z.idx in res["passed"]:
                                hits["passed"] += 1
                                if forb:
                                    out["viol"].append(dict(cfg=cfg.name, x=str(x), msg="passed forbidden Z"))
                        # backward search from bot(Z) points (exhaustive per point)
                        B = Back(cfg, tau, delta, h)
                        bfound = 0
                        for j in range(81):
                            q = Z.pt(Fr(j, 80) - HALF, -HALF)
                            if q[1] > h:
                                continue
                            o = []
                            B.nodes = 0
                            B.preds(q, delta, 0, [], o, limit=3)
                            if o:
                                bfound += 1
                        if forb and bfound:
                            out["viol"].append(dict(cfg=cfg.name, msg="backward found bot path", n=bfound))
                        case = dict(cfg=cfg.name, forb=bool(forb), Zs=bool(zs), P2f=bool(p2f),
                                    slack=float(slack), hits=hits, back_bot_points_reached=bfound,
                                    Z=Z.to_json(), v=float(Z.v))
                        out["cases"].append(case)
                        out["n_cases"] += 1
                        if forb:
                            out["n_forb"] += 1
                            out["side_hits"] += hits["side"]
                            out["bot_hits_forb"] += hits["bot"]
                            out["pass_forb"] += hits["passed"]
                        else:
                            out["ctrl_cases"] += 1
                            out["ctrl_bot_hits"] += hits["bot"] + bfound
                        if rss_mb() > 700:
                            raise SystemExit(4)
out["elapsed"] = time.time() - T0
with open(out_path("side_contacts_out.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, default=str)
print(json.dumps({k: v for k, v in out.items() if k != "cases"}, indent=1, default=str))
print("cases with side hits on forbidden Z:", sum(1 for c in out["cases"] if c["forb"] and c["hits"]["side"] > 0))
