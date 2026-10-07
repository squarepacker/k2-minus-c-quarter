# -*- coding: utf-8 -*-
"""Scaled analogue of Theorem 4.13 (P2 for the squares of Z(s)), parts (b)-(d).

alpha(s) = c0 s^{-1/2}, F_s = flow F(alpha(s), delta, s+2) without rule R1 (a superset of
the contacts of the flow with R1; see tracer.py).
Squares Z of the Z(s)-analogue: a(Z) < b and some y* in Phi(Z) cap [(1-eps)s, s] with
dist(y*,Z) >= w0, where (analogue of the constant w0 of Section 4.1; the exact corollary
constant alpha^2/(2 cos alpha) is used in place of 0.50001 alpha^2, since alpha may exceed 1e-3)
   w0 = alpha^2/(2 cos alpha) * (s + 1 + b) + delta + (sin b + 1 - cos b) + 1e-12.
Also 'P2-forbidden' squares: some ramp point y with dist(y,Z) >= K*(v+sin a) + delta + f(a(Z)).
Statement tested (analogue of Theorem 4.13(b),(c)): no F_s trajectory point on closed bot(Z);
paths meet Z only in their termination point (E, D, W) on the boundary minus bot.
Two independent searches: forward dense sampling, and exhaustive backward predecessor search
from sample points of bot(Z).

Usage:  python consequence.py <budget_s> <prefix> [quick|full] [c0=<value>]
  writes out/<prefix>_summary.json.  Without c0=..., c0 runs over 3/10, 45/100, 6/10 in turn
  (until the budget is used up).
"""
import sys, time, json, random, math
from fractions import Fraction as Fr
import mpmath as mp
from tracer import (Sq, Cfg, trace, rot, HALF, mpf, atan2t, Kcoef, dist_int,
                    maxdist_interval, verify_contact, disjoint)
from common import out_path, rss_mb, place, place_min_gap, build_column, ok_with

T0 = time.time()
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 600
OUT = sys.argv[2] if len(sys.argv) > 2 else "consequence"
QUICK = len(sys.argv) > 3 and sys.argv[3] == "quick"


def tl():
    return BUDGET - (time.time() - T0)


def f_of(a):
    return mp.sin(a) + 1 - mp.cos(a)


def phi_intervals(Z):
    """open intervals of Phi(Z) = {0 < chord < 1}"""
    if Z.t == 0:
        return []
    v, sa, ca = Z.v, Z.sa, Z.ca
    return [(v, v + sa * ca), (v + ca + sa - sa * ca, v + ca + sa)]


def best_offint_point(lo, hi, w0):
    """sup of dist(y,Z) over the open interval (lo,hi) restricted to dist >= w0;
    returns (best_dist_float, witness y rational) or None"""
    lo, hi = Fr(lo), Fr(hi)
    best = None
    n0 = math.floor(lo) - 1
    n1 = math.floor(hi) + 1
    for n in range(n0, n1 + 1):
        # target region [n + w0, n + 1 - w0]; best point = n + 1/2 if inside
        cand = [Fr(n) + HALF]
        for y in cand:
            pass
        A = max(lo, Fr(n))
        B = min(hi, Fr(n + 1))
        if A >= B:
            continue
        # point in (A,B) maximizing dist: closest to n+1/2
        mid = Fr(n) + HALF
        if A < mid < B:
            y = mid
        elif mid <= A:
            y = A + (B - A) / 10 ** 6
        else:
            y = B - (B - A) / 10 ** 6
        dd = dist_int(y)
        if mpf(dd) >= w0 and (best is None or dd > best[0]):
            best = (dd, y)
    return best


def classify_Z(Z, s, eps, w0, b, alpha, delta):
    aZ = Z.a_mp()
    zs = False
    wit = None
    if aZ < b:
        for (lo, hi) in phi_intervals(Z):
            A = max(lo, (1 - eps) * s)
            B = min(hi, Fr(s))
            if A >= B:
                continue
            r = best_offint_point(A, B, w0)
            if r is not None:
                zs = True
                wit = r
                break
    # P2-forbidden: some ramp point with dist >= K*(v + sa) + delta + f(a)
    wmax = Kcoef(alpha) * mpf(Z.v + Z.sa) + mpf(delta) + mpf(Z.sa + 1 - Z.ca)
    md = Fr(0)
    if Z.t != 0:
        for (lo, hi) in Z.ramps_closure():
            dd, yy = maxdist_interval(lo, hi)
            md = max(md, dd)
    p2f = (Z.t != 0) and mpf(md) >= wmax
    slack = mpf(md) - wmax
    return zs, wit, p2f, slack


# ------------------------------------------------------------ backward predecessor search
def ray_interval(p, d, S):
    """param interval {t : p + t d in S} (any direction), or None"""
    rx, ry = p[0] - S.cx, p[1] - S.cy
    xi0 = rx * S.c + ry * S.s
    eta0 = -rx * S.s + ry * S.c
    du = d[0] * S.c + d[1] * S.s
    dn = -d[0] * S.s + d[1] * S.c
    L, H = None, None
    for (v0, dv) in ((xi0, du), (eta0, dn)):
        if dv == 0:
            if abs(v0) > HALF:
                return None
            continue
        a = (-HALF - v0) / dv
        bb = (HALF - v0) / dv
        lo_, hi_ = (a, bb) if a <= bb else (bb, a)
        L = lo_ if L is None else max(L, lo_)
        H = hi_ if H is None else min(H, hi_)
    if L is None:
        return (Fr(-10 ** 9), Fr(10 ** 9))
    if L > H:
        return None
    return (L, H)


def reverse_clear(cfg, q, d, t_end, allow_end_idx):
    """True iff the reverse open segment q - t d, 0 < t < t_end meets no square, and at t_end
    only square allow_end_idx (if not None) is met; also no wall crossing."""
    nd = (-d[0], -d[1])
    for S in cfg.candidates(q[0], q[1], nd[0], nd[1], t_end):
        r = ray_interval(q, nd, S)
        if r is None:
            continue
        L, H = r
        if H <= 0:
            continue
        if L >= t_end:
            if L == t_end and S.idx != allow_end_idx:
                return False
            continue
        # interval meets (0, t_end)
        return False
    end = (q[0] - t_end * d[0], q[1] - t_end * d[1])
    if not (0 <= end[0] <= cfg.k):
        return False
    return True


class Back:
    def __init__(self, cfg, tau, delta, h):
        self.cfg = cfg
        self.tau = Fr(tau)
        self.delta = Fr(delta)
        self.h = Fr(h)
        self.pass_sq = [S for S in cfg.sq if abs(S.t) < self.tau]
        self.nodes = 0

    def preds(self, q, g_rem, depth, chain, out, limit=200):
        """q: a point that must be a first contact (on the trajectory end of a gap).
        Enumerate all predecessor chains with total gap <= g_rem."""
        self.nodes += 1
        if self.nodes > 200000 or len(out) >= limit:
            return
        # floor source: vertical ray from (q_x, 0)
        if q[1] <= g_rem and 0 < q[0] < self.cfg.k:
            if q[1] == 0 or reverse_clear(self.cfg, q, (Fr(0), Fr(1)), q[1], None):
                out.append((q[0], list(reversed(chain)), q[1]))
        # square sources
        qf = (float(q[0]), float(q[1]))
        gf = float(g_rem)
        for X in self.pass_sq:
            if X.idx in chain[-1:]:
                pass
            nx, ny = -float(X.s), float(X.c)
            ux, uy = float(X.c), float(X.s)
            rx, ry = qf[0] - float(X.cx), qf[1] - float(X.cy)
            tn = rx * nx + ry * ny - 0.5
            tu = rx * ux + ry * uy
            if tn < -1e-9 or tn > gf + 1e-9 or abs(tu) > 0.5 + 1e-9:
                continue
            # exact
            n = X.n()
            erx, ery = q[0] - X.cx, q[1] - X.cy
            t = erx * n[0] + ery * n[1] - HALF
            xi = erx * X.c + ery * X.s
            if not (0 < t <= g_rem) or not (abs(xi) < HALF):
                continue
            r = (q[0] - t * n[0], q[1] - t * n[1])     # exit point on relint top(X)
            if r[1] >= self.h:
                continue        # path would end by H at the exit
            if not reverse_clear(self.cfg, q, n, t, X.idx):
                continue
            e = (r[0] - n[0], r[1] - n[1])               # entry point on relint bot(X)
            self.preds(e, g_rem - t, depth + 1, chain + [X.idx], out, limit)


# ------------------------------------------------------------ configurations
def make_cfg(c0, s, delta, b, eps, seed, alpha, tau, w0):
    r = random.Random(seed)
    k = Fr(2 * s + 8)
    sq = []
    ncol = max(2, (int(k) - 2) // 2)
    styles = ["maxtilt", "maxtilt_neg", "alt", "rand", "axis", "late"]
    col_tops = []
    for ci in range(ncol):
        x0 = Fr(2 + 2 * ci) + Fr(r.randint(0, 999), 4000)
        if x0 > k - 2:
            break
        st = styles[(ci + seed) % len(styles)]
        m = max(1, s + r.randint(-3, 1))
        specs = []
        for i in range(m):
            if st == "maxtilt":
                t = tau * Fr(999, 1000)
            elif st == "maxtilt_neg":
                t = -tau * Fr(999, 1000)
            elif st == "alt":
                t = tau * Fr(999, 1000) * (1 if i % 2 else -1)
            elif st == "axis":
                t = Fr(0)
            else:
                t = tau * Fr(r.randint(-999, 999), 1000)
            xi = Fr(r.randint(-30, 30), 100) if st in ("rand", "late") else Fr(0)
            specs.append((t, xi))
        if st == "late":
            gaps = [Fr(1, 10 ** 12)] * (m - 1) + [delta * Fr(9, 10)]
        elif st == "alt":
            gaps = [delta * Fr(1, 10 ** 6) / (m + 1)] * m
        else:
            gaps = [delta * Fr(r.randint(1, 900), 1000) / (m + 1)] * m
        c1, s1 = rot(specs[0][0])
        gaps[0] = gaps[0] + (HALF + abs(specs[0][1])) * abs(s1)
        col = build_column(x0, specs, gaps, k, existing=sq, gcap=Fr(1, 2))
        if col is None:
            continue
        sq += col[0]
        col_tops.append((col[1], col[2]))
    # forbidden Z candidates near column tops and in free space
    zcount = 0
    for (p, d) in col_tops:
        for trial in range(6):
            aZf = r.choice([Fr(1, 2), Fr(99, 100), Fr(1, 10)])
            tZ = Fr(math.tan(float(b) * float(aZf) / 2)).limit_denominator(10 ** 9) * r.choice([1, -1])
            if tZ == 0:
                continue
            cz, sz = rot(tZ)
            n0 = math.floor(p[1]) + r.choice([0, 0, 1, -1])
            frc = r.choice([Fr(mp.nstr(w0, 12)) + Fr(1, 10 ** 6), HALF - abs(sz) / 2,
                            1 - Fr(mp.nstr(w0, 12)) - abs(sz) * cz - Fr(1, 10 ** 6),
                            Fr(r.randint(0, 999), 1000)])
            v = Fr(n0) + frc
            cy = v + (cz + abs(sz)) / 2
            cx = p[0] + Fr(r.randint(-60, 60), 100)
            S = Sq(cx, cy, tZ, tag="Zc")
            if ok_with(S, sq, k):
                sq.append(S)
                zcount += 1
    # forbidden-type squares BESIDE column tops (so that paths can touch their sides)
    for (p, d) in col_tops:
        for side in (-1, 1):
            for trial in range(4):
                aZf = r.choice([Fr(1, 2), Fr(99, 100)])
                tZ = Fr(math.tan(float(b) * float(aZf) / 2)).limit_denominator(10 ** 9) * r.choice([1, -1])
                cz, sz = rot(tZ)
                frc = r.choice([HALF - abs(sz) / 2, Fr(mp.nstr(w0, 12)) + Fr(1, 10 ** 6), Fr(r.randint(300, 700), 1000)])
                n0 = math.floor(p[1]) - r.choice([0, 1])
                v = Fr(n0) + frc
                cy = v + (cz + abs(sz)) / 2
                gapx = Fr(r.choice([1, 10, 100]), 10 ** 4)
                cx = p[0] + side * (1 + gapx + abs(sz))
                S = Sq(cx, cy, tZ, tag="Zside")
                if ok_with(S, sq, k):
                    sq.append(S)
                    break
    # random squares (stepping stones / obstacles) of all kinds
    for trial in range(60):
        t = r.choice([Fr(0), tau * Fr(r.randint(-999, 999), 1000),
                      Fr(math.tan(float(b) * r.random() / 2)).limit_denominator(10 ** 9),
                      Fr(r.randint(-40, 40), 100)])
        S = Sq(Fr(r.randint(600, int(k) * 1000 - 600), 1000), Fr(r.randint(600, (s + 3) * 1000), 1000), t,
               tag="rnd")
        if ok_with(S, sq, k):
            sq.append(S)
    return Cfg(k, sq, "conseq_c%s_s%d_d%s_b%s_seed%d" % (c0, s, delta, b, seed))


def main():
    summ = dict(runs=[], n_cfg=0, n_Zs=0, n_P2f=0, n_traces=0, n_contacts_withZ=0, viol=[],
                kinds_withZ={}, back_q=0, back_found_forb=0, back_sanity=dict(checked=0, matched=0),
                min_slack_P2f=None, min_margin_Zs=None, anomalies=[])
    c0s = [Fr(3, 10), Fr(45, 100), Fr(6, 10)]
    if len(sys.argv) > 4 and sys.argv[4].startswith("c0="):
        c0s = [Fr(sys.argv[4][3:])]
    ss = [4, 9, 16, 25] if not QUICK else [9]
    deltas = [Fr(1, 1000), Fr(1, 100), Fr(3, 100)] if not QUICK else [Fr(1, 100)]
    bs = [Fr(1, 1000), Fr(1, 100)]
    eps = Fr(3, 10)
    seed = 0
    for c0 in c0s:
        for s in ss:
            for delta in deltas:
                for b in bs:
                    if tl() < 30:
                        break
                    seed += 1
                    alpha_true = mp.mpf(c0.numerator) / c0.denominator / mp.sqrt(s)
                    tau = Fr(mp.nstr(mp.tan(alpha_true / 2), 15)).limit_denominator(10 ** 9)
                    alpha = atan2t(tau)
                    w0 = Kcoef(alpha) * (s + 1 + mpf(b)) + mpf(delta) + f_of(mpf(b)) + mp.mpf("1e-12")
                    if w0 >= mp.mpf("0.5"):
                        summ["runs"].append(dict(c0=str(c0), s=s, delta=str(delta), b=str(b),
                                                 skipped="w0>=1/2", w0=float(w0)))
                        continue
                    for rep in range(2 if not QUICK else 1):
                        if tl() < 30:
                            break
                        cfg = make_cfg(c0, s, delta, b, eps, seed * 10 + rep, alpha, tau, w0)
                        ok, msg = cfg.validate()
                        if not ok:
                            summ["anomalies"].append("invalid " + cfg.name + " " + msg)
                            continue
                        summ["n_cfg"] += 1
                        h = Fr(s + 2)
                        Zs, P2f = set(), set()
                        for Z in cfg.sq:
                            zs, wit, p2f, slack = classify_Z(Z, s, eps, w0, b, alpha, delta)
                            if zs:
                                Zs.add(Z.idx)
                                # margin of the Z(s)-analogue: dist(y*) - w(s,Z) with q_y<=s+1+b
                                wsZ = Kcoef(alpha) * (s + 1 + mpf(b)) + mpf(delta) + f_of(Z.a_mp())
                                mg = mpf(wit[0]) - wsZ
                                if summ["min_margin_Zs"] is None or mg < summ["min_margin_Zs"][0]:
                                    summ["min_margin_Zs"] = (float(mg), cfg.name, Z.to_json())
                            if p2f:
                                P2f.add(Z.idx)
                                if summ["min_slack_P2f"] is None or slack < summ["min_slack_P2f"][0]:
                                    summ["min_slack_P2f"] = (float(slack), cfg.name, Z.to_json())
                        summ["n_Zs"] += len(Zs)
                        summ["n_P2f"] += len(P2f)
                        forb = Zs | P2f
                        # forward dense sampling
                        nx = int(cfg.k) * (60 if not QUICK else 20)
                        rr = random.Random(seed)
                        xs = [Fr(j, 60) + Fr(rr.randint(0, 10 ** 6), 10 ** 8) for j in range(1, int(cfg.k) * 60)]
                        if QUICK:
                            xs = xs[::3]
                        found_reach = []
                        for x in xs:
                            if tl() < 20:
                                break
                            if rss_mb() > 700:
                                raise SystemExit(4)
                            try:
                                res = trace(cfg, x, delta, tau, h)
                            except RuntimeError as e:
                                summ["anomalies"].append("trace err %s x=%s %s" % (cfg.name, x, e))
                                continue
                            summ["n_traces"] += 1
                            for i, rec in enumerate(res["contacts"]):
                                if rec["Z"] in forb:
                                    summ["n_contacts_withZ"] += 1
                                    kk = rec["kind"] + ("/" + rec["flag"] if rec["flag"] else "")
                                    summ["kinds_withZ"][kk] = summ["kinds_withZ"].get(kk, 0) + 1
                                    last = (i == len(res["contacts"]) - 1)
                                    onbot = rec["kind"] in ("bot", "vbot") and rec["eta"] == -HALF
                                    if onbot or not last or res["term"] not in ("E", "D", "W"):
                                        v = dict(cfg=cfg.name, x=str(x), Z=cfg.sq[rec["Z"]].to_json(),
                                                 kind=rec["kind"], flag=rec["flag"], term=res["term"],
                                                 inZs=rec["Z"] in Zs, q=[str(rec["q"][0]), str(rec["q"][1])])
                                        summ["viol"].append(v)
                            if any(pidx in forb for pidx in res["passed"]):
                                summ["viol"].append(dict(cfg=cfg.name, x=str(x), msg="passes forbidden Z"))
                        # backward exhaustive predecessor search from bot(Z) samples
                        B = Back(cfg, tau, delta, h)
                        for zi in sorted(forb):
                            if tl() < 20:
                                break
                            Z = cfg.sq[zi]
                            pts = [Z.pt(Fr(j, 40) - HALF, -HALF) for j in range(41)]
                            pts += [Z.pt(Fr(rr.randint(-10 ** 6, 10 ** 6), 2 * 10 ** 6), -HALF) for _ in range(10)]
                            for q in pts:
                                if q[1] > h:
                                    continue
                                summ["back_q"] += 1
                                out = []
                                B.nodes = 0
                                B.preds(q, delta, 0, [], out)
                                if out:
                                    summ["back_found_forb"] += 1
                                    if len(summ["viol"]) < 100:
                                        summ["viol"].append(dict(cfg=cfg.name, Z=Z.to_json(), q=[str(q[0]), str(q[1])],
                                                                 msg="backward search found path",
                                                                 x=str(out[0][0])))
                        # side contacts with forbidden squares: backward search from side points, then
                        # forward confirmation that the path terminates there (E/D/W) and never enters Z
                        for zi in sorted(forb):
                            if tl() < 20:
                                break
                            Z = cfg.sq[zi]
                            for xi_s in (-HALF, HALF):
                                for j in range(1, 16):
                                    q = Z.pt(xi_s, Fr(j, 16) - HALF)
                                    if q[1] > h:
                                        continue
                                    out = []
                                    B.nodes = 0
                                    B.preds(q, delta, 0, [], out, limit=3)
                                    for (xx, chain, gtot) in out[:2]:
                                        res = trace(cfg, xx, delta, tau, h)
                                        summ["side_checked"] = summ.get("side_checked", 0) + 1
                                        cz_ = [rec for rec in res["contacts"] if rec["Z"] == zi]
                                        if not cz_:
                                            summ["side_notreached"] = summ.get("side_notreached", 0) + 1
                                            if len(summ["anomalies"]) < 50:
                                                summ["anomalies"].append("side: backward found x=%s but forward did not reach Z (%s, term %s)" % (xx, cfg.name, res["term"]))
                                            continue
                                        rec = cz_[0]
                                        okk = (rec["q"] == q and rec["kind"] == "side" and rec is res["contacts"][-1]
                                               and res["term"] in ("E", "D", "W")
                                               and not any(pp == zi for pp in res["passed"]))
                                        if okk:
                                            summ["side_ok"] = summ.get("side_ok", 0) + 1
                                            kk = "side/" + res["term"]
                                            summ["kinds_withZ"][kk] = summ["kinds_withZ"].get(kk, 0) + 1
                                        else:
                                            summ["viol"].append(dict(cfg=cfg.name, x=str(xx), Z=Z.to_json(),
                                                                     msg="side contact not terminal/not side",
                                                                     kind=rec["kind"], term=res["term"]))
                        # sanity of the backward search: reproduce forward contacts on non-forbidden squares
                        chk = 0
                        for x in xs[::25]:
                            if tl() < 20 or chk > 40:
                                break
                            res = trace(cfg, x, delta, tau, h)
                            for rec in res["contacts"]:
                                if rec["kind"] == "bot" and rec["Z"] not in forb:
                                    out = []
                                    B.nodes = 0
                                    B.preds(rec["q"], delta, 0, [], out)
                                    summ["back_sanity"]["checked"] += 1
                                    if any(o[0] == x and o[1] == rec["passed"] for o in out):
                                        summ["back_sanity"]["matched"] += 1
                                    chk += 1
                        summ["runs"].append(dict(cfg=cfg.name, c0=str(c0), s=s, delta=str(delta), b=str(b),
                                                 alpha=float(alpha), w0=float(w0), nsq=len(cfg.sq),
                                                 nZs=len(Zs), nP2f=len(P2f), traces=len(xs),
                                                 elapsed=time.time() - T0))
                        print("[%6.1fs] %s nsq=%d Zs=%d P2f=%d viol=%d backq=%d" %
                              (time.time() - T0, cfg.name, len(cfg.sq), len(Zs), len(P2f), len(summ["viol"]),
                               summ["back_q"]), flush=True)
                        with open(out_path(OUT + "_summary.json"), "w", encoding="utf-8") as fh:
                            json.dump(summ, fh, indent=1, default=str)
    summ["elapsed"] = time.time() - T0
    summ["rss_mb"] = rss_mb()
    with open(out_path(OUT + "_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summ, fh, indent=1, default=str)
    print("done", summ["n_cfg"], summ["n_Zs"], summ["n_P2f"], len(summ["viol"]), summ["back_q"],
          summ["back_found_forb"], summ["back_sanity"])


if __name__ == "__main__":
    main()
