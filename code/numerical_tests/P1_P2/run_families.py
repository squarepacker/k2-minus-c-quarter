# -*- coding: utf-8 -*-
"""Driver: adversarial configurations for Theorem 4.2 (P1), Corollary 4.4 and
Lemma 4.8 / Theorem 4.11 (P2), traced in the flow without rule R1 (a superset
of the contacts of the flow with R1; see tracer.py).

Usage:  python run_families.py <budget_seconds> <prefix> [quick|full] [families]

  budget_seconds  wall-clock budget; the driver stops starting new work when it is used up
  prefix          output prefix; files out/<prefix>_summary.json and out/<prefix>_viol.jsonl
  quick           reduced parameter sets and family sizes (any other word, e.g. "full", = all)
  families        comma-separated subset of exact,floorZ,overhang,valley,pile,tower
                  (default: all six), or "fuzz" for the random-cluster fuzzer.

The fuzzer's random seed is sum(map(ord, prefix)); the recorded fuzz run used the
prefix "fuzz" (seed 463).
"""
import sys, time, json, random, math
from fractions import Fraction as Fr
import mpmath as mp
from tracer import (Sq, Cfg, trace, verify_contact, analyse_contact, rot, HALF,
                    disjoint, mpf, atan2t)
from common import out_path, rss_mb, place, ok_with, place_min_gap, build_column

T0 = time.time()
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 600
OUT = sys.argv[2] if len(sys.argv) > 2 else "run"
QUICK = len(sys.argv) > 3 and sys.argv[3] == "quick"
FAMS = set(sys.argv[4].split(",")) if len(sys.argv) > 4 else \
    {"exact", "floorZ", "overhang", "valley", "pile", "tower"}


def time_left():
    return BUDGET - (time.time() - T0)


def frac(xf, den=10 ** 9):
    return Fr(round(xf * den), den)


# ------------------------------------------------------------ bookkeeping
class Stats:
    def __init__(self):
        self.n_contacts = 0
        self.n_bot = 0
        self.n_vbot = 0
        self.n_traces = 0
        self.n_cfg = 0
        self.fail = []
        self.worst = {}
        self.kinds = {}
        self.terms = {}
        self.flags = {}
        self.anom = []
        self.aux_viol = []
        self.aux_ineq_false = 0
        self.fam = {}

    def upd(self, key, val, info):
        if val is None:
            return
        v = float(val)
        if key not in self.worst or v > self.worst[key][0]:
            self.worst[key] = (v, info)


ST = Stats()
OUTF = open(out_path(OUT + "_viol.jsonl"), "w", encoding="utf-8")


def info_of(cfg, rec, an, fam, par):
    Z = cfg.sq[rec["Z"]]
    return dict(fam=fam, cfg=cfg.name, par=par, x=str(rec["x"]), m=an["m"],
                qy=float(an["qy"]), g=float(an["g"]), Dtilt=float(an["D"]), Gvert=float(an["G"]),
                kind=rec["kind"], flag=rec["flag"], Z=Z.to_json(),
                aZ=float(Z.a_mp()), passed_t=[str(cfg.sq[i].t) for i in rec["passed"]][:6],
                maxdist=(float(an["maxdist"]) if an.get("onbot") else None))


def run_cfg(cfg, xs, tauT, delta, h, fam, par, check_verify=True, extra_tag=""):
    ok, msg = cfg.validate()
    if not ok:
        ST.anom.append(dict(fam=fam, cfg=cfg.name, msg="invalid cfg: " + msg))
        return []
    ST.n_cfg += 1
    fs = ST.fam.setdefault(fam, dict(cfg=0, traces=0, contacts=0, bot=0, worst_r_alphaT=0.0,
                                     worst_r_g=0.0, worst_r_eff=0.0, worst_tilt_eff=0.0,
                                     worst_gap=0.0))
    fs["cfg"] += 1
    recs_all = []
    for x in xs:
        if time_left() < 5:
            break
        if rss_mb() > 700:
            ST.anom.append(dict(msg="RSS>700MB abort"))
            raise SystemExit(4)
        x = Fr(x)
        if not (0 < x < cfg.k):
            continue
        try:
            res = trace(cfg, x, delta, tauT, h)
        except RuntimeError as e:
            ST.anom.append(dict(fam=fam, cfg=cfg.name, x=str(x), msg=str(e)))
            continue
        ST.n_traces += 1
        fs["traces"] += 1
        ST.terms[res["term"]] = ST.terms.get(res["term"], 0) + 1
        for rec in res["contacts"]:
            ST.n_contacts += 1
            fs["contacts"] += 1
            ST.kinds[rec["kind"]] = ST.kinds.get(rec["kind"], 0) + 1
            if rec["flag"]:
                ST.flags[rec["flag"]] = ST.flags.get(rec["flag"], 0) + 1
            if rec["kind"] in ("top", "inside"):
                ST.anom.append(dict(fam=fam, cfg=cfg.name, x=str(x), msg="contact kind " + rec["kind"]))
            if check_verify:
                fl = verify_contact(cfg, rec)
                if fl:
                    ST.fail.append(dict(fam=fam, cfg=cfg.name, x=str(x), fails=fl))
                    OUTF.write(json.dumps(dict(type="verify_fail", fam=fam, cfg=cfg.name,
                                               x=str(x), fails=fl)) + "\n")
            an = analyse_contact(cfg, rec, tauT, delta)
            inf = None
            bad = []
            if not an["id_qy"]:
                bad.append("id_qy")
            if not an["cor_low_T"]:
                bad.append("cor_low_T")
            if not an["cor_low_eff"]:
                bad.append("cor_low_eff")
            if not an["cor_up"]:
                bad.append("cor_up")
            if an["cor_050001"] is False and rec["flag"] != "D":
                bad.append("cor_050001")
            if an["onbot"]:
                ST.n_bot += 1
                fs["bot"] += 1
                if rec["kind"] == "vbot":
                    ST.n_vbot += 1
                if an["r_alphaT"] >= 1:
                    bad.append("P2_alphaT")
                if an["r_g"] > 1:
                    bad.append("P2_g")
                if an["r_eff"] > 1:
                    bad.append("P2_eff")
                if an["r_050001"] is not None and an["r_050001"] >= 1:
                    bad.append("P2_050001")
                if not an["contain_ok"]:
                    bad.append("contain")
                if not an["aux_ineq_holds"]:
                    ST.aux_ineq_false += 1
                if an["r_aux"] >= 1:
                    inf = info_of(cfg, rec, an, fam, par)
                    if len(ST.aux_viol) < 400:
                        ST.aux_viol.append(dict(info=inf, r_aux=float(an["r_aux"]),
                                                r_alphaT=float(an["r_alphaT"])))
            if bad:
                inf = inf or info_of(cfg, rec, an, fam, par)
                ST.fail.append(dict(type="claim_fail", bad=bad, info=inf))
                OUTF.write(json.dumps(dict(type="claim_fail", bad=bad, info=inf)) + "\n")
            # worst ratios
            if an["onbot"]:
                for key in ("r_alphaT", "r_g", "r_eff", "r_050001"):
                    if an[key] is not None and (key not in ST.worst or float(an[key]) > ST.worst[key][0]):
                        ST.upd(key, an[key], info_of(cfg, rec, an, fam, par))
                fs["worst_r_alphaT"] = max(fs["worst_r_alphaT"], float(an["r_alphaT"]))
                fs["worst_r_g"] = max(fs["worst_r_g"], float(an["r_g"]))
                fs["worst_r_eff"] = max(fs["worst_r_eff"], float(an["r_eff"]))
            for key in ("ratio_tilt_T", "ratio_tilt_eff", "ratio_gap", "ratio_gap_delta"):
                v = an[key]
                if key == "ratio_gap_delta" and rec["flag"] == "D":
                    key = "ratio_gap_delta_Dtie"
                if key not in ST.worst or float(v) > ST.worst[key][0]:
                    ST.upd(key, v, info_of(cfg, rec, an, fam, par))
            fs["worst_tilt_eff"] = max(fs["worst_tilt_eff"], float(an["ratio_tilt_eff"]))
            fs["worst_gap"] = max(fs["worst_gap"], float(an["ratio_gap"]))
        recs_all.append(res)
    return recs_all


def affine_targets(cfg, results, tauT, delta, h, zsel=None):
    """From traced results, find pairs of x in the same branch hitting relint bot(Z)
    and solve the affine map for xi = +-1/2 (vertex) and +-(1/2 - 1e-12)."""
    branches = {}
    for res in results:
        for rec in res["contacts"]:
            if rec["kind"] != "bot":
                continue
            if zsel is not None and rec["Z"] not in zsel:
                continue
            key = (tuple(rec["passed"]), rec["Z"])
            branches.setdefault(key, []).append((rec["x"], rec["xi"]))
    xs = []
    for key, lst in branches.items():
        if len(lst) < 2:
            continue
        lst = sorted(set(lst))
        (x1, xi1), (x2, xi2) = lst[0], lst[-1]
        if xi1 == xi2:
            continue
        for target in (-HALF, HALF, -HALF + Fr(1, 10 ** 12), HALF - Fr(1, 10 ** 12)):
            xt = x1 + (target - xi1) * (x2 - x1) / (xi2 - xi1)
            if 0 < xt < cfg.k:
                xs.append(xt)
    return xs


def run_with_targets(cfg, xs, tauT, delta, h, fam, par):
    r1 = run_cfg(cfg, xs, tauT, delta, h, fam, par)
    xt = affine_targets(cfg, r1, tauT, delta, h)
    if xt:
        run_cfg(cfg, xt, tauT, delta, h, fam + "+vert", par)
    return r1


# ------------------------------------------------------------ parameter sets
PARS = [
    # name, tau (alpha_T = 2 atan tau), delta
    ("a0.02_d3e-3", Fr(1, 100), Fr(3, 1000)),
    ("a0.05_d1e-2", Fr(1, 40), Fr(1, 100)),
    ("a0.1_d3e-2", Fr(1, 20), Fr(3, 100)),
    ("a0.2_d5e-2", Fr(1, 10), Fr(5, 100)),
    ("a0.2_d2e-3", Fr(1, 10), Fr(2, 1000)),
    ("a0.1_d1e-3", Fr(1, 20), Fr(1, 1000)),
    ("a1e-3_d1e-5", Fr(1, 2000), Fr(1, 10 ** 5)),
    ("theory_a1.5e-6_d1e-5", Fr(3, 4000000), Fr(1, 10 ** 5)),
]
if QUICK:
    PARS = PARS[1:2] + PARS[-1:]

rng = random.Random(20261006)


def tfr(tau, f):
    """rational strictly-inside tilt param: f*tau (f<1)"""
    return Fr(tau) * Fr(f).limit_denominator(10 ** 6)


def fam_tower(pname, tau, delta, m, mode, tZlist, xiZlist, gapfrac, seed):
    """mode: same / alt / rand.  Builds a column + target Z, traces dense x + vertex targets."""
    k = Fr(2 * m + 12)
    x0 = Fr(m + 6) + Fr(rng.randint(0, 10 ** 6), 10 ** 7)
    r = random.Random(seed)
    specs = []
    for i in range(m):
        if mode == "same":
            t = tfr(tau, Fr(999, 1000))
            xi = Fr(0)
        elif mode == "same_neg":
            t = -tfr(tau, Fr(999, 1000))
            xi = Fr(r.randint(-40, 40), 100)
        elif mode == "alt":
            t = tfr(tau, Fr(999, 1000)) * (1 if i % 2 == 0 else -1)
            xi = Fr(0)
        else:
            t = tfr(tau, Fr(r.randint(-999, 999), 1000))
            xi = Fr(r.randint(-45, 45), 100)
        specs.append((t, xi))
    gb = Fr(delta) * Fr(gapfrac).limit_denominator(10 ** 6) / (m + 1)
    gb = max(gb, Fr(1, 10 ** 14))
    gaps = [gb] * m
    if specs:
        # first square: make its lowest vertex >= 0 by a floor gap
        t1, xi1 = specs[0]
        c1, s1 = rot(t1)
        need = (HALF + (xi1 if s1 > 0 else -xi1)) * abs(s1)
        gaps[0] = max(gb, need + Fr(1, 10 ** 12) if need > 0 else gb)
    col = build_column(x0, specs, gaps, k, gcap=Fr(1, 2))
    if col is None:
        ST.anom.append(dict(fam="tower", msg="build failed %s %s" % (mode, m)))
        return
    sq, p, d, used = col
    for tZ in tZlist:
        for xiZ in xiZlist:
            if time_left() < 10 or time.time() - T_PAR > PAR_BUDGET:
                return
            S, gZ = place_min_gap(sq, p, d, tZ, xiZ, gb, k, Fr(1, 2))
            if S is None:
                continue
            S.tag = "Z"
            cfg = Cfg(k, sq + [S], "tower_%s_m%d_tZ%s_xi%s" % (mode, m, tZ, xiZ))
            # dense x around x0
            xs = [x0 + Fr(j, 20) - Fr(1, 2) for j in range(21)]
            xs += [x0 + Fr(r.randint(-10 ** 6, 10 ** 6), 2 * 10 ** 6) for _ in range(6)]
            run_with_targets(cfg, xs, tau, delta, cfg.k / 2 - 1, "tower_" + mode, pname)


def fam_floorZ(pname, tau, delta):
    k = Fr(12)
    for tZ in [Fr(0), tfr(tau, Fr(1, 2)), tfr(tau, Fr(999, 1000)), Fr(tau), Fr(1, 50), Fr(1, 10),
               Fr(1, 5), Fr(2, 5), -Fr(1, 10), -Fr(2, 5), Fr(41421, 100000), -Fr(41421, 100000)]:
        c, s = rot(tZ)
        for h0 in [Fr(0), Fr(delta) / 3, Fr(delta) * Fr(999, 1000), Fr(delta) * 2]:
            cy = h0 + (c + abs(s)) / 2
            S = Sq(Fr(6), cy, tZ, tag="Z")
            cfg = Cfg(k, [S], "floorZ_t%s_h%s" % (tZ, h0))
            xs = [Fr(6) + Fr(j, 30) - Fr(1, 1) for j in range(61)]
            # exact lowest vertex x and bot endpoints
            for v in (S.pt(-HALF, -HALF), S.pt(HALF, -HALF)):
                xs += [v[0], v[0] + Fr(1, 10 ** 12), v[0] - Fr(1, 10 ** 12)]
            run_cfg(cfg, xs, tau, delta, k / 2 - 1, "floorZ", pname)


def fam_overhang(pname, tau, delta, m, tZ, g1frac, sgn):
    """tower of m squares (tilt ~ +-tau, same sign) then a large-tilt Z overhanging the
    edge of the last tower square; path near that edge hits the upper vertex of bot(Z)."""
    k = Fr(2 * m + 14)
    x0 = Fr(m + 7)
    tt = tfr(tau, Fr(999, 1000)) * sgn
    specs = [(tt, Fr(0))] * m
    gb = max(Fr(delta) / (10 ** 6 * (m + 1)), Fr(1, 10 ** 14))
    gaps = [gb] * m
    if m:
        c1, s1 = rot(tt)
        gaps[0] = abs(s1) / 2 + Fr(1, 10 ** 12)
    col = build_column(x0, specs, gaps, k)
    if col is None:
        return
    sq, p, d, used = col
    # path near the edge of the last square: shift x so that exit is close to the
    # left (for phi_Z>0) / right (phi_Z<0) end of top(X_m).
    tZ = Fr(tZ)
    if m:
        L = sq[-1]
        endpt = L.pt(-HALF if tZ > 0 else HALF, HALF)
        # exit point ~ endpt + small inward along u
        eps = Fr(1, 10 ** 6)
        uu = L.u()
        sign_in = 1 if tZ > 0 else -1
        ex = (endpt[0] + sign_in * eps * uu[0], endpt[1] + sign_in * eps * uu[1])
        dd = L.n()
    else:
        X1 = Sq(x0, HALF, 0, tag="X1")
        sq = [X1]
        endpt = X1.pt(-HALF if tZ > 0 else HALF, HALF)
        eps = Fr(1, 10 ** 6)
        sign_in = 1 if tZ > 0 else -1
        ex = (endpt[0] + sign_in * eps, endpt[1])
        dd = (Fr(0), Fr(1))
    g1 = max(Fr(delta) * Fr(g1frac).limit_denominator(10 ** 6), Fr(1, 10 ** 13))
    xiZ = (HALF - Fr(1, 10 ** 9)) if tZ > 0 else (-HALF + Fr(1, 10 ** 9))
    S, gZ = place_min_gap(sq, ex, dd, tZ, xiZ, g1, k, Fr(1, 2))
    if S is None:
        return
    S.tag = "Z"
    cfg = Cfg(k, sq + [S], "overhang_m%d_tZ%s_s%d" % (m, tZ, sgn))
    # x of the designed path: invert through the column (trace two nearby x and solve)
    # simplest: design x by back-substitution is complex; use dense sampling + affine targets
    base = x0 if m else x0
    xs = [base + Fr(j, 200) - Fr(1, 2) for j in range(201)]
    run_with_targets(cfg, xs, tau, delta, cfg.k / 2 - 1, "overhang", pname)


def fam_valley(pname, tau, delta):
    """two squares tilted -th (left) and +th (right) forming a V; Z above."""
    k = Fr(12)
    th = tfr(tau, Fr(1, 2))
    c, s = rot(th)
    for epsf in (Fr(1, 10), Fr(1, 2)):
        J = (Fr(6), Fr(1) + abs(s) + Fr(delta) / 4)
        eps = Fr(delta) * Fr(th) * epsf
        X1 = Sq(*_corner_center(J[0] - eps, J[1], -th, "tr"), -th, tag="X1")
        X2 = Sq(*_corner_center(J[0] + eps, J[1], th, "tl"), th, tag="X2")
        for tZ in (Fr(0), th / 3, -th / 3, Fr(1, 7)):
            for hz in (Fr(delta) / 3, Fr(delta) * 2 / 3):
                cz, sz = rot(tZ)
                # Z with bottom near height J_y + hz centred above J
                S = place((J[0], J[1]), (Fr(0), Fr(1)), hz, tZ, Fr(0))
                S.tag = "Z"
                cfg = Cfg(k, [X1, X2, S], "valley_e%s_tZ%s_h%s" % (epsf, tZ, hz))
                ok, msg = cfg.validate()
                if not ok:
                    continue
                xs = [Fr(6) + Fr(j, 100) - Fr(3, 2) for j in range(301)]
                res = run_cfg(cfg, xs, tau, delta, k / 2 - 1, "valley", pname)
                # exact merge points: pairs of branches (passed X1) and (passed X2) hitting Z
                br = {}
                for rr in res:
                    for rec in rr["contacts"]:
                        if rec["kind"] == "bot" and rec["Z"] == 2:
                            br.setdefault(tuple(rec["passed"]), []).append((rec["x"], rec["xi"]))
                keys = [kk for kk in br if len(br[kk]) >= 2]
                xm = []
                for i1 in range(len(keys)):
                    for i2 in range(i1 + 1, len(keys)):
                        A = sorted(set(br[keys[i1]]))
                        B = sorted(set(br[keys[i2]]))
                        (a1, ax1), (a2, ax2) = A[0], A[-1]
                        (b1, bx1), (b2, bx2) = B[0], B[-1]
                        if ax1 == ax2 or bx1 == bx2:
                            continue
                        lo = max(min(ax1, ax2), min(bx1, bx2))
                        hi = min(max(ax1, ax2), max(bx1, bx2))
                        if lo >= hi:
                            continue
                        for xi_t in (lo + (hi - lo) / 3, (lo + hi) / 2):
                            xm.append(a1 + (xi_t - ax1) * (a2 - a1) / (ax2 - ax1))
                            xm.append(b1 + (xi_t - bx1) * (b2 - b1) / (bx2 - bx1))
                if xm:
                    rr2 = run_cfg(cfg, xm, tau, delta, k / 2 - 1, "valley+merge", pname)
                    # count exact merges (same q on Z from different branches)
                    qs = {}
                    for rr in rr2:
                        for rec in rr["contacts"]:
                            if rec["Z"] == 2 and rec["kind"] == "bot":
                                qs.setdefault(rec["q"], set()).add(tuple(rec["passed"]))
                    nm = sum(1 for q in qs if len(qs[q]) >= 2)
                    ST.flags["exact_merge_points"] = ST.flags.get("exact_merge_points", 0) + nm


def _corner_center(x, y, t, which):
    c, s = rot(t)
    u = (c, s)
    n = (-s, c)
    if which == "tr":   # corner = c + u/2 + n/2
        return x - u[0] / 2 - n[0] / 2, y - u[1] / 2 - n[1] / 2
    else:               # tl: corner = c - u/2 + n/2
        return x + u[0] / 2 - n[0] / 2, y + u[1] / 2 - n[1] / 2


def fam_pile(pname, tau, delta, ncol, m, nobs, seed):
    r = random.Random(seed)
    k = Fr(max(2 * m + 8, 2 * ncol + 6))
    sq = []
    for cidx in range(ncol):
        x0 = Fr(2 + 2 * cidx) + Fr(r.randint(0, 999), 2000)
        specs = []
        for i in range(m):
            t = tfr(tau, Fr(r.randint(-999, 999), 1000))
            if r.random() < 0.2:
                t = Fr(0)
            specs.append((t, Fr(r.randint(-30, 30), 100)))
        gb = Fr(delta) * Fr(r.randint(1, 900), 1000) / (m + 1)
        gaps = [gb] * m
        c1, s1 = rot(specs[0][0])
        need = (HALF + abs(specs[0][1])) * abs(s1)
        gaps[0] = gb + need
        col = build_column(x0, specs, gaps, k, existing=sq, gcap=Fr(1, 2))
        if col is None:
            continue
        sq += col[0]
    # obstacles at random positions / fractional heights / random tilts (incl. large)
    tries = 0
    while nobs > 0 and tries < 400:
        tries += 1
        t = r.choice([Fr(0), tfr(tau, Fr(r.randint(-999, 999), 1000)), Fr(r.randint(-41, 41), 100)])
        S = Sq(Fr(r.randint(1000, int(k) * 1000 - 1000), 1000), Fr(r.randint(1000, int(k) * 500), 1000), t,
               tag="obs")
        if ok_with(S, sq, k):
            sq.append(S)
            nobs -= 1
    cfg = Cfg(k, sq, "pile_s%d" % seed)
    xs = [Fr(j, 50) + Fr(r.randint(0, 999), 10 ** 6) for j in range(1, int(k) * 50)]
    res = run_cfg(cfg, xs, tau, delta, k / 2 - 1, "pile", pname)
    xt = affine_targets(cfg, res, tau, delta, k / 2 - 1)
    if xt:
        run_cfg(cfg, xt[:400], tau, delta, k / 2 - 1, "pile+vert", pname)
    # reflected configuration = ceiling flow of the original
    cr = cfg.reflected()
    run_cfg(cr, xs[::3], tau, delta, k / 2 - 1, "pile_refl", pname)


def fam_exact_h_and_Dtie(pname, tau, delta, m):
    """designed contact on bot(Z); rerun with h = q_y exactly and with delta' = g exactly."""
    k = Fr(2 * m + 12)
    x0 = Fr(m + 6)
    tt = tfr(tau, Fr(9, 10))
    specs = [(tt if i % 3 else -tt, Fr(0)) for i in range(m)]
    gb = Fr(delta) / (4 * (m + 1))
    gaps = [gb] * m
    c1, s1 = rot(specs[0][0])
    gaps[0] = gb + abs(s1) / 2
    col = build_column(x0, specs, gaps, k)
    if col is None:
        return
    sq, p, d, used = col
    for tZ in (tt / 2, Fr(1, 5), -tt):
        S, gZ = place_min_gap(sq, p, d, tZ, Fr(1, 7), gb, k, Fr(1, 2))
        if S is None:
            continue
        S.tag = "Z"
        cfg = Cfg(k, sq + [S], "exact_m%d_tZ%s" % (m, tZ))
        res = trace(cfg, x0, delta, tau, k / 2 - 1)
        hit = [rec for rec in res["contacts"] if rec["Z"] == len(sq)]
        if not hit:
            continue
        rec = hit[0]
        hq = rec["q"][1]
        gq = rec["g"]
        # (1) height cap exactly at the contact height: contact must be processed
        run_cfg(cfg, [x0, x0 + Fr(1, 10 ** 6), x0 - Fr(1, 10 ** 6)], tau, delta, hq, "exact_h", pname)
        r2 = trace(cfg, x0, delta, tau, hq)
        lastc = r2["contacts"][-1] if r2["contacts"] else None
        if lastc is None or lastc["Z"] != len(sq):
            ST.anom.append(dict(fam="exact_h", msg="contact at h not processed", cfg=cfg.name))
        # (2) death tie: delta' = g at the contact
        if gq > 0:
            run_cfg(cfg, [x0], tau, gq, k / 2 - 1, "Dtie", pname)
            r3 = trace(cfg, x0, gq, tau, k / 2 - 1)
            if r3["term"] != "D" or not r3["contacts"] or r3["contacts"][-1]["flag"] != "D":
                ST.anom.append(dict(fam="Dtie", msg="D tie not produced", cfg=cfg.name,
                                    term=r3["term"]))


def fam_fuzz(nconf, seed):
    """random small clusters, tilts over the whole range, alpha_T up to 2 atan(0.4) (~43.6 deg),
    delta in {0.01,0.05,0.2}; squares stacked loosely so that paths pass several of them."""
    r = random.Random(seed)
    for ic in range(nconf):
        if time_left() < 15:
            return
        k = Fr(8)
        tau = r.choice([Fr(1, 100), Fr(1, 20), Fr(1, 5), Fr(2, 5)])
        delta = r.choice([Fr(1, 100), Fr(1, 20), Fr(1, 5)])
        sq = []
        # loose columns: place squares above each other with random small gaps (< delta)
        for col in range(r.randint(1, 3)):
            x0 = Fr(r.randint(1500, 6500), 1000)
            p = (x0, Fr(0))
            d = (Fr(0), Fr(1))
            for lvl in range(r.randint(1, 3)):
                if r.random() < 0.7:
                    t = tau * Fr(r.randint(-999, 999), 1000)
                else:
                    t = Fr(r.randint(-410, 410), 1000)
                xi = Fr(r.randint(-45, 45), 100)
                c_, s_ = rot(t)
                gmin = (HALF + (xi if s_ > 0 else -xi)) * abs(s_) if lvl == 0 else Fr(0)
                S, g = place_min_gap(sq, p, d, t, xi, gmin + delta * Fr(r.randint(1, 900), 1000),
                                     k, Fr(1, 2))
                if S is None:
                    break
                sq.append(S)
                q = S.pt(xi, -HALF)
                nn = S.n()
                p = (q[0] + nn[0], q[1] + nn[1])
                d = nn
                if abs(t) >= tau:
                    break
        for _ in range(r.randint(0, 6)):
            t = Fr(r.randint(-414, 414), 1000)
            S = Sq(Fr(r.randint(700, 7300), 1000), Fr(r.randint(700, 3300), 1000), t, tag="fz")
            if ok_with(S, sq, k):
                sq.append(S)
        if not sq:
            continue
        cfg = Cfg(k, sq, "fuzz_%d_%d" % (seed, ic))
        xs = [Fr(j, 40) + Fr(r.randint(0, 999), 10 ** 6) for j in range(1, 320)]
        res = run_cfg(cfg, xs, tau, delta, k / 2 - 1, "fuzz", "tau%s_d%s" % (tau, delta))
        xt = affine_targets(cfg, res, tau, delta, k / 2 - 1)
        if xt:
            run_cfg(cfg, xt[:200], tau, delta, k / 2 - 1, "fuzz+vert", "tau%s_d%s" % (tau, delta))


# ------------------------------------------------------------ main
T_PAR = time.time()
PAR_BUDGET = 1e9


def main():
    global T_PAR, PAR_BUDGET
    if FAMS == {"fuzz"}:
        fam_fuzz(10 ** 6, seed=sum(map(ord, OUT)))
        print("[%7.1fs] fuzz done; cfg=%d traces=%d contacts=%d bot=%d vbot=%d fails=%d anom=%d" %
              (time.time() - T0, ST.n_cfg, ST.n_traces, ST.n_contacts, ST.n_bot, ST.n_vbot,
               len(ST.fail), len(ST.anom)), flush=True)
        return
    for ip, (pname, tau, delta) in enumerate(PARS):
        if time_left() < 20:
            break
        T_PAR = time.time()
        PAR_BUDGET = time_left() / (len(PARS) - ip)
        dl = lambda: time_left() < 20 or time.time() - T_PAR > PAR_BUDGET
        print("[%7.1fs] params %s (budget %.0fs)" % (time.time() - T0, pname, PAR_BUDGET), flush=True)
        if "exact" in FAMS:
            for m in (2, 9):
                fam_exact_h_and_Dtie(pname, tau, delta, m)
        if "floorZ" in FAMS:
            fam_floorZ(pname, tau, delta)
        for m in (([0, 1, 4, 12] if not QUICK else [0, 4]) if "overhang" in FAMS else []):
            for tZ in (Fr(1, 100), Fr(1, 20), Fr(1, 10), Fr(1, 5), Fr(3, 10), Fr(2, 5), -Fr(1, 10), -Fr(3, 10)):
                for g1f in (Fr(1, 10 ** 6), Fr(1, 2)):
                    if dl():
                        break
                    fam_overhang(pname, tau, delta, m, tZ, g1f, 1 if tZ > 0 else -1)
                    fam_overhang(pname, tau, delta, m, tZ, g1f, -1 if tZ > 0 else 1)
        print("[%7.1fs]   targeted families done; contacts=%d bot=%d fails=%d" %
              (time.time() - T0, ST.n_contacts, ST.n_bot, len(ST.fail)), flush=True)
        if not dl() and "valley" in FAMS:
            fam_valley(pname, tau, delta)
        if not dl() and "pile" in FAMS:
            fam_pile(pname, tau, delta, ncol=5, m=10, nobs=25, seed=100 * len(pname))
        mlist = [1, 4, 16] if not QUICK else [3, 8]
        if "tower" not in FAMS:
            mlist = []
        if pname.startswith("a1e-3"):
            mlist = [8, 40] if not QUICK else [8]
        tZs = [Fr(0), tfr(tau, Fr(1, 3)), -tfr(tau, Fr(999, 1000)), tfr(tau, Fr(999, 1000)), Fr(tau),
               Fr(1, 10), -Fr(3, 10), Fr(41, 100)]
        xiZs = [Fr(0), Fr(2, 5)]
        combos = []
        for mode in ("same", "alt", "rand", "same_neg"):
            for m in mlist[::-1]:
                combos.append((m, mode))
        # balanced order: interleave long and short towers, all modes early
        combos.sort(key=lambda z: (["same", "alt", "rand", "same_neg"].index(z[1]) + mlist[::-1].index(z[0])) % 4)
        for (m, mode) in combos:
            if dl():
                break
            fam_tower(pname, tau, delta, m, mode, tZs[:6] if m >= 16 else tZs, xiZs[:1] if m >= 16 else xiZs,
                      Fr(1, 10 ** 6), seed=m * 7 + len(mode))
            if dl():
                break
            fam_tower(pname, tau, delta, m, mode, tZs[:3], xiZs[:1], Fr(9, 10), seed=m * 11 + len(mode))
            print("[%7.1fs]   tower m=%d %s done; contacts=%d bot=%d fails=%d" %
                  (time.time() - T0, m, mode, ST.n_contacts, ST.n_bot, len(ST.fail)), flush=True)
        for seed in (range(1, 3) if "pile" in FAMS else []):
            if dl():
                break
            fam_pile(pname, tau, delta, ncol=5, m=10, nobs=25, seed=seed + 100 * len(pname))
        print("[%7.1fs] %s done; cfg=%d traces=%d contacts=%d bot=%d vbot=%d fails=%d anom=%d" %
              (time.time() - T0, pname, ST.n_cfg, ST.n_traces, ST.n_contacts, ST.n_bot, ST.n_vbot,
               len(ST.fail), len(ST.anom)), flush=True)
        dump()
    dump()


def dump():
    summ = dict(elapsed=time.time() - T0, n_cfg=ST.n_cfg, n_traces=ST.n_traces,
                n_contacts=ST.n_contacts, n_bot=ST.n_bot, n_vbot=ST.n_vbot,
                kinds=ST.kinds, terms=ST.terms, flags=ST.flags,
                n_fail=len(ST.fail), fails=ST.fail[:50], anomalies=ST.anom[:50],
                n_anom=len(ST.anom),
                aux_ineq_false_contacts=ST.aux_ineq_false,
                n_aux_viol=len(ST.aux_viol),
                aux_viol_examples=sorted(ST.aux_viol, key=lambda z: -z["r_aux"])[:5],
                worst={k_: dict(val=v[0], info=v[1]) for k_, v in ST.worst.items()},
                families=ST.fam, rss_mb=rss_mb())
    with open(out_path(OUT + "_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summ, f, indent=1, default=str)


if __name__ == "__main__":
    try:
        main()
    finally:
        dump()
        OUTF.close()
        print("done elapsed %.1fs rss %.1f MB" % (time.time() - T0, rss_mb()))
