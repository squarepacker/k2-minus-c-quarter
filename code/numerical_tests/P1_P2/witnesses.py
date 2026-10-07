# -*- coding: utf-8 -*-
"""Explicit near-extremal configurations (exact coordinates) for Corollary 4.4,
Lemma 4.8 and Theorem 4.11(a).  They show how close the ratios of the tests come to 1.
E1  P2, lower end of the upper ramp: a square Z with large a(Z) overhanging an
    axis-parallel floor square (alpha_T = 2 atan(3/4000000) ~ 1.5e-6, delta = 1e-5)
E2  P2, upper end of the lower ramp: Z next to the floor with tiny tilt, gap ~ delta
    (same alpha_T and delta)
E3  tilt side of Corollary 4.4 and the 0.50001 form of Theorem 4.11(a) after a long
    leaning tower at (almost) maximal tilt, alpha ~ 1e-3, tiny delta, Z overhanging the
    edge of the last square
E4  as E3 with alternating tilts (gaps forced by disjointness); keys "E3_alt_..."
Each case is traced at the exact vertex x (contact at an endpoint of bot(Z)) and,
for E1 and E2, also 1e-18 inside (contact in the relative interior of bot(Z)).

Usage:  python witnesses.py          (writes out/witnesses_out.json; about 1 s)
"""
import json, time, sys
from fractions import Fraction as Fr
import mpmath as mp
from tracer import Sq, Cfg, trace, verify_contact, analyse_contact, rot, HALF, mpf, atan2t
from common import out_path, place, place_min_gap, build_column

T0 = time.time()
res = {}


def report(name, cfg, x, tau, delta, h=None, zidx=None):
    ok, msg = cfg.validate()
    h = cfg.k / 2 - 1 if h is None else h
    r = trace(cfg, x, delta, tau, h)
    out = dict(valid=ok, msg=msg, term=r["term"], x=str(x), nsq=len(cfg.sq))
    recs = [c for c in r["contacts"] if (zidx is None or c["Z"] == zidx)]
    if not recs:
        out["note"] = "target not reached"
        return out
    rc = recs[-1]
    an = analyse_contact(cfg, rc, tau, delta)
    out.update(kind=rc["kind"], m=an["m"], qy=mp.nstr(mpf(an["qy"]), 20), g=mp.nstr(mpf(an["g"]), 12),
               Dtilt=mp.nstr(mpf(an["D"]), 15), Gvert=mp.nstr(mpf(an["G"]), 12),
               verify=verify_contact(cfg, rc),
               ratio_tilt_T=mp.nstr(an["ratio_tilt_T"], 12), ratio_tilt_eff=mp.nstr(an["ratio_tilt_eff"], 12),
               cor_ok=bool(an["cor_low_T"] and an["cor_up"]))
    if an["onbot"]:
        out.update(maxdist=mp.nstr(mpf(an["maxdist"]), 20), argy=mp.nstr(mpf(an["argy"]), 20),
                   w_alphaT=mp.nstr(an["w_alphaT"], 20), r_alphaT=mp.nstr(an["r_alphaT"], 15),
                   one_minus_r_alphaT=mp.nstr(1 - an["r_alphaT"], 6),
                   r_g=mp.nstr(an["r_g"], 15), r_050001=(mp.nstr(an["r_050001"], 15) if an["r_050001"] is not None else None),
                   r_aux=mp.nstr(an["r_aux"], 15), contain_ok=an["contain_ok"],
                   Z=cfg.sq[rc["Z"]].to_json(), aZ=mp.nstr(cfg.sq[rc["Z"]].a_mp(), 12))
    return out


# ---------------- E1
tau = Fr(3, 4000000)
delta = Fr(1, 10 ** 5)
for tZ in (Fr(1, 5), Fr(1, 10), Fr(1, 100)):
    k = Fr(8)
    X1 = Sq(Fr(4), HALF, 0, tag="X1")
    c, s = rot(tZ)
    g1 = Fr(1, 10 ** 15)
    eps = g1 / (2 * (s / c))
    xx = Fr(4) - HALF + eps
    q = (xx, 1 + g1)
    Z = Sq(q[0] - s / 2 - HALF * c, q[1] + c / 2 - HALF * s, tZ, tag="Z")
    cfg = Cfg(k, [X1, Z], "E1")
    res["E1_vertex_tZ" + str(tZ)] = report("E1", cfg, xx, tau, delta, zidx=1)
    res["E1_relint_tZ" + str(tZ)] = report("E1", cfg, xx - Fr(1, 10 ** 18), tau, delta, zidx=1)

# ---------------- E2
for tZ in (Fr(1, 10 ** 7), Fr(1, 10 ** 5)):
    k = Fr(8)
    c, s = rot(tZ)
    h0 = delta - Fr(1, 10 ** 15)
    Z = Sq(Fr(4), h0 + (c + s) / 2, tZ, tag="Z")
    cfg = Cfg(k, [Z], "E2")
    low = Z.pt(-HALF, -HALF)    # lowest vertex for t>0
    res["E2_lowvertex_tZ" + str(tZ)] = report("E2", cfg, low[0], tau, delta, zidx=0)
    res["E2_relint_tZ" + str(tZ)] = report("E2", cfg, low[0] + Fr(1, 10 ** 18), tau, delta, zidx=0)

# ---------------- E3 / E4
for (mode, m, tau3, delta3) in (("same", 60, Fr(1, 2000), Fr(1, 10 ** 9)),
                                ("same", 25, Fr(1, 20), Fr(1, 10 ** 6)),
                                ("alt", 12, Fr(1, 2000), Fr(1, 10 ** 4))):
    k = Fr(2 * m + 14)
    x0 = Fr(m + 7)
    tt = tau3 * Fr(999999, 1000000)
    specs = [(tt if (mode == "same" or i % 2 == 0) else -tt, Fr(0)) for i in range(m)]
    gb = Fr(1, 10 ** 15)
    gaps = [gb] * m
    c1, s1 = rot(tt)
    gaps[0] = s1 / 2 + gb
    col = build_column(x0, specs, gaps, k, gcap=Fr(1, 2))
    if col is None:
        res["E3_%s_m%d" % (mode, m)] = "build failed"
        continue
    sq, p, d, used = col
    L = sq[-1]
    for tZ in (tt / 2, Fr(1, 10)):
        # exit near the left end of top(L) (Z with phi>0 overhangs to the left)
        endpt = L.pt(-HALF, HALF)
        e_in = Fr(1, 10 ** 9)
        ex = (endpt[0] + e_in * L.u()[0], endpt[1] + e_in * L.u()[1])
        S, gZ = place_min_gap(sq, ex, L.n(), tZ, HALF - Fr(1, 10 ** 12), gb, k, Fr(1, 2))
        if S is None:
            res["E3_%s_m%d_tZ%s" % (mode, m, tZ)] = "Z placement failed"
            continue
        S.tag = "Z"
        cfg = Cfg(k, sq + [S], "E3")
        # find x hitting the exit point ex: invert the column map using two traces
        xa, xb = x0, x0 + Fr(1, 10 ** 6)
        ra = trace(Cfg(k, sq, "col"), xa, Fr(1), tau3, k / 2 - 1)
        rb = trace(Cfg(k, sq, "col"), xb, Fr(1), tau3, k / 2 - 1)
        # positions along u of the last exit
        def last_exit(rr):
            st = rr["contacts"][-1] if rr["contacts"] else None
            return None
        # use contacts with last square: its entry xi
        ca_ = [c_ for c_ in ra["contacts"] if c_["Z"] == len(sq) - 1]
        cb_ = [c_ for c_ in rb["contacts"] if c_["Z"] == len(sq) - 1]
        if not ca_ or not cb_:
            res["E3_%s_m%d_tZ%s" % (mode, m, tZ)] = "column not traversed"
            continue
        xia, xib = ca_[0]["xi"], cb_[0]["xi"]
        target_xi = -HALF + e_in
        xt = xa + (target_xi - xia) * (xb - xa) / (xib - xia)
        key = "E3_%s_m%d_tau%s_tZ%s" % (mode, m, tau3, tZ)
        res[key] = report("E3", cfg, xt, tau3, delta3, zidx=len(sq))
        res[key]["sum_gaps_column"] = mp.nstr(mpf(sum(used, Fr(0))), 6)

res["elapsed"] = time.time() - T0
with open(out_path("witnesses_out.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, indent=1, default=str)
print(json.dumps(res, indent=1, default=str)[:12000])
