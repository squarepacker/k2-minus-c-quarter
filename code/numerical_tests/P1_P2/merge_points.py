# -*- coding: utf-8 -*-
"""Exact merge points (the situation in which rule R1 decides).
Valley: X1 tilted -th (leans right), X2 tilted +th (leans left), lowest vertices near the floor next
to the junction; paths exiting near the facing top corners cross; Z above the crossing.
Exact merge points q on bot(Z) reached by a path through X1 and a path through X2 are computed
from the affine branch maps, traced, and both contacts are checked with Corollary 4.4 and
Lemma 4.8 (Theorem 4.11(a) applies to every entry candidate reached alive, before R1 is applied);
the R1 winner by lexicographic (g, x) is recorded.  Also sweeps Z tilt and height.

Usage:  python merge_points.py        (writes out/merge_points_out.json; about 1 s)
"""
import json, time
from fractions import Fraction as Fr
import mpmath as mp
from tracer import Sq, Cfg, trace, verify_contact, analyse_contact, rot, HALF, mpf
from common import out_path

T0 = time.time()
out = dict(cases=[], n_merge_points=0, n_merge_contacts_checked=0, fails=[], worst_r_alphaT=0.0, worst_r_g=0.0)
for (tau, delta) in ((Fr(1, 20), Fr(1, 100)), (Fr(1, 40), Fr(3, 1000)), (Fr(1, 2000), Fr(1, 10 ** 5))):
    for thf in (Fr(1, 2), Fr(9, 10)):
        th = min(tau * thf, delta / 8 * thf)
        c, s = rot(th)          # s = sin(angle) > 0, angle = 2 atan(th) < alpha_T
        k = Fr(10)
        Jx = Fr(5)
        eps2 = th * th / 16
        eps = s + eps2
        efl = delta / 100       # floor gap at the lowest vertices
        # X1: phi=-th, lowest vertex = bottom-right (xi=+1/2, eta=-1/2) at (Jx-eps, efl)
        c1, s1 = rot(-th)
        u1 = (c1, s1); n1 = (-s1, c1)
        X1c = (Jx - eps - u1[0] / 2 + n1[0] / 2, efl - u1[1] / 2 + n1[1] / 2)
        X1 = Sq(X1c[0], X1c[1], -th, tag="X1")
        c2, s2 = rot(th)
        u2 = (c2, s2); n2 = (-s2, c2)
        X2c = (Jx + eps + u2[0] / 2 + n2[0] / 2, efl + u2[1] / 2 + n2[1] / 2)
        X2 = Sq(X2c[0], X2c[1], th, tag="X2")
        ycorner = X1.pt(HALF, HALF)[1]
        ycross = ycorner + eps2 / (s / c)
        for tZ in (Fr(0), th / 4, -th / 4, th / 2):
            for dz in (Fr(1, 4), Fr(1, 2), Fr(3, 4)):
                cz, sz = rot(tZ)
                yZ = ycorner + th * (1 + dz) * Fr(11, 10) + abs(sz)
                # Z centred above Jx with its lowest point at yZ
                Z = Sq(Jx, yZ + (cz + abs(sz)) / 2, tZ, tag="Z")
                cfg = Cfg(k, [X1, X2, Z], "merge_tau%s_th%s_tZ%s_dz%s" % (tau, thf, tZ, dz))
                ok, msg = cfg.validate()
                if not ok:
                    out["cases"].append(dict(cfg=cfg.name, invalid=msg))
                    continue
                h = k / 2 - 1
                # branch samples: x just left / right of the junction under X1 / X2 lowest vertices
                br = {}
                xs = [Jx - eps - Fr(j, 4000) * delta for j in range(1, 40)] + \
                     [Jx + eps + Fr(j, 4000) * delta for j in range(1, 40)]
                for x in xs:
                    res = trace(cfg, x, delta, tau, h)
                    for rec in res["contacts"]:
                        if rec["Z"] == 2 and rec["kind"] == "bot":
                            br.setdefault(tuple(rec["passed"]), []).append((x, rec["xi"]))
                keys = sorted(br)
                nm = 0
                if (0,) in br and (1,) in br and len(br[(0,)]) >= 2 and len(br[(1,)]) >= 2:
                    A = sorted(set(br[(0,)])); B = sorted(set(br[(1,)]))
                    (a1, ax1), (a2, ax2) = A[0], A[-1]
                    (b1, bx1), (b2, bx2) = B[0], B[-1]
                    lo = max(min(ax1, ax2), min(bx1, bx2))
                    hi = min(max(ax1, ax2), max(bx1, bx2))
                    if lo < hi:
                        for frac_ in (Fr(1, 7), Fr(1, 2), Fr(6, 7)):
                            xi_t = lo + (hi - lo) * frac_
                            xa = a1 + (xi_t - ax1) * (a2 - a1) / (ax2 - ax1)
                            xb = b1 + (xi_t - bx1) * (b2 - b1) / (bx2 - bx1)
                            ra = trace(cfg, xa, delta, tau, h)
                            rb = trace(cfg, xb, delta, tau, h)
                            ca = [r_ for r_ in ra["contacts"] if r_["Z"] == 2]
                            cb = [r_ for r_ in rb["contacts"] if r_["Z"] == 2]
                            if not ca or not cb or ca[0]["q"] != cb[0]["q"]:
                                out["fails"].append(dict(cfg=cfg.name, msg="merge point not reproduced"))
                                continue
                            nm += 1
                            out["n_merge_points"] += 1
                            winner = min((ca[0]["g"], xa), (cb[0]["g"], xb))
                            for rec in (ca[0], cb[0]):
                                fl = verify_contact(cfg, rec)
                                an = analyse_contact(cfg, rec, tau, delta)
                                out["n_merge_contacts_checked"] += 1
                                bad = fl + [kk for kk in ("cor_low_T", "cor_up", "contain_ok") if not an[kk]]
                                if an["r_alphaT"] >= 1 or an["r_g"] > 1:
                                    bad.append("P2")
                                if bad:
                                    out["fails"].append(dict(cfg=cfg.name, x=str(rec["x"]), bad=bad))
                                out["worst_r_alphaT"] = max(out["worst_r_alphaT"], float(an["r_alphaT"]))
                                out["worst_r_g"] = max(out["worst_r_g"], float(an["r_g"]))
                            out["cases"].append(dict(cfg=cfg.name, q=[str(ca[0]["q"][0]), str(ca[0]["q"][1])],
                                                     xa=str(xa), xb=str(xb), ga=float(ca[0]["g"]), gb=float(cb[0]["g"]),
                                                     winner_x=str(winner[1]),
                                                     dir_angle_diff=float(2 * mp.atan(mpf(th)) * 2)))
                if nm == 0:
                    out["cases"].append(dict(cfg=cfg.name, note="no overlap", branches=[list(kk) for kk in keys]))
out["elapsed"] = time.time() - T0
with open(out_path("merge_points_out.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, default=str)
print(json.dumps(dict(n_merge_points=out["n_merge_points"], checked=out["n_merge_contacts_checked"],
                      fails=out["fails"][:5], worst_r_alphaT=out["worst_r_alphaT"], worst_r_g=out["worst_r_g"],
                      n_cases=len(out["cases"]),
                      invalid=sum(1 for c_ in out["cases"] if "invalid" in c_),
                      nooverlap=sum(1 for c_ in out["cases"] if "note" in c_),
                      elapsed=out["elapsed"]), indent=1))
