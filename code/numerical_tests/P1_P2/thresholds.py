# -*- coding: utf-8 -*-
"""Thresholds related to Lemma 4.6 (the auxiliary inequality) and Theorem 4.11(b):
 (1) a* : exact threshold of the auxiliary inequality
        w(s,Z) <= 0.50001 c0^2 (1+(q_y-s)/s) + delta + 1.0001 a(Z)
     <=>  f(a) := sin a + 1 - cos a <= 1.0001 a      (first terms and delta cancel exactly).
 (2) a_delta : threshold above which the *conclusion* of P2 written with the auxiliary
     width (dist < 0.50001 alpha^2 q_y + delta + 1.0001 a) can fail:  f(a) - 1.0001 a > delta
     (+ slack that can be made arbitrarily small).  Explicit exact configuration.
     (Theorem 4.11(b) asserts the auxiliary form only for a(Z) <= a*.)
 (3) numerical values used in Corollary 4.4 (sec x - 1 <= 0.50001 x^2 for x <= 1e-3, and the
     constant 1/(2 cos x) of the tested version), Lemma 4.6(c) and Lemma 4.12(b).

Usage:  python thresholds.py          (writes out/thresholds_out.json; a few seconds)
"""
import json, time
from fractions import Fraction as Fr
import mpmath as mp
from tracer import Sq, Cfg, trace, verify_contact, analyse_contact, rot, HALF, mpf, atan2t
from common import out_path

mp.mp.dps = 60
T0 = time.time()
out = {}

f = lambda a: mp.sin(a) + 1 - mp.cos(a)

# (1) a*
astar = mp.findroot(lambda a: f(a) - mp.mpf("1.0001") * a, mp.mpf("2.0001e-4"))
out["a_star"] = mp.nstr(astar, 30)
out["a_star_series_check"] = mp.nstr(mp.mpf("2e-4") + mp.mpf("4e-8") / 3, 20)
# monotonicity of f(a)/a on (0, pi/4]: h(a) = a cos a + a sin a - sin a - 1 + cos a >= 0 ; h' = a(cos a - sin a)
grid = [mp.pi / 4 * j / 20000 for j in range(1, 20001)]
out["min_h_on_grid"] = mp.nstr(min(a * mp.cos(a) + a * mp.sin(a) - mp.sin(a) - 1 + mp.cos(a) for a in grid), 10)
# unique root: sign of f(a)-1.0001a on a grid
sg = [(f(a) - mp.mpf("1.0001") * a) > 0 for a in grid]
out["sign_changes_on_grid"] = sum(1 for i in range(1, len(sg)) if sg[i] != sg[i - 1])
# exact-rational bracketing with t = tan(a/2) rational: a1 < a* < a2
t_star = mp.tan(astar / 2)
den = 10 ** 15
t1 = Fr(int(mp.floor(t_star * den)), den)
t2 = t1 + Fr(1, den)
a1, a2 = atan2t(t1), atan2t(t2)
c1, s1 = rot(t1)
c2, s2 = rot(t2)
out["bracket"] = dict(t1=str(t1), t2=str(t2), a1=mp.nstr(a1, 25), a2=mp.nstr(a2, 25),
                      aux_holds_at_a1=bool(mpf(s1 + 1 - c1) <= mp.mpf("1.0001") * a1),
                      aux_holds_at_a2=bool(mpf(s2 + 1 - c2) <= mp.mpf("1.0001") * a2),
                      margin_a1=mp.nstr(mp.mpf("1.0001") * a1 - mpf(s1 + 1 - c1), 5),
                      margin_a2=mp.nstr(mp.mpf("1.0001") * a2 - mpf(s2 + 1 - c2), 5))
out["ratio_at"] = {str(a): mp.nstr(f(mp.mpf(a)) / mp.mpf(a), 12) for a in
                   ("1e-4", "2e-4", "2.0001e-4", "2.00013e-4", "2.00014e-4", "3e-4", "1e-3", "1e-2")}
out["f_pi4"] = mp.nstr(f(mp.pi / 4), 15)
out["1.0001_pi4"] = mp.nstr(mp.mpf("1.0001") * mp.pi / 4, 15)

# (3) other numerical values
x = mp.mpf("1e-3")
out["sec_minus1_over_x2_at_1e-3"] = mp.nstr((1 / mp.cos(x) - 1) / x ** 2, 12)
out["1_over_2cos_1e-3"] = mp.nstr(1 / (2 * mp.cos(x)), 12)
# largest x with 1/(2 cos x) <= 0.50001, i.e. cos x >= 1/1.00002 (constant of the tested version)
out["C1_exact_threshold_arccos(1/1.00002)"] = mp.nstr(mp.acos(1 / mp.mpf("1.00002")), 12)
out["f(2e-4)/2e-4"] = mp.nstr(f(mp.mpf("2e-4")) / mp.mpf("2e-4"), 12)
out["f(3e-4)/3e-4"] = mp.nstr(f(mp.mpf("3e-4")) / mp.mpf("3e-4"), 12)
# b + b^2/2 <= 1.0001 b  <=> b <= 2e-4 ; and f(b) <= b + b^2/2 check on grid
out["f_le_b_plus_b2over2_grid"] = all(f(b) <= b + b * b / 2 for b in [mp.mpf(j) / 10 ** 7 for j in range(1, 3001)])
# cos t + sin t <= 1 + t on grid
out["cos+sin<=1+t_grid"] = all(mp.cos(t) + mp.sin(t) <= 1 + t for t in [mp.pi / 4 * j / 2000 for j in range(0, 2001)])

# (2) a_delta: root of f(a) - 1.0001 a = delta
ad = {}
for dl in ("1e-3", "1e-5", "1e-7", "1e-9"):
    d = mp.mpf(dl)
    r = mp.findroot(lambda a: f(a) - mp.mpf("1.0001") * a - d, mp.sqrt(2 * d) + mp.mpf("2e-4"))
    ad[dl] = mp.nstr(r, 20)
out["a_delta"] = ad

# explicit exact configuration violating the auxiliary-width conclusion:
# m = 1, X1 axis-parallel on the floor [x1-1/2, x1+1/2] x [0,1]; Z tilted (phi>0) with its
# upper bot vertex at (x, 1+g1), x = x1 - 1/2 + eps, bot(Z) overhanging the left edge of X1.
# The path from (x,0) passes X1 vertically and touches the upper bot vertex of Z (E contact; also
# a relint contact just left of it).  delta = 1e-5, alpha_T = 2 atan(1/2000) (alpha <= 1e-3).
def build(tZ, g1, eps, delta, tau):
    k = Fr(8)
    x1 = Fr(4)
    X1 = Sq(x1, HALF, 0, tag="X1")
    xx = x1 - HALF + eps
    q = (xx, 1 + g1)
    c, s = rot(tZ)
    # Z with bot point at q at xi = +1/2 (upper end for phi > 0): centre = q + n/2 - xi u
    xi = HALF
    Z = Sq(q[0] - s / 2 - xi * c, q[1] + c / 2 - xi * s, tZ, tag="Z")
    cfg = Cfg(k, [X1, Z], "auxviol_t%s" % tZ)
    return cfg, xx


res_cfg = []
delta = Fr(1, 10 ** 5)
tau = Fr(1, 2000)
for tZ in (Fr(1, 400), Fr(1, 200), Fr(1, 100), Fr(1, 40), Fr(1, 10)):
    c, s = rot(tZ)
    g1 = Fr(1, 10 ** 9)
    eps = g1 / (2 * (s / c))          # bot(Z) passes above X1's top-left corner
    cfg, xx = build(tZ, g1, eps, delta, tau)
    ok, msg = cfg.validate()
    r = trace(cfg, xx, delta, tau, cfg.k / 2 - 1)
    rec = [c_ for c_ in r["contacts"] if c_["Z"] == 1]
    entry = dict(tZ=str(tZ), aZ=mp.nstr(atan2t(tZ), 12), valid=ok, msg=msg, term=r["term"])
    if rec:
        rc = rec[0]
        an = analyse_contact(cfg, rc, tau, delta)
        entry.update(kind=rc["kind"], q=[str(rc["q"][0]), str(rc["q"][1])], m=an["m"],
                     g=str(rc["g"]), maxdist=float(an["maxdist"]), argy=str(an["argy"]),
                     w_alphaT=mp.nstr(an["w_alphaT"], 12), w_aux=mp.nstr(an["w_aux"], 12),
                     r_alphaT=mp.nstr(an["r_alphaT"], 12), r_aux=mp.nstr(an["r_aux"], 12),
                     aux_conclusion_violated=bool(an["r_aux"] >= 1),
                     alphaT_holds=bool(an["r_alphaT"] < 1),
                     verify=verify_contact(cfg, rc),
                     Z=cfg.sq[1].to_json(), X1=cfg.sq[0].to_json(), x=str(xx))
        # also a relint contact slightly left of the vertex: x' = xx - 1e-12
        r2 = trace(cfg, xx - Fr(1, 10 ** 12), delta, tau, cfg.k / 2 - 1)
        rec2 = [c_ for c_ in r2["contacts"] if c_["Z"] == 1]
        if rec2:
            an2 = analyse_contact(cfg, rec2[0], tau, delta)
            entry.update(relint_kind=rec2[0]["kind"], relint_r_aux=mp.nstr(an2["r_aux"], 12),
                         relint_r_alphaT=mp.nstr(an2["r_alphaT"], 12))
    res_cfg.append(entry)
out["aux_violation_configs"] = res_cfg

# minimal a(Z) for which this construction violates the aux conclusion at delta=1e-5:
# bisection on t = tan(a/2) rational, with g1 = 1e-12
def viol(tZ):
    c, s = rot(tZ)
    g1 = Fr(1, 10 ** 12)
    eps = g1 / (2 * (s / c))
    cfg, xx = build(tZ, g1, eps, delta, tau)
    if not cfg.validate()[0]:
        return None
    r = trace(cfg, xx, delta, tau, cfg.k / 2 - 1)
    rec = [c_ for c_ in r["contacts"] if c_["Z"] == 1]
    if not rec:
        return None
    an = analyse_contact(cfg, rec[0], tau, delta)
    return an["r_aux"] >= 1, an


lo, hi = Fr(1, 10 ** 4), Fr(1, 100)
assert viol(hi)[0] and not viol(lo)[0]
for _ in range(60):
    mid = (lo + hi) / 2
    mid = mid.limit_denominator(10 ** 14)
    v = viol(mid)
    if v is None:
        break
    if v[0]:
        hi = mid
    else:
        lo = mid
out["aux_conclusion_construction_threshold"] = dict(t_lo=str(lo), t_hi=str(hi), a_lo=mp.nstr(atan2t(lo), 15),
                                                    a_hi=mp.nstr(atan2t(hi), 15),
                                                    a_delta_theory=ad["1e-5"])
out["elapsed"] = time.time() - T0
with open(out_path("thresholds_out.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, default=str)
print(json.dumps(out, indent=1, default=str)[:6000])
