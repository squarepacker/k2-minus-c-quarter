# -*- coding: utf-8 -*-
"""
Exact tracer and checks for Section 4 of the manuscript (P1 and P2).

Written from scratch for these tests; it does not import code from the other
test folders.  It traces the flow F(alpha_T, delta, h) of Section 3 WITHOUT
rule R1 (every path that reaches an entry candidate traverses the square if
a(Y) < alpha_T).  A path of a flow with R1 is an initial part of the path of
the same floor point in the flow without R1 (R1 only terminates paths), so
every first contact of the flow with R1 is also a first contact checked here.

Geometry: closed unit squares with rational centre and rational
t = tan(phi/2), phi in (-pi/4, pi/4]  ->  cos phi, sin phi rational.
u = (cos, sin), n = (-sin, cos).  bot = c - n/2 + xi*u, |xi| <= 1/2.
All positions, gaps, contact points are exact Fractions.
Only alpha (threshold) and a(Z) enter through mpmath (60 digits).

Checks implemented here (manuscript numbering):
  verify_contact  : Theorem 4.2 (height identity (4.3)) recomputed from the
                    recorded gaps, and Lemma 4.1 (F2), (F4) re-checked with an
                    independent segment-clipping routine.
  analyse_contact : Corollary 4.4 (inequality chain (4.5)), Lemma 4.8
                    (inclusions (4.6), (4.7) and the bound on dist(y, Z)),
                    Theorem 4.11(a) (the form with 0.50001 alpha^2 q_y) and the
                    auxiliary form of Theorem 4.11(b).
"""
from fractions import Fraction as Fr
import math
import mpmath as mp

mp.mp.dps = 60
HALF = Fr(1, 2)


def rot(t):
    t = Fr(t)
    den = 1 + t * t
    return (1 - t * t) / den, 2 * t / den


def mpf(q):
    q = Fr(q)
    return mp.mpf(q.numerator) / q.denominator


def atan2t(t):
    """angle 2*atan(|t|) as mpf"""
    return 2 * mp.atan(abs(mpf(t)))


class Sq:
    __slots__ = ("cx", "cy", "t", "c", "s", "tag", "idx", "verts", "bb", "bbf",
                 "sa", "ca", "v")

    def __init__(self, cx, cy, t, tag=""):
        self.cx = Fr(cx)
        self.cy = Fr(cy)
        self.t = Fr(t)
        tt = self.t
        # |phi| < pi/4  <=>  t^2 + 2|t| - 1 < 0  (rational t never hits tan(pi/8))
        if not (tt * tt + 2 * abs(tt) - 1 < 0):
            raise ValueError("tilt out of (-pi/4, pi/4]")
        self.c, self.s = rot(tt)
        self.tag = tag
        self.idx = None
        vs = []
        for xi in (-HALF, HALF):
            for eta in (-HALF, HALF):
                vs.append(self.pt(xi, eta))
        self.verts = vs
        xs = [p[0] for p in vs]
        ys = [p[1] for p in vs]
        self.bb = (min(xs), max(xs), min(ys), max(ys))
        self.bbf = tuple(float(z) for z in self.bb)
        self.sa = abs(self.s)          # sin a(S)
        self.ca = self.c               # cos a(S)
        self.v = self.cy - (self.ca + self.sa) / 2   # lowest height

    def pt(self, xi, eta):
        return (self.cx + xi * self.c - eta * self.s,
                self.cy + xi * self.s + eta * self.c)

    def n(self):
        return (-self.s, self.c)

    def u(self):
        return (self.c, self.s)

    def a_mp(self):
        return atan2t(self.t)

    def ramps_closure(self):
        """closure of ramp set R(S): [v, v+sin a] and [v+cos a, v+cos a+sin a]"""
        v = self.v
        return [(v, v + self.sa), (v + self.ca, v + self.ca + self.sa)]

    def reflect(self, k):
        return Sq(self.cx, Fr(k) - self.cy, -self.t, tag=self.tag + "_refl")

    def to_json(self):
        return [str(self.cx), str(self.cy), str(self.t), self.tag]


# ---------------------------------------------------------------- disjointness
def _proj(verts, ax):
    vals = [p[0] * ax[0] + p[1] * ax[1] for p in verts]
    return min(vals), max(vals)


def disjoint(A, B):
    """closed convex polygons: strictly separated on one of the edge normals"""
    if A.bbf[1] < B.bbf[0] - 1e-9 or B.bbf[1] < A.bbf[0] - 1e-9 or \
       A.bbf[3] < B.bbf[2] - 1e-9 or B.bbf[3] < A.bbf[2] - 1e-9:
        return True
    for ax in (A.u(), A.n(), B.u(), B.n()):
        a0, a1 = _proj(A.verts, ax)
        b0, b1 = _proj(B.verts, ax)
        if a1 < b0 or b1 < a0:
            return True
    return False


def min_sep(A, B):
    """max over the 4 axes of the separation gap (exact); >0 iff disjoint"""
    best = None
    for ax in (A.u(), A.n(), B.u(), B.n()):
        a0, a1 = _proj(A.verts, ax)
        b0, b1 = _proj(B.verts, ax)
        sep = max(b0 - a1, a0 - b1)
        if best is None or sep > best:
            best = sep
    return best


# ---------------------------------------------------------------- ray contact
def first_contact(px, py, dx, dy, S):
    """First contact of ray p + t d (t>=0) with closed square S.
    Returns (t, kind, xi, eta) or None.  kind in
    'bot' (relint bot), 'vbot' (endpoint of bot), 'vtop' (top vertex),
    'side' (relint of a side), 'top' (relint top: should be impossible),
    'inside' (start strictly inside: invalid)."""
    rx = px - S.cx
    ry = py - S.cy
    xi0 = rx * S.c + ry * S.s
    eta0 = -rx * S.s + ry * S.c
    du = dx * S.c + dy * S.s
    dn = -dx * S.s + dy * S.c
    if dn <= 0:
        raise RuntimeError("d.n_S <= 0: direction outside allowed cone")
    lo = (-HALF - eta0) / dn
    hi = (HALF - eta0) / dn
    if du > 0:
        lo2 = (-HALF - xi0) / du
        hi2 = (HALF - xi0) / du
        L = max(lo, lo2)
        H = min(hi, hi2)
    elif du < 0:
        lo2 = (HALF - xi0) / du
        hi2 = (-HALF - xi0) / du
        L = max(lo, lo2)
        H = min(hi, hi2)
    else:
        if abs(xi0) > HALF:
            return None
        L, H = lo, hi
    if L > H or H < 0:
        return None
    if L < 0:
        if H > 0:
            return (Fr(0), "inside", xi0, eta0)
        tc = Fr(0)
    else:
        tc = L
    xi = xi0 + tc * du
    eta = eta0 + tc * dn
    onb = eta == -HALF
    ont = eta == HALF
    ons = abs(xi) == HALF
    if onb and not ons:
        kind = "bot"
    elif onb and ons:
        kind = "vbot"
    elif ont and ons:
        kind = "vtop"
    elif ons:
        kind = "side"
    elif ont:
        kind = "top"
    else:
        kind = "inside"
    return (tc, kind, xi, eta)


# ---------------------------------------------------------------- configuration
class Cfg:
    def __init__(self, k, squares, name=""):
        self.k = Fr(k)
        self.sq = list(squares)
        self.name = name
        for i, S in enumerate(self.sq):
            S.idx = i
        self.grid = {}
        for S in self.sq:
            x0, x1, y0, y1 = S.bbf
            for gx in range(math.floor(x0), math.floor(x1) + 1):
                for gy in range(math.floor(y0), math.floor(y1) + 1):
                    self.grid.setdefault((gx, gy), []).append(S)

    def validate(self):
        k = self.k
        for S in self.sq:
            for p in S.verts:
                if not (0 <= p[0] <= k and 0 <= p[1] <= k):
                    return False, "container: square %d (%s)" % (S.idx, S.tag)
        n = len(self.sq)
        for i in range(n):
            A = self.sq[i]
            for j in range(i + 1, n):
                if not disjoint(A, self.sq[j]):
                    return False, "overlap %d %d (%s,%s)" % (i, j, A.tag, self.sq[j].tag)
        return True, "ok"

    def candidates(self, px, py, dx, dy, tmax):
        fx, fy = float(px), float(py)
        ex = fx + float(tmax) * float(dx)
        ey = fy + float(tmax) * float(dy)
        m = 1e-7
        x0, x1 = min(fx, ex) - m, max(fx, ex) + m
        y0, y1 = min(fy, ey) - m, max(fy, ey) + m
        seen = set()
        out = []
        for gx in range(math.floor(x0), math.floor(x1) + 1):
            for gy in range(math.floor(y0), math.floor(y1) + 1):
                for S in self.grid.get((gx, gy), ()):
                    if S.idx in seen:
                        continue
                    seen.add(S.idx)
                    b = S.bbf
                    if b[1] < x0 or b[0] > x1 or b[3] < y0 or b[2] > y1:
                        continue
                    out.append(S)
        return out

    def reflected(self):
        return Cfg(self.k, [S.reflect(self.k) for S in self.sq], self.name + "_refl")


# ---------------------------------------------------------------- tracer
def trace(cfg, x, delta, tauT, h, maxsteps=100000):
    """R1-free flow F(alpha_T, delta, h) with alpha_T = 2 atan(tauT).
    Priority D > W > contact > H; entry only at relint bot; T if a(Y) >= alpha_T
    (<=> |t_Y| >= tauT, exact).  Returns dict(term=..., contacts=[...]).
    Every first contact (also D-/W-ties) is recorded with the full history."""
    k = cfg.k
    delta = Fr(delta)
    tauT = Fr(tauT)
    h = Fr(h)
    px, py = Fr(x), Fr(0)
    dx, dy = Fr(0), Fr(1)
    X = None
    g = Fr(0)
    passed = []
    gaps = []
    starts = [(px, py)]
    contacts = []
    term = None
    steps = 0
    while True:
        steps += 1
        if steps > maxsteps:
            term = "MAXSTEPS"
            break
        tD = delta - g
        tH = (h - py) / dy
        if dx > 0:
            tW = (k - px) / dx
        elif dx < 0:
            tW = -px / dx
        else:
            tW = None
        tlim = min(tD, tH) if tW is None else min(tD, tH, tW)
        if tlim < 0:
            tlim = Fr(0)
        best = None
        for S in cfg.candidates(px, py, dx, dy, tlim):
            if X is not None and S.idx == X:
                continue
            r = first_contact(px, py, dx, dy, S)
            if r is None:
                continue
            t, kind, xi, eta = r
            if kind == "inside":
                raise RuntimeError("start inside square %d" % S.idx)
            if t > tlim:
                continue
            if best is None or t < best[0]:
                best = (t, kind, S, xi, eta)
            elif t == best[0]:
                raise RuntimeError("two squares at same contact parameter")
        tS = best[0] if best is not None else None
        cands = [tD, tH]
        if tW is not None:
            cands.append(tW)
        if tS is not None:
            cands.append(tS)
        tstar = min(cands)
        if tS is not None and tS == tstar:
            t, kind, S, xi, eta = best
            q = (px + t * dx, py + t * dy)
            flag = "D" if tD == tstar else ("W" if (tW is not None and tW == tstar) else "")
            contacts.append(dict(x=Fr(x), passed=list(passed), gaps=list(gaps) + [t],
                                 q=q, Z=S.idx, kind=kind, xi=xi, eta=eta,
                                 g=g + t, flag=flag, p=(px, py), d=(dx, dy), X=X,
                                 starts=list(starts)))
        if tD == tstar:
            term = "D"
            break
        if tW is not None and tW == tstar:
            term = "W"
            break
        if tS is not None and tS == tstar:
            t, kind, S, xi, eta = best
            g = g + t
            gaps.append(t)
            q = (px + t * dx, py + t * dy)
            if kind == "bot":
                if abs(S.t) >= tauT:
                    term = "T"
                    break
                passed.append(S.idx)
                nn = S.n()
                px, py = q[0] + nn[0], q[1] + nn[1]
                dx, dy = nn
                X = S.idx
                starts.append((px, py))
                if py >= h:
                    term = "Hexit"
                    break
                continue
            term = "E"
            break
        term = "H"
        break
    return dict(term=term, contacts=contacts, passed=passed, g=g, x=Fr(x))


# ---------------------------------------------------------------- independent verifier
def _edges(S):
    # CCW order of vertices: (xi,eta) = (-,-), (+,-), (+,+), (-,+)
    P = [S.pt(-HALF, -HALF), S.pt(HALF, -HALF), S.pt(HALF, HALF), S.pt(-HALF, HALF)]
    return [(P[i], P[(i + 1) % 4]) for i in range(4)]


def clip_segment(P, Q, S):
    """Cyrus-Beck in world frame from vertices: returns [t0,t1] subset [0,1] with
    P + t(Q-P) in closed S, or None.  Independent of first_contact()."""
    t0, t1 = Fr(0), Fr(1)
    Dx, Dy = Q[0] - P[0], Q[1] - P[1]
    for (A, B) in _edges(S):
        ex, ey = B[0] - A[0], B[1] - A[1]
        # inside: cross(e, X - A) >= 0
        f0 = ex * (P[1] - A[1]) - ey * (P[0] - A[0])
        f1 = ex * Dy - ey * Dx          # derivative in t
        if f1 == 0:
            if f0 < 0:
                return None
            continue
        tc = -f0 / f1
        if f1 > 0:
            if tc > t0:
                t0 = tc
        else:
            if tc < t1:
                t1 = tc
        if t0 > t1:
            return None
    return (t0, t1)


def verify_contact(cfg, rec):
    """Independent re-verification of one recorded contact:
    (1) P1 identity (4.3) from recorded gaps and passed squares,
    (2) q on the boundary of Z at recorded (xi,eta), and classification by world-frame test,
    (3) every gap segment touches no square except its own endpoints' squares,
    (4) every passage segment lies in its square and enters through relint bot.
    Returns list of failure strings (empty = ok)."""
    fails = []
    x = rec["x"]
    ns = [(Fr(0), Fr(1))] + [cfg.sq[i].n() for i in rec["passed"]]
    qx, qy = Fr(x), Fr(0)
    for gi, nv in zip(rec["gaps"], ns):
        qx += gi * nv[0]
        qy += gi * nv[1]
    for nv in ns[1:]:
        qx += nv[0]
        qy += nv[1]
    if (qx, qy) != rec["q"]:
        fails.append("P1 identity")
    if sum(rec["gaps"], Fr(0)) != rec["g"]:
        fails.append("P1 g-state")
    if len(rec["gaps"]) != len(rec["passed"]) + 1:
        fails.append("gap count")
    Z = cfg.sq[rec["Z"]]
    if Z.pt(rec["xi"], rec["eta"]) != rec["q"]:
        fails.append("q not at recorded frame point")
    # rebuild polyline
    pts = [(Fr(x), Fr(0))]
    cur = (Fr(x), Fr(0))
    seg_owner = [None]
    for i, sidx in enumerate(rec["passed"]):
        gi = rec["gaps"][i]
        dv = ns[i]
        qq = (cur[0] + gi * dv[0], cur[1] + gi * dv[1])
        S = cfg.sq[sidx]
        # entry must be relint bot of S: world test via bot endpoints
        A = S.pt(-HALF, -HALF)
        B = S.pt(HALF, -HALF)
        cr = (B[0] - A[0]) * (qq[1] - A[1]) - (B[1] - A[1]) * (qq[0] - A[0])
        lam = ((qq[0] - A[0]) * (B[0] - A[0]) + (qq[1] - A[1]) * (B[1] - A[1]))
        if cr != 0 or not (0 < lam < 1):
            fails.append("entry %d not relint bot" % i)
        # gap segment check
        fails += _gap_free(cfg, cur, qq, prev=(rec["passed"][i - 1] if i > 0 else None), nxt=sidx)
        ex = (qq[0] + S.n()[0], qq[1] + S.n()[1])
        # exit in top relint
        C = S.pt(-HALF, HALF)
        D = S.pt(HALF, HALF)
        cr2 = (D[0] - C[0]) * (ex[1] - C[1]) - (D[1] - C[1]) * (ex[0] - C[0])
        lam2 = ((ex[0] - C[0]) * (D[0] - C[0]) + (ex[1] - C[1]) * (D[1] - C[1]))
        if cr2 != 0 or not (0 < lam2 < 1):
            fails.append("exit %d not relint top" % i)
        cur = ex
    # final gap segment
    gi = rec["gaps"][-1]
    dv = ns[-1]
    qq = (cur[0] + gi * dv[0], cur[1] + gi * dv[1])
    if qq != rec["q"]:
        fails.append("final point mismatch")
    fails += _gap_free(cfg, cur, qq, prev=(rec["passed"][-1] if rec["passed"] else None),
                       nxt=rec["Z"])
    return fails


def _gap_free(cfg, P, Q, prev, nxt):
    fails = []
    if P == Q:
        # zero-length gap (floor-touching square): only allowed at floor start
        if prev is not None:
            fails.append("zero gap after square")
        return fails
    tl = Q
    for S in cfg.candidates(P[0], P[1], Q[0] - P[0], Q[1] - P[1], Fr(1)):
        r = clip_segment(P, Q, S)
        if r is None:
            continue
        if S.idx == nxt:
            if r != (Fr(1), Fr(1)):
                fails.append("gap hits target %d before q: %s" % (S.idx, (float(r[0]), float(r[1]))))
        elif prev is not None and S.idx == prev:
            if r != (Fr(0), Fr(0)):
                fails.append("gap re-enters previous %d" % S.idx)
        else:
            fails.append("gap meets square %d at %s" % (S.idx, (float(r[0]), float(r[1]))))
    return fails


# ---------------------------------------------------------------- P1 / P2 quantities
def dist_int(y):
    y = Fr(y)
    f = y - math.floor(y)
    return min(f, 1 - f)


def maxdist_interval(lo, hi):
    """max over y in [lo,hi] of dist(y, Z); exact"""
    lo, hi = Fr(lo), Fr(hi)
    # half-integer inside?
    hmin = math.ceil(lo - HALF)
    if hmin + HALF <= hi:
        return HALF, hmin + HALF
    d1, d2 = dist_int(lo), dist_int(hi)
    return (d1, lo) if d1 >= d2 else (d2, hi)


def alpha_of(tau):
    return atan2t(tau)


def Kcoef(alpha):
    """alpha^2/(2 cos alpha) as mpf.

    This is the constant of the version of Corollary 4.4 / Lemma 4.8 that was
    tested.  The manuscript states these results with the coefficient
    sec(alpha) - 1 = (1 - cos alpha)/cos alpha, which is not larger
    (sec(alpha) - 1 = alpha^2/(2 cos alpha) * (1 - alpha^2/12 + O(alpha^4)));
    so the inequalities checked with this coefficient are slightly weaker
    than the manuscript's forms."""
    return alpha * alpha / (2 * mp.cos(alpha))


def analyse_contact(cfg, rec, tauT, delta, c0sq_over=None):
    """P1 corollary and P2 quantities for one contact.  Returns dict of mp floats."""
    sq = cfg.sq
    m = len(rec["passed"])
    q = rec["q"]
    qy = q[1]
    g = rec["g"]
    phis_t = [sq[i].t for i in rec["passed"]]
    cosv = [sq[i].c for i in rec["passed"]]
    D = sum((1 - c for c in cosv), Fr(0))            # tilt loss, exact
    G = sum((gi * c for gi, c in zip(rec["gaps"], [Fr(1)] + cosv)), Fr(0))   # vertical gap
    out = dict(m=m, qy=qy, g=g, D=D, G=G)
    # exact identity check of q_y - m = -D + G
    out["id_qy"] = (qy - m == -D + G)
    alT = alpha_of(tauT)
    out["alphaT"] = alT
    qy_mp = mpf(qy)
    # alpha_eff = max tilt among passed (the smallest alpha for which the corollary hypothesis
    # holds in the closed form); 0 if none
    if phis_t:
        tmax = max(abs(t) for t in phis_t)
        aeff = atan2t(tmax)
    else:
        aeff = mp.mpf(0)
    out["aeff"] = aeff
    KT = Kcoef(alT) * qy_mp
    Keff = Kcoef(aeff) * qy_mp if aeff > 0 else mp.mpf(0)
    out["KT"] = KT
    out["Keff"] = Keff
    lhs = mpf(qy - m)
    # corollary: -K <= q_y - m <= g
    out["cor_low_T"] = (lhs >= -KT)
    out["cor_low_eff"] = (D == 0) or (mpf(D) <= Keff * (1 + mp.mpf(10) ** -50))
    out["cor_up"] = (qy - m <= g)
    out["ratio_tilt_T"] = (mpf(D) / KT) if KT > 0 else mp.mpf(0)
    out["ratio_tilt_eff"] = (mpf(D) / Keff) if Keff > 0 else mp.mpf(0)
    out["ratio_gap"] = (mpf(qy - m) / mpf(g)) if g > 0 else mp.mpf(0)
    out["ratio_gap_delta"] = mpf(qy - m) / mpf(delta)
    # numerical form of Corollary 4.4 (valid for alpha <= 1e-3): q_y - m in [-0.50001 alpha^2 q_y, delta)
    if alT <= mp.mpf("1e-3"):
        out["cor_050001"] = (lhs >= -mp.mpf("0.50001") * alT ** 2 * qy_mp) and (qy - m < delta)
    else:
        out["cor_050001"] = None
    # P2 only for contacts on closed bot(Z)
    out["onbot"] = rec["kind"] in ("bot", "vbot") and rec["eta"] == -HALF
    if out["onbot"]:
        Z = sq[rec["Z"]]
        sa, ca = Z.sa, Z.ca
        aZ = Z.a_mp()
        out["aZ"] = aZ
        R = Z.ramps_closure()
        md = HALF * 0
        argy = None
        for (lo, hi) in R:
            dd, yy = maxdist_interval(lo, hi)
            if dd > md:
                md, argy = dd, yy
        if Z.t == 0:
            md = Fr(0)   # R(Z) empty for a = 0
        out["maxdist"] = md
        out["argy"] = argy
        f_a = mpf(sa + (1 - ca))
        # widths compared with max dist(y, Z) over the closure of R(Z) (Lemma 4.8):
        #   w_alphaT = K_T + delta + f(a)   (K with the flow threshold alpha_T; strict bound)
        #   w_g      = K_T + g + f(a)       (cumulative gap g <= delta in place of delta)
        #   w_eff    = K_eff + g + f(a)     (alpha_eff = largest inclination actually traversed)
        wP = KT + mpf(delta) + f_a
        wg = KT + mpf(g) + f_a
        weff = Keff + mpf(g) + f_a
        out["w_alphaT"] = wP
        out["w_g"] = wg
        out["w_eff"] = weff
        def sdiv(a_, b_):
            if b_ == 0:
                return mp.mpf(0) if a_ == 0 else mp.inf
            return a_ / b_
        out["r_alphaT"] = sdiv(mpf(md), wP)
        out["r_g"] = sdiv(mpf(md), wg)
        out["r_eff"] = sdiv(mpf(md), weff)
        # Theorem 4.11(a) form (alpha <= 1e-3): w = 0.50001 alpha^2 q_y + delta + f(a(Z))
        if alT <= mp.mpf("1e-3"):
            wD = mp.mpf("0.50001") * alT ** 2 * qy_mp + mpf(delta) + f_a
            out["r_050001"] = mpf(md) / wD
        else:
            out["r_050001"] = None
        # auxiliary form of Theorem 4.11(b): 0.50001 alpha^2 q_y + delta + 1.0001 a(Z)
        # (the manuscript asserts it only for a(Z) <= a*, Lemma 4.6; it is recorded, not required)
        waux = mp.mpf("0.50001") * alT ** 2 * qy_mp + mpf(delta) + mp.mpf("1.0001") * aZ
        out["w_aux"] = waux
        out["r_aux"] = mpf(md) / waux
        out["aux_ineq_holds"] = bool(f_a <= mp.mpf("1.0001") * aZ)
        # interval inclusions (4.6), (4.7) of Lemma 4.8 with K = Keff and g (strongest form)
        v = Z.v
        Kf = Keff
        low_ok = (mpf(v) >= m - Kf - mpf(sa) - mp.mpf(10) ** -50) and (v + sa <= m + g + sa)
        up_lo = mpf(v + ca) >= (m + 1 - Kf - mpf(sa) - mpf(1 - ca) - mp.mpf(10) ** -50)
        up_hi = (v + ca + sa) <= (m + 1 + g + sa - (1 - ca))
        out["contain_ok"] = bool(low_ok and up_lo and up_hi)
    return out
