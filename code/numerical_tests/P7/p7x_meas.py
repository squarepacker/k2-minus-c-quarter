"""p7x_meas.py -- quantities of Section 3.8 and Section 9 (P7) of the manuscript for the scaled analogue, computed
exactly: the scale flows F_s (Definition 3.21), T* (Definition 3.22), Z(s), N(s), Lambda_0, the wall term, the line
quantities omega, E, Y_b, W_s, V_s, and the width test of Theorem 4.11 (P2) with w(s, Z) of Lemma 9.3(e).

Scaling (see README.md):
  alpha(s) = c0 s^-1/2, alpha_max = alpha(y0), beta^(t) = c1 t^-3/4, beta(y) = beta^(d(y)),
  b := beta^((1-eps) y0)  (replaces 1.0001-1 = 1e-4: max tilt of a square of Z(s)),  h := 1 + b  (replaces 1.0001)
  P2 width (exact form):  w(s,Z) = (sec alpha(s) - 1) q_y + delta + sin a + 1 - cos a      (replaces 0.50001 alpha^2 q_y ...)
  w0 := sup_{s in [y0,y1]} (sec alpha(s) - 1)(s + h) + delta + sin b + 1 - cos b + 1e-12   (rounded up to a rational)
  V_s = [(1-eps)s - h, s + h],  |V_s| = eps s + 2h   (b_W := 2h replaces 2.0002)
"""
from fractions import Fraction as Fr
import math
import mpmath
from p7x_core import mpf_fr, HALF

mpmath.mp.dps = 60


def rat_up(x, den=10 ** 30):
    q = Fr(int(mpmath.ceil(x * den)), den)
    assert mpf_fr(q) >= x
    return q


def rat_near(x, den=10 ** 40):
    return Fr(int(mpmath.nint(x * den)), den)


class Params:
    def __init__(self, k, delta, c0, c1, y0, eps, omega0):
        self.k = Fr(k); self.delta = Fr(delta); self.c0 = Fr(c0); self.c1 = Fr(c1)
        self.y0 = Fr(y0); self.eps = Fr(eps); self.omega0 = Fr(omega0)
        self.y1 = self.k / 2 - 3
        self.hmax = self.k / 2 - 1
        self.alpha_max = self.alpha(self.y0)
        self.b = mpf_fr(self.c1) * mpf_fr((1 - self.eps) * self.y0) ** mpmath.mpf(-0.75)
        self.h = 1 + self.b
        # sup_s (sec alpha(s) - 1)(s+h): check decreasing on a grid, then take value at y0
        def f(s):
            return (mpmath.sec(mpf_fr(self.c0) / mpmath.sqrt(s)) - 1) * (s + self.h)
        y0m, y1m = mpf_fr(self.y0), mpf_fr(self.y1)
        grid = [y0m + (y1m - y0m) * i / 400 for i in range(401)]
        vals = [f(s) for s in grid]
        self.sup_drift = max(vals)
        self.decreasing = all(vals[i] >= vals[i + 1] for i in range(len(vals) - 1))
        w0 = self.sup_drift + mpf_fr(self.delta) + mpmath.sin(self.b) + 1 - mpmath.cos(self.b) + mpmath.mpf('1e-12')
        self.w0 = rat_up(w0)
        self.h_rat_up = rat_up(self.h)

    def alpha(self, s):
        return mpf_fr(self.c0) / mpmath.sqrt(mpf_fr(Fr(s)))

    def beta_hat(self, t):
        return mpf_fr(self.c1) * mpf_fr(Fr(t)) ** mpmath.mpf(-0.75)

    def Vs(self, s):
        s = Fr(s)
        return ((1 - self.eps) * s - self.h, s + self.h)      # mpf endpoints

    def Vlen(self, s):
        return mpf_fr(self.eps * Fr(s)) + 2 * self.h

    def wall(self, s):
        return 2 * (mpf_fr(Fr(s)) + 2) * mpmath.tan(self.alpha(s))

    def desc(self):
        return dict(k=float(self.k), delta=float(self.delta), c0=float(self.c0), c1=float(self.c1),
                    y0=float(self.y0), eps=float(self.eps), omega0=float(self.omega0), y1=float(self.y1),
                    alpha_max=float(self.alpha_max), b=float(self.b), w0=float(self.w0),
                    sup_drift=float(self.sup_drift), drift_decreasing=self.decreasing)


# ------------------------------------------------------------------ intervals of x
def merge(iv):
    iv = sorted((a, b) for (a, b) in iv if a < b)
    out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return out


def measure(iv):
    return sum((b - a) for a, b in merge(iv))


def symdiff_measure(A, B):
    A = merge(A); B = merge(B)
    pts = sorted(set([p for ab in A + B for p in ab]))
    tot = Fr(0)
    def inside(L, x):
        return any(a < x < b for a, b in L)
    for i in range(len(pts) - 1):
        m = (pts[i] + pts[i + 1]) / 2
        if inside(A, m) != inside(B, m):
            tot += pts[i + 1] - pts[i]
    return tot


def split_affine(xlo, xhi, funcs, val):
    """split (xlo,xhi) at points where any affine f(x) = val."""
    pts = {xlo, xhi}
    for f in funcs:
        if f[1] != 0:
            r = (val - f[0]) / f[1]
            if xlo < r < xhi:
                pts.add(r)
    pts = sorted(pts)
    return list(zip(pts[:-1], pts[1:]))


# ------------------------------------------------------------------ F_s pieces (direct truncation walk)
class FsPiece:
    __slots__ = ('xlo', 'xhi', 'typ', 'way', 'seg', 'dirs', 'tau', 'tsq')

    def __init__(self, xlo, xhi, typ, way, seg, dirs, tau, tsq=None):
        self.xlo = xlo; self.xhi = xhi; self.typ = typ; self.way = way; self.seg = seg
        self.dirs = dirs; self.tau = tau; self.tsq = tsq


def trunc_map(flow, P, s):
    tan_as = mpmath.tan(P.alpha(s))
    tm = {}; near = []
    for sid, S in flow.sq.items():
        d = mpf_fr(S.tan_a) - tan_as
        if abs(d) < mpmath.mpf(10) ** -40:
            near.append(sid)
        tm[sid] = d >= 0
    return tm, near


def fs_pieces(flow, P, s, tm):
    """F_s = F^0 truncated at the first of: T_s (R1-won entry with a >= alpha(s)), H_s (height s+2)."""
    H = Fr(s) + 2
    out = []
    for r in flow.recs:
        entset = {idx: sid for (sid, idx) in r.ent}
        last = len(r.way) - 1
        stack = [(r.xlo, r.xhi, 0)]
        while stack:
            xlo, xhi, i = stack.pop()
            w0, w1 = r.way[i], r.way[i + 1]
            he = (w1[2], w1[3])
            for (a, b) in split_affine(xlo, xhi, [he], H):
                m = (a + b) / 2
                hm = he[0] + he[1] * m
                if r.seg[i] == -1:
                    if hm > H:          # H_s inside this gap
                        d = r.dirs[i]
                        # point = w0 + ((H - w0y)/dy) d
                        t0 = (H - w0[2]) / d[1]; t1 = -w0[3] / d[1]
                        pt = (w0[0] + d[0] * t0, w0[1] + d[0] * t1, w0[2] + d[1] * t0, w0[3] + d[1] * t1)
                        out.append(FsPiece(a, b, 'H', r.way[:i + 1] + (pt,), r.seg[:i + 1], r.dirs[:i + 1], (H, Fr(0))))
                        continue
                    j = i + 1
                    if j in entset:
                        sid = entset[j]
                        if tm[sid]:
                            out.append(FsPiece(a, b, 'T', r.way[:j + 1], r.seg[:j], r.dirs[:j], (w1[2], w1[3]), sid))
                            continue
                        if j == last:     # T_max that does not truncate F_s -- impossible
                            out.append(FsPiece(a, b, 'BAD_T', r.way, r.seg, r.dirs, (w1[2], w1[3]), sid))
                            continue
                        stack.append((a, b, j))
                        continue
                    if j == last:
                        typ = r.typ
                        if typ == 'H':
                            out.append(FsPiece(a, b, 'H', r.way, r.seg, r.dirs, (H, Fr(0))))
                        else:
                            out.append(FsPiece(a, b, typ, r.way, r.seg, r.dirs, (w1[2], w1[3]), r.tsq))
                        continue
                    stack.append((a, b, j))
                else:                    # passage
                    if hm >= H:
                        out.append(FsPiece(a, b, 'H', r.way[:i + 2], r.seg[:i + 1], r.dirs[:i + 1], (H, Fr(0))))
                        continue
                    j = i + 1
                    if j == last:
                        out.append(FsPiece(a, b, 'BAD_end', r.way, r.seg, r.dirs, (w1[2], w1[3])))
                        continue
                    stack.append((a, b, j))
    return out


def recs_as_pieces(flow):
    """F^0 records in the FsPiece format (tau = termination height, H -> h_max)."""
    out = []
    for r in flow.recs:
        w = r.way[-1]
        tau = (flow.hF, Fr(0)) if r.typ == 'H' else (w[2], w[3])
        out.append(FsPiece(r.xlo, r.xhi, r.typ, r.way, r.seg, r.dirs, tau, r.tsq))
    return out


# ------------------------------------------------------------------ T* and the sets of Lemma 9.6
def tstar_sets(flow, P, s, tm):
    """returns (the set {T* <= s} by the formula of Definition 3.22, the set of x having an R1-passed entry
    with a(Z_j) >= alpha(s) and eta_j <= s+2, as in Lemma 9.6(b),(c)) as interval lists."""
    s = Fr(s); H = s + 2
    A = []; B = []
    for r in flow.recs:
        for (sid, idx) in r.ent:
            S = flow.sq[sid]
            eta = (r.way[idx][2], r.way[idx][3])
            # T*: s_j = max(y0, eta-2, (c0/a)^2) <= s
            ok_tilt = S.a > 0 and (mpf_fr(P.c0) / S.a) ** 2 <= mpf_fr(s)
            for (a, b) in split_affine(r.xlo, r.xhi, [eta], H):
                m = (a + b) / 2
                if eta[0] + eta[1] * m <= H:
                    if ok_tilt and P.y0 <= s:
                        A.append((a, b))
                    if tm[sid]:
                        B.append((a, b))
    return merge(A), merge(B)


# ------------------------------------------------------------------ per-height bookkeeping
def height_stats(pieces, y, flow, squares, want_gap_check=True):
    y = Fr(y)
    k = flow.k
    loss = {t: Fr(0) for t in ('T', 'D', 'E', 'M', 'W', 'H')}
    mu = {}
    Fgap = Fr(0)
    gapints = []
    degen = 0
    for p in pieces:
        funcs = [p.tau] + [(w[2], w[3]) for w in p.way]
        for (a, b) in split_affine(p.xlo, p.xhi, funcs, y):
            m = (a + b) / 2
            L = b - a
            tv = p.tau[0] + p.tau[1] * m
            if tv < y:
                loss[p.typ if p.typ in loss else 'H'] += L
                continue
            if tv == y:
                degen += 1
            found = False
            for i in range(len(p.way) - 1):
                h0 = p.way[i][2] + p.way[i][3] * m
                h1 = p.way[i + 1][2] + p.way[i + 1][3] * m
                if h0 == y or h1 == y:
                    degen += 1
                if h0 < y < h1:
                    found = True
                    if p.seg[i] == -1:
                        d = p.dirs[i]
                        w = p.way[i]
                        r = d[0] / d[1]
                        X0 = w[0] + (y - w[2]) * r
                        X1 = w[1] - w[3] * r
                        xa, xb = X0 + X1 * a, X0 + X1 * b
                        if xa > xb:
                            xa, xb = xb, xa
                        Fgap += L
                        gapints.append((xa, xb, L))
                    else:
                        mu[p.seg[i]] = mu.get(p.seg[i], Fr(0)) + L
                    break
            if not found:
                degen += 1
    chords = {S.id: S.chord(y) for S in squares if S.ymin < y < S.ymax}
    omega = k - sum(chords.values())
    Sh = Fr(0); Ov = Fr(0)
    for S in squares:
        c = chords.get(S.id, Fr(0)); mz = mu.get(S.id, Fr(0))
        if c > mz:
            Sh += c - mz
        elif mz > c:
            Ov += mz - c
    # Ov_gap: integral of (sum of densities - 1)_+
    ev = []
    for (xa, xb, L) in gapints:
        if xb > xa:
            rho = L / (xb - xa)
            ev.append((xa, rho)); ev.append((xb, -rho))
        else:
            degen += 1
    ev.sort(key=lambda z: z[0])
    Ovg = Fr(0); cur = Fr(0); last = None
    for (xp, dr) in ev:
        if last is not None and cur > 1:
            Ovg += (cur - 1) * (xp - last)
        cur += dr; last = xp
    Ltot = sum(loss.values())
    musum = sum(mu.values())
    res = dict(loss=loss, mu=mu, Fgap=Fgap, omega=omega, Sh=Sh, Ov=Ov, Ovgap=Ovg, degen=degen,
               cons=k - (Ltot + musum + Fgap), ident=Sh - (Ltot + Fgap - omega + Ov))
    if want_gap_check:
        bad = 0
        act = [S for S in squares if S.ymin < y < S.ymax]
        for (xa, xb, L) in gapints[:400]:
            xm = (xa + xb) / 2
            if any(S.contains_open(xm, y) for S in act):
                bad += 1
        res['gap_in_square'] = bad
    return res


# ------------------------------------------------------------------ lines, H, kinds, Z(s)
class Lines:
    def __init__(self, squares, P):
        self.sq = squares; self.P = P
        self.k = P.k
        hs = set()
        for S in squares:
            for v in S.V:
                hs.add(v[1])
        self.vh = sorted(hs)
        # Phi / c=1 breakpoints (rational)
        eb = set(self.vh)
        for S in squares:
            if S.s != 0:
                sca = abs(S.s) * S.c
                eb.add(S.ymin + sca); eb.add(S.ymax - sca)
        self.eb = sorted(eb)
        self.ystar = {}
        for S in squares:
            if S.s != 0:
                self.ystar[S.id] = (mpf_fr(P.c1) / S.a) ** (mpmath.mpf(4) / 3)
        self.byy = sorted(squares, key=lambda S: S.ymin)
        self.ymins = [S.ymin for S in self.byy]
        self._ch = {}
        self._kind = {}

    def active(self, y):
        import bisect
        lo = bisect.bisect_left(self.ymins, y - Fr(3, 2))
        hi = bisect.bisect_right(self.ymins, y)
        return [S for S in self.byy[lo:hi] if S.ymax >= y]

    def chords(self, y):
        if y not in self._ch:
            self._ch[y] = [(S, S.chord(y)) for S in self.active(y)]
        return self._ch[y]

    def omega(self, y):
        return self.k - sum((c for (S, c) in self.chords(y)), Fr(0))

    def E(self, y):
        tot = Fr(0)
        for (S, c) in self.chords(y):
            if c > 1:
                tot += c - 1
        return tot

    def short_sum(self, y):
        tot = Fr(0)
        for (S, c) in self.chords(y):
            if 0 < c < 1:
                tot += c
        return tot

    def in_H(self, y):
        P = self.P
        d = y if y <= self.k / 2 else self.k - y
        if not (P.y0 <= d <= (1 - P.eps) * P.y1):
            return False
        f = y - math.floor(y)
        return P.w0 <= f <= 1 - P.w0

    def kind(self, y):
        """kind of line y in H: 'i','ii','iii' (exact; tilt comparisons in 60-digit arithmetic)."""
        if y in self._kind:
            return self._kind[y]
        P = self.P
        if self.omega(y) >= P.omega0:
            self._kind[y] = 'i'
            return 'i'
        d = y if y <= self.k / 2 else self.k - y
        be = mpf_fr(P.c1) * mpf_fr(d) ** mpmath.mpf(-0.75)
        res = 'iii'
        for (S, c) in self.chords(y):
            if S.a >= be:
                res = 'ii'; break
        self._kind[y] = res
        return res

    def in_Yb(self, y):
        return 0 < y < self.k / 2 and self.in_H(y) and self.kind(y) == 'iii'

    def in_Yt(self, y):
        return self.k / 2 < y < self.k and self.in_H(y) and self.kind(y) == 'iii'

    def Zset_ceiling(self, s):
        """tilde Z(s) = {S : Phi(S) cap Y_t cap (k - W_s) nonempty}, by exact tests on a fine candidate set."""
        P = self.P
        lo, hi = self.k - Fr(s), self.k - (1 - P.eps) * Fr(s)
        bp = sorted(set(self.k - b for b in self.breakpoints((1 - P.eps) * Fr(s), Fr(s))) |
                    {S.ymin for S in self.sq} | {S.ymax for S in self.sq} | {lo, hi})
        out = {}
        for S in self.sq:
            if S.s == 0:
                continue
            sca = abs(S.s) * S.c
            for (p, q) in ((S.ymin, S.ymin + sca), (S.ymax - sca, S.ymax)):
                inner = sorted(z for z in bp if p < z < q)
                pl = [p] + inner + [q]
                cand = inner + [(pl[i] + pl[i + 1]) / 2 for i in range(len(pl) - 1)]
                if any(lo <= z <= hi and self.in_Yt(z) for z in cand):
                    out[S.id] = True; break
        return set(out)

    def breakpoints(self, lo, hi):
        """candidate breakpoints of Y_b within [lo,hi] (rational)."""
        P = self.P
        pts = {lo, hi, P.y0, (1 - P.eps) * P.y1, self.k / 2}
        for m in range(int(math.floor(lo)) - 1, int(math.ceil(hi)) + 2):
            pts.add(m + P.w0); pts.add(m + 1 - P.w0)
        for S in self.sq:
            pts.add(S.ymin); pts.add(S.ymax)
        for v in self.ystar.values():
            if v < 10 ** 6:
                pts.add(rat_near(v))
        # omega = omega0 crossings
        vh = [v for v in self.vh if lo - 2 <= v <= hi + 2]
        for i in range(len(vh) - 1):
            a, b = vh[i], vh[i + 1]
            y1 = a + (b - a) / 3; y2 = a + 2 * (b - a) / 3
            o1, o2 = self.omega(y1), self.omega(y2)
            if o1 != o2:
                sl = (o2 - o1) / (y2 - y1)
                r = y1 + (P.omega0 - o1) / sl
                if a < r < b:
                    pts.add(r)
        return sorted(p for p in pts if lo <= p <= hi)

    def Yb_W(self, s):
        """Y_b cap W_s as list of (lo,hi) closed-ish intervals plus isolated points (exact membership tests)."""
        P = self.P
        lo, hi = (1 - P.eps) * Fr(s), Fr(s)
        bp = self.breakpoints(lo, hi)
        ivs = []; pts = []
        memo = {}
        def mem(y):
            if y not in memo:
                memo[y] = self.in_Yb(y)
            return memo[y]
        for i in range(len(bp) - 1):
            a, b = bp[i], bp[i + 1]
            if mem((a + b) / 2):
                if ivs and ivs[-1][1] == a:
                    ivs[-1] = (ivs[-1][0], b)
                else:
                    ivs.append((a, b))
        for p in bp:
            if mem(p) and not any(a <= p <= b for a, b in ivs):
                pts.append(p)
        return ivs, pts, bp

    def Zset(self, s, ivs, pts, bp):
        """Z(s): squares S with Phi(S) cap Y_b cap W_s nonempty.  returns dict sid -> witness y."""
        out = {}
        for S in self.sq:
            if S.s == 0:
                continue
            sca = abs(S.s) * S.c
            for (p, q) in ((S.ymin, S.ymin + sca), (S.ymax - sca, S.ymax)):
                inner = [z for z in bp if p < z < q] + [z for z in pts if p < z < q]
                tests = sorted(set(inner))
                pl = [p] + tests + [q]
                cand = tests + [(pl[i] + pl[i + 1]) / 2 for i in range(len(pl) - 1)]
                hit = None
                lo_w, hi_w = (1 - self.P.eps) * Fr(s), Fr(s)
                for z in cand:
                    if lo_w <= z <= hi_w and self.in_Yb(z):
                        hit = z; break
                if hit is not None:
                    out[S.id] = hit; break
        return out

    def integral_I(self, ivs):
        """int over ivs of (1 - omega0 - E(y))_+ dy (exact; E piecewise linear)."""
        P = self.P
        tot = Fr(0)
        for (a, b) in ivs:
            pts = sorted(set([a, b] + [z for z in self.eb if a < z < b]))
            for i in range(len(pts) - 1):
                u, v = pts[i], pts[i + 1]
                y1 = u + (v - u) / 3; y2 = u + 2 * (v - u) / 3
                f1 = 1 - P.omega0 - self.E(y1); f2 = 1 - P.omega0 - self.E(y2)
                sl = (f2 - f1) / (y2 - y1)
                fu = f1 + sl * (u - y1); fv = f1 + sl * (v - y1)
                if fu >= 0 and fv >= 0:
                    tot += (fu + fv) / 2 * (v - u)
                elif fu > 0 or fv > 0:
                    r = u + (0 - fu) / sl
                    if fu > 0:
                        tot += fu / 2 * (r - u)
                    else:
                        tot += fv / 2 * (v - r)
        return tot

    def J_Z(self, Z, ivs):
        """sum over Z of int over Phi(Z) cap ivs of c_Z."""
        tot = Fr(0)
        for sid in Z:
            S = next(T for T in self.sq if T.id == sid)
            sca = abs(S.s) * S.c
            for (p, q, lower) in ((S.ymin, S.ymin + sca, True), (S.ymax - sca, S.ymax, False)):
                for (a, b) in ivs:
                    u, v = max(a, p), min(b, q)
                    if u < v:
                        if lower:
                            tot += ((u - p) + (v - p)) / 2 * (v - u) / sca
                        else:
                            tot += ((q - u) + (q - v)) / 2 * (v - u) / sca
        return tot


# ------------------------------------------------------------------ P2 test on all live F_s contacts
def ramp_maxdist(S):
    """max over the ramp set R(S) of dist(y, Z) (exact); None if a=0."""
    if S.s == 0:
        return None
    sa = abs(S.s)
    best = Fr(0)
    for (p, q) in ((S.ymin, S.ymin + sa), (S.ymax - sa, S.ymax)):
        m = math.ceil(p - HALF)          # smallest integer m with m + 1/2 >= p
        if m + HALF <= q:
            return HALF
        for z in (p, q):
            d = min(z - math.floor(z), math.ceil(z) - z)
            best = max(best, d)
    return best


def p2_test(flow, P, s, tm):
    """for every pre-R1 entry candidate that is alive in F_s with q_y <= s+2: maxdist(R(Y)) <= w(s,Y,q_y)."""
    H = Fr(s) + 2
    sec1 = mpmath.sec(P.alpha(s)) - 1
    worst = None; viol = []; ncheck = 0
    live = {}    # sid -> measure of live F_s contact
    for cd in flow.cands:
        if any(tm[sid] for (sid, idx) in cd.ent):
            continue
        qy = (cd.q[2], cd.q[3])
        for (a, b) in split_affine(cd.xlo, cd.xhi, [qy], H):
            m = (a + b) / 2
            if qy[0] + qy[1] * m > H:
                continue
            live[cd.sid] = live.get(cd.sid, Fr(0)) + (b - a)
            S = flow.sq[cd.sid]
            md = ramp_maxdist(S)
            if md is None:
                continue
            ncheck += 1
            qmin = min(qy[0] + qy[1] * a, qy[0] + qy[1] * b)
            sa = mpmath.sin(S.a); ca = mpmath.cos(S.a)
            w = sec1 * mpf_fr(qmin) + mpf_fr(P.delta) + sa + 1 - ca
            ratio = mpf_fr(md) / w
            if worst is None or ratio > worst[0]:
                worst = (ratio, cd.sid, float(qmin), float(md), float(w))
            if mpf_fr(md) > w:
                viol.append((cd.sid, float(qmin), float(md), float(w), float(a), float(b)))
    return dict(worst=worst, viol=viol, ncheck=ncheck, live=live)
