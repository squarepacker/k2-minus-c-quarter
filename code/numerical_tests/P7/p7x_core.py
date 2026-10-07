"""p7x_core.py -- independent exact (rational) flow tracer for the test of Section 9 (P7) of the manuscript.

Written from scratch for this test from the definitions of Section 3 of the manuscript; no other code is
imported, run or copied.

Geometry
  square: closed unit square, phase phi in (-pi/4, pi/4) given by rational t = tan(phi/2):
          cos = (1-t^2)/(1+t^2), sin = 2t/(1+t^2) (exact rationals);  u = (cos, sin), n = (-sin, cos)
          bot = c - n/2 + tau*u (|tau| <= 1/2), top = c + n/2 + tau*u, right = c + u/2 + tau*n, left = c - u/2 + tau*n
  (exactly 45 degrees is not representable with rationals; such squares never occur here.)
Flow (Definitions 3.9, 3.10, 3.14 and 3.20 of the manuscript: the master flow F^0)
  state (p, d, X, g); ray step with t_S, t_W, t_D = delta - g, t_H; priority D > W > contact > H.
  contact on relint bot(Y) -> entry candidate; any other first contact (side relint, vertex) -> E.
  R1: squares processed in increasing centre height; at each entry point among live arrivals the
      lexicographically smallest (g, x) wins, the others are M-terminated.
  R2: winner with a(Y) >= alpha_F -> T at q; else pass along n_Y to q + n_Y; exit height >= h_F -> H.
Implementation: beams = open x-intervals with affine (in x) state; the first-hit map is computed as
  the lower envelope of the facing edges in the ray frame (breakpoints = projected vertices), and the
  event choice by splitting at crossings of the affine event times.  Measure-zero x are dropped.
"""
from fractions import Fraction as Fr
import mpmath

mpmath.mp.dps = 60
HALF = Fr(1, 2)


def mpf_fr(q):
    return mpmath.mpf(q.numerator) / q.denominator


class Sq:
    __slots__ = ('id', 'cx', 'cy', 't', 'c', 's', 'u', 'n', 'V', 'E', 'ymin', 'ymax',
                 'xmin', 'xmax', 'tan_a', 'a', 'af', 'Vf', 'order')

    def __init__(self, sid, cx, cy, t):
        t = Fr(t); cx = Fr(cx); cy = Fr(cy)
        den = 1 + t * t
        c = (1 - t * t) / den
        s = 2 * t / den
        if not (c > 0 and abs(s) < c):
            raise ValueError('phase must be in (-pi/4, pi/4)')
        self.id = sid; self.cx = cx; self.cy = cy; self.t = t; self.c = c; self.s = s
        self.u = (c, s); self.n = (-s, c)
        hx, hy = c * HALF, s * HALF          # u/2
        nx, ny = -s * HALF, c * HALF         # n/2
        BL = (cx - hx - nx, cy - hy - ny)
        BR = (cx + hx - nx, cy + hy - ny)
        TR = (cx + hx + nx, cy + hy + ny)
        TL = (cx - hx + nx, cy - hy + ny)
        self.V = (BL, BR, TR, TL)
        # edges: name -> (P, Q, outward normal)
        self.E = {'bot': (BL, BR, (s, -c)), 'right': (BR, TR, (c, s)),
                  'top': (TL, TR, (-s, c)), 'left': (BL, TL, (-c, -s))}
        ys = [v[1] for v in self.V]; xs = [v[0] for v in self.V]
        self.ymin = min(ys); self.ymax = max(ys); self.xmin = min(xs); self.xmax = max(xs)
        self.tan_a = abs(s) / c
        self.a = mpmath.atan(mpf_fr(self.tan_a))
        self.af = float(self.a)
        self.Vf = tuple((float(v[0]), float(v[1])) for v in self.V)
        self.order = None

    def chord(self, y):
        """length of the horizontal chord of the closed square at height y (exact)."""
        if y < self.ymin or y > self.ymax:
            return Fr(0)
        xs = []
        for (P, Q, _) in self.E.values():
            if P[1] == Q[1]:
                if P[1] == y:
                    xs += [P[0], Q[0]]
            elif min(P[1], Q[1]) <= y <= max(P[1], Q[1]):
                xs.append(P[0] + (y - P[1]) * (Q[0] - P[0]) / (Q[1] - P[1]))
        return (max(xs) - min(xs)) if xs else Fr(0)

    def contains_open(self, X, Y):
        dx = X - self.cx; dy = Y - self.cy
        xi = dx * self.c + dy * self.s
        ze = -dx * self.s + dy * self.c
        return abs(xi) < HALF and abs(ze) < HALF

    def bottom_vertex_height(self):
        return self.ymin


def sat_disjoint(A, B):
    """closed squares disjoint <=> strict separation along one of the 4 edge normals (exact)."""
    for (ax, ay) in (A.u, A.n, B.u, B.n):
        pa = [v[0] * ax + v[1] * ay for v in A.V]
        pb = [v[0] * ax + v[1] * ay for v in B.V]
        if max(pa) < min(pb) or max(pb) < min(pa):
            return True
    return False


def validate(squares, k):
    """exact check: inside [0,k]^2 and pairwise disjoint (closed).  Returns list of problems."""
    prob = []
    for S in squares:
        for v in S.V:
            if not (0 <= v[0] <= k and 0 <= v[1] <= k):
                prob.append(('outside', S.id)); break
    srt = sorted(squares, key=lambda S: S.xmin)
    for i, A in enumerate(srt):
        for B in srt[i + 1:]:
            if B.xmin > A.xmax:
                break
            if B.ymin > A.ymax or A.ymin > B.ymax:
                continue
            if not sat_disjoint(A, B):
                prob.append(('overlap', A.id, B.id))
    return prob


# ---------------- affine helpers (functions of the floor coordinate x) ----------------
# scalar affine: (a, b) = a + b x ; point affine: (ax, bx, ay, by)

def pev(P, x):
    return (P[0] + P[1] * x, P[2] + P[3] * x)


def lin_cross(f, g):
    """u where f(u) = g(u) for linear f=(a,b): a + b u ; None if parallel."""
    if f[1] == g[1]:
        return None
    return (g[0] - f[0]) / (f[1] - g[1])


class Beam:
    __slots__ = ('xlo', 'xhi', 'P', 'd', 'src', 'g', 'way', 'seg', 'dirs', 'ent')

    def __init__(self, xlo, xhi, P, d, src, g, way, seg, dirs, ent):
        self.xlo = xlo; self.xhi = xhi; self.P = P; self.d = d; self.src = src; self.g = g
        self.way = way; self.seg = seg; self.dirs = dirs; self.ent = ent


class Rec:
    """final piece of the master flow F^0.
    way: waypoints (affine points), way[0]=(x,0), last = termination point
    seg[i]: label of segment way[i]->way[i+1]: -1 gap, sid passage inside square sid
    dirs[i]: unit direction of segment i
    ent: tuple of (sid, index into way of the entry point) for R1-won entries (incl. a final T)
    typ: D W E M T H ; tsq: square of the E/M/T contact"""
    __slots__ = ('xlo', 'xhi', 'typ', 'way', 'seg', 'dirs', 'ent', 'tsq')

    def __init__(self, xlo, xhi, typ, way, seg, dirs, ent, tsq=None):
        self.xlo = xlo; self.xhi = xhi; self.typ = typ; self.way = way; self.seg = seg
        self.dirs = dirs; self.ent = ent; self.tsq = tsq


class Cand:
    """entry candidate (live contact with relint bot(Y), before R1)."""
    __slots__ = ('sid', 'xlo', 'xhi', 'q', 'g', 'way', 'seg', 'dirs', 'ent')

    def __init__(self, sid, xlo, xhi, q, g, way, seg, dirs, ent):
        self.sid = sid; self.xlo = xlo; self.xhi = xhi; self.q = q; self.g = g
        self.way = way; self.seg = seg; self.dirs = dirs; self.ent = ent


class Flow:
    """master flow F(alpha_F, delta, h_F) for a closed configuration in [0,k]^2."""

    def __init__(self, squares, k, delta, alpha_F, hF):
        self.sq = {S.id: S for S in squares}
        self.k = Fr(k); self.delta = Fr(delta); self.hF = Fr(hF)
        self.alpha_F = mpmath.mpf(alpha_F)
        tanF = mpmath.tan(self.alpha_F)
        self.passes = {}
        self.near = []
        for S in squares:
            diff = mpf_fr(S.tan_a) - tanF
            if abs(diff) < mpmath.mpf(10) ** -40:
                self.near.append(S.id)
            self.passes[S.id] = diff < 0          # a(S) < alpha_F
        self.order = sorted(squares, key=lambda S: (S.cy, S.id))
        for i, S in enumerate(self.order):
            S.order = i
        self.recs = []
        self.cands = []           # all entry candidates (pre-R1) -- for the P2 test
        self.pending = {S.id: [] for S in squares}
        self.cur = -1
        self.anom = []
        self.nray = 0

    # ---------------------------------------------------------------- ray step
    def ray_step(self, bm):
        self.nray += 1
        dx, dy = bm.d
        ex, ey = dy, -dx
        P = bm.P
        U0 = P[0] * ex + P[2] * ey; U1 = P[1] * ex + P[3] * ey
        T0 = P[0] * dx + P[2] * dy; T1 = P[1] * dx + P[3] * dy
        if T1 != 0 or U1 == 0:
            self.anom.append(('frame', T1, U1)); return
        uA = U0 + U1 * bm.xlo; uB = U0 + U1 * bm.xhi
        ulo, uhi = (uA, uB) if uA < uB else (uB, uA)
        g0, g1 = bm.g
        gu = g1 / U1; gc = g0 - gu * U0
        fD = (self.delta - gc, -gu)
        pxu = P[1] / U1; pxc = P[0] - pxu * U0
        pyu = P[3] / U1; pyc = P[2] - pyu * U0
        if dx > 0:
            fW = ((self.k - pxc) / dx, -pxu / dx)
        elif dx < 0:
            fW = ((0 - pxc) / dx, -pxu / dx)
        else:
            fW = None
        fH = ((self.hF - pyc) / dy, -pyu / dy)
        reachD = max(fD[0] + fD[1] * ulo, fD[0] + fD[1] * uhi)
        reachH = max(fH[0] + fH[1] * ulo, fH[0] + fH[1] * uhi)
        reach = min(reachD, reachH)
        # float prefilter
        dxf, dyf, exf, eyf = float(dx), float(dy), float(ex), float(ey)
        ulf, uhf, T0f, rf = float(ulo), float(uhi), float(T0), float(reach)
        edges = []
        for S in self.order:
            if S.id == bm.src:
                continue
            us = [vx * exf + vy * eyf for (vx, vy) in S.Vf]
            if max(us) < ulf - 1e-9 or min(us) > uhf + 1e-9:
                continue
            ts = [vx * dxf + vy * dyf for (vx, vy) in S.Vf]
            if min(ts) > T0f + rf + 1e-9 or max(ts) < T0f - 1e-9:
                continue
            fr = [(v[0] * ex + v[1] * ey, v[0] * dx + v[1] * dy) for v in S.V]
            VV = dict(zip(('BL', 'BR', 'TR', 'TL'), fr))
            for name, (a, b) in (('bot', ('BL', 'BR')), ('right', ('BR', 'TR')),
                                 ('top', ('TL', 'TR')), ('left', ('BL', 'TL'))):
                nu = S.E[name][2]
                if nu[0] * dx + nu[1] * dy >= 0:
                    continue
                (u1, t1), (u2, t2) = VV[a], VV[b]
                if u1 == u2:
                    continue
                if u1 > u2:
                    u1, t1, u2, t2 = u2, t2, u1, t1
                if u2 <= ulo or u1 >= uhi:
                    continue
                sl = (t2 - t1) / (u2 - u1)
                lf = (t1 - sl * u1, sl)              # t(u) = lf0 + lf1 u
                um = (max(u1, ulo) + min(u2, uhi)) / 2
                tm = lf[0] + lf[1] * um
                if bm.src >= 0:
                    if tm <= T0:
                        if tm > T0 - 1:
                            self.anom.append(('inside_src_column', bm.src, S.id))
                        continue
                else:
                    if tm < T0:
                        self.anom.append(('below_floor', S.id)); continue
                edges.append((u1, u2, lf, S.id, name))
        bps = {ulo, uhi}
        for (u1, u2, lf, sid, name) in edges:
            if ulo < u1 < uhi:
                bps.add(u1)
            if ulo < u2 < uhi:
                bps.add(u2)
        bps = sorted(bps)
        outs = []
        for i in range(len(bps) - 1):
            ul, ur = bps[i], bps[i + 1]
            um = (ul + ur) / 2
            best = None
            for (u1, u2, lf, sid, name) in edges:
                if u1 <= ul and u2 >= ur:
                    tv = lf[0] + lf[1] * um
                    if best is None or tv < best[0]:
                        best = (tv, lf, sid, name)
                    elif tv == best[0] and sid != best[2]:
                        self.anom.append(('tie_two_squares', sid, best[2]))
            funcs = [('D', fD, None, None)]
            if fW is not None:
                funcs.append(('W', fW, None, None))
            if best is not None:
                lf = best[1]
                funcs.append(('S', (lf[0] - T0, lf[1]), best[2], best[3]))
            funcs.append(('H', fH, None, None))
            cps = {ul, ur}
            for a in range(len(funcs)):
                for b in range(a + 1, len(funcs)):
                    c = lin_cross(funcs[a][1], funcs[b][1])
                    if c is not None and ul < c < ur:
                        cps.add(c)
            cps = sorted(cps)
            for j in range(len(cps) - 1):
                a_, b_ = cps[j], cps[j + 1]
                m = (a_ + b_) / 2
                bestf = None
                for (nm, f, sid, ename) in funcs:   # list order = priority D > W > S > H
                    v = f[0] + f[1] * m
                    if bestf is None or v < bestf[0]:
                        bestf = (v, nm, f, sid, ename)
                if bestf[0] < 0:
                    self.anom.append(('negative_time', bestf[1]))
                key = (bestf[1], bestf[3], bestf[4])
                if outs and outs[-1][0] == key and outs[-1][2] == a_:
                    outs[-1][2] = b_
                else:
                    outs.append([key, a_, b_, bestf[2]])
        for (key, a_, b_, f) in outs:
            nm, sid, ename = key
            xa = (a_ - U0) / U1; xb = (b_ - U0) / U1
            if xa > xb:
                xa, xb = xb, xa
            # event time as affine function of x: f0 + f1*(U0 + U1 x)
            tx = (f[0] + f[1] * U0, f[1] * U1)
            pt = (P[0] + dx * tx[0], P[1] + dx * tx[1], P[2] + dy * tx[0], P[3] + dy * tx[1])
            way = bm.way + (pt,)
            seg = bm.seg + (-1,)
            dirs = bm.dirs + (bm.d,)
            if nm == 'S':
                if ename == 'bot':
                    gnew = (bm.g[0] + tx[0], bm.g[1] + tx[1])
                    S = self.sq[sid]
                    if S.order <= self.cur:
                        self.anom.append(('dag_violation', bm.src, sid))
                    cd = Cand(sid, xa, xb, pt, gnew, way, seg, dirs, bm.ent)
                    self.pending[sid].append(cd)
                    self.cands.append(cd)
                else:
                    self.recs.append(Rec(xa, xb, 'E', way, seg, dirs, bm.ent, sid))
            else:
                self.recs.append(Rec(xa, xb, nm, way, seg, dirs, bm.ent, None))

    # ---------------------------------------------------------------- R1 at Y
    def resolve(self, Y):
        cands = self.pending[Y.id]
        items = []
        for cd in cands:
            q = cd.q
            A = (q[0] - Y.cx) * Y.c + (q[2] - Y.cy) * Y.s
            B = q[1] * Y.c + q[3] * Y.s
            if B == 0:
                self.anom.append(('degenerate_bot_map', Y.id)); continue
            ta, tb = A + B * cd.xlo, A + B * cd.xhi
            lo, hi = (ta, tb) if ta < tb else (tb, ta)
            if lo < -HALF or hi > HALF:
                self.anom.append(('entry_outside_bot', Y.id, float(lo), float(hi)))
            # x(tau) = (tau - A)/B -> (-A/B, 1/B); g(tau) = g0 + g1 x(tau)
            xl = (-A / B, 1 / B)
            gl = (cd.g[0] + cd.g[1] * xl[0], cd.g[1] * xl[1])
            items.append((lo, hi, xl, gl, cd))
        crit = set()
        for it in items:
            crit.add(it[0]); crit.add(it[1])
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                a, b = items[i], items[j]
                lo = max(a[0], b[0]); hi = min(a[1], b[1])
                if lo >= hi:
                    continue
                for (f, g) in ((a[3], b[3]), (a[2], b[2])):
                    c = lin_cross(f, g)
                    if c is not None and lo < c < hi:
                        crit.add(c)
        crit = sorted(crit)
        res = {}   # id(cd) -> list of [role, tlo, thi]
        for i in range(len(crit) - 1):
            a_, b_ = crit[i], crit[i + 1]
            m = (a_ + b_) / 2
            cov = [it for it in items if it[0] <= a_ and it[1] >= b_]
            if not cov:
                continue
            keyed = []
            for it in cov:
                gv = it[3][0] + it[3][1] * m
                xv = it[2][0] + it[2][1] * m
                keyed.append(((gv, xv), it))
            keyed.sort(key=lambda z: z[0])
            if len(keyed) > 1 and keyed[0][0] == keyed[1][0]:
                self.anom.append(('identical_arrivals', Y.id))
            for r, (kv, it) in enumerate(keyed):
                role = 'W' if r == 0 else 'M'
                lst = res.setdefault(id(it[4]), [])
                if lst and lst[-1][0] == role and lst[-1][2] == a_:
                    lst[-1][2] = b_
                else:
                    lst.append([role, a_, b_])
        cdmap = {id(it[4]): it for it in items}
        for key, lst in res.items():
            lo, hi, xl, gl, cd = cdmap[key]
            for role, ta, tb in lst:
                xa = xl[0] + xl[1] * ta; xb = xl[0] + xl[1] * tb
                if xa > xb:
                    xa, xb = xb, xa
                if role == 'M':
                    self.recs.append(Rec(xa, xb, 'M', cd.way, cd.seg, cd.dirs, cd.ent, Y.id))
                    continue
                idx = len(cd.way) - 1
                ent = cd.ent + ((Y.id, idx),)
                if not self.passes[Y.id]:
                    self.recs.append(Rec(xa, xb, 'T', cd.way, cd.seg, cd.dirs, ent, Y.id))
                    continue
                q = cd.q
                ex = (q[0] + Y.n[0], q[1], q[2] + Y.n[1], q[3])     # exit point q + n_Y
                way = cd.way + (ex,)
                seg = cd.seg + (Y.id,)
                dirs = cd.dirs + (Y.n,)
                # split at exit height = hF
                pieces = [(xa, xb)]
                if ex[3] != 0:
                    xr = (self.hF - ex[2]) / ex[3]
                    if xa < xr < xb:
                        pieces = [(xa, xr), (xr, xb)]
                for (pa, pb) in pieces:
                    hm = ex[2] + ex[3] * ((pa + pb) / 2)
                    if hm >= self.hF:
                        self.recs.append(Rec(pa, pb, 'H', way, seg, dirs, ent, None))
                    else:
                        self.ray_step(Beam(pa, pb, ex, Y.n, Y.id, cd.g, way, seg, dirs, ent))

    def run(self):
        start = (Fr(0), Fr(1), Fr(0), Fr(0))
        self.ray_step(Beam(Fr(0), self.k, start, (Fr(0), Fr(1)), -1, (Fr(0), Fr(0)),
                           (start,), (), (), ()))
        for i, Y in enumerate(self.order):
            self.cur = i
            if self.pending[Y.id]:
                self.resolve(Y)
        tot = sum((r.xhi - r.xlo) for r in self.recs)
        if tot != self.k:
            self.anom.append(('measure_sum', float(tot)))
        return self


# ---------------------------------------------------------------- independent point tracer (no R1)
def clip(p, d, S):
    """parameter interval of the ray p + t d inside the closed square S (exact), or None."""
    dx0 = p[0] - S.cx; dy0 = p[1] - S.cy
    lo, hi = None, None
    for (ax, ay) in (S.u, S.n):
        o = dx0 * ax + dy0 * ay
        v = d[0] * ax + d[1] * ay
        if v == 0:
            if abs(o) > HALF:
                return None
            continue
        t1 = (-HALF - o) / v; t2 = (HALF - o) / v
        if t1 > t2:
            t1, t2 = t2, t1
        lo = t1 if lo is None or t1 > lo else lo
        hi = t2 if hi is None or t2 < hi else hi
    if lo is None or lo > hi:
        return None
    return lo, hi


def trace_point(x, squares, k, delta, alpha_F, hF, max_steps=10000):
    """trace one path without merging. Returns (type, list of passed square ids, termination point, g, contacts)."""
    tanF = mpmath.tan(mpmath.mpf(alpha_F))
    k = Fr(k); delta = Fr(delta); hF = Fr(hF)
    p = (Fr(x), Fr(0)); d = (Fr(0), Fr(1)); src = None; g = Fr(0)
    passed = []; contacts = []
    for _ in range(max_steps):
        tS, Ysel = None, None
        for S in squares:
            if S is src:
                continue
            r = clip(p, d, S)
            if r is None or r[1] < 0:
                continue
            te = r[0] if r[0] > 0 else Fr(0)
            if tS is None or te < tS:
                tS, Ysel = te, S
        tD = delta - g
        if d[0] > 0:
            tW = (k - p[0]) / d[0]
        elif d[0] < 0:
            tW = (0 - p[0]) / d[0]
        else:
            tW = None
        tH = (hF - p[1]) / d[1]
        cand = [('D', tD), ('W', tW), ('S', tS), ('H', tH)]
        tmin = min(v for (_, v) in cand if v is not None)
        for nm, v in cand:
            if v is not None and v == tmin:
                ev = nm; break
        if ev == 'D':
            return ('D', passed, (p[0] + tD * d[0], p[1] + tD * d[1]), delta, contacts)
        if ev == 'W':
            return ('W', passed, (p[0] + tW * d[0], p[1] + tW * d[1]), g + tW, contacts)
        if ev == 'H':
            return ('H', passed, (p[0] + tH * d[0], p[1] + tH * d[1]), g + tH, contacts)
        q = (p[0] + tS * d[0], p[1] + tS * d[1]); g = g + tS
        Y = Ysel
        xi = (q[0] - Y.cx) * Y.c + (q[1] - Y.cy) * Y.s
        ze = -(q[0] - Y.cx) * Y.s + (q[1] - Y.cy) * Y.c
        contacts.append((Y.id, q, g))
        if ze == -HALF and abs(xi) < HALF:
            if mpf_fr(Y.tan_a) >= tanF:
                return ('T', passed, q, g, contacts)
            passed.append(Y.id)
            p = (q[0] + Y.n[0], q[1] + Y.n[1]); d = Y.n; src = Y
            if p[1] >= hF:
                return ('H', passed, p, g, contacts)
        else:
            return ('E', passed, q, g, contacts)
    return ('LOOP', passed, p, g, contacts)
