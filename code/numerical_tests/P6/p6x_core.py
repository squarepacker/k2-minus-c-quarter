# p6x_core.py -- independent EXACT (rational) implementation of the master flow
# F^0 = F(alpha_max, delta, h_max) of Section 3 of the manuscript (Definition 3.9: ray steps and the
# priority D > W > contact > H; Definition 3.10: rule R2 and passing; Definition 3.14: rule R1;
# Definition 3.20: master flow; M is decided before T at an entry candidate) and of its ceiling
# version (Definition 3.24, by reflection).  It is used by p6x_checks.py to test Section 8
# (Lemmas 8.5, 8.6, 8.7, 8.9 and Step 1 of the proof of Theorem 8.3).
# Written from scratch for this test; it does not import or run any other code.
#
# Exactness: every square has a Pythagorean rotation (cos, sin rational), centres rational,
# delta, k, h rational, tilt threshold given by a rational tangent Ta (alpha = atan Ta).
# Then all path positions / gaps / breakpoints are rational and every comparison is exact.
from fractions import Fraction as Fr
from collections import defaultdict
import math
import mpmath as mp

mp.mp.dps = 40
HALF = Fr(1, 2)


def mpf(q):
    return mp.mpf(q.numerator) / q.denominator


def pyth(t):
    t = Fr(t)
    d = 1 + t * t
    return ((1 - t * t) / d, 2 * t / d)


def rot_for_angle(phi, maxden=10 ** 5):
    """rational unit vector at angle ~phi (|phi| < pi)."""
    t = Fr(math.tan(phi / 2)).limit_denominator(maxden)
    return pyth(t)


def canon(u):
    ux, uy = u
    for c in ((ux, uy), (-uy, ux), (-ux, -uy), (uy, -ux)):
        if c[0] > abs(c[1]):
            return c
    raise ValueError("exact 45 deg impossible for rationals")


class Sq:
    __slots__ = ('id', 'c', 'u', 'n', 'V', 'E', 'tan_abs', 'phi', 'bb', 'cy')

    def __init__(self, id, c, u):
        u = canon(u)
        self.id = id
        self.c = (Fr(c[0]), Fr(c[1]))
        self.u = u
        self.n = (-u[1], u[0])
        cx, cy = self.c
        ux, uy = u
        nx, ny = self.n
        h = HALF
        self.V = [(cx - h * nx - h * ux, cy - h * ny - h * uy), (cx - h * nx + h * ux, cy - h * ny + h * uy),
                  (cx + h * nx + h * ux, cy + h * ny + h * uy), (cx + h * nx - h * ux, cy + h * ny - h * uy)]
        cn = nx * cx + ny * cy
        cu = ux * cx + uy * cy
        self.E = [('bot', (-nx, -ny), -cn + h), ('right', (ux, uy), cu + h),
                  ('top', (nx, ny), cn + h), ('left', (-ux, -uy), -cu + h)]
        self.tan_abs = abs(uy / ux)
        self.phi = mp.atan2(mpf(uy), mpf(ux))
        xs = [float(v[0]) for v in self.V]
        ys = [float(v[1]) for v in self.V]
        self.bb = (min(xs), min(ys), max(xs), max(ys))
        self.cy = cy


def reflect_squares(sqs, k):
    out = []
    for S in sqs:
        out.append(Sq(S.id, (S.c[0], k - S.c[1]), (S.u[0], -S.u[1])))
    return out


# ---------------- exact polygon tools ----------------
def _axes(P):
    m = len(P)
    for i in range(m):
        a, b = P[i], P[(i + 1) % m]
        ex, ey = b[0] - a[0], b[1] - a[1]
        if ex != 0 or ey != 0:
            yield (-ey, ex)


def _proj(P, ax):
    vals = [p[0] * ax[0] + p[1] * ax[1] for p in P]
    return min(vals), max(vals)


def closed_disjoint(A, B):
    for ax in list(_axes(A)) + list(_axes(B)):
        a0, a1 = _proj(A, ax)
        b0, b1 = _proj(B, ax)
        if a1 < b0 or b1 < a0:
            return True
    return False


def interiors_disjoint(A, B):
    for ax in list(_axes(A)) + list(_axes(B)):
        a0, a1 = _proj(A, ax)
        b0, b1 = _proj(B, ax)
        if a1 <= b0 or b1 <= a0:
            return True
    return False


def sat_gap(A, B):
    """float: max over axes of normalised separation (>0 separated, <0 penetration depth)."""
    best = -1e300
    for ax in list(_axes(A)) + list(_axes(B)):
        L = math.hypot(float(ax[0]), float(ax[1]))
        a0, a1 = _proj(A, ax)
        b0, b1 = _proj(B, ax)
        gap = max(float(b0 - a1), float(a0 - b1)) / L
        if gap > best:
            best = gap
    return best


def poly_area2(P):
    s = Fr(0)
    m = len(P)
    for i in range(m):
        s += P[i][0] * P[(i + 1) % m][1] - P[(i + 1) % m][0] * P[i][1]
    return s


def dedup(P):
    out = []
    for p in P:
        if not out or out[-1] != p:
            out.append(p)
    if len(out) > 1 and out[0] == out[-1]:
        out.pop()
    return out


def bbox(P):
    xs = [float(p[0]) for p in P]
    ys = [float(p[1]) for p in P]
    return (min(xs), min(ys), max(xs), max(ys))


def bb_hit(a, b, m=1e-9):
    return not (a[2] < b[0] - m or b[2] < a[0] - m or a[3] < b[1] - m or b[3] < a[1] - m)


# ---------------- ray / square hit (exact) ----------------
def hit(Z, p, d):
    """first-hit parameter of ray p+t d (t>=0) on closed square Z, and edge name.
    returns None (miss), or (t, edge, degenerate_flag)."""
    t_in = None
    e_in = None
    t_out = None
    degen = False
    for (nm, nu, ce) in Z.E:
        nd = nu[0] * d[0] + nu[1] * d[1]
        val = ce - (nu[0] * p[0] + nu[1] * p[1])
        if nd == 0:
            if val < 0:
                return None
            if val == 0:
                degen = True
            continue
        t = val / nd
        if nd < 0:
            if t_in is None or t > t_in:
                t_in, e_in = t, nm
            elif t == t_in:
                degen = True
        else:
            if t_out is None or t < t_out:
                t_out = t
    if t_in > t_out:
        return None
    if t_out < 0:
        return None
    if t_in == t_out:
        degen = True
    if t_in < 0:
        raise RuntimeError("ray start strictly inside a square")
    return (t_in, e_in, degen)


def lin_root(a, b):
    if a[1] == b[1]:
        return None
    return (b[0] - a[0]) / (a[1] - b[1])


class Beam:
    __slots__ = ('xa', 'xb', 'P0', 'P1', 'd', 'last', 'g', 'wex', 'hist')

    def __init__(self, xa, xb, P0, P1, d, last, g, wex, hist):
        self.xa, self.xb, self.P0, self.P1, self.d = xa, xb, P0, P1, d
        self.last, self.g, self.wex, self.hist = last, g, wex, hist


class Flow:
    def __init__(self, sqs, k, delta, Ta, h, fam='floor'):
        self.sqs = sqs
        self.k = Fr(k)
        self.delta = Fr(delta)
        self.Ta = Fr(Ta)
        self.h = Fr(h)
        self.fam = fam
        self.byid = {S.id: S for S in sqs}
        self.anom = []          # anomalies (should stay empty)
        self.degen_mid = 0
        self.variant = set()    # negative controls only: 'no_R1', 'T_before_M', 'no_R2'

    # ---- one ray step of a beam: returns list of (xa,xb,kind,Zid,edge)
    def ray_step(self, B):
        k, delta, h = self.k, self.delta, self.h
        xa, xb, P0, P1, d = B.xa, B.xb, B.P0, B.P1, B.d
        dx, dy = d
        g0, g1 = B.g
        pa = (P0[0] + xa * P1[0], P0[1] + xa * P1[1])
        pb = (P0[0] + xb * P1[0], P0[1] + xb * P1[1])
        tmax = float(delta - min(g0 + g1 * xa, g0 + g1 * xb))
        pts = [pa, pb]
        fx = [float(p[0]) for p in pts] + [float(p[0]) + tmax * float(dx) for p in pts]
        fy = [float(p[1]) for p in pts] + [float(p[1]) + tmax * float(dy) for p in pts]
        bb = (min(fx), min(fy), max(fx), max(fy))
        rel = [Z for Z in self.sqs if Z.id != B.last and bb_hit(bb, Z.bb, 1e-7)]
        tD = (delta - g0, -g1)
        tH = ((h - P0[1]) / dy, -P1[1] / dy)
        tW = None
        if dx > 0:
            tW = ((k - P0[0]) / dx, -P1[0] / dx)
        elif dx < 0:
            tW = ((-P0[0]) / dx, -P1[0] / dx)
        base = [tD, tH, (Fr(0), Fr(0))] + ([tW] if tW is not None else [])
        bps = set()

        def add(r):
            if r is not None and xa < r < xb:
                bps.add(r)
        for i in range(len(base)):
            for j in range(i + 1, len(base)):
                add(lin_root(base[i], base[j]))
        cP1d = P1[0] * dy - P1[1] * dx
        for Z in rel:
            for v in Z.V:
                c0 = (v[0] - P0[0]) * dy - (v[1] - P0[1]) * dx
                if cP1d != 0:
                    add(c0 / cP1d)
            for (nm, nu, ce) in Z.E:
                nd = nu[0] * dx + nu[1] * dy
                np0 = nu[0] * P0[0] + nu[1] * P0[1]
                np1 = nu[0] * P1[0] + nu[1] * P1[1]
                if nd == 0:
                    if np1 != 0:
                        add((ce - np0) / np1)
                else:
                    te = ((ce - np0) / nd, -np1 / nd)
                    for b in base:
                        add(lin_root(te, b))
        cuts = [xa] + sorted(bps) + [xb]
        out = []
        for i in range(len(cuts) - 1):
            a, b = cuts[i], cuts[i + 1]
            xm = (a + b) / 2
            p = (P0[0] + xm * P1[0], P0[1] + xm * P1[1])
            best = None
            for Z in rel:
                r = hit(Z, p, d)
                if r is None:
                    continue
                if r[2]:
                    self.degen_mid += 1
                if best is None or r[0] < best[0]:
                    best = (r[0], Z.id, r[1])
                elif r[0] == best[0]:
                    self.anom.append(('tie_two_squares', float(xm)))
            tDv = tD[0] + tD[1] * xm
            tHv = tH[0] + tH[1] * xm
            tWv = (tW[0] + tW[1] * xm) if tW is not None else None
            cand = [tDv, tHv]
            if tWv is not None:
                cand.append(tWv)
            if best is not None:
                cand.append(best[0])
            tstar = min(cand)
            if tDv == tstar:
                o = ('D', None, None)
            elif tWv is not None and tWv == tstar:
                o = ('W', None, None)
            elif best is not None and best[0] == tstar:
                if best[2] == 'bot':
                    o = ('C', best[1], 'bot')
                else:
                    if best[2] == 'top':
                        self.anom.append(('top_first_contact', float(xm)))
                    o = ('E', best[1], best[2])
            else:
                o = ('H', None, None)
            if out and out[-1][2:] == o and out[-1][1] == a:
                out[-1] = (out[-1][0], b) + o
            else:
                out.append((a, b) + o)
        return out

    def term(self, kind, xa, xb, hist, extra=None):
        self.terms.append((kind, xa, xb, hist, extra))

    def advance(self, B):
        """ray step of beam B, dispatch outcomes"""
        for (a, b, kind, Zid, ed) in self.ray_step(B):
            if kind == 'E' and 'no_R2' in self.variant and ed in ('left', 'right') \
                    and self.byid[Zid].tan_abs >= self.Ta:
                w0, w1 = B.wex
                wa, wb = sorted((w0 + w1 * a, w0 + w1 * b))
                self.pairs[(B.last, Zid)].append((a, b, wa, wb))
                self.term('T', a, b, B.hist + (Zid,), Zid)
                continue
            if kind != 'C':
                self.term(kind, a, b, B.hist, Zid)
                continue
            Y = self.byid[Zid]
            if Zid in self.processed:
                self.anom.append(('DAG_violation_processed', B.last, Zid))
            if B.last >= 0 and not (Y.cy > self.byid[B.last].cy):
                self.anom.append(('DAG_center_order', B.last, Zid))
            # contact parameter on bot edge
            nm, nu, ce = Y.E[0]
            d = B.d
            nd = nu[0] * d[0] + nu[1] * d[1]
            np0 = nu[0] * B.P0[0] + nu[1] * B.P0[1]
            np1 = nu[0] * B.P1[0] + nu[1] * B.P1[1]
            tS = ((ce - np0) / nd, -np1 / nd)
            Q0 = (B.P0[0] + tS[0] * d[0], B.P0[1] + tS[0] * d[1])
            Q1 = (B.P1[0] + tS[1] * d[0], B.P1[1] + tS[1] * d[1])
            g = (B.g[0] + tS[0], B.g[1] + tS[1])
            m = (Y.c[0] - HALF * Y.n[0], Y.c[1] - HALF * Y.n[1])
            tau = ((Q0[0] - m[0]) * Y.u[0] + (Q0[1] - m[1]) * Y.u[1], Q1[0] * Y.u[0] + Q1[1] * Y.u[1])
            if abs(tau[1]) < 1:
                self.anom.append(('branch_slope_lt_1', float(tau[1])))
            self.cands[Zid].append(dict(xa=a, xb=b, tau=tau, g=g, Q0=Q0, Q1=Q1, hist=B.hist,
                                        last=B.last, wex=B.wex))

    def resolve(self, Y, cl):
        items = []
        for c in cl:
            t0, t1 = c['tau']
            ta, tb = t0 + t1 * c['xa'], t0 + t1 * c['xb']
            lo, hi = (ta, tb) if ta < tb else (tb, ta)
            xf = (-t0 / t1, 1 / t1)
            gf = (c['g'][0] + c['g'][1] * xf[0], c['g'][1] * xf[1])
            items.append([lo, hi, xf, gf, c])
        items.sort(key=lambda it: it[0])
        bps = set()
        for it in items:
            bps.add(it[0]); bps.add(it[1])
        n = len(items)
        for i in range(n):
            for j in range(i + 1, n):
                if items[j][0] >= items[i][1]:
                    break
                lo, hi = max(items[i][0], items[j][0]), min(items[i][1], items[j][1])
                for f, gg in ((items[i][3], items[j][3]), (items[i][2], items[j][2])):
                    r = lin_root(f, gg)
                    if r is not None and lo < r < hi:
                        bps.add(r)
        pts = sorted(bps)
        won = defaultdict(list)
        lost = defaultdict(list)
        start = 0
        active = []
        for s in range(len(pts) - 1):
            a, b = pts[s], pts[s + 1]
            m = (a + b) / 2
            while start < n and items[start][0] < m:
                active.append(start); start += 1
            active = [i for i in active if items[i][1] > m]
            if not active:
                continue
            best = None
            for i in active:
                lo, hi, xf, gf, c = items[i]
                key = (gf[0] + gf[1] * m, xf[0] + xf[1] * m)
                if 'R1_leftmost' in self.variant:
                    key = (xf[0] + xf[1] * m,)
                if best is None or key < best[0]:
                    best = (key, i)
                elif key == best[0]:
                    self.anom.append(('R1_exact_tie', Y.id))
            if len(active) > 1:
                lasts = [items[i][4]['last'] for i in active]
                dirs = set()
                for L_ in lasts:
                    dirs.add((Fr(0), Fr(1)) if L_ < 0 else self.byid[L_].n)
                if len(set(lasts)) < len(lasts) or len(dirs) < len(lasts):
                    self.anom.append(('merge_same_direction', Y.id, lasts))
                self.merge_meas += float(b - a) * (len(active) - 1)
                if Y.tan_abs >= self.Ta:
                    self.merge_at_T += float(b - a) * (len(active) - 1)
            allwin = ('no_R1' in self.variant) or ('T_before_M' in self.variant and Y.tan_abs >= self.Ta)
            for i in active:
                L = won[i] if (i == best[1] or allwin) else lost[i]
                if L and L[-1][1] == a:
                    L[-1] = (L[-1][0], b)
                else:
                    L.append((a, b))
        # record winners per square (tau pieces) for pointwise R1 verification
        for i, segs in won.items():
            for (a, b) in segs:
                self.winners[Y.id].append((a, b, items[i][2], items[i][3]))
        for i, segs in lost.items():
            xf = items[i][2]
            for (a, b) in segs:
                xa_, xb_ = sorted((xf[0] + xf[1] * a, xf[0] + xf[1] * b))
                self.term('M', xa_, xb_, items[i][4]['hist'], Y.id)
        isT = Y.tan_abs >= self.Ta
        for i, segs in won.items():
            c = items[i][4]
            xf = items[i][2]
            hist2 = c['hist'] + (Y.id,)
            for (a, b) in segs:
                xa_, xb_ = sorted((xf[0] + xf[1] * a, xf[0] + xf[1] * b))
                w0, w1 = c['wex']
                wa, wb = sorted((w0 + w1 * xa_, w0 + w1 * xb_))
                self.pairs[(c['last'], Y.id)].append((xa_, xb_, wa, wb))
                self.entries.append((Y.id, xa_, xb_, hist2))
                if isT:
                    self.term('T', xa_, xb_, hist2, Y.id)
                    continue
                P0 = (c['Q0'][0] + Y.n[0], c['Q0'][1] + Y.n[1])
                P1 = c['Q1']
                # H check on exit
                y0, y1 = P0[1], P1[1]
                cuts = [xa_, xb_]
                if y1 != 0:
                    r = (self.h - y0) / y1
                    if xa_ < r < xb_:
                        cuts = [xa_, r, xb_]
                for j in range(len(cuts) - 1):
                    a2, b2 = cuts[j], cuts[j + 1]
                    xm = (a2 + b2) / 2
                    if y0 + y1 * xm >= self.h:
                        self.term('H', a2, b2, hist2, Y.id)
                    else:
                        self.advance(Beam(a2, b2, P0, P1, Y.n, Y.id, c['g'], c['tau'], hist2))

    def run(self):
        self.cands = defaultdict(list)
        self.terms = []
        self.pairs = defaultdict(list)
        self.entries = []
        self.winners = defaultdict(list)
        self.processed = set()
        self.merge_meas = 0.0     # tau-measure (on bot edges) of merge overlaps
        self.merge_at_T = 0.0     # ... at tilt terminators
        B = Beam(Fr(0), self.k, (Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(0), Fr(1)), -1,
                 (Fr(0), Fr(0)), (Fr(0), Fr(1)), ())
        self.advance(B)
        for Y in sorted(self.sqs, key=lambda S: (S.cy, S.id)):
            self.processed.add(Y.id)
            cl = self.cands.pop(Y.id, None)
            if cl:
                self.resolve(Y, cl)
        if self.cands:
            self.anom.append(('unprocessed_cands', list(self.cands.keys())))
        tot = sum(b - a for (_, a, b, _, _) in self.terms)
        self.conservation = (tot == self.k)
        return self


# ---------------- pointwise tracer (no R1), exact ----------------
def trace_point(sqs, x, k, delta, Ta, h):
    """returns (kind, entries list [(Yid, q, g)], end info). Ignores R1 (no merging)."""
    p = (Fr(x), Fr(0))
    d = (Fr(0), Fr(1))
    last = -1
    g = Fr(0)
    ents = []
    byid = {S.id: S for S in sqs}
    for _ in range(10000):
        best = None
        for Z in sqs:
            if Z.id == last:
                continue
            r = hit(Z, p, d)
            if r is None:
                continue
            if best is None or r[0] < best[0]:
                best = (r[0], Z.id, r[1], r[2])
        tD = delta - g
        tH = (h - p[1]) / d[1]
        tW = None
        if d[0] > 0:
            tW = (k - p[0]) / d[0]
        elif d[0] < 0:
            tW = -p[0] / d[0]
        cand = [tD, tH] + ([tW] if tW is not None else []) + ([best[0]] if best else [])
        ts = min(cand)
        if tD == ts:
            return ('D', ents, None)
        if tW is not None and tW == ts:
            return ('W', ents, None)
        if best is not None and best[0] == ts:
            Y = byid[best[1]]
            if best[2] != 'bot' or best[3]:
                return ('E', ents, best[1])
            q = (p[0] + ts * d[0], p[1] + ts * d[1])
            g = g + ts
            ents.append((Y.id, q, g))
            if Y.tan_abs >= Ta:
                return ('T', ents, Y.id)
            p = (q[0] + Y.n[0], q[1] + Y.n[1])
            d = Y.n
            last = Y.id
            if p[1] >= h:
                return ('H', ents, Y.id)
            continue
        return ('H', ents, None)
    raise RuntimeError('loop')


# ---------------- theta ----------------
def theta(phi1, phi2):
    q = mp.pi / 2
    dd = abs(phi1 - phi2) % q
    return min(dd, q - dd)
