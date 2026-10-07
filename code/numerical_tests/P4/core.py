"""Exact (rational) implementation of the path flows of Section 3, for the tests of Section 6
(Theorem 6.9, P4*: the shadow inequality).

Written from scratch for these tests; it does not import code from the other test folders.

Two independent tracers:
  * Tracer      : beam tracer.  x-intervals with affine data; ray steps done in the source frame
                  (lower envelope of each square = front edges); R1 resolved per square in centre-height
                  order by a kinetic lexicographic sweep in the bottom-side coordinate sigma.
  * point_trace : one floor point at a time; first contact by the slab method on all 4 edges;
                  R1 resolved by an explicit backward search of all live arrivals (arrivals()).
"""
from fractions import Fraction as Fr
import ctypes
import os

ZERO, ONE, HALF = Fr(0), Fr(1), Fr(1, 2)
PRIO = {'D': 0, 'W': 1, 'C': 2, 'H': 3}


# ---------------------------------------------------------------- output folder
def out_path(name):
    """path of an output file in the folder out/ next to the scripts (created on first use)"""
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)


# ---------------------------------------------------------------- memory self-check
class _PMC(ctypes.Structure):
    _fields_ = [('cb', ctypes.c_ulong), ('PageFaultCount', ctypes.c_ulong),
                ('PeakWorkingSetSize', ctypes.c_size_t), ('WorkingSetSize', ctypes.c_size_t),
                ('QuotaPeakPagedPoolUsage', ctypes.c_size_t), ('QuotaPagedPoolUsage', ctypes.c_size_t),
                ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t), ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                ('PagefileUsage', ctypes.c_size_t), ('PeakPagefileUsage', ctypes.c_size_t)]


def rss_mb():
    try:
        k32 = ctypes.windll.kernel32
        ps = ctypes.windll.psapi
        k32.GetCurrentProcess.restype = ctypes.c_void_p
        ps.GetProcessMemoryInfo.argtypes = [ctypes.c_void_p, ctypes.POINTER(_PMC), ctypes.c_ulong]
        p = _PMC()
        p.cb = ctypes.sizeof(_PMC)
        ps.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(p), p.cb)
        return p.WorkingSetSize / 2 ** 20
    except Exception:
        return -1.0


# ---------------------------------------------------------------- geometry
def norm_t(t):
    """tan(phi/2) normalised so that phi in (-pi/4, pi/4]  (rational t never hits +-pi/4 exactly)."""
    t = Fr(t)
    for _ in range(16):
        if t > 0 and (t + 1) ** 2 > 2:          # phi > pi/4  -> phi - pi/2
            t = (t - 1) / (1 + t)
        elif t < 0 and (1 - t) ** 2 >= 2:       # phi <= -pi/4 -> phi + pi/2
            t = (t + 1) / (1 - t)
        else:
            return t
    raise ValueError('norm_t')


class Sq:
    """closed unit square; phase phi = 2 atan t in (-pi/4, pi/4]; u=(cos,sin), n=(-sin,cos)."""
    __slots__ = ('idx', 'c', 't', 'cs', 'sn', 'u', 'n', 'V', 'E', 'xmin', 'xmax', 'ymin', 'ymax')

    def __init__(self, cx, cy, t, idx=-1):
        t = norm_t(t)
        self.idx = idx
        self.t = t
        d = 1 + t * t
        cs, sn = (1 - t * t) / d, 2 * t / d
        self.cs, self.sn = cs, sn
        cx, cy = Fr(cx), Fr(cy)
        self.c = (cx, cy)
        u = (cs, sn)
        n = (-sn, cs)
        self.u, self.n = u, n

        def P(a, b):
            return (cx + a * u[0] + b * n[0], cy + a * u[1] + b * n[1])
        BL, BR, TR, TL = P(-HALF, -HALF), P(HALF, -HALF), P(HALF, HALF), P(-HALF, HALF)
        self.V = (BL, BR, TR, TL)
        # (A, B, outward normal, name); |B-A| = 1 exactly
        self.E = ((BL, BR, (-n[0], -n[1]), 'bottom'),
                  (BR, TR, u, 'right'),
                  (TL, TR, n, 'top'),
                  (BL, TL, (-u[0], -u[1]), 'left'))
        xs = [v[0] for v in self.V]
        ys = [v[1] for v in self.V]
        self.xmin, self.xmax, self.ymin, self.ymax = min(xs), max(xs), min(ys), max(ys)

    def chord(self, y):
        """closed horizontal chord [xl, xr] of the square at height y, or None."""
        if y < self.ymin or y > self.ymax:
            return None
        pts = []
        for (A, B, m, nm) in self.E:
            ya, yb = A[1], B[1]
            if ya == yb:
                if ya == y:
                    pts += [A[0], B[0]]
                continue
            lo, hi = (ya, yb) if ya < yb else (yb, ya)
            if lo <= y <= hi:
                pts.append(A[0] + (y - ya) * (B[0] - A[0]) / (yb - ya))
        return (min(pts), max(pts))

    def vchord(self, x):
        """closed vertical chord [yl, yh] at abscissa x, or None."""
        if x < self.xmin or x > self.xmax:
            return None
        pts = []
        for (A, B, m, nm) in self.E:
            xa, xb = A[0], B[0]
            if xa == xb:
                if xa == x:
                    pts += [A[1], B[1]]
                continue
            lo, hi = (xa, xb) if xa < xb else (xb, xa)
            if lo <= x <= hi:
                pts.append(A[1] + (x - xa) * (B[1] - A[1]) / (xb - xa))
        return (min(pts), max(pts))


def disjoint(A, B):
    if A.xmax < B.xmin or B.xmax < A.xmin or A.ymax < B.ymin or B.ymax < A.ymin:
        return True
    for ax in (A.u, A.n, B.u, B.n):
        pa = [ax[0] * v[0] + ax[1] * v[1] for v in A.V]
        pb = [ax[0] * v[0] + ax[1] * v[1] for v in B.V]
        if max(pa) < min(pb) or max(pb) < min(pa):
            return True
    return False


def validate(k, sqs):
    k = Fr(k)
    for S in sqs:
        if S.xmin < 0 or S.ymin < 0 or S.xmax > k or S.ymax > k:
            return False, 'outside container: square %d' % S.idx
    for i in range(len(sqs)):
        for j in range(i + 1, len(sqs)):
            if not disjoint(sqs[i], sqs[j]):
                return False, 'overlap %d %d' % (i, j)
    return True, ''


class Cfg:
    def __init__(self, k, squares, delta, tF, hF, name=''):
        self.k = Fr(k)
        self.sq = [Sq(cx, cy, t, i) for i, (cx, cy, t) in enumerate(squares)]
        self.delta = Fr(delta)
        self.tF = Fr(tF)
        self.hF = Fr(hF)
        self.name = name


# ---------------------------------------------------------------- affine helpers (f = (c0, c1): c0 + c1 x)
def ev(f, x):
    return f[0] + f[1] * x


def where(f, a, b, v, op):
    """open sub-interval of (a,b) on which f(x) op v (op in < > <= >=), up to endpoints; None if empty."""
    c0 = f[0] - v
    c1 = f[1]
    if c1 == 0:
        ok = {'<': c0 < 0, '>': c0 > 0, '<=': c0 <= 0, '>=': c0 >= 0}[op]
        return (a, b) if ok else None
    r = -c0 / c1
    if op in ('<', '<='):
        lo, hi = (a, min(b, r)) if c1 > 0 else (max(a, r), b)
    else:
        lo, hi = (max(a, r), b) if c1 > 0 else (a, min(b, r))
    return (lo, hi) if lo < hi else None


def pt_plus(P, t, d):
    """affine point P(x) + t(x) * d  (d constant vector)."""
    return ((P[0][0] + t[0] * d[0], P[0][1] + t[1] * d[0]),
            (P[1][0] + t[0] * d[1], P[1][1] + t[1] * d[1]))


def dotc(v, P):
    return (v[0] * P[0][0] + v[1] * P[1][0], v[0] * P[0][1] + v[1] * P[1][1])


class Piece:
    __slots__ = ('lo', 'hi', 'kind', 'tau', 'hist', '_segs', '_evp')

    def __init__(self, lo, hi, kind, tau, hist):
        self.lo, self.hi, self.kind, self.tau, self.hist = lo, hi, kind, tau, hist
        self._segs = None
        self._evp = None


def envelope(cands, a, b):
    """lower envelope on (a,b) of affine candidates (kind, info, (c0,c1)); ties by PRIO."""
    out = []
    cur = a
    while cur < b:
        best, bkey = None, None
        for c in cands:
            key = (c[2][0] + c[2][1] * cur, c[2][1], PRIO[c[0]])
            if bkey is None or key < bkey:
                best, bkey = c, key
        nxt = b
        bs = best[2][1]
        for c in cands:
            if c is best:
                continue
            cs = c[2][1]
            if cs < bs:
                xc = (c[2][0] - best[2][0]) / (bs - cs)
                if cur < xc < nxt:
                    nxt = xc
        if out and out[-1][2] is best and out[-1][1] == cur:
            out[-1] = (out[-1][0], nxt, best)
        else:
            out.append((cur, nxt, best))
        cur = nxt
    return out


# ---------------------------------------------------------------- beam tracer (master flow F^0)
class Tracer:
    def __init__(self, cfg, r1=True):
        self.cfg = cfg
        self.r1 = r1
        n = len(cfg.sq)
        self.order = sorted(range(n), key=lambda j: (cfg.sq[j].c[1], j))
        self.pos = {j: i for i, j in enumerate(self.order)}
        self.pending = [[] for _ in range(n)]
        self.processed = [False] * n
        self.final = []
        self.anom = []
        self.minlam = None
        self.minmu = None
        self.merge_meas = ZERO
        self.nentry_multi = 0

    def _anom(self, msg):
        if len(self.anom) < 50:
            self.anom.append(msg)

    def run(self):
        k = self.cfg.k
        self.ray_step(-1, ZERO, k, ((ZERO, ONE), (ZERO, ZERO)), (ZERO, ZERO), [])
        for j in self.order:
            self.processed[j] = True
            arr = self.pending[j]
            self.pending[j] = None
            if arr:
                self.resolve(j, arr)
        return self

    def ray_step(self, src, lo, hi, P, G, hist):
        cfg = self.cfg
        if src < 0:
            u, d = (ONE, ZERO), (ZERO, ONE)
        else:
            S = cfg.sq[src]
            u, d = S.u, S.n
        s_aff = dotc(u, P)
        h_aff = dotc(d, P)
        if h_aff[1] != 0:
            self._anom('start points not on a line orthogonal to d (src %d)' % src)
        h0 = h_aff[0]
        lam = s_aff[1]
        al = abs(lam)
        if self.minlam is None or al < self.minlam:
            self.minlam = al
        if al < 1:
            self._anom('|lambda|<1 at src %d: %s' % (src, lam))
        s1, s2 = ev(s_aff, lo), ev(s_aff, hi)
        s_lo, s_hi = min(s1, s2), max(s1, s2)
        edges = []
        for Y in cfg.sq:
            if Y.idx == src:
                continue
            ss = [u[0] * v[0] + u[1] * v[1] for v in Y.V]
            if max(ss) <= s_lo or min(ss) >= s_hi:
                continue
            for (A, B, m, nm) in Y.E:
                md = m[0] * d[0] + m[1] * d[1]
                if md >= 0:
                    continue
                sa = u[0] * A[0] + u[1] * A[1]
                sb = u[0] * B[0] + u[1] * B[1]
                ha = d[0] * A[0] + d[1] * A[1]
                hb = d[0] * B[0] + d[1] * B[1]
                ke = (hb - ha) / (sb - sa)
                t_x = (ha + (s_aff[0] - sa) * ke - h0, lam * ke)
                xa, xb = (sa - s_aff[0]) / lam, (sb - s_aff[0]) / lam
                xl, xr = max(min(xa, xb), lo), min(max(xa, xb), hi)
                if xl >= xr:
                    continue
                tm = ev(t_x, (xl + xr) / 2)
                if tm < 0:
                    if ev(t_x, xl) > 0 or ev(t_x, xr) > 0:
                        self._anom('square %d straddles start line of src %d' % (Y.idx, src))
                    continue
                if src >= 0 and (ev(t_x, xl) <= 0 or ev(t_x, xr) <= 0):
                    self._anom('square %d touches start line of src %d' % (Y.idx, src))
                edges.append((xl, xr, t_x, Y.idx, nm))
        c0 = [('D', None, (cfg.delta - G[0], -G[1]))]
        if d[0] < 0:
            c0.append(('W', None, (-P[0][0] / d[0], -P[0][1] / d[0])))
        elif d[0] > 0:
            c0.append(('W', None, ((cfg.k - P[0][0]) / d[0], -P[0][1] / d[0])))
        c0.append(('H', None, ((cfg.hF - P[1][0]) / d[1], -P[1][1] / d[1])))
        bps = sorted(set([lo, hi] + [e[0] for e in edges] + [e[1] for e in edges]))
        for i in range(len(bps) - 1):
            a, b = bps[i], bps[i + 1]
            act = list(c0) + [('C', (e[3], e[4]), e[2]) for e in edges if e[0] <= a and e[1] >= b]
            for (a2, b2, w) in envelope(act, a, b):
                self.emit(src, a2, b2, P, G, d, hist, w)

    def emit(self, src, a, b, P, G, d, hist, w):
        kind, info, t = w
        E = pt_plus(P, t, d)
        xm = (a + b) / 2
        if ev(t, xm) < 0:
            self._anom('negative step length')
        if kind == 'D':
            self.final.append(Piece(a, b, 'D', E[1], hist + [('g', src, P, t, d, 'D', None)]))
        elif kind == 'W':
            wx = ev(E[0], xm)
            if wx != 0 and wx != self.cfg.k:
                self._anom('W point not on wall')
            self.final.append(Piece(a, b, 'W', E[1], hist + [('g', src, P, t, d, 'W', None)]))
        elif kind == 'H':
            self.final.append(Piece(a, b, 'H', (self.cfg.hF, ZERO), hist + [('g', src, P, t, d, 'H', None)]))
        else:
            j, nm = info
            g2 = (G[0] + t[0], G[1] + t[1])
            if ev(g2, xm) >= self.cfg.delta:
                self._anom('contact with g >= delta')
            if nm == 'bottom':
                if self.processed[j]:
                    self._anom('DAG violation: entry into already processed square %d from src %d' % (j, src))
                    self.final.append(Piece(a, b, 'X', E[1], hist + [('g', src, P, t, d, 'B', j)]))
                    return
                self.pending[j].append((a, b, E, g2, hist + [('g', src, P, t, d, 'B', j)]))
            else:
                self.final.append(Piece(a, b, 'E', E[1], hist + [('g', src, P, t, d, 'E', j)]))

    def resolve(self, j, arr):
        Y = self.cfg.sq[j]
        recs = []
        for (a, b, q, g, h) in arr:
            sig = (Y.u[0] * (q[0][0] - Y.c[0]) + Y.u[1] * (q[1][0] - Y.c[1]),
                   Y.u[0] * q[0][1] + Y.u[1] * q[1][1])
            s1 = sig[1]
            am = abs(s1)
            if self.minmu is None or am < self.minmu:
                self.minmu = am
            if am < 1:
                self._anom('|mu|<1 at square %d' % j)
            sa, sb = ev(sig, a), ev(sig, b)
            if sa > sb:
                sa, sb = sb, sa
            if sa < -HALF or sb > HALF:
                self._anom('entry outside bottom side of %d' % j)
            xs = (-sig[0] / s1, 1 / s1)
            gs = (g[0] + g[1] * xs[0], g[1] * xs[1])
            recs.append((sa, sb, xs, gs, (a, b, q, g, h)))
        if (not self.r1) or len(recs) == 1:
            for r in recs:
                self.winner(j, r[4][0], r[4][1], r[4])
            return
        bps = sorted(set([r[0] for r in recs] + [r[1] for r in recs]))
        for i in range(len(bps) - 1):
            A_, B_ = bps[i], bps[i + 1]
            cov = [r for r in recs if r[0] <= A_ and r[1] >= B_]
            if not cov:
                continue
            if len(cov) == 1:
                self._sig_part(j, cov[0], A_, B_, True)
                continue
            self.nentry_multi += 1
            cur = A_
            while cur < B_:
                cc = cur
                best = min(cov, key=lambda r: (ev(r[3], cc), r[3][1], ev(r[2], cc), r[2][1]))
                nxt = B_
                for r in cov:
                    if r is best:
                        continue
                    if r[3][1] < best[3][1]:
                        xc = (r[3][0] - best[3][0]) / (best[3][1] - r[3][1])
                        if cur < xc < nxt:
                            nxt = xc
                    elif r[3] == best[3]:
                        # identical g: x-lines must not cross inside (same floor point twice)
                        if r[2][1] != best[2][1]:
                            xc = (r[2][0] - best[2][0]) / (best[2][1] - r[2][1])
                            if cur < xc < nxt:
                                self._anom('x-lines cross inside an R1 cell at %d' % j)
                for r in cov:
                    self._sig_part(j, r, cur, nxt, r is best)
                cur = nxt

    def _sig_part(self, j, r, sa, sb, win):
        xa, xb = ev(r[2], sa), ev(r[2], sb)
        if xa > xb:
            xa, xb = xb, xa
        a0, b0, q, g, h = r[4]
        if xa < a0 or xb > b0:
            self._anom('sigma part maps outside arrival piece')
        if win:
            self.winner(j, xa, xb, (xa, xb, q, g, h))
        else:
            self.merge_meas += xb - xa
            self.final.append(Piece(xa, xb, 'M', q[1], h + [('e', j, q, 'M')]))

    def winner(self, j, a, b, rec):
        cfg = self.cfg
        Y = cfg.sq[j]
        _, _, q, g, h = rec
        if abs(Y.t) >= cfg.tF:
            self.final.append(Piece(a, b, 'T', q[1], h + [('e', j, q, 'T')]))
            return
        ex = ((q[0][0] + Y.n[0], q[0][1]), (q[1][0] + Y.n[1], q[1][1]))
        h2 = h + [('e', j, q, 'P'), ('p', j, q, ex)]
        r = where(ex[1], a, b, cfg.hF, '>=')
        if r:
            self.final.append(Piece(r[0], r[1], 'H', (cfg.hF, ZERO), h2))
        r2 = where(ex[1], a, b, cfg.hF, '<')
        if r2:
            self.ray_step(j, r2[0], r2[1], ex, g, h2)


def check_partition(pieces, k):
    ps = sorted(pieces, key=lambda p: p.lo)
    tot = ZERO
    cur = ZERO
    bad = 0
    for p in ps:
        if p.lo != cur:
            bad += 1
        if p.hi <= p.lo:
            bad += 1
        tot += p.hi - p.lo
        cur = p.hi
    if cur != k:
        bad += 1
    return bad == 0 and tot == k, tot


# ---------------------------------------------------------------- scale flow F_s by truncation of F^0
def truncate(cfg, pieces, ts, hs):
    out = []
    for pc in pieces:
        _trunc(cfg, pc.lo, pc.hi, pc.hist, ts, hs, out)
    return out


def _trunc(cfg, lo, hi, hist, ts, hs, out):
    nh = []
    i = 0
    while i < len(hist):
        st = hist[i]
        if st[0] == 'g':
            _, src, P, t, d, ev_, Y = st
            Ey = (P[1][0] + t[0] * d[1], P[1][1] + t[1] * d[1])
            r = where(Ey, lo, hi, hs, '>')
            if r:
                tH = ((hs - P[1][0]) / d[1], -P[1][1] / d[1])
                out.append(Piece(r[0], r[1], 'H', (hs, ZERO), nh + [('g', src, P, tH, d, 'Hs', None)]))
            r2 = where(Ey, lo, hi, hs, '<=')
            if not r2:
                return
            lo, hi = r2
            nh = nh + [st]
            if ev_ in ('D', 'W', 'E', 'H', 'X'):
                out.append(Piece(lo, hi, ev_, (hs, ZERO) if ev_ == 'H' else Ey, nh))
                return
            i += 1
            continue
        if st[0] == 'e':
            _, j, q, res = st
            if res == 'M':
                out.append(Piece(lo, hi, 'M', q[1], nh + [st]))
                return
            if abs(cfg.sq[j].t) >= ts:
                out.append(Piece(lo, hi, 'T', q[1], nh + [('e', j, q, 'Ts')]))
                return
            if res != 'P':
                raise RuntimeError('truncate: winner with |t|<ts but F0 says %s' % res)
            nh = nh + [st]
            i += 1
            continue
        if st[0] == 'p':
            _, j, q, ex = st
            r = where(ex[1], lo, hi, hs, '>=')
            if r:
                out.append(Piece(r[0], r[1], 'H', (hs, ZERO), nh + [st]))
            r2 = where(ex[1], lo, hi, hs, '<')
            if not r2:
                return
            lo, hi = r2
            nh = nh + [st]
            i += 1
            continue
    raise RuntimeError('truncate fell off the history')


# ---------------------------------------------------------------- per-piece geometry caches
def piece_segments(pc, cfg):
    if pc._segs is None:
        segs = []
        for st in pc.hist:
            if st[0] == 'g':
                _, src, P, t, d, ev_, Y = st
                E = pt_plus(P, t, d)
                segs.append(('g', src, P, E, d[0] / d[1]))
            elif st[0] == 'p':
                _, j, q, ex = st
                n = cfg.sq[j].n
                segs.append(('p', j, q, ex, n[0] / n[1]))
        pc._segs = segs
    return pc._segs


def piece_events(pc):
    """distinct event points along the path: (affine point, label)."""
    if pc._evp is None:
        pts = []
        first = True
        for st in pc.hist:
            if st[0] == 'g':
                _, src, P, t, d, ev_, Y = st
                if first:
                    pts.append((P, ('floor',)))
                    first = False
                E = pt_plus(P, t, d)
                lab = {'D': ('D',), 'W': ('W',), 'H': ('H',), 'Hs': ('H',), 'B': ('bot', Y),
                       'E': ('E', Y)}[ev_]
                pts.append((E, lab))
            elif st[0] == 'p':
                _, j, q, ex = st
                pts.append((ex, ('top', j)))
        pc._evp = pts
    return pc._evp


def exc_set(pieces, hF):
    """heights in (0,hF) at which some event-point height is constant on some piece -> {y: set(labels)}."""
    out = {}
    for pc in pieces:
        for (pt, lab) in piece_events(pc):
            hy = pt[1]
            if hy[1] == 0 and 0 < hy[0] < hF:
                out.setdefault(hy[0], set()).add(lab)
    return out


# ---------------------------------------------------------------- horizontal-line bookkeeping
def analyse(cfg, pieces, y):
    """exact bookkeeping on the line l_y for one flow.  Returns dict of Fractions + anomaly list."""
    k = cfg.k
    L = {'T': ZERO, 'D': ZERO, 'E': ZERO, 'M': ZERO, 'W': ZERO, 'X': ZERO}
    Btau = ZERO
    Bev = ZERO
    mu = {}
    Fgap = ZERO
    gimgs = []
    pimgs = {}
    an = []
    min_gslope = None
    min_pslope = None
    for pc in pieces:
        lo, hi = pc.lo, pc.hi
        if pc.kind == 'H':
            A = (lo, hi)
        else:
            I = where(pc.tau, lo, hi, y, '<')
            if I:
                L[pc.kind] += I[1] - I[0]
            if pc.tau[1] == 0 and pc.tau[0] == y:
                Btau += hi - lo
                continue
            A = where(pc.tau, lo, hi, y, '>')
        if not A:
            continue
        alo, ahi = A
        alive = ahi - alo
        cov = ZERO
        for (typ, Z, S, E, r) in piece_segments(pc, cfg):
            J = where(S[1], alo, ahi, y, '<')
            if not J:
                continue
            J = where(E[1], J[0], J[1], y, '>')
            if not J:
                continue
            m = J[1] - J[0]
            cov += m
            X = (S[0][0] + (y - S[1][0]) * r, S[0][1] - S[1][1] * r)
            x1, x2 = ev(X, J[0]), ev(X, J[1])
            il, ih = (x1, x2) if x1 < x2 else (x2, x1)
            sl = abs(X[1])
            if typ == 'g':
                Fgap += m
                gimgs.append((il, ih, 1 / sl))
                if min_gslope is None or sl < min_gslope:
                    min_gslope = sl
            else:
                mu[Z] = mu.get(Z, ZERO) + m
                pimgs.setdefault(Z, []).append((il, ih))
                if min_pslope is None or sl < min_pslope:
                    min_pslope = sl
        nconst = 0
        for (pt, lab) in piece_events(pc):
            if pt[1][1] == 0 and pt[1][0] == y:
                nconst += 1
        if nconst > 1:
            an.append('two constant event points at y on one piece')
        if nconst == 1:
            Bev += alive
        if cov + (alive if nconst == 1 else ZERO) != alive:
            an.append('cover mismatch on piece [%s,%s]: cov=%s alive=%s nconst=%d'
                      % (float(lo), float(hi), float(cov), float(alive), nconst))
    # chords
    chords = {}
    sumc = ZERO
    for Z in cfg.sq:
        ch = Z.chord(y)
        if ch is not None:
            chords[Z.idx] = ch
            sumc += ch[1] - ch[0]
    # sweep: density of live gap crossings vs waste
    drho, dch = {}, {}
    for (il, ih, dn) in gimgs:
        drho[il] = drho.get(il, ZERO) + dn
        drho[ih] = drho.get(ih, ZERO) - dn
        if il < 0 or ih > k:
            an.append('gap image outside (0,k)')
    for (xl, xr) in chords.values():
        dch[xl] = dch.get(xl, 0) + 1
        dch[xr] = dch.get(xr, 0) - 1
    Ps = sorted(set([ZERO, k]) | set(drho) | set(dch))
    rho = ZERO
    chc = 0
    Ovg = U = Wm = intr = ZERO
    in_chord = ZERO
    maxrho = ZERO
    for i in range(len(Ps) - 1):
        p = Ps[i]
        rho += drho.get(p, ZERO)
        chc += dch.get(p, 0)
        ln = Ps[i + 1] - p
        if p < 0 or Ps[i + 1] > k:
            continue
        if chc > 0:
            if rho != 0:
                in_chord += rho * ln
        else:
            Wm += ln
            intr += rho * ln
            if rho > maxrho:
                maxrho = rho
            if rho > 1:
                Ovg += (rho - 1) * ln
            else:
                U += (1 - rho) * ln
    if in_chord != 0:
        an.append('gap crossings inside a chord interior: %s' % float(in_chord))
    # pass images: inside own chord, pairwise disjoint (R1 injectivity), slope >= 1
    overlap_pass = ZERO
    outside_chord = 0
    for Z, lst in pimgs.items():
        ch = chords.get(Z)
        lst.sort()
        for (il, ih) in lst:
            if ch is None or il < ch[0] or ih > ch[1]:
                outside_chord += 1
        for i in range(len(lst) - 1):
            if lst[i][1] > lst[i + 1][0]:
                overlap_pass += min(lst[i][1], lst[i + 1][1]) - lst[i + 1][0]
    if outside_chord:
        an.append('pass image outside its chord (%d)' % outside_chord)
    Sh = Ov = ZERO
    for Z in cfg.sq:
        c = (chords[Z.idx][1] - chords[Z.idx][0]) if Z.idx in chords else ZERO
        m = mu.get(Z.idx, ZERO)
        if c > m:
            Sh += c - m
        else:
            Ov += m - c
    Ltot = sum(L.values())
    B = Btau + Bev
    LT = L['T']
    Lam = L['D'] + L['E'] + L['M'] + L['W'] + L['X'] + Ovg
    res = dict(y=y, L=L, Ltot=Ltot, Btau=Btau, Bev=Bev, B=B, summu=sum(mu.values(), ZERO), Fgap=Fgap,
               sumc=sumc, omega=Wm, Ovg=Ovg, U=U, intrho=intr, Sh=Sh, Ov=Ov, LT=LT, Lam=Lam,
               margin=LT + Lam - Sh, maxrho=maxrho, overlap_pass=overlap_pass,
               min_gslope=min_gslope, min_pslope=min_pslope, an=an)
    return res


def identity_checks(r, k):
    """exact identities that must hold at every y (general form with B)."""
    bad = []
    if r['Ltot'] + r['summu'] + r['Fgap'] + r['B'] != k:
        bad.append('id1(k=L+sum mu+Fgap+B)')
    if r['sumc'] + r['omega'] != k:
        bad.append('id2(k=sum c+omega)')
    if r['intrho'] != r['Fgap']:
        bad.append('int rho != Fgap')
    if r['Sh'] - (r['Ltot'] + r['Fgap'] - r['omega'] + r['Ov'] + r['B']) != 0:
        bad.append('general shadow identity')
    if r['Fgap'] - r['omega'] != r['Ovg'] - r['U']:
        bad.append('(e) Fgap-omega = Ovg-U')
    if r['margin'] != r['U'] - r['Ov'] - r['B']:
        bad.append('margin != U - Ov - B')
    if r['B'] < 0 or r['U'] < 0:
        bad.append('negative B or U')
    return bad


# ---------------------------------------------------------------- independent point tracer
def first_hit(cfg, p, d, exclude=None):
    """slab method on all squares: first contact of ray p + t d (t >= 0).
    returns (t, j, edge-name or 'vertex', point) or None."""
    best = None
    for Y in cfg.sq:
        if Y.idx == exclude:
            continue
        tin = None
        tout = None
        arg = []
        miss = False
        for (A, B, m, nm) in Y.E:
            md = m[0] * d[0] + m[1] * d[1]
            mp = m[0] * (p[0] - A[0]) + m[1] * (p[1] - A[1])
            if md == 0:
                if mp > 0:
                    miss = True
                    break
                continue
            tb = -mp / md
            if md < 0:
                if tin is None or tb > tin:
                    tin, arg = tb, [nm]
                elif tb == tin:
                    arg.append(nm)
            else:
                if tout is None or tb < tout:
                    tout = tb
        if miss or tin is None or tout is None:
            continue
        if tin > tout or tout < 0:
            continue
        if tin < 0:
            raise RuntimeError('ray starts inside square %d' % Y.idx)
        if best is None or tin < best[0]:
            best = (tin, Y.idx, arg, tout)
    if best is None:
        return None
    t, j, arg, tout = best
    Y = cfg.sq[j]
    z = (p[0] + t * d[0], p[1] + t * d[1])
    if len(arg) == 1 and t < tout:
        nm = arg[0]
        for (A, B, m, n2) in Y.E:
            if n2 == nm:
                w = (z[0] - A[0]) * (B[0] - A[0]) + (z[1] - A[1]) * (B[1] - A[1])
                if 0 < w < 1:
                    return (t, j, nm, z)
    return (t, j, 'vertex', z)


def arrivals(cfg, j, q, memo):
    """all (g, x) of floor points whose master-flow path reaches the entry candidate (j, q) alive.
    Backward search over possible last sources; R1 at earlier squares resolved recursively."""
    key = (j, q)
    if key in memo:
        return memo[key]
    res = []
    if q[1] <= cfg.hF:
        if q[1] < cfg.delta and 0 < q[0] < cfg.k:
            ok = True
            if q[1] > 0:
                h = first_hit(cfg, q, (ZERO, -ONE), exclude=j)
                ok = (h is None) or (h[0] > q[1])
            if ok:
                res.append((q[1], q[0]))
        for X in cfg.sq:
            if X.idx == j or abs(X.t) >= cfg.tF:
                continue
            d = X.n
            if (q[0] - X.c[0]) * d[0] + (q[1] - X.c[1]) * d[1] <= HALF:
                continue
            h = first_hit(cfg, q, (-d[0], -d[1]), exclude=j)
            if h is None or h[1] != X.idx or h[2] != 'top':
                continue
            tb = h[0]
            if tb <= 0 or tb >= cfg.delta:
                continue
            e = h[3]
            if e[1] >= cfg.hF:
                continue
            qX = (e[0] - d[0], e[1] - d[1])
            arrX = arrivals(cfg, X.idx, qX, memo)
            if not arrX:
                continue
            gw, xw = min(arrX)
            if gw + tb < cfg.delta:
                res.append((gw + tb, xw))
    memo[key] = res
    return res


def point_trace(cfg, x, memo, r1=True):
    """returns (kind, tau, tuple of passed squares)."""
    x = Fr(x)
    k, delta, hF, tF = cfg.k, cfg.delta, cfg.hF, cfg.tF
    p = (x, ZERO)
    d = (ZERO, ONE)
    src = None
    g = ZERO
    seq = []
    for _ in range(100000):
        hit = first_hit(cfg, p, d, exclude=src)
        tS = hit[0] if hit else None
        tW = None
        if d[0] < 0:
            tW = -p[0] / d[0]
        elif d[0] > 0:
            tW = (k - p[0]) / d[0]
        tD = delta - g
        tH = (hF - p[1]) / d[1]
        ts = min(v for v in (tS, tW, tD, tH) if v is not None)
        if tD == ts:
            return ('D', p[1] + tD * d[1], tuple(seq))
        if tW is not None and tW == ts:
            return ('W', p[1] + tW * d[1], tuple(seq))
        if tS is not None and tS == ts:
            t, j, nm, q = hit
            g = g + t
            if nm != 'bottom':
                return ('E', q[1], tuple(seq))
            Y = cfg.sq[j]
            if r1:
                arr = arrivals(cfg, j, q, memo)
                if (g, x) not in arr:
                    return ('ANOM', q[1], tuple(seq))
                if min(arr) != (g, x):
                    return ('M', q[1], tuple(seq))
            if abs(Y.t) >= tF:
                return ('T', q[1], tuple(seq))
            seq.append(j)
            ex = (q[0] + Y.n[0], q[1] + Y.n[1])
            if ex[1] >= hF:
                return ('H', hF, tuple(seq))
            p, d, src = ex, Y.n, j
            continue
        return ('H', hF, tuple(seq))
    raise RuntimeError('point_trace loop')


def piece_passed(pc):
    return tuple(st[1] for st in pc.hist if st[0] == 'p')
