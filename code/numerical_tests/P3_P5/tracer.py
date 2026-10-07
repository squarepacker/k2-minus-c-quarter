"""Exact beam tracer of the flow F^0 of Section 3, for the tests of Sections 5 (P3) and 7 (P5).
Written from scratch from the definitions; it does not import code from the other test folders.

All squares have rational-parametrized phases: s rational, cos = (1-s^2)/(1+s^2), sin = 2s/(1+s^2),
phi = 2 atan(s) in (-pi/4, pi/4]  <=>  s in (-tan(pi/8), tan(pi/8)].
Centers, delta, k rational -> every computation is exact (fractions.Fraction).

Master flow F^0 (bottom flow): priority D > W > contact > H; contact on relint bot(Y) -> entry candidate,
other contact -> E. R1: squares processed by center height; at each entry point the live arrival with
lexicographically smallest (g, x) wins, others M. R2: winner with |s_Y| >= sF (i.e. a(Y) >= alpha_F) -> T.
"""
from fractions import Fraction as Fr
import math, time, os

MUT = {'no_r1': False, 'contact_before_death': False}   # mutation switches (tests of the checker only)


def out_path(name):
    """path of an output file in the folder out/ next to the scripts (created on first use)"""
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)

class TimeUp(Exception):
    pass

class Aff:
    __slots__ = ('c0', 'c1')
    def __init__(self, c0, c1=Fr(0)):
        self.c0 = Fr(c0); self.c1 = Fr(c1)
    def __call__(self, x): return self.c0 + self.c1 * x
    def __add__(self, o):
        if isinstance(o, Aff): return Aff(self.c0 + o.c0, self.c1 + o.c1)
        return Aff(self.c0 + o, self.c1)
    def __sub__(self, o):
        if isinstance(o, Aff): return Aff(self.c0 - o.c0, self.c1 - o.c1)
        return Aff(self.c0 - o, self.c1)
    def __rsub__(self, o): return Aff(o - self.c0, -self.c1)
    def scale(self, a): return Aff(self.c0 * a, self.c1 * a)
    def compose(self, a0, a1):  # f(a0 + a1*y)
        return Aff(self.c0 + self.c1 * a0, self.c1 * a1)
    def root(self):
        if self.c1 == 0: return None
        return -self.c0 / self.c1

class Sq:
    def __init__(self, idx, cx, cy, s):
        self.idx = idx
        s = Fr(s); D = 1 + s * s
        self.s = s; self.cs = (1 - s * s) / D; self.sn = 2 * s / D
        self.c = (Fr(cx), Fr(cy))
        self.u = (self.cs, self.sn); self.n = (-self.sn, self.cs)
        self.phi = 2 * math.atan(float(s))
        h = Fr(1, 2)
        u, n, c = self.u, self.n, self.c
        # vertices: BL, BR, TR, TL ; edges: 0 bottom (BL-BR), 1 right (BR-TR), 2 top (TR-TL), 3 left (TL-BL)
        self.V = [(c[0] + a * h * u[0] + b * h * n[0], c[1] + a * h * u[1] + b * h * n[1])
                  for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))]

FLOOR = None

def dot(a, b): return a[0] * b[0] + a[1] * b[1]

def disjoint(A, B):
    for ax in (A.u, A.n, B.u, B.n):
        pa = [dot(v, ax) for v in A.V]; pb = [dot(v, ax) for v in B.V]
        if min(pb) > max(pa) or min(pa) > max(pb): return True
    return False

def check_config(sqs, k):
    for S in sqs:
        for v in S.V:
            if not (0 <= v[0] <= k and 0 <= v[1] <= k): return f'square {S.idx} outside container'
    for i in range(len(sqs)):
        for j in range(i + 1, len(sqs)):
            if not disjoint(sqs[i], sqs[j]): return f'squares {i},{j} not disjoint'
    return None

class Rec:
    def __init__(self):
        self.gaps = []      # (src, xa, xb, px, py, d, tstar, lam)
        self.D = Fr(0); self.W = Fr(0); self.H = Fr(0); self.T = Fr(0)
        self.E = {}         # (src, Y) -> measure
        self.M = []         # (Y, winner_src, loser_src, measure, sigma_len)
        self.Dcells_area = Fr(0)   # sum |gamma(C)| over D paths
        self.anom = []      # anomalies (lemma-level violations / internal)
        self.minlam = None
        self.entries = 0

def frame_of(sqs, src):
    if src < 0: return (Fr(1), Fr(0)), (Fr(0), Fr(1))
    S = sqs[src]; return S.u, S.n

def trace(sqs, k, delta, h, sF, tlimit=600.0, hist_cells=True):
    """returns Rec. sF: tilt-termination threshold on |s| (a(Y) >= alpha_F  <=> |s_Y| >= sF)."""
    t_start = time.time()
    k = Fr(k); delta = Fr(delta); h = Fr(h); sF = Fr(sF)
    rec = Rec()
    order = sorted(range(len(sqs)), key=lambda i: (sqs[i].c[1], i))
    arrivals = {i: [] for i in range(len(sqs))}

    def ray_step(src, xa, xb, px, py, g, lam, hist):
        if time.time() - t_start > tlimit: raise TimeUp()
        u, d = frame_of(sqs, src)
        # height cut: portion with p_y >= h is complete (H)
        r = (py - h).root() if py.c1 != 0 else None
        pieces = []
        cuts = sorted(set([xa, xb] + ([r] if r is not None and xa < r < xb else [])))
        for a, b in zip(cuts[:-1], cuts[1:]):
            m = (a + b) / 2
            if py(m) >= h: rec.H += b - a
            else: pieces.append((a, b))
        for (a, b) in pieces:
            _ray_piece(src, a, b, px, py, g, lam, hist, u, d)

    def _ray_piece(src, xa, xb, px, py, g, lam, hist, u, d):
        eta0 = px.c0 * d[0] + py.c0 * d[1]
        assert px.c1 * d[0] + py.c1 * d[1] == 0
        xi = Aff(px.c0 * u[0] + py.c0 * u[1], px.c1 * u[0] + py.c1 * u[1])
        assert xi.c1 == lam and lam != 0
        if rec.minlam is None or abs(lam) < rec.minlam: rec.minlam = abs(lam)
        lo_x, hi_x = xa, xb
        xis = sorted([xi(xa), xi(xb)])
        cand = []
        bps = set([xa, xb])
        for Y in sqs:
            if src >= 0 and Y.idx == src: continue
            fv = [(dot(v, u), dot(v, d)) for v in Y.V]
            ximin = min(f[0] for f in fv); ximax = max(f[0] for f in fv)
            if ximax <= xis[0] or ximin >= xis[1]: continue
            if max(f[1] for f in fv) < eta0: continue
            cand.append((Y, fv))
            for f in fv:
                xv = (f[0] - xi.c0) / xi.c1
                if xa < xv < xb: bps.add(xv)
        bps = sorted(bps)
        for a, b in zip(bps[:-1], bps[1:]):
            m = (a + b) / 2; xm = xi(m)
            best = None
            for (Y, fv) in cand:
                cr = []
                for e in range(4):
                    f1 = fv[e]; f2 = fv[(e + 1) % 4]
                    if f1[0] == f2[0]: continue
                    lo_, hi_ = (f1[0], f2[0]) if f1[0] < f2[0] else (f2[0], f1[0])
                    if lo_ < xm < hi_:
                        et = f1[1] + (xm - f1[0]) * (f2[1] - f1[1]) / (f2[0] - f1[0])
                        cr.append((et, e, f1, f2))
                if not cr: continue
                assert len(cr) == 2
                cr.sort(key=lambda z: z[0])
                lo, hi = cr[0][0], cr[1][0]
                if hi < eta0: continue
                if lo < eta0:
                    rec.anom.append(('start_inside', src, Y.idx)); continue
                if best is None or lo < best[0]:
                    best = (lo, Y, cr[0][1], cr[0][2], cr[0][3])
                elif lo == best[0]:
                    rec.anom.append(('tie_hit', Y.idx, best[1].idx))
            funcs = []
            tD = Aff(delta) - g
            funcs.append(('D', tD))
            if d[0] < 0: funcs.append(('W', px.scale(-1 / d[0])))
            elif d[0] > 0: funcs.append(('W', (Aff(k) - px).scale(1 / d[0])))
            tS = None
            if best is not None:
                lo, Y, e, f1, f2 = best
                slope = (f2[1] - f1[1]) / (f2[0] - f1[0])
                # eta_edge(xi(x)) - eta0
                tS = Aff(f1[1] - eta0 + slope * (xi.c0 - f1[0]), slope * xi.c1)
                funcs.append(('C', tS))
            # (H handled at height cut + after pass) still include tH for safety
            funcs.append(('H', (Aff(h) - py).scale(1 / d[1])))
            sub = set([a, b])
            for i in range(len(funcs)):
                for j in range(i + 1, len(funcs)):
                    rt = (funcs[i][1] - funcs[j][1]).root()
                    if rt is not None and a < rt < b: sub.add(rt)
            sub = sorted(sub)
            pr = {'D': 0, 'W': 1, 'C': 2, 'H': 3} if not MUT['contact_before_death'] else {'C': 0, 'D': 1, 'W': 2, 'H': 3}
            for sa, sb in zip(sub[:-1], sub[1:]):
                sm = (sa + sb) / 2
                vals = [(f(sm), pr[nm], nm, f) for nm, f in funcs]
                mv = min(v[0] for v in vals)
                ev = min([v for v in vals if v[0] == mv], key=lambda v: v[1])
                nm, tf = ev[2], ev[3]
                if tf(sa) < 0 or tf(sb) < 0:
                    rec.anom.append(('neg_t', nm, src, float(tf(sa)), float(tf(sb))))
                rec.gaps.append((src, sa, sb, px, py, d, tf, lam, g))
                if nm == 'D':
                    rec.D += sb - sa
                    if hist_cells:
                        cells = hist + [(lam, tf)]
                        for (lm, tt) in cells:
                            rec.Dcells_area += abs(lm) * (tt(sa) + tt(sb)) / 2 * (sb - sa)
                elif nm == 'W':
                    rec.W += sb - sa
                elif nm == 'H':
                    rec.H += sb - sa
                else:
                    lo, Y, e, f1, f2 = best
                    if e == 0:
                        qx = px + tS.scale(d[0]); qy = py + tS.scale(d[1])
                        sig = Aff((qx.c0 - Y.c[0]) * Y.u[0] + (qy.c0 - Y.c[1]) * Y.u[1],
                                  qx.c1 * Y.u[0] + qy.c1 * Y.u[1])
                        garr = g + tS
                        if src >= 0 and not (Y.c[1] > sqs[src].c[1]):
                            rec.anom.append(('DAG', src, Y.idx))
                        arrivals[Y.idx].append(dict(src=src, xa=sa, xb=sb, qx=qx, qy=qy, sig=sig, g=garr,
                                                    lam=lam, hist=hist + [(lam, tS)]))
                    elif e == 2:
                        rec.anom.append(('top_contact', src, Y.idx))
                    else:
                        key = (src, Y.idx)
                        rec.E[key] = rec.E.get(key, Fr(0)) + (sb - sa)

    def resolve(Y):
        A = arrivals[Y.idx]
        if not A: return
        # convert to sigma parametrization
        items = []
        for i, ar in enumerate(A):
            sg = ar['sig']; mu = sg.c1
            s1, s2 = sg(ar['xa']), sg(ar['xb'])
            # x(sigma) = (sigma - c0)/mu
            xs = Aff(-sg.c0 / mu, 1 / mu)
            gs = ar['g'].compose(xs.c0, xs.c1)
            items.append((min(s1, s2), max(s1, s2), xs, gs, i))
        bps = set()
        for it in items: bps.add(it[0]); bps.add(it[1])
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                a = items[i]; b = items[j]
                lo = max(a[0], b[0]); hi = min(a[1], b[1])
                if lo >= hi: continue
                for f in ((a[3] - b[3]), (a[2] - b[2])):
                    rt = f.root()
                    if rt is not None and lo < rt < hi: bps.add(rt)
        bps = sorted(bps)
        win = {}
        for sa, sb in zip(bps[:-1], bps[1:]):
            sm = (sa + sb) / 2
            act = [it for it in items if it[0] <= sa and sb <= it[1]]
            if not act: continue
            act.sort(key=lambda it: (it[3](sm), it[2](sm)))
            w = act[0]
            if MUT['no_r1']:
                for it in act: win.setdefault(it[4], []).append((sa, sb))
                continue
            win.setdefault(w[4], []).append((sa, sb))
            for lo_it in act[1:]:
                arw = A[w[4]]; arl = A[lo_it[4]]
                if arw['src'] == arl['src']:
                    rec.anom.append(('same_source_merge', Y.idx, arw['src']))
                if arw['g'] is not None and lo_it[3](sm) == w[3](sm) and lo_it[2](sm) == w[2](sm):
                    rec.anom.append(('identical_arrival', Y.idx))
                meas = abs(lo_it[2](sb) - lo_it[2](sa))
                rec.M.append((Y.idx, arw['src'], arl['src'], meas, sb - sa))
        # winners -> pass or T
        tilt = abs(Y.s) >= sF
        for i, ivs in win.items():
            ar = A[i]
            xs = items[i][2]
            # merge contiguous sigma intervals
            ivs.sort(); merged = []
            for a, b in ivs:
                if merged and merged[-1][1] == a: merged[-1] = (merged[-1][0], b)
                else: merged.append((a, b))
            for a, b in merged:
                x1, x2 = sorted([xs(a), xs(b)])
                if tilt:
                    rec.T += x2 - x1; continue
                rec.entries += 1
                npx = ar['qx'] + Y.n[0]; npy = ar['qy'] + Y.n[1]
                nlam = ar['sig'].c1
                ray_step(Y.idx, x1, x2, npx, npy, ar['g'], nlam, ar['hist'])
        arrivals[Y.idx] = []

    ray_step(-1, Fr(0), k, Aff(0, 1), Aff(0), Aff(0), Fr(1), [])
    for i in order:
        resolve(sqs[i])
    rec.k = k
    rec.total = rec.D + rec.W + rec.H + rec.T + sum(rec.E.values(), Fr(0)) + sum((m[3] for m in rec.M), Fr(0))
    return rec
