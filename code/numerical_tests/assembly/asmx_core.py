# -*- coding: utf-8 -*-
"""asmx_core.py -- independent end-to-end numerical test of the assembled path-flow chain of the manuscript
(Theorems 9.1 and 9.2 and Lemma 9.6 of Section 9, and Steps 6-12 of the proof of Theorem 10.3).

Written from scratch for this test on 2026-10-06.  No other tracer is imported or run.
Float64 piecewise-affine tracer: the floor interval (0,k) is split into pieces on which
the combinatorial type is constant and every point / gap / time is affine in x.

Conventions (Definitions 3.9, 3.10, 3.14, 3.20 and 3.21 of the manuscript):
  phase angle phi in (-pi/4, pi/4]; u=(cos,sin), n=(-sin,cos); bot = c - n/2 + t u.
  state (p, d=n_X, X, g); ray step t* = min(t_S, t_W, t_D, t_H); priority D > W > contact > H.
  contact in relint bot(Y) -> entry candidate; side/vertex -> E.
  R1 (master flow only): squares processed by increasing centre height; among the alive
  candidates reaching the same point q, the lexicographic minimum of (g, x) wins, others M.
  R2: winner with a(Y) >= alpha_F -> T at q; else pass [q, q+n_Y], exit -> H if exit_y >= h_F.
  F_s = truncation of F^0 at the first R1-passed entry with a >= alpha(s) (T_s), at height
  s+2 (H_s: gap point with t_Hs < t*, or exit >= s+2), or an earlier F^0 termination.
"""
import math, heapq, time
import numpy as np

PI4 = math.pi / 4.0


def norm_phi(phi):
    q = math.pi / 2.0
    p = (phi + PI4) % q - PI4
    if p <= -PI4 + 1e-15:
        p = PI4
    return p


class Sq:
    __slots__ = ('i', 'cx', 'cy', 'phi', 'a', 'co', 'si', 'u', 'n', 'V', 'v', 'yhi',
                 'xlo', 'xhi', 'sa', 'ca', 'ext', 'nrm')

    def __init__(self, i, cx, cy, phi):
        phi = norm_phi(phi)
        self.i = i; self.cx = cx; self.cy = cy; self.phi = phi; self.a = abs(phi)
        co = math.cos(phi); si = math.sin(phi)
        self.co = co; self.si = si
        self.u = (co, si); self.n = (-si, co)
        V = []
        for (p, q) in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5)):
            V.append((cx + p * co - q * si, cy + p * si + q * co))
        self.V = V   # V0 bot-left, V1 bot-right, V2 top-right, V3 top-left (local frame)
        ys = [t[1] for t in V]; xs = [t[0] for t in V]
        self.v = min(ys); self.yhi = max(ys); self.xlo = min(xs); self.xhi = max(xs)
        self.sa = math.sin(self.a); self.ca = math.cos(self.a)
        self.ext = self.yhi - self.v
        # outward normals of edges (V0V1)=bot, (V1V2)=right, (V2V3)=top, (V3V0)=left
        self.nrm = ((si, -co), (co, si), (-si, co), (-co, -si))


def sat_gap(A, B):
    """max over the 4 face normals of the projection gap (positive => disjoint with margin)."""
    best = -1e300
    for S in (A, B):
        for ax in (S.u, S.n):
            pa = [vx * ax[0] + vy * ax[1] for (vx, vy) in A.V]
            pb = [vx * ax[0] + vy * ax[1] for (vx, vy) in B.V]
            g = max(min(pb) - max(pa), min(pa) - max(pb))
            if g > best:
                best = g
    return best


class Packing:
    def __init__(self, k, sqs):
        self.k = k
        self.S = [Sq(i, *t) for i, t in enumerate(sqs)]
        self.N = len(self.S)
        self.W = k * k - self.N
        self.grid = {}
        for s in self.S:
            for gx in range(int(math.floor(s.xlo)), int(math.floor(s.xhi)) + 1):
                for gy in range(int(math.floor(s.v)), int(math.floor(s.yhi)) + 1):
                    self.grid.setdefault((gx, gy), []).append(s.i)
        self.v_arr = np.array([s.v for s in self.S])
        self.sa_arr = np.array([s.sa for s in self.S])
        self.ca_arr = np.array([s.ca for s in self.S])
        self.ext_arr = np.array([s.ext for s in self.S])
        self.a_arr = np.array([s.a for s in self.S])
        self.yhi_arr = np.array([s.yhi for s in self.S])

    def query(self, x0, x1, y0, y1):
        out = set()
        for gx in range(int(math.floor(x0)) - 1, int(math.floor(x1)) + 1):
            for gy in range(int(math.floor(y0)) - 1, int(math.floor(y1)) + 1):
                l = self.grid.get((gx, gy))
                if l:
                    out.update(l)
        return out

    def validate(self, tol=1e-9):
        k = self.k
        for s in self.S:
            for (vx, vy) in s.V:
                if vx < -1e-12 or vx > k + 1e-12 or vy < -1e-12 or vy > k + 1e-12:
                    return False, ('outside', s.i)
        for s in self.S:
            for j in self.query(s.xlo, s.xhi, s.v, s.yhi):
                if j <= s.i:
                    continue
                t = self.S[j]
                if t.xlo > s.xhi + 1e-6 or t.xhi < s.xlo - 1e-6 or t.v > s.yhi + 1e-6 or t.yhi < s.v - 1e-6:
                    continue
                if sat_gap(s, t) <= tol:
                    return False, ('overlap', s.i, j)
        return True, None

    def reflect(self):
        return Packing(self.k, [(s.cx, self.k - s.cy, -s.phi) for s in self.S])

    def chords(self, y):
        t = y - self.v_arr
        sa = self.sa_arr; ca = self.ca_arr; ext = self.ext_arr
        tilt = self.a_arr > 1e-15
        sc = np.where(tilt, sa * ca, 1.0)
        c = np.where(t < sa, t / sc, np.where(t <= ca, 1.0 / ca, (ext - t) / sc))
        c = np.where(tilt, c, 1.0)
        c = np.where((t < 0) | (t > ext), 0.0, c)
        return c


# ------------------------------------------------------------------------------------
#  Flow (master flow F^0 and its reflection), piecewise affine in the floor coordinate x
# ------------------------------------------------------------------------------------
class Flow:
    """recs: (x0,x1,typ,node,T0x,T0y,T1x,T1y,sq); termination point = T0 + x*T1.
    node: (parent,'R',(P0x,P0y,P1x,P1y,dx,dy,t0,t1))  ray segment p(x)+t d, t in [0,t0+t1 x]
          (parent,'E',(j,Q0x,Q0y,Q1x,Q1y,passed))    R1-passed entry at q(x) into square j
    """

    def __init__(self, P, delta, alphaF, hF, orderB=False):
        self.P = P; self.k = P.k; self.delta = delta; self.alphaF = alphaF; self.hF = hF
        self.orderB = orderB      # negative control: "T before M" at T_max squares
        self.recs = []
        self.cands = {}
        self.heap = []
        self.inheap = set()
        self.done = set()
        self.anom = []
        self.dropped = 0.0
        self.nR1 = 0
        self.nR1multi = 0

    def run(self, tlimit=60.0):
        t0 = time.time()
        self._ray((0.0, float(self.k), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1, None))
        while self.heap:
            cy, j = heapq.heappop(self.heap)
            self.inheap.discard(j)
            self._square(j)
            if time.time() - t0 > tlimit:
                raise TimeoutError('flow time limit')
        self.paths = [node_path(r[3]) for r in self.recs]

    # -- one ray step for a piece --
    def _ray(self, pc):
        xa, xb, P0x, P0y, P1x, P1y, g0, g1, src, node = pc
        if xb - xa < 1e-13:
            self.dropped += max(0.0, xb - xa)
            return
        S = self.P.S; k = self.k; delta = self.delta; hF = self.hF
        if src < 0:
            ux, uy, dx, dy = 1.0, 0.0, 0.0, 1.0; cyX = -1e9
        else:
            X = S[src]; ux, uy = X.u; dx, dy = X.n; cyX = X.cy
        wA = P0x * ux + P0y * uy; wB = P1x * ux + P1y * uy
        hA = P0x * dx + P0y * dy; hB = P1x * dx + P1y * dy
        if wB <= 1e-12:
            self.anom.append(('wB<=0', src, wB)); return
        Wlo = wA + wB * xa; Whi = wA + wB * xb
        ga = g0 + g1 * xa; gb = g0 + g1 * xb
        reach = delta - min(ga, gb) + 1e-9
        pax = P0x + P1x * xa; pay = P0y + P1y * xa; pbx = P0x + P1x * xb; pby = P0y + P1y * xb
        xs4 = (pax, pbx, pax + reach * dx, pbx + reach * dx)
        ys4 = (pay, pby, pay + reach * dy, pby + reach * dy)
        cand = self.P.query(min(xs4), max(xs4), min(ys4), max(ys4))
        hmax_start = max(hA + hB * xa, hA + hB * xb)
        edges = []
        for j in cand:
            if j == src:
                continue
            Y = S[j]
            fw = [vx * ux + vy * uy for (vx, vy) in Y.V]
            if max(fw) <= Wlo or min(fw) >= Whi:
                continue
            fh = [vx * dx + vy * dy for (vx, vy) in Y.V]
            if min(fh) > hmax_start + reach + 1e-9:
                continue
            if max(fh) < min(hA + hB * xa, hA + hB * xb) - 1e-9:
                continue
            for e in range(4):
                nx, ny = Y.nrm[e]
                if nx * dx + ny * dy >= -1e-15:
                    continue
                a = e; b = (e + 1) % 4
                w1, h1, w2, h2 = fw[a], fh[a], fw[b], fh[b]
                if w1 > w2:
                    w1, h1, w2, h2 = w2, h2, w1, h1
                if w2 - w1 < 1e-12 or w2 <= Wlo or w1 >= Whi:
                    continue
                edges.append((w1, h1, w2, h2, j, e))
        bps = {Wlo, Whi}
        for ed in edges:
            if Wlo < ed[0] < Whi:
                bps.add(ed[0])
            if Wlo < ed[2] < Whi:
                bps.add(ed[2])
        bps = sorted(bps)
        segs = []
        for i in range(len(bps) - 1):
            b0, b1 = bps[i], bps[i + 1]
            if b1 - b0 < 1e-13:
                continue
            wm = 0.5 * (b0 + b1); xm = (wm - wA) / wB; hs = hA + hB * xm
            best = None; bt = 1e300
            for ed in edges:
                if ed[0] < wm < ed[2]:
                    he = ed[1] + (ed[3] - ed[1]) * (wm - ed[0]) / (ed[2] - ed[0])
                    t = he - hs
                    if t >= -1e-9 and t < bt:
                        bt = t; best = ed
            x0 = (b0 - wA) / wB; x1 = (b1 - wA) / wB
            if segs and segs[-1][2] is best:
                segs[-1][1] = x1
            else:
                segs.append([x0, x1, best])
        if not segs:
            self.dropped += xb - xa
            return
        segs[0][0] = xa; segs[-1][1] = xb
        for (x0, x1, ed) in segs:
            if x1 - x0 < 1e-14:
                self.dropped += max(0.0, x1 - x0); continue
            F = [(0, delta - g0, -g1)]
            if dx > 1e-15:
                F.append((1, (k - P0x) / dx, -P1x / dx))
            elif dx < -1e-15:
                F.append((1, -P0x / dx, -P1x / dx))
            if ed is not None:
                w1, h1, w2, h2, j, e = ed
                sl = (h2 - h1) / (w2 - w1)
                F.append((2, h1 + sl * (wA - w1) - hA, sl * wB - hB))
            F.append((3, (hF - P0y) / dy, -P1y / dy))
            xs = [x0, x1]
            for i in range(len(F)):
                for i2 in range(i + 1, len(F)):
                    d1 = F[i][2] - F[i2][2]
                    if abs(d1) > 1e-300:
                        xc = (F[i2][1] - F[i][1]) / d1
                        if x0 < xc < x1:
                            xs.append(xc)
            xs.sort()
            out = []
            for i in range(len(xs) - 1):
                u0, u1 = xs[i], xs[i + 1]
                if u1 - u0 < 1e-14:
                    if out:
                        out[-1][1] = u1
                    continue
                xm = 0.5 * (u0 + u1)
                vals = [f[1] + f[2] * xm for f in F]
                vmin = min(vals)
                tol = 1e-12 * max(1.0, abs(vmin))
                ch = None
                for idx in range(len(F)):
                    if vals[idx] <= vmin + tol:
                        ch = idx; break
                if out and out[-1][2] == ch:
                    out[-1][1] = u1
                else:
                    out.append([u0, u1, ch])
            if out:
                out[0][0] = x0; out[-1][1] = x1
            for (u0, u1, ch) in out:
                f = F[ch]
                self._emit(u0, u1, P0x, P0y, P1x, P1y, g0, g1, node, dx, dy, f[0], f[1], f[2], ed, cyX)

    def _emit(self, xa, xb, P0x, P0y, P1x, P1y, g0, g1, node, dx, dy, prio, t0, t1, ed, cyX):
        rn = (node, 'R', (P0x, P0y, P1x, P1y, dx, dy, t0, t1))
        T0x = P0x + t0 * dx; T0y = P0y + t0 * dy; T1x = P1x + t1 * dx; T1y = P1y + t1 * dy
        if prio == 0:
            self.recs.append((xa, xb, 'D', rn, T0x, T0y, T1x, T1y, -1))
        elif prio == 1:
            self.recs.append((xa, xb, 'W', rn, T0x, T0y, T1x, T1y, -1))
        elif prio == 3:
            self.recs.append((xa, xb, 'H', rn, T0x, T0y, T1x, T1y, -1))
        else:
            j, e = ed[4], ed[5]
            if e == 0:
                Y = self.P.S[j]
                if Y.cy <= cyX + 1e-12:
                    self.anom.append(('DAG', j, Y.cy, cyX))
                if j in self.done:
                    self.anom.append(('late-candidate', j))
                    self.recs.append((xa, xb, 'E', rn, T0x, T0y, T1x, T1y, j))
                    return
                self.cands.setdefault(j, []).append((xa, xb, T0x, T0y, T1x, T1y, g0 + t0, g1 + t1, rn))
                if j not in self.inheap:
                    heapq.heappush(self.heap, (Y.cy, j)); self.inheap.add(j)
            else:
                self.recs.append((xa, xb, 'E', rn, T0x, T0y, T1x, T1y, j))

    def _square(self, j):
        self.done.add(j)
        C = self.cands.pop(j, [])
        if not C:
            return
        Y = self.P.S[j]
        ux, uy = Y.u; nx, ny = Y.n
        self.nR1 += 1
        wins = []; loses = []
        if len(C) == 1:
            c = C[0]
            wins.append((c, c[0], c[1]))
        else:
            self.nR1multi += 1
            its = []
            for ci, c in enumerate(C):
                xa, xb, Q0x, Q0y, Q1x, Q1y, gq0, gq1, rn = c
                xi0 = (Q0x - Y.cx) * ux + (Q0y - Y.cy) * uy; xi1 = Q1x * ux + Q1y * uy
                if xi1 <= 1e-12:
                    self.anom.append(('xi1<=0', j, xi1))
                    xi1 = 1e-12
                B = gq1 / xi1; A = gq0 - B * xi0
                its.append([xi0 + xi1 * xa, xi0 + xi1 * xb, A, B, xa, ci, xi0, xi1])
            bps = set()
            for it in its:
                bps.add(it[0]); bps.add(it[1])
            for i in range(len(its)):
                for i2 in range(i + 1, len(its)):
                    lo = max(its[i][0], its[i2][0]); hi = min(its[i][1], its[i2][1])
                    if hi > lo:
                        dB = its[i][3] - its[i2][3]
                        if abs(dB) > 1e-300:
                            xc = (its[i2][2] - its[i][2]) / dB
                            if lo < xc < hi:
                                bps.add(xc)
            bps = sorted(bps)
            lab = {it[5]: [] for it in its}
            for i in range(len(bps) - 1):
                b0, b1 = bps[i], bps[i + 1]
                if b1 - b0 < 1e-14:
                    continue
                m = 0.5 * (b0 + b1)
                cover = [it for it in its if it[0] < m < it[1]]
                if not cover:
                    continue
                if len(cover) == 1:
                    lab[cover[0][5]].append((b0, b1, True)); continue
                gv = [it[2] + it[3] * m for it in cover]
                gmin = min(gv)
                tie = [it for it, g in zip(cover, gv) if g <= gmin + 1e-12]
                w = min(tie, key=lambda it: it[4])
                for it in cover:
                    lab[it[5]].append((b0, b1, it is w))
            for it in its:
                c = C[it[5]]; xi0, xi1 = it[6], it[7]
                ivs = []
                for (b0, b1, wflag) in lab[it[5]]:
                    x0 = (b0 - xi0) / xi1; x1 = (b1 - xi0) / xi1
                    if ivs and ivs[-1][2] == wflag:
                        ivs[-1][1] = x1
                    else:
                        ivs.append([x0, x1, wflag])
                if not ivs:
                    continue
                ivs[0][0] = max(ivs[0][0], c[0]); ivs[-1][1] = min(ivs[-1][1], c[1])
                if abs(ivs[0][0] - c[0]) < 1e-11:
                    ivs[0][0] = c[0]
                if abs(ivs[-1][1] - c[1]) < 1e-11:
                    ivs[-1][1] = c[1]
                for (x0, x1, wflag) in ivs:
                    if x1 - x0 <= 0:
                        continue
                    if wflag:
                        wins.append((c, x0, x1))
                    else:
                        loses.append((c, x0, x1))
        for (c, x0, x1) in loses:
            xa, xb, Q0x, Q0y, Q1x, Q1y, gq0, gq1, rn = c
            if self.orderB and Y.a >= self.alphaF:
                ln = (rn, 'L', (j, Q0x, Q0y, Q1x, Q1y, False))   # T-terminated loser (not R1-passed)
                self.recs.append((x0, x1, 'T', ln, Q0x, Q0y, Q1x, Q1y, j))
                continue
            self.recs.append((x0, x1, 'M', rn, Q0x, Q0y, Q1x, Q1y, j))
        for (c, x0, x1) in wins:
            xa, xb, Q0x, Q0y, Q1x, Q1y, gq0, gq1, rn = c
            if Y.a >= self.alphaF:
                en = (rn, 'E', (j, Q0x, Q0y, Q1x, Q1y, False))
                self.recs.append((x0, x1, 'T', en, Q0x, Q0y, Q1x, Q1y, j))
                continue
            en = (rn, 'E', (j, Q0x, Q0y, Q1x, Q1y, True))
            E0x = Q0x + nx; E0y = Q0y + ny
            for (u0, u1, pos) in split_affine(x0, x1, E0y - self.hF, Q1y, ge=True):
                if pos:
                    self.recs.append((u0, u1, 'H', en, E0x, E0y, Q1x, Q1y, -1))
                else:
                    self._ray((u0, u1, E0x, E0y, Q1x, Q1y, gq0, gq1, j, en))


def node_path(node):
    out = []
    while node is not None:
        out.append(node); node = node[0]
    out.reverse()
    return out


def split_affine(a, b, c0, c1, ge=False):
    """split [a,b] by the sign of c0+c1 x; returns [(u,v,pos)], pos = (>0) or (>=0 if ge)."""
    if abs(c1) < 1e-15:
        val = c0
        pos = (val >= 0) if ge else (val > 0)
        return [(a, b, pos)]
    xr = -c0 / c1
    if xr <= a or xr >= b:
        xm = 0.5 * (a + b)
        val = c0 + c1 * xm
        pos = (val >= 0) if ge else (val > 0)
        return [(a, b, pos)]
    if c1 > 0:
        return [(a, xr, False), (xr, b, True)]
    return [(a, xr, True), (xr, b, False)]


# ------------------------------------------------------------------------------------
#  Truncation walk: F^0 (hcap=None) or F_s (alpha_s, hcap=s+2)
# ------------------------------------------------------------------------------------
TYP = {'D': 0, 'W': 1, 'E': 2, 'M': 3, 'T': 4, 'H': 5}


def walk(flow, alpha_s, hcap):
    """returns segs (list), terms (list), contacts (list).
    segs : (x0,x1,kind,sq,A0x,A0y,A1x,A1y,Dx,Dy,t0,t1)   kind 0=gap, 1=pass
    terms: (x0,x1,typ,T0x,T0y,T1x,T1y,sq)
    contacts: (sq,x0,x1) alive bottom contacts (R1 winners and M losers) before truncation."""
    S = flow.P.S
    segs = []; terms = []; contacts = []
    for rec, path in zip(flow.recs, flow.paths):
        x0, x1, typ, node, T0x, T0y, T1x, T1y, sq = rec
        iv = [(x0, x1)]
        for nd in path:
            if not iv:
                break
            kind = nd[1]; d = nd[2]
            if kind == 'R':
                P0x, P0y, P1x, P1y, dx, dy, t0, t1 = d
                if hcap is None:
                    for (a, b) in iv:
                        segs.append((a, b, 0, -1, P0x, P0y, P1x, P1y, dx, dy, t0, t1))
                else:
                    th0 = (hcap - P0y) / dy; th1 = -P1y / dy
                    niv = []
                    for (a, b) in iv:
                        for (u, v, pos) in split_affine(a, b, t0 - th0, t1 - th1):
                            if v - u <= 0:
                                continue
                            if pos:
                                segs.append((u, v, 0, -1, P0x, P0y, P1x, P1y, dx, dy, th0, th1))
                                terms.append((u, v, 'H', P0x + th0 * dx, P0y + th0 * dy,
                                              P1x + th1 * dx, P1y + th1 * dy, -1))
                            else:
                                segs.append((u, v, 0, -1, P0x, P0y, P1x, P1y, dx, dy, t0, t1))
                                niv.append((u, v))
                    iv = niv
            else:
                j, Q0x, Q0y, Q1x, Q1y, passed = d
                Y = S[j]
                for (a, b) in iv:
                    contacts.append((j, a, b))
                if Y.a >= alpha_s:
                    for (a, b) in iv:
                        terms.append((a, b, 'T', Q0x, Q0y, Q1x, Q1y, j))
                    iv = []
                else:
                    if not passed:
                        raise RuntimeError('T-entry with a < alpha_s')
                    nx, ny = Y.n
                    for (a, b) in iv:
                        segs.append((a, b, 1, j, Q0x, Q0y, Q1x, Q1y, nx, ny, 1.0, 0.0))
                    if hcap is not None:
                        E0y = Q0y + ny
                        niv = []
                        for (a, b) in iv:
                            for (u, v, pos) in split_affine(a, b, E0y - hcap, Q1y, ge=True):
                                if v - u <= 0:
                                    continue
                                if pos:
                                    terms.append((u, v, 'H', Q0x + nx, E0y, Q1x, Q1y, -1))
                                else:
                                    niv.append((u, v))
                        iv = niv
        for (a, b) in iv:
            if typ == 'T':
                raise RuntimeError('T rec survived walk')
            terms.append((a, b, typ, T0x, T0y, T1x, T1y, sq))
            if typ == 'M':
                contacts.append((sq, a, b))
    return segs, terms, contacts


class LineEval:
    """vectorised line quantities on l_y for a flow given by (segs, terms)."""

    def __init__(self, P, segs, terms):
        self.P = P; self.k = P.k
        if segs:
            A = np.array([s[:2] + s[2:4] + s[4:] for s in segs], dtype=float)
        else:
            A = np.zeros((0, 12))
        self.x0, self.x1, self.kind, self.sq = A[:, 0], A[:, 1], A[:, 2].astype(int), A[:, 3].astype(int)
        self.A0x, self.A0y, self.A1x, self.A1y = A[:, 4], A[:, 5], A[:, 6], A[:, 7]
        self.Dx, self.Dy, self.t0, self.t1 = A[:, 8], A[:, 9], A[:, 10], A[:, 11]
        self.B0y = self.A0y + self.t0 * self.Dy; self.B1y = self.A1y + self.t1 * self.Dy
        T = np.array([(t[0], t[1], TYP[t[2]], t[3], t[4], t[5], t[6], t[7]) for t in terms], dtype=float) \
            if terms else np.zeros((0, 8))
        self.tx0, self.tx1, self.ttyp = T[:, 0], T[:, 1], T[:, 2].astype(int)
        self.T0y, self.T1y = T[:, 4], T[:, 6]
        self.tsq = T[:, 7].astype(int)

    @staticmethod
    def _below(x0, x1, c0, c1, y):
        """measure of {x in [x0,x1]: c0 + c1 x < y} (vectorised)."""
        lo = x0.copy(); hi = x1.copy()
        pos = c1 > 1e-15; neg = c1 < -1e-15; zer = ~(pos | neg)
        with np.errstate(divide='ignore', invalid='ignore'):
            r = (y - c0) / np.where(zer, 1.0, c1)
        hi = np.where(pos, np.minimum(hi, r), hi)
        lo = np.where(neg, np.maximum(lo, r), lo)
        m = np.maximum(hi - lo, 0.0)
        m = np.where(zer, np.where(c0 < y, x1 - x0, 0.0), m)
        return m

    def losses(self, y):
        m = self._below(self.tx0, self.tx1, self.T0y, self.T1y, y)
        L = np.bincount(self.ttyp, weights=m, minlength=6)
        return L  # index by TYP

    def at(self, y, want_ov=True):
        k = self.k
        lo = self.x0.copy(); hi = self.x1.copy()
        with np.errstate(divide='ignore', invalid='ignore'):
            # a(x) < y
            c1 = self.A1y; c0 = self.A0y
            pos = c1 > 1e-15; neg = c1 < -1e-15; zer = ~(pos | neg)
            r = (y - c0) / np.where(zer, 1.0, c1)
            hi = np.where(pos, np.minimum(hi, r), hi)
            lo = np.where(neg, np.maximum(lo, r), lo)
            ok = ~(zer & (c0 >= y))
            # b(x) > y
            c1 = self.B1y; c0 = self.B0y
            pos = c1 > 1e-15; neg = c1 < -1e-15; zer = ~(pos | neg)
            r = (y - c0) / np.where(zer, 1.0, c1)
            lo = np.where(pos, np.maximum(lo, r), lo)
            hi = np.where(neg, np.minimum(hi, r), hi)
            ok &= ~(zer & (c0 <= y))
        m = np.where(ok, np.maximum(hi - lo, 0.0), 0.0)
        isp = self.kind == 1
        mu = np.bincount(self.sq[isp], weights=m[isp], minlength=self.P.N) if self.P.N else np.zeros(0)
        gap = (~isp) & (m > 0)
        Fgap = float(m[gap].sum())
        ov = 0.0
        if want_ov and gap.any():
            r_ = self.Dx[gap] / self.Dy[gap]
            X0 = self.A0x[gap] + (y - self.A0y[gap]) * r_
            sig = self.A1x[gap] - self.A1y[gap] * r_
            la = X0 + sig * lo[gap]; lb = X0 + sig * hi[gap]
            l = np.minimum(la, lb); rr = np.maximum(la, lb)
            dens = 1.0 / np.maximum(np.abs(sig), 1e-300)
            pts = np.concatenate([l, rr]); dd = np.concatenate([dens, -dens])
            o = np.argsort(pts, kind='mergesort')
            pts = pts[o]; dd = dd[o]
            rho = np.cumsum(dd)
            widths = np.diff(pts)
            ov = float(np.sum(np.maximum(rho[:-1] - 1.0, 0.0) * widths))
        L = self.losses(y)
        c = self.P.chords(y)
        Sh = float(np.sum(np.maximum(c - mu, 0.0)))
        Ov = float(np.sum(np.maximum(mu - c, 0.0)))
        Lsum = float(L[:5].sum())   # D,W,E,M,T (H = completion, not a loss)
        Bd = k - Lsum - float(mu.sum()) - Fgap
        return dict(L=L, mu=mu, Fgap=Fgap, Ovgap=ov, Sh=Sh, Ov=Ov, Bd=Bd, c=c)

    def int_losses(self, ya, yb, types):
        """exact  int_{ya}^{yb} sum_{typ in types} L_typ(y) dy."""
        sel = np.isin(self.ttyp, types)
        x0 = self.tx0[sel]; x1 = self.tx1[sel]; c0 = self.T0y[sel]; c1 = self.T1y[sel]
        tot = 0.0
        for i in range(len(x0)):
            tot += _int_clamp(x0[i], x1[i], c0[i], c1[i], ya, yb)
        return tot


def _int_clamp(x0, x1, c0, c1, ya, yb):
    """int_{x0}^{x1} (yb - clamp(c0+c1 x, ya, yb)) dx  (= int_{ya}^{yb} 1[tau(x)<y] dy dx)."""
    if x1 <= x0:
        return 0.0
    pts = [x0, x1]
    if abs(c1) > 1e-15:
        for yy in (ya, yb):
            xr = (yy - c0) / c1
            if x0 < xr < x1:
                pts.append(xr)
    pts.sort()
    tot = 0.0
    for i in range(len(pts) - 1):
        a, b = pts[i], pts[i + 1]
        if b <= a:
            continue
        tm = c0 + c1 * 0.5 * (a + b)
        if tm <= ya:
            tot += (yb - ya) * (b - a)
        elif tm >= yb:
            pass
        else:
            ta = c0 + c1 * a; tb = c0 + c1 * b
            tot += (b - a) * (yb - 0.5 * (ta + tb))
    return tot


# ------------------------------------------------------------------------------------
#  T* from the master flow (R1-passed entries only; deaths and E contacts excluded)
# ------------------------------------------------------------------------------------
def tstar_pieces(flow, c0, y0, y1):
    """list of (x0,x1,kind,v0,v1): T*(x) = v0 if kind==0 (const) else v0+v1 x (affine);
    only the part of dom T* (T* <= y1)."""
    S = flow.P.S
    out = []
    for rec, path in zip(flow.recs, flow.paths):
        x0, x1 = rec[0], rec[1]
        ents = []
        for nd in path:
            if nd[1] == 'E':
                j, Q0x, Q0y, Q1x, Q1y, passed = nd[2]
                a = S[j].a
                if a <= 0:
                    continue
                cj = max(y0, (c0 / a) ** 2)
                ents.append((cj, Q0y - 2.0, Q1y))
        if not ents:
            continue
        consts = [e[0] for e in ents] + [y1]
        pts = {x0, x1}
        for (cj, e0, e1) in ents:
            if abs(e1) > 1e-15:
                for cc in consts:
                    xr = (cc - e0) / e1
                    if x0 < xr < x1:
                        pts.add(xr)
        for i in range(len(ents)):
            for i2 in range(i + 1, len(ents)):
                d1 = ents[i][2] - ents[i2][2]
                if abs(d1) > 1e-15:
                    xr = (ents[i2][1] - ents[i][1]) / d1
                    if x0 < xr < x1:
                        pts.add(xr)
        pts = sorted(pts)
        for i in range(len(pts) - 1):
            a, b = pts[i], pts[i + 1]
            if b - a <= 0:
                continue
            xm = 0.5 * (a + b)
            best = None; bv = 1e300
            for (cj, e0, e1) in ents:
                ev = e0 + e1 * xm
                if ev > cj:
                    val = ev; kd = (1, e0, e1)
                else:
                    val = cj; kd = (0, cj, 0.0)
                if val < bv:
                    bv = val; best = kd
            if bv <= y1:
                out.append((a, b, best[0], best[1], best[2]))
    return out


def tstar_stats(tp, c0, y1):
    intA = 0.0; dom = 0.0
    for (a, b, kd, v0, v1) in tp:
        dom += b - a
        if kd == 0 or abs(v1) < 1e-15:
            intA += c0 * (v0 + (v1 * 0.5 * (a + b) if kd else 0.0)) ** -0.5 * (b - a)
        else:
            ta = v0 + v1 * a; tb = v0 + v1 * b
            intA += 2.0 * c0 / v1 * (math.sqrt(tb) - math.sqrt(ta))
    return intA, dom


def tstar_le(tp, s):
    """measure of {T* <= s} and the interval list."""
    tot = 0.0; ivs = []
    for (a, b, kd, v0, v1) in tp:
        if kd == 0 or abs(v1) < 1e-15:
            val = v0 if kd == 0 else v0 + v1 * 0.5 * (a + b)
            if val <= s:
                tot += b - a; ivs.append((a, b))
        else:
            xr = (s - v0) / v1
            if v1 > 0:
                lo, hi = a, min(b, xr)
            else:
                lo, hi = max(a, xr), b
            if hi > lo:
                tot += hi - lo; ivs.append((lo, hi))
    return tot, ivs


def union_measure(ivs):
    if not ivs:
        return 0.0
    ivs = sorted(ivs)
    tot = 0.0; ca, cb = ivs[0]
    for (a, b) in ivs[1:]:
        if a > cb:
            tot += cb - ca; ca, cb = a, b
        else:
            cb = max(cb, b)
    tot += cb - ca
    return tot


def inter_measure(A, B):
    """measure of (union A) cap (union B)."""
    def merged(ivs):
        ivs = sorted(ivs); out = []
        for (a, b) in ivs:
            if out and a <= out[-1][1]:
                out[-1][1] = max(out[-1][1], b)
            else:
                out.append([a, b])
        return out
    A = merged(A); B = merged(B)
    i = j = 0; tot = 0.0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0]); hi = min(A[i][1], B[j][1])
        if hi > lo:
            tot += hi - lo
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return tot


# ------------------------------------------------------------------------------------
#  Gauss-Legendre helper
# ------------------------------------------------------------------------------------
_GLX, _GLW = np.polynomial.legendre.leggauss(10)


def gl(f, a, b):
    if b <= a:
        return 0.0
    xm = 0.5 * (a + b); xr = 0.5 * (b - a)
    xs = xm + xr * _GLX
    return float(xr * np.sum(_GLW * f(xs)))
