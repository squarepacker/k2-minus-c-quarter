# -*- coding: utf-8 -*-
"""asmx_chain.py -- chain evaluation: scaled-analogue parameters, line structure, per-family flow data,
and the four links of the assembled chain evaluated with MEASURED quantities:
  link 1: Theorem 9.1 (per-scale count) and the steps of its proof, Theorem 9.2 (lower bound for |Z(s)|);
  link 2: Lemma 9.6(c) ({T_s-terminated} = {T* <= s}) and Step 8 of the proof of Theorem 10.3 (first use of Tonelli);
  link 3: Step 10 of that proof (second use of Tonelli, the constant invF) and the integrated per-family chain;
  link 4: Steps 3-5 and 12 of that proof ((1 - omega0) l(H) <= the four contributions), floor and ceiling.
Written from scratch for this test (2026-10-06)."""
import math, time
import numpy as np
from asmx_core import (Packing, Flow, walk, LineEval, tstar_pieces, tstar_stats, tstar_le,
                       union_measure, inter_measure, gl, TYP, PI4)

TD, TW, TE, TM, TT, TH = 0, 1, 2, 3, 4, 5


class Params:
    def __init__(self, k, delta, c0, y0, c1, eps, omega0):
        self.k = k; self.delta = delta; self.c0 = c0; self.y0 = y0
        self.c1 = c1; self.eps = eps; self.omega0 = omega0
        self.y1 = k / 2.0 - 3.0; self.hmax = k / 2.0 - 1.0
        self.amax = c0 / math.sqrt(y0)
        self.bstar = c1 * ((1 - eps) * y0) ** -0.75
        bs = min(self.bstar, PI4)
        self.h = math.cos(bs) + math.sin(bs)        # replaces 1.0001 (max vertical extent of Z)
        self.bW = 2.0 * self.h                      # replaces b_W = 2.0002
        ss = np.linspace(y0, self.y1, 4001)
        q = (1.0 / np.cos(c0 / np.sqrt(ss)) - 1.0) * (ss + self.h)   # (sec alpha(s) - 1) q_y
        self.wq = float(q.max())
        self.w0 = self.wq * (1 + 1e-9) + delta + math.sin(bs) + (1 - math.cos(bs)) + 1e-12
        self.h0 = 1.0 - 2.0 * self.w0
        e = eps
        self.invF = 4 * (1 - e) ** 0.75 * (1 - (1 - e) ** 0.75) / (3 * e * (1 + self.bW / (e * y0)))
        self.top = (1 - eps) * self.y1

    def ok(self):
        msgs = []
        if not (self.w0 < 0.5): msgs.append('w0>=1/2')
        if not (self.y0 + 1 <= self.top): msgs.append('C7')
        if not (self.delta <= 0.5 * math.cos(2 * self.amax)): msgs.append('C6')
        if not (self.amax < PI4): msgs.append('amax')
        if not ((1 - self.eps) * self.y0 - self.h > 0): msgs.append('Vs<0')
        return msgs

    def alpha(self, s):
        return self.c0 / np.sqrt(s)

    def dalpha(self, s):          # |alpha'(s)|
        return 0.5 * self.c0 * s ** -1.5

    def betahat(self, t):
        return self.c1 * t ** -0.75

    def d(self):
        return dict(k=self.k, delta=self.delta, c0=self.c0, y0=self.y0, c1=self.c1, eps=self.eps,
                    omega0=self.omega0, y1=self.y1, amax=self.amax, bstar=self.bstar, h=self.h,
                    bW=self.bW, w0=self.w0, h0=self.h0, invF=self.invF)


def ell(ya, yb, k):
    """int_{ya}^{yb} d(y)^{-3/4} dy, with [ya,yb] on one side of k/2."""
    if yb <= ya:
        return 0.0
    if yb <= k / 2.0 + 1e-12:
        return 4.0 * (yb ** 0.25 - ya ** 0.25)
    return 4.0 * ((k - ya) ** 0.25 - (k - yb) ** 0.25)


def int_lin_ell(ya, yb, va, vb, k):
    """int of the linear function (va at ya, vb at yb) against d^{-3/4} dy."""
    if yb <= ya:
        return 0.0
    def f(y):
        lam = (y - ya) / (yb - ya)
        dd = np.minimum(y, k - y)
        return (va + (vb - va) * lam) * dd ** -0.75
    return gl(f, ya, yb)


class LineStruct:
    def __init__(self, P, prm):
        self.P = P; self.prm = prm
        k = P.k; y0 = prm.y0; w0 = prm.w0; om = prm.omega0; c1 = prm.c1; top = prm.top
        B = [0.0, float(k), k / 2.0, y0, k - y0, top, k - top]
        for m in range(k + 1):
            B += [m + w0, m + 1 - w0]
        for s in P.S:
            B += [s.v, s.v + s.sa, s.v + s.ca, s.yhi]
            if s.a > 1e-15:
                B += [s.v + s.sa * s.ca, s.yhi - s.sa * s.ca]
                ds = (c1 / s.a) ** (4.0 / 3.0)
                B += [ds, k - ds]
        B = np.unique(np.clip(np.array(B, dtype=float), 0.0, float(k)))
        rows = []
        self.l35_min = 1e9
        for i in range(len(B) - 1):
            ya, yb = float(B[i]), float(B[i + 1])
            if yb - ya < 1e-13:
                continue
            vals = []
            for y in (ya + 0.25 * (yb - ya), ya + 0.75 * (yb - ya)):
                c = P.chords(y)
                omg = k - float(c.sum())
                E = float(np.sum(np.where(c >= 1.0, c - 1.0, 0.0)))
                fs = float(np.sum(np.where((c > 0) & (c < 1.0), c, 0.0)))
                vals.append((omg, E, fs))
            def lin(idx, y):
                a = vals[0][idx]; b = vals[1][idx]
                y1_ = ya + 0.25 * (yb - ya); y2_ = ya + 0.75 * (yb - ya)
                return a + (b - a) * (y - y1_) / (y2_ - y1_)
            cuts = {ya, yb}
            for (idx, thr) in ((0, om), (1, 1.0 - om)):
                a = lin(idx, ya); b = lin(idx, yb)
                if (a - thr) * (b - thr) < 0:
                    yc = ya + (thr - a) / (b - a) * (yb - ya)
                    if ya < yc < yb:
                        cuts.add(yc)
            cuts = sorted(cuts)
            for j in range(len(cuts) - 1):
                u, v = cuts[j], cuts[j + 1]
                if v - u < 1e-14:
                    continue
                ym = 0.5 * (u + v)
                dd = min(ym, k - ym)
                fr = ym - math.floor(ym)
                inH = (y0 <= dd <= top) and (w0 <= fr <= 1 - w0)
                omu, omv, omm = lin(0, u), lin(0, v), lin(0, ym)
                Eu, Ev, Em = max(lin(1, u), 0.0), max(lin(1, v), 0.0), max(lin(1, ym), 0.0)
                fsu, fsv = max(lin(2, u), 0.0), max(lin(2, v), 0.0)
                kind = 0
                if inH:
                    if omm >= om:
                        kind = 1
                    else:
                        meets = (P.v_arr <= ym) & (ym <= P.yhi_arr)
                        maxa = float(P.a_arr[meets].max()) if meets.any() else 0.0
                        beta = c1 * dd ** -0.75
                        kind = 2 if maxa >= beta else 3
                if kind == 3:
                    fsm = max(lin(2, ym), 0.0)
                    self.l35_min = min(self.l35_min, fsm - (1.0 - omm - Em))
                fu = max(1.0 - om - Eu, 0.0); fv = max(1.0 - om - Ev, 0.0)
                rows.append((u, v, kind, 0 if ym < k / 2.0 else 1, omu, omv, Eu, Ev, fu, fv, fsu, fsv))
        self.rows = rows
        R = np.array(rows) if rows else np.zeros((0, 12))
        self.R = R
        lH = lYi = lYii = lYiii = lYb = lYt = KE = fb = ft = fsb = fst = 0.0
        omYi = lenYi = 0.0
        Yb = []
        for r in rows:
            u, v, kind, side, omu, omv, Eu, Ev, fu, fv, fsu, fsv = r
            if kind == 0:
                continue
            le = ell(u, v, k)
            lH += le
            if kind == 1:
                lYi += le; lenYi += v - u; omYi += 0.5 * (omu + omv) * (v - u)
            elif kind == 2:
                lYii += le
            else:
                lYiii += le
                KE += int_lin_ell(u, v, Eu, Ev, k)
                fi = int_lin_ell(u, v, fu, fv, k)
                fsi = int_lin_ell(u, v, fsu, fsv, k)
                if side == 0:
                    lYb += le; fb += fi; fsb += fsi
                    Yb.append((u, v, fu, fv, fsu, fsv))
                else:
                    lYt += le; ft += fi; fst += fsi
        self.lH, self.lYi, self.lYii, self.lYiii = lH, lYi, lYii, lYiii
        self.lYb, self.lYt, self.KE = lYb, lYt, KE
        self.fb, self.ft, self.fsb, self.fst = fb, ft, fsb, fst
        self.omYi, self.lenYi = omYi, lenYi
        self.Yb = Yb
        # G_Z = Phi(Z) cap Y_b for tilted Z
        self.GZ = {}
        if Yb:
            ya_arr = np.array([t[0] for t in Yb]); yb_arr = np.array([t[1] for t in Yb])
            for s in P.S:
                if s.a <= 1e-15:
                    continue
                r1 = s.sa * s.ca
                G = []
                for (pa, pb) in ((s.v, s.v + r1), (s.yhi - r1, s.yhi)):
                    lo = np.maximum(ya_arr, pa); hi = np.minimum(yb_arr, pb)
                    ok = hi - lo > 1e-12
                    for a_, b_ in zip(lo[ok], hi[ok]):
                        G.append((float(a_), float(b_)))
                if G:
                    self.GZ[s.i] = G

    # --- scale-dependent pieces ---
    def Zset(self, s):
        lo = (1 - self.prm.eps) * s
        out = []
        for j, G in self.GZ.items():
            for (a, b) in G:
                if a < s and b > lo:
                    out.append(j); break
        return out

    def int_f_window(self, s, star=False):
        lo = (1 - self.prm.eps) * s; hi = s
        tot = 0.0
        for (u, v, fu, fv, fsu, fsv) in self.Yb:
            a = max(u, lo); b = min(v, hi)
            if b <= a:
                continue
            A, Bv = (fsu, fsv) if star else (fu, fv)
            fa = A + (Bv - A) * (a - u) / (v - u); fb_ = A + (Bv - A) * (b - u) / (v - u)
            tot += 0.5 * (fa + fb_) * (b - a)
        return tot

    def Gb(self, s, star=False):
        prm = self.prm
        return ((1 - prm.eps) * s) ** 0.75 / prm.c1 * self.int_f_window(s, star)

    def short_chord_mass(self, s, Z):
        """sum_Z int_{Phi(Z) cap Y_b cap W_s} c_Z  (the middle term of Theorem 9.2(a))."""
        lo = (1 - self.prm.eps) * s
        tot = 0.0
        for j in Z:
            sq = self.P.S[j]
            r1 = sq.sa * sq.ca
            for (a, b) in self.GZ[j]:
                a2 = max(a, lo); b2 = min(b, s)
                if b2 <= a2:
                    continue
                # c_Z linear on ramps: lower ramp c = (y - v)/r1 ; upper c = (yhi - y)/r1
                if b2 <= sq.v + r1 + 1e-12:
                    ca_ = (a2 - sq.v) / r1; cb_ = (b2 - sq.v) / r1
                else:
                    ca_ = (sq.yhi - a2) / r1; cb_ = (sq.yhi - b2) / r1
                tot += 0.5 * (ca_ + cb_) * (b2 - a2)
        return tot

    def Zint(self, wfun):
        """int_{y0}^{y1} |Z(s)| wfun(s) ds, exact up to Gauss-Legendre on each s-interval."""
        prm = self.prm; e = prm.eps
        tot = 0.0
        for j, G in self.GZ.items():
            ivs = []
            for (a, b) in G:
                lo = max(a, prm.y0); hi = min(b / (1 - e), prm.y1)
                if hi > lo:
                    ivs.append((lo, hi))
            ivs.sort(); mer = []
            for (a, b) in ivs:
                if mer and a <= mer[-1][1]:
                    mer[-1][1] = max(mer[-1][1], b)
                else:
                    mer.append([a, b])
            for (a, b) in mer:
                tot += gl(wfun, a, b)
        return tot

    def Mb(self, wfun, star=False):
        """int_{y0}^{y1} G_b(s) wfun(s) ds with breakpoints where window ends cross Y_b ends."""
        prm = self.prm; e = prm.eps
        pts = {prm.y0, prm.y1}
        for (u, v, *_r) in self.Yb:
            for yv in (u, v):
                for sv in (yv, yv / (1 - e)):
                    if prm.y0 < sv < prm.y1:
                        pts.add(sv)
        pts = sorted(pts)
        tot = 0.0
        for i in range(len(pts) - 1):
            a, b = pts[i], pts[i + 1]
            f = lambda ss: np.array([self.Gb(float(x), star) for x in ss]) * wfun(ss)
            tot += gl(f, a, b)
        return tot

    def Mb_swap(self, star=False):
        """Tonelli-swapped exact value: c0(1-e)^{3/4}/(2c1) int_Yb f(y) [int_y^{y/(1-e)} s^{-3/4}/(e s+bW) ds] dy."""
        prm = self.prm; e = prm.eps
        inner = lambda y: gl(lambda s: s ** -0.75 / (e * s + prm.bW), y, y / (1 - e))
        tot = 0.0
        for (u, v, fu, fv, fsu, fsv) in self.Yb:
            A, Bv = (fsu, fsv) if star else (fu, fv)
            f = lambda yy: np.array([(A + (Bv - A) * (y - u) / (v - u)) * inner(float(y)) for y in yy])
            tot += gl(f, u, v)
        return prm.c0 * (1 - e) ** 0.75 / (2 * prm.c1) * tot


# ------------------------------------------------------------------------------------
#  per-family flow data (independent of c1, eps, omega0)
# ------------------------------------------------------------------------------------
class FamilyFlows:
    def __init__(self, P, k, delta, c0, y0, s_grid, tlimit=60.0, nov=500, orderB=False):
        t0 = time.time()
        self.P = P; self.k = k; self.delta = delta; self.c0 = c0; self.y0 = y0
        self.y1 = k / 2.0 - 3.0; self.hmax = k / 2.0 - 1.0
        self.amax = c0 / math.sqrt(y0)
        F0 = Flow(P, delta, self.amax, self.hmax, orderB=orderB)
        F0.run(tlimit)
        self.F0 = F0
        segs, terms, contacts = walk(F0, self.amax, None)
        self.LE0 = LineEval(P, segs, terms)
        meas = np.zeros(6)
        for t in terms:
            meas[TYP[t[2]]] += t[1] - t[0]
        self.meas0 = meas
        self.total_measure = float(meas.sum())
        # sup_y Ov_gap^{F0}
        ys = np.linspace(1e-3, self.hmax - 1e-3, nov) + 1e-7 * np.pi
        ovs = []; worst_bd = 0.0; worst_ov = 0.0; worst_p4 = 1e9
        for y in ys:
            r = self.LE0.at(float(y))
            ovs.append(r['Ovgap'])
            worst_bd = max(worst_bd, abs(r['Bd']))
            worst_ov = max(worst_ov, r['Ov'])
            L = r['L']
            worst_p4 = min(worst_p4, L[TT] + L[TD] + L[TE] + L[TM] + L[TW] + r['Ovgap'] - r['Sh'])
        self.supOv0 = float(max(ovs)) if ovs else 0.0
        self.F0_bd = worst_bd; self.F0_ov = worst_ov; self.F0_p4min = worst_p4
        self.Lam0 = float(meas[TD] + meas[TM] + meas[TE]) + self.supOv0
        # T*
        self.tp = tstar_pieces(F0, c0, y0, self.y1)
        self.intA, self.dom = tstar_stats(self.tp, c0, self.y1)
        # scale flows
        self.S = []
        for s in s_grid:
            a_s = c0 / math.sqrt(s)
            sg, tm, ct = walk(F0, a_s, s + 2.0)
            Tset = [(t[0], t[1]) for t in tm if t[2] == 'T']
            N = sum(b - a for (a, b) in Tset)
            n_s, Tstar_iv = tstar_le(self.tp, s)
            inter = inter_measure(Tset, Tstar_iv)
            symd = N + n_s - 2 * inter
            LE = LineEval(P, sg, tm)
            LW = float(LE.losses(s + 2.0)[TW])
            Lall = LE.losses(s + 2.0)
            contacted = set(c[0] for c in ct if c[2] - c[1] > 1e-13)
            self.S.append(dict(s=float(s), N=float(N), n=float(n_s), symd=float(symd), LW=LW,
                               LE=LE, contacted=contacted,
                               LD=float(Lall[TD]), LEE=float(Lall[TE]), LM=float(Lall[TM])))
        self.time = time.time() - t0
        self.anom = list(F0.anom)
        self.dropped = F0.dropped
        self.npieces = len(F0.recs)
        self.nR1multi = F0.nR1multi


# ------------------------------------------------------------------------------------
#  links for one family and one parameter variant
# ------------------------------------------------------------------------------------
def family_links(FF, LS, prm, nypts=24):
    k = prm.k; c0 = prm.c0; e = prm.eps
    out = dict()
    alpha = prm.alpha; dal = prm.dalpha
    wfun = lambda ss: dal(ss) / (e * ss + prm.bW)
    # ---- link 1 (P7, Theorem 9.1) per s ----
    per = []
    for rec in FF.S:
        s = rec['s']; N = rec['N']
        Z = LS.Zset(s)
        nZ = len(Z)
        Vlo = (1 - e) * s - prm.h; Vhi = s + prm.h
        Vlen = Vhi - Vlo
        wallb = 2 * (s + 2) * math.tan(alpha(s))
        Lstar_paper = FF.Lam0 + wallb
        Lstar_meas = FF.Lam0 + rec['LW']
        rhs_paper = Vlen * (N + Lstar_paper)
        rhs_meas = Vlen * (N + Lstar_meas)
        d = dict(s=s, nZ=nZ, N=N, n=rec['n'], symd=rec['symd'], LW=rec['LW'], wallb=wallb,
                 r_paper=nZ / rhs_paper if rhs_paper > 0 else (0.0 if nZ == 0 else math.inf),
                 r_meas=nZ / rhs_meas if rhs_meas > 0 else (0.0 if nZ == 0 else math.inf))
        # tighter (diagnostic) variants: drop b_W ; drop wall
        rb = e * s * (N + Lstar_meas)
        d['r_nobW'] = nZ / rb if rb > 0 else (0.0 if nZ == 0 else math.inf)
        # per-Z conditions
        p2v = 0; outV = 0; tiltv = 0
        bh = prm.betahat((1 - e) * s)
        for j in Z:
            sq = LS.P.S[j]
            if j in rec['contacted']:
                p2v += 1
            if sq.v < Vlo - 1e-12 or sq.yhi > Vhi + 1e-12:
                outV += 1
            if not (sq.a < bh):
                tiltv += 1
        d.update(p2v=p2v, outV=outV, tiltv=tiltv)
        # lower bound (Theorem 9.2(a)) at this s
        If = LS.int_f_window(s)
        mid = LS.short_chord_mass(s, Z) if nZ else 0.0
        d['If'] = If; d['mid'] = mid; d['Gb'] = LS.Gb(s)
        d['lb1'] = If / mid if mid > 0 else (0.0 if If <= 1e-15 else math.inf)
        d['lb2'] = mid / (nZ * bh) if nZ else 0.0
        d['lbZ'] = d['Gb'] / nZ if nZ else (0.0 if d['Gb'] <= 1e-15 else math.inf)
        # sharp form and pointwise checks (only when Z(s) nonempty, plus a few others)
        if nZ or rec is FF.S[0]:
            LE = rec['LE']
            ya = max(Vlo, 1e-9); yb = min(Vhi, s + 2.0 - 1e-9)
            Lint = LE.int_losses(ya, yb, [TT, TD, TE, TM, TW])
            ys = ya + (np.arange(nypts) + 0.5) / nypts * (yb - ya) + 1e-7 * math.e
            ovint = 0.0; shmin = 1e9; p4min = 1e9; bdmax = 0.0; ovmax = 0.0; ltmax = 0.0; lammax = 0.0
            Zarr = np.array(Z, dtype=int)
            shint = 0.0
            for y in ys:
                r = LE.at(float(y))
                ovint += r['Ovgap'] * (yb - ya) / nypts
                Lv = r['L']
                lam = Lv[TD] + Lv[TE] + Lv[TM] + Lv[TW] + r['Ovgap']
                lammax = max(lammax, lam)
                ltmax = max(ltmax, Lv[TT] - N)
                zc = float(r['c'][Zarr].sum()) if nZ else 0.0
                shmin = min(shmin, r['Sh'] - zc)
                p4min = min(p4min, Lv[TT] + lam - r['Sh'])
                bdmax = max(bdmax, abs(r['Bd'])); ovmax = max(ovmax, r['Ov'])
                shint += r['Sh'] * (yb - ya) / nypts
            sharp = Lint + ovint
            d.update(sharp=sharp, r_sharp=nZ / sharp if sharp > 0 else (0.0 if nZ == 0 else math.inf),
                     step3min=shmin, p4min=p4min, bdmax=bdmax, ovmax=ovmax, LTminusN=ltmax,
                     lam_sup_s=lammax, lam_bound=FF.Lam0 + rec['LW'],
                     r_shadow=nZ / shint if shint > 0 else (0.0 if nZ == 0 else math.inf))
        per.append(d)
    out['per_s'] = per
    # ---- link 2: inclusion and Tonelli I ----
    out['symd_max'] = max(d['symd'] for d in per)
    out['N_vs_n_max'] = max(abs(d['N'] - d['n']) for d in per)
    sg = [d['s'] for d in per]; Ns = [d['N'] for d in per]
    mono = all(Ns[i] <= Ns[i + 1] + 1e-12 for i in range(len(Ns) - 1))
    out['N_monotone'] = mono
    # grid upper Riemann sum of int N |alpha'| (N nondecreasing)
    up = 0.0
    for i in range(len(sg) - 1):
        up += Ns[i + 1] * (alpha(sg[i]) - alpha(sg[i + 1]))
    exact_n = FF.intA - alpha(prm.y1) * FF.dom      # = int_{y0}^{y1} n(s)|alpha'(s)| ds exactly
    out['T1_exact_lhs'] = exact_n; out['T1_grid_upper'] = up; out['intA'] = FF.intA; out['dom'] = FF.dom
    out['T1_ratio_exact'] = exact_n / FF.intA if FF.intA > 0 else 0.0
    out['T1_ratio_grid'] = up / FF.intA if FF.intA > 0 else (0.0 if up <= 0 else math.inf)
    # ---- link 3 (Tonelli II) and the integrated per-family chain ----
    LHS3 = prm.c0 * prm.invF / (2 * prm.c1) * LS.fb
    Mb = LS.Mb(wfun); Mb_sw = LS.Mb_swap()
    Zint = LS.Zint(wfun)
    LHS3s = prm.c0 * prm.invF / (2 * prm.c1) * LS.fsb
    Mbs = LS.Mb(wfun, star=True)
    Lam_int = FF.Lam0 * c0 * (prm.y0 ** -0.5 - prm.y1 ** -0.5)
    OmegaW = gl(lambda ss: 2 * (ss + 2) * np.tan(alpha(ss)) * dal(ss), prm.y0, prm.y1)
    # measured wall: trapezoid on the s-grid of L_W^{F_s}(s+2)|alpha'(s)|
    wv = [d['LW'] * dal(d['s']) for d in per]
    Wall_meas = 0.0
    for i in range(len(sg) - 1):
        Wall_meas += 0.5 * (wv[i] + wv[i + 1]) * (sg[i + 1] - sg[i])
    P7int_paper = exact_n + Lam_int + OmegaW
    P7int_meas = exact_n + Lam_int + Wall_meas
    RHS_paper = FF.intA + float(alpha(prm.y0)) * FF.Lam0 + OmegaW
    out.update(LHS3=LHS3, Mb=Mb, Mb_swap=Mb_sw, Zint=Zint, LHS3star=LHS3s, Mbstar=Mbs,
               Lam_int=Lam_int, OmegaW=OmegaW, Wall_meas=Wall_meas,
               P7int_paper=P7int_paper, P7int_meas=P7int_meas, RHS_paper=RHS_paper)
    def rat(a, b):
        return a / b if b > 0 else (0.0 if a <= 1e-15 else math.inf)
    out['T2_ratio'] = rat(LHS3, Mb)
    out['T2_swap_rel'] = abs(Mb - Mb_sw) / max(Mb, 1e-300) if Mb > 0 else 0.0
    out['T2s_ratio'] = rat(LHS3s, Mbs)
    out['LBint_ratio'] = rat(Mb, Zint)
    out['LBint_ratio_star'] = rat(Mbs, Zint)
    out['P7int_ratio_paper'] = rat(Zint, P7int_paper)
    out['P7int_ratio_meas'] = rat(Zint, P7int_meas)
    out['chain_ratio_paper'] = rat(LHS3, RHS_paper)
    out['chain_ratio_meas'] = rat(LHS3, P7int_meas)
    out['chain_ratio_star_meas'] = rat(LHS3s, P7int_meas)
    out['Lam0'] = FF.Lam0; out['supOv0'] = FF.supOv0
    out['meas0'] = [float(x) for x in FF.meas0]
    return out


def assembly(LSb, LSt_direct, fam_b, fam_t, FFb, FFt, prm, W):
    """link 4: (1-omega0) l(H) <= four contributions; floor/ceiling split."""
    om = prm.omega0
    lH = LSb.lH
    K = 2 * prm.c1 / (prm.c0 * prm.invF)
    Ki_paper = W / (om * prm.y0 ** 0.75)
    Ki_meas = LSb.lYi
    Kii = LSb.lYii
    KE = LSb.KE
    Aprime_W = FFb.intA + FFt.intA
    lam = FFb.Lam0 + FFt.Lam0
    Kcol_paper = K * (Aprime_W + float(prm.alpha(prm.y0)) * lam + fam_b['OmegaW'] + fam_t['OmegaW'])
    Kcol_meas = K * (fam_b['P7int_meas'] + fam_t['P7int_meas'])
    LHS = (1 - om) * lH
    r = dict(LHS4=LHS, Ki_paper=Ki_paper, Ki_meas=Ki_meas, Kii=Kii, KE=KE,
             Kcol_paper=Kcol_paper, Kcol_meas=Kcol_meas,
             ratio4_paper=LHS / (Ki_paper + Kii + KE + Kcol_paper) if LHS > 0 else 0.0,
             ratio4_sharp=LHS / (Ki_meas + Kii + KE + Kcol_meas) if LHS > 0 else 0.0)
    # sub-links
    r['kind_i_ratio'] = Ki_meas / Ki_paper if Ki_paper > 0 else 0.0
    r['iii_split_ratio'] = (1 - om) * LSb.lYiii / (LSb.fb + LSb.ft + KE) if (1 - om) * LSb.lYiii > 0 else 0.0
    r['floor_ratio_meas'] = LSb.fb / (K * fam_b['P7int_meas']) if LSb.fb > 0 else 0.0
    r['ceil_ratio_meas'] = LSt_direct.fb / (K * fam_t['P7int_meas']) if LSt_direct.fb > 0 else 0.0
    r['floor_ratio_paper'] = LSb.fb / (K * fam_b['RHS_paper']) if LSb.fb > 0 else 0.0
    r['ceil_ratio_paper'] = LSt_direct.fb / (K * fam_t['RHS_paper']) if LSt_direct.fb > 0 else 0.0
    r['Yt_reflect_diff'] = abs(LSb.ft - LSt_direct.fb)
    r['lH_reflect_diff'] = abs(LSb.lH - LSt_direct.lH)
    # measured constants
    r['A_ii'] = 2 * prm.c1 * Kii / W
    r['Aprime_meas'] = Aprime_W / W
    r['CLam_meas'] = lam * prm.delta / W
    r['Gamma_meas'] = KE / W
    r['lYiii'] = LSb.lYiii; r['fb'] = LSb.fb; r['ft'] = LSb.ft; r['lYi'] = LSb.lYi; r['lYii'] = Kii
    r['lH'] = lH
    return r
