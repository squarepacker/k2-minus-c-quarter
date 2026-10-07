"""Configuration families (exact rational squares) for the tests of Sections 5 and 7."""
from fractions import Fraction as Fr
import math, random
from tracer import Sq, disjoint, dot

def s_of(phi, den=10**12):
    return Fr(math.tan(phi / 2)).limit_denominator(den)

H = Fr(1, 2)

def sq_vertex_at(idx, s, which, P):
    """square with phase param s whose vertex `which` (BL,BR,TR,TL) is at P."""
    tmp = Sq(0, 0, 0, s)
    a, b = {'BL': (-1, -1), 'BR': (1, -1), 'TR': (1, 1), 'TL': (-1, 1)}[which]
    c = (P[0] - a * H * tmp.u[0] - b * H * tmp.n[0], P[1] - a * H * tmp.u[1] - b * H * tmp.n[1])
    return Sq(idx, c[0], c[1], s)

def sq_lowest_at(idx, s, P):
    return sq_vertex_at(idx, s, 'BL' if s >= 0 else 'BR', P)

def vtx(S, which):
    return S.V[{'BL': 0, 'BR': 1, 'TR': 2, 'TL': 3}[which]]

def reindex(sqs):
    return [Sq(i, S.c[0], S.c[1], S.s) for i, S in enumerate(sqs)]

def valley(theta, delta, gap_frac=Fr(1, 2), eps=None, xc=Fr(2), base=False):
    """L (phi=-theta/2) and R (phi=+theta/2); top corners w (L.TR), v (R.TL) at horizontal distance gap."""
    sL = s_of(-theta / 2); sR = s_of(theta / 2)
    if eps is None: eps = Fr(delta) / 4
    y0 = eps + (Fr(11, 10) if base else 0)
    L = sq_lowest_at(0, sL, (xc, y0))
    w = vtx(L, 'TR')
    tth = math.tan(theta)
    gap = Fr(gap_frac) * Fr(delta) * Fr(tth).limit_denominator(10**15)
    R = sq_vertex_at(1, sR, 'TL', (w[0] + gap, w[1]))
    sqs = [L, R]
    if base:
        sqs = [Sq(0, Fr(3, 2) - Fr(1, 100), Fr(1, 2), 0), Sq(1, Fr(5, 2) + Fr(1, 100), Fr(1, 2), 0)] + sqs
        # base squares with tops at 1 ; valley squares lowest at 1.1+eps -> gap 0.1+eps (live only if delta large)
    return reindex(sqs)

def place_on_bottom(idx, sY, P, tau0):
    """square with phase sY whose bottom edge passes through P at local coordinate tau0."""
    tmp = Sq(0, 0, 0, sY)
    c = (P[0] + H * tmp.n[0] - tau0 * tmp.u[0], P[1] + H * tmp.n[1] - tau0 * tmp.u[1])
    return Sq(idx, c[0], c[1], sY)

def merge_cfg(theta, delta, phiY, tau0, a_frac=Fr(1, 2), b_frac=Fr(1, 2), gap_frac=Fr(1, 2), eps=None):
    sqs = valley(theta, delta, gap_frac=gap_frac, eps=eps)
    L, R = sqs[0], sqs[1]
    v = vtx(R, 'TL')
    cth = dot(R.u, L.u); sth = L.u[0] * R.u[1] - L.u[1] * R.u[0]
    we = Fr(delta) * sth / (cth * cth)
    P = (v[0] + a_frac * we * R.u[0] + b_frac * Fr(delta) * R.n[0], v[1] + a_frac * we * R.u[1] + b_frac * Fr(delta) * R.n[1])
    Y = place_on_bottom(2, s_of(phiY), P, tau0)
    return reindex(sqs + [Y])

# ---------------- float drop placement (then exactified) ----------------
def fverts(cx, cy, phi):
    u = (math.cos(phi), math.sin(phi)); n = (-math.sin(phi), math.cos(phi))
    return [(cx + a * .5 * u[0] + b * .5 * n[0], cy + a * .5 * u[1] + b * .5 * n[1]) for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))], u, n

def fdisj(A, B):
    VA, uA, nA = A; VB, uB, nB = B
    for ax in (uA, nA, uB, nB):
        pa = [v[0] * ax[0] + v[1] * ax[1] for v in VA]; pb = [v[0] * ax[0] + v[1] * ax[1] for v in VB]
        if min(pb) > max(pa) + 1e-13 or min(pa) > max(pb) + 1e-13: return True
    return False

def drop(placed_f, cx, phi, k):
    V, u, n = fverts(cx, 0, phi)
    ymin = min(v[1] for v in V)
    cy_floor = -ymin
    xs = [v[0] for v in V]
    if min(xs) < 0 or max(xs) > k: return None
    hi = k + 2.0
    def free(cy):
        S = fverts(cx, cy, phi)
        return all(fdisj(S, P) for P in placed_f)
    cy = hi; last_free = hi
    while cy > cy_floor:
        c2 = max(cy - 0.02, cy_floor)
        if free(c2): last_free = c2; cy = c2
        else:
            lo_, hi_ = c2, last_free
            for _ in range(60):
                m = (lo_ + hi_) / 2
                if free(m): hi_ = m
                else: lo_ = m
            return hi_
        if c2 == cy_floor: return cy_floor
    return cy_floor

def jam(n, k, delta, tilt_sampler, seed, gapmax_factor=1.5, rowcap=None):
    rnd = random.Random(seed)
    placed_f = []; sqs = []
    tries = 0
    while len(sqs) < n and tries < 20 * n:
        tries += 1
        phi = tilt_sampler(rnd)
        s = s_of(phi, 10**9); phi = 2 * math.atan(float(s))
        cx = rnd.uniform(0.75, k - 0.75)
        cy = drop(placed_f, cx, phi, k)
        if cy is None: continue
        gap = 1e-9 + rnd.random() * gapmax_factor * float(delta)
        cyx = Fr(cy + gap).limit_denominator(10**12)
        cxx = Fr(cx).limit_denominator(10**9)
        S = Sq(len(sqs), cxx, cyx, s)
        if any(v[1] < 0 or v[1] > k or v[0] < 0 or v[0] > k for v in S.V): continue
        if rowcap is not None and float(cyx) > rowcap: continue
        if not all(disjoint(S, T) for T in sqs):
            cyx += Fr(1, 10**10)
            S = Sq(len(sqs), cxx, cyx, s)
            if not all(disjoint(S, T) for T in sqs): continue
        sqs.append(S)
        placed_f.append(fverts(float(cxx), float(cyx), phi))
    return sqs

def _free(placed_f, S):
    cx = sum(v[0] for v in S[0]) / 4; cy = sum(v[1] for v in S[0]) / 4
    for P in placed_f:
        px = sum(v[0] for v in P[0]) / 4; py = sum(v[1] for v in P[0]) / 4
        if abs(px - cx) > 1.5 or abs(py - cy) > 1.5: continue
        if not fdisj(S, P): return False
    return True

def rowjam(k, delta, nrows, tilt_sampler, seed, gapmax=1.5):
    """rows packed side by side (slide left) and stacked (drop), all near-contacts with random gaps in [0, gapmax*delta]."""
    rnd = random.Random(seed)
    placed_f = []; sqs = []
    d = float(delta)
    for r in range(nrows):
        xr = 0.0
        while True:
            phi = tilt_sampler(rnd)
            s = s_of(phi, 10**9); phi = 2 * math.atan(float(s))
            V0, _, _ = fverts(0, 0, phi)
            hw = max(abs(v[0]) for v in V0); hh = max(abs(v[1]) for v in V0)
            cx = xr + hw + 0.3
            if cx + hw > k: break
            # vertical: drop onto floor/row below
            if r == 0:
                cy = hh + rnd.random() * 0.5 * d
            else:
                hi = r * 1.2 + 2.0
                cy = hi
                if not _free(placed_f, fverts(cx, cy, phi)): break
                lo = hh
                # step down until collision
                step = 0.02; c = hi
                while c - step > lo and _free(placed_f, fverts(cx, c - step, phi)): c -= step
                a_, b_ = max(lo, c - step), c
                if _free(placed_f, fverts(cx, a_, phi)): cy = a_
                else:
                    for _ in range(55):
                        m = (a_ + b_) / 2
                        if _free(placed_f, fverts(cx, m, phi)): b_ = m
                        else: a_ = m
                    cy = b_
                cy += 1e-11 + rnd.random() ** 4 * gapmax * d
            # horizontal: slide left until contact
            lo_x = hw
            if _free(placed_f, fverts(lo_x, cy, phi)):
                cxs = lo_x
            else:
                a_, b_ = lo_x, cx
                if not _free(placed_f, fverts(cx, cy, phi)):
                    xr = cx + 0.05; continue
                for _ in range(55):
                    m = (a_ + b_) / 2
                    if _free(placed_f, fverts(m, cy, phi)): b_ = m
                    else: a_ = m
                cxs = b_ + 1e-11 + rnd.random() ** 4 * gapmax * d
            S = Sq(len(sqs), Fr(cxs).limit_denominator(10**12), Fr(cy).limit_denominator(10**12), s)
            ok = all(0 <= v[0] <= k and 0 <= v[1] <= k for v in S.V)
            if ok:
                near = [T for T in sqs if abs(float(T.c[0]) - cxs) < 1.5 and abs(float(T.c[1]) - cy) < 1.5]
                if not all(disjoint(S, T) for T in near):
                    S = Sq(len(sqs), Fr(cxs + 1e-9).limit_denominator(10**12), Fr(cy + 1e-9).limit_denominator(10**12), s)
                    ok = all(disjoint(S, T) for T in near)
            if ok:
                sqs.append(S); placed_f.append(fverts(float(S.c[0]), float(S.c[1]), phi))
                xr = float(max(v[0] for v in S.V))
            else:
                xr = cx
    return sqs

def regime_cluster(k, delta, amax, rnd, nrow=5, ytries=6, ybig=0.3):
    """exact in-regime cluster: a row of squares with alternating tiny tilts (|phi| < amax), valley pairs placed corner to
    corner at the top with gap in (0, delta*tan(theta)), peak pairs at the bottom; then Y squares whose bottom edge passes
    through the double rectangle of a valley (tilt tiny or up to ybig)."""
    delta = Fr(delta)
    sqs = []
    x = Fr(1)
    sign = rnd.choice([1, -1])
    prev = None
    for i in range(nrow):
        phi = sign * rnd.uniform(0.05, 1.0) * amax; sign = -sign
        s = s_of(phi, 10**15)
        eps = delta * Fr(rnd.randint(1, 1000), 10000)
        if prev is None:
            S = sq_lowest_at(0, s, (x, eps))
        else:
            P = prev
            if P.s < 0 and s > 0:      # valley: L=P, R=new ; top corners
                w = vtx(P, 'TR')
                L_, R_ = P, Sq(0, 0, 0, s)
                sth = L_.u[0] * R_.u[1] - L_.u[1] * R_.u[0]; cth = dot(L_.u, R_.u)
                g = delta * sth / cth * Fr(rnd.randint(1, 999), 1000)
                S = sq_vertex_at(0, s, 'TL', (w[0] + g, w[1] + delta * Fr(rnd.randint(-5, 5), 10**8)))
            else:                       # peak (or same sign): bottom corners with a small gap
                b = vtx(P, 'BR')
                S = sq_vertex_at(0, s, 'BL', (b[0] + delta * Fr(rnd.randint(1, 1000), 1000), b[1] + delta * Fr(rnd.randint(-100, 100), 10**4)))
                if min(v[1] for v in S.V) < 0:
                    S = Sq(0, S.c[0], S.c[1] - min(v[1] for v in S.V) + eps, S.s)
        S = Sq(len(sqs), S.c[0], S.c[1], S.s)
        if all(disjoint(S, T) for T in sqs) and all(0 <= v[0] <= k and 0 <= v[1] <= k for v in S.V):
            sqs.append(S); prev = S
        else:
            prev = S if not sqs else sqs[-1]
    # Y squares over valleys
    for _ in range(ytries):
        if len(sqs) < 2: break
        i = rnd.randrange(len(sqs) - 1)
        L_, R_ = sqs[i], sqs[i + 1]
        if not (L_.s < 0 < R_.s): continue
        v = vtx(R_, 'TL')
        sth = L_.u[0] * R_.u[1] - L_.u[1] * R_.u[0]; cth = dot(L_.u, R_.u)
        we = delta * sth / (cth * cth)
        a = Fr(rnd.randint(0, 1000), 1000); bb = Fr(rnd.randint(1, 999), 1000)
        P = (v[0] + a * we * R_.u[0] + bb * delta * R_.n[0], v[1] + a * we * R_.u[1] + bb * delta * R_.n[1])
        phY = rnd.choice([rnd.uniform(-amax, amax), rnd.uniform(-ybig, ybig), float(R_.phi), float(L_.phi)])
        tau0 = Fr(rnd.randint(-499, 499), 1000)
        Y = place_on_bottom(len(sqs), s_of(phY, 10**15), P, tau0)
        if all(disjoint(Y, T) for T in sqs) and all(0 <= vv[0] <= k and 0 <= vv[1] <= k for vv in Y.V):
            sqs.append(Y)
    return sqs

def rows_cfg(k, delta, tilts_rows, seed, spacing=None):
    """rows of squares: row r uses tilt pattern tilts_rows[r] (cycled), dropped left to right."""
    rnd = random.Random(seed)
    placed_f = []; sqs = []
    for pat in tilts_rows:
        x = 0.62
        i = 0
        while x < k - 0.62:
            phi = pat[i % len(pat)]; i += 1
            s = s_of(phi, 10**9); phi = 2 * math.atan(float(s))
            cy = drop(placed_f, x, phi, k)
            if cy is not None:
                gap = 1e-9 + rnd.random() * 1.5 * float(delta)
                S = Sq(len(sqs), Fr(x).limit_denominator(10**9), Fr(cy + gap).limit_denominator(10**12), s)
                if all(0 <= v[0] <= k and 0 <= v[1] <= k for v in S.V) and all(disjoint(S, T) for T in sqs):
                    sqs.append(S); placed_f.append(fverts(float(S.c[0]), float(S.c[1]), phi))
            x += (spacing or 1.0) + abs(math.sin(phi)) + 1e-3 + rnd.random() * 2 * float(delta)
    return sqs

def mechA(k, delta, phi, ncol, seed):
    """staggered tilted columns: each square dropped at x offset alternating, all tilt phi."""
    rnd = random.Random(seed)
    placed_f = []; sqs = []
    s = s_of(phi, 10**9); phif = 2 * math.atan(float(s))
    for c in range(ncol):
        x0 = 0.9 + 1.5 * c
        for j in range(int(k / 2)):
            x = x0 + (0.25 if j % 2 else 0.0)
            cy = drop(placed_f, x, phif, k)
            if cy is None: continue
            gap = 1e-9 + rnd.random() * 1.5 * float(delta)
            S = Sq(len(sqs), Fr(x).limit_denominator(10**9), Fr(cy + gap).limit_denominator(10**12), s)
            if all(0 <= v[0] <= k and 0 <= v[1] <= k for v in S.V) and all(disjoint(S, T) for T in sqs):
                sqs.append(S); placed_f.append(fverts(float(S.c[0]), float(S.c[1]), phif))
    return sqs

def strip_cfg(k, delta, phi, b, seed):
    """1 x b rigid block tilted by phi (L-shape strip) resting with lowest corner near floor, plus drops."""
    s = s_of(phi, 10**9)
    tmp = Sq(0, 0, 0, s)
    g = Fr(1, 10**7)
    base = sq_lowest_at(0, s, (Fr(1), Fr(delta) / 3))
    sqs = [base]
    for j in range(1, b):
        c = (base.c[0] + j * (1 + g) * tmp.u[0], base.c[1] + j * (1 + g) * tmp.u[1])
        sqs.append(Sq(j, c[0], c[1], s))
    sqs = [S for S in sqs if all(0 <= v[0] <= k and 0 <= v[1] <= k for v in S.V)]
    rnd = random.Random(seed)
    placed_f = [fverts(float(S.c[0]), float(S.c[1]), 2 * math.atan(float(S.s))) for S in sqs]
    for _ in range(4 * k):
        ph = rnd.choice([0.0, 0.0, phi, -phi, rnd.uniform(-0.05, 0.05)])
        ss = s_of(ph, 10**9); ph = 2 * math.atan(float(ss))
        x = rnd.uniform(0.7, k - 0.7)
        cy = drop(placed_f, x, ph, k)
        if cy is None or cy > k / 2 + 0.5: continue
        S = Sq(len(sqs), Fr(x).limit_denominator(10**9), Fr(cy + 1e-9 + rnd.random() * 1.5 * float(delta)).limit_denominator(10**12), ss)
        if all(0 <= v[0] <= k and 0 <= v[1] <= k for v in S.V) and all(disjoint(S, T) for T in sqs):
            sqs.append(S); placed_f.append(fverts(float(S.c[0]), float(S.c[1]), ph))
    return reindex(sqs)

