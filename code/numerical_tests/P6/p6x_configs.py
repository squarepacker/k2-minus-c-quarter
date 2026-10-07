# p6x_configs.py -- configuration generators (exact rational squares, Pythagorean rotations).
import math
import random
CFG = dict(rotden=10 ** 4, posden=10 ** 7)
from fractions import Fraction as Fr
from p6x_core import Sq, rot_for_angle, closed_disjoint, canon


def offsets(u):
    ux, uy = float(u[0]), float(u[1])
    nx, ny = -uy, ux
    return [(-0.5 * nx - 0.5 * ux, -0.5 * ny - 0.5 * uy), (-0.5 * nx + 0.5 * ux, -0.5 * ny + 0.5 * uy),
            (0.5 * nx + 0.5 * ux, 0.5 * ny + 0.5 * uy), (0.5 * nx - 0.5 * ux, 0.5 * ny - 0.5 * uy)]


def _vert_range(poly, x):
    """y-range of convex polygon (float verts) on vertical line x, or None"""
    ys = []
    m = len(poly)
    for i in range(m):
        (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % m]
        if (x1 - x) * (x2 - x) <= 0 and x1 != x2:
            t = (x - x1) / (x2 - x1)
            ys.append(y1 + t * (y2 - y1))
        elif x1 == x2 == x:
            ys += [y1, y2]
    if not ys:
        return None
    return min(ys), max(ys)


def rest_height(placed_polys, off, cx):
    """lowest cy so that polygon (cx+off, cy+off) is above floor and above every placed polygon
    it overlaps horizontally (Tetris drop from +inf)."""
    cy = -min(o[1] for o in off)
    Z0 = [(cx + o[0], o[1]) for o in off]
    zx0, zx1 = min(p[0] for p in Z0), max(p[0] for p in Z0)
    for P in placed_polys:
        px0, px1 = min(p[0] for p in P), max(p[0] for p in P)
        if px1 < zx0 or zx1 < px0:
            continue
        for v in P:
            r = _vert_range(Z0, v[0])
            if r is not None:
                cy = max(cy, v[1] - r[0])
        for w in Z0:
            r = _vert_range(P, w[0])
            if r is not None:
                cy = max(cy, r[1] - w[1])
    return cy


class Builder:
    def __init__(self, k, maxden=None):
        maxden = maxden or CFG['posden']
        self.k = k
        self.sqs = []
        self.polys = []
        self.maxden = maxden

    def try_place(self, cx, u, gap, ymax=None, cy=None):
        off = offsets(u)
        if cx + min(o[0] for o in off) < 0 or cx + max(o[0] for o in off) > self.k:
            return None
        if cy is None:
            cy = rest_height(self.polys, off, cx) + gap
        if ymax is not None and cy + max(o[1] for o in off) > ymax:
            return None
        if cy + max(o[1] for o in off) > self.k:
            return None
        for attempt in range(8):
            c = (Fr(cx).limit_denominator(self.maxden), Fr(cy).limit_denominator(self.maxden))
            S = Sq(len(self.sqs), c, u)
            ok = all(0 <= v[0] <= self.k and 0 <= v[1] <= self.k for v in S.V)
            if ok:
                for T in self.sqs:
                    if T.bb[2] < S.bb[0] - 1e-6 or S.bb[2] < T.bb[0] - 1e-6 or T.bb[3] < S.bb[1] - 1e-6 or S.bb[3] < T.bb[1] - 1e-6:
                        continue
                    if not closed_disjoint(S.V, T.V):
                        ok = False
                        break
            if ok:
                self.sqs.append(S)
                self.polys.append([(float(v[0]), float(v[1])) for v in S.V])
                return S
            cy += max(gap, 1e-9) * (2 ** attempt) + 1e-12
            if ymax is not None and cy + max(o[1] for o in off) > ymax:
                return None
        return None


def rand_tilt(rng, amax, mix):
    r = rng.random()
    if r < mix[0]:
        return 0.0
    r -= mix[0]
    if r < mix[1]:
        return rng.uniform(-1.6 * amax, 1.6 * amax)
    r -= mix[1]
    if r < mix[2]:
        return rng.uniform(-0.6, 0.6)
    s = rng.choice((-1, 1))
    return s * (math.pi / 4 + rng.uniform(-0.03, 0.03))


def rand_gap(rng, delta):
    r = rng.random()
    if r < 0.35:
        return delta * rng.uniform(0.0005, 0.05)
    if r < 0.75:
        return delta * rng.uniform(0.05, 0.6)
    if r < 0.92:
        return delta * rng.uniform(0.6, 1.0)
    return delta * rng.uniform(1.0, 3.0)


def pile(rng, k, n, delta, amax, mix=(0.15, 0.6, 0.15, 0.1), ymax=None, lattice=False, maxden=None):
    B = Builder(k, maxden)
    fails = 0
    i = 0
    while len(B.sqs) < n and fails < 6 * n:
        if lattice:
            col = i % (k - 1) if k > 1 else 0
            cx = 0.55 + col * (k - 1.1) / max(1, k - 2) + rng.uniform(-0.04, 0.04)
            i += 1
        else:
            cx = rng.uniform(0.5, k - 0.5)
        phi = rand_tilt(rng, amax, mix)
        u = rot_for_angle(phi, CFG['rotden'])
        if rng.random() < 0.12:        # snap to a wall (vertex/side touching x=0 or x=k)
            hw = max(abs(o[0]) for o in offsets(u))
            cx = hw if rng.random() < 0.5 else k - hw
            cx = float(Fr(cx).limit_denominator(10 ** 12))
            if cx - hw < 0:
                cx += 1e-12
            if cx + hw > k:
                cx -= 1e-12
        if B.try_place(cx, u, rand_gap(rng, delta), ymax) is None:
            fails += 1
    return B


def two_sided(rng, k, n, delta, amax, mix=(0.15, 0.6, 0.15, 0.1), lattice=False):
    """bottom pile up to k/2-0.2 and an independent top pile (built reflected)."""
    lim = k / 2 - 0.2
    B1 = pile(rng, k, n, delta, amax, mix, ymax=lim, lattice=lattice)
    B2 = pile(rng, k, n, delta, amax, mix, ymax=lim, lattice=lattice)
    sqs = list(B1.sqs)
    for S in B2.sqs:
        sqs.append(Sq(len(sqs), (S.c[0], k - S.c[1]), (S.u[0], -S.u[1])))
    return sqs


def columns(rng, k, delta, amax, tilt, alternate=False, height=None):
    """columns of stacked tilted squares, each column tilt +-tilt (or alternating); kind 'colA'."""
    B = Builder(k)
    ymax = height if height is not None else k / 2 - 0.2
    ncol = k - 1
    for j in range(ncol):
        cx = 0.6 + j * (k - 1.2) / max(1, ncol - 1)
        for r in range(int(ymax) + 1):
            s = (1 if (r % 2 == 0 or not alternate) else -1) * (1 if j % 2 == 0 else -1)
            u = rot_for_angle(s * tilt * rng.uniform(0.7, 1.0), CFG['rotden'])
            B.try_place(cx + rng.uniform(-0.02, 0.02), u, delta * rng.uniform(0.01, 0.4), ymax)
    return B.sqs


def valleys(rng, k, delta, a, rows):
    """zigzag rows (tilts -a,+a alternating => V valleys) alternating with axis rows (fan/merge)."""
    B = Builder(k)
    ymax = k / 2 - 0.2
    for r in range(rows):
        if r % 2 == 0:
            n = k - 1
            for j in range(n):
                cx = 0.62 + j * (k - 1.24) / max(1, n - 1)
                s = -1 if j % 2 == 0 else 1
                u = rot_for_angle(s * a, CFG['rotden'])
                B.try_place(cx, u, delta * rng.uniform(0.001, 0.05), ymax)
        else:
            n = k - 1
            for j in range(n):
                cx = 0.5 + 0.5 + j * (k - 2.0) / max(1, n - 1) + rng.uniform(-0.05, 0.05)
                u = rot_for_angle(rng.uniform(-0.3, 0.3) * a, CFG['rotden'])
                B.try_place(cx, u, delta * rng.uniform(0.001, 0.1), ymax)
    return B.sqs


def valleys2(rng, k, delta, amax, layers=3, ytilt=None):
    """explicit merge machines: pairs X1 (phi=-a), X2 (phi=+a) with top corners eps apart at the
    same height, cap Y deposited on both corners plus gap g<delta; eps < 2 g tan a so the two exit
    bands cross on bot(Y) (positive-measure R1 merges). Y tilt: 0, small, or >= alpha (T)."""
    B = Builder(k)
    ymax = k / 2 - 0.2
    for L in range(layers):
        x = 0.05 + rng.uniform(0, 0.3)
        while x < k - 2.4:
            a = rng.uniform(0.4, 0.97) * amax
            g = delta * rng.uniform(0.3, 0.95)
            eps = rng.uniform(0.1, 0.95) * 2 * g * math.tan(a)
            w = math.sin(a) + math.cos(a)
            cx1 = x + w / 2
            cx2 = cx1 + w + eps
            u1 = rot_for_angle(-a, CFG['rotden'])
            u2 = rot_for_angle(a, CFG['rotden'])
            off1, off2 = offsets(u1), offsets(u2)
            base = max(rest_height(B.polys, off1, cx1), rest_height(B.polys, off2, cx2)) + delta * rng.uniform(0.001, 0.05)
            S1 = B.try_place(cx1, u1, 0, ymax, cy=base)
            S2 = B.try_place(cx2, u2, 0, ymax, cy=base) if S1 is not None else None
            if S2 is not None:
                r = rng.random()
                if ytilt is not None:
                    phiY = ytilt
                elif r < 0.4:
                    phiY = 0.0
                elif r < 0.7:
                    phiY = rng.uniform(-0.9, 0.9) * amax
                elif r < 0.9:
                    phiY = rng.choice((-1, 1)) * rng.uniform(amax, 0.6)
                else:
                    phiY = rng.choice((-1, 1)) * (math.pi / 4 + rng.uniform(-0.02, 0.02))
                uY = rot_for_angle(phiY, CFG['rotden'])
                B.try_place((cx1 + cx2) / 2 + rng.uniform(-0.1, 0.1), uY, g, ymax)
            x = cx2 + w / 2 + rng.uniform(0.05, 0.6)
    for _ in range(2 * k):
        u = rot_for_angle(rand_tilt(rng, amax, (0.2, 0.6, 0.15, 0.05)), CFG['rotden'])
        B.try_place(rng.uniform(0.5, k - 0.5), u, rand_gap(rng, delta) * 0.3, ymax)
    return B.sqs


def bridge(rng, k, delta, amax):
    """attack on Lemmas 8.6/8.7: X bridging two supports A,B over a hole (so I_e(X,Y) gets holes: middle of
    bot(X) unreachable / paths with large g die), Y above X with random overlap shift and tilt,
    and a side square Z (often steeply tilted) next to Y whose corner dips toward the X-Y gap."""
    B = Builder(k)
    ymax = k / 2 - 0.2
    x = rng.uniform(0.0, 0.3)
    while x < k - 2.3:
        tA = rng.uniform(-1.5, 1.5) * amax
        tB = rng.uniform(-1.5, 1.5) * amax
        hole = rng.uniform(0.02, 0.85)
        cA = x + 0.55
        cB = cA + 1.05 + hole
        SA = B.try_place(cA, rot_for_angle(tA, CFG['rotden']), delta * rng.uniform(0.0005, 0.6), ymax)
        SB = B.try_place(cB, rot_for_angle(tB, CFG['rotden']), delta * rng.uniform(0.0005, 0.6), ymax)
        if SA is not None and SB is not None:
            phX = rng.uniform(-0.98, 0.98) * amax
            cX = (cA + cB) / 2 + rng.uniform(-0.2, 0.2)
            SX = B.try_place(cX, rot_for_angle(phX, CFG['rotden']), delta * rng.uniform(0.0005, 0.7), ymax)
            if SX is not None:
                r = rng.random()
                if r < 0.3:
                    phY = rng.uniform(-0.98, 0.98) * amax
                elif r < 0.6:
                    phY = rng.choice((-1, 1)) * rng.uniform(amax, 0.7)
                elif r < 0.8:
                    phY = rng.choice((-1, 1)) * (math.pi / 4 + rng.uniform(-0.03, 0.03))
                else:
                    phY = 0.0
                sh = rng.choice((-1, 1)) * rng.choice((rng.uniform(0, 0.5), rng.uniform(0.85, 0.999)))
                SY = B.try_place(float(SX.c[0]) + sh, rot_for_angle(phY, CFG['rotden']),
                                 delta * rng.uniform(0.0005, 0.97), ymax)
                if SY is not None and rng.random() < 0.8:
                    side = -1 if sh > 0 else 1
                    phZ = rng.choice((rng.uniform(0.2, 0.7), math.pi / 4 + rng.uniform(-0.03, 0.03),
                                      rng.uniform(-1, 1) * amax)) * rng.choice((-1, 1))
                    uZ = rot_for_angle(phZ, CFG['rotden'])
                    wZ = (abs(math.cos(phZ)) + abs(math.sin(phZ))) / 2
                    wY = (abs(float(SY.u[0])) + abs(float(SY.u[1]))) / 2
                    cZ = float(SY.c[0]) + side * (wY + wZ + rng.uniform(0.001, 0.05))
                    B.try_place(cZ, uZ, delta * rng.uniform(0.0005, 0.5), ymax)
        x = cB + 0.6 + rng.uniform(0.05, 0.4)
    for _ in range(2 * k):
        u = rot_for_angle(rand_tilt(rng, amax, (0.2, 0.6, 0.15, 0.05)), CFG['rotden'])
        B.try_place(rng.uniform(0.5, k - 0.5), u, rand_gap(rng, delta) * 0.5, ymax)
    return B.sqs


def fan(rng, k, delta, amax, rows):
    B = Builder(k)
    ymax = k / 2 - 0.2
    for r in range(rows):
        n = k - 1
        for j in range(n):
            cx = 0.6 + j * (k - 1.2) / max(1, n - 1) + rng.uniform(-0.03, 0.03)
            phi = (j - (n - 1) / 2) / max(1, (n - 1) / 2) * 1.5 * amax * (1 if r % 2 == 0 else -1)
            u = rot_for_angle(phi, CFG['rotden'])
            B.try_place(cx, u, delta * rng.uniform(0.001, 0.3), ymax)
    return B.sqs


def lstrip(rng, k, delta, amax):
    """L-shape: a horizontal strip of slightly tilted squares on the floor and a vertical strip
    along the left wall, plus random fill."""
    B = Builder(k)
    ymax = k / 2 - 0.2
    for j in range(k - 1):
        u = rot_for_angle(rng.choice((1, -1)) * rng.uniform(0.3, 1.4) * amax, CFG['rotden'])
        B.try_place(1.0 + j * (k - 1.6) / max(1, k - 2), u, delta * rng.uniform(0.0, 0.3), ymax)
    for r in range(int(ymax)):
        u = rot_for_angle(rng.uniform(-1.4, 1.4) * amax, CFG['rotden'])
        B.try_place(0.52, u, delta * rng.uniform(0.0, 0.3), ymax)
    for _ in range(3 * k):
        u = rot_for_angle(rand_tilt(rng, amax, (0.2, 0.6, 0.15, 0.05)), CFG['rotden'])
        B.try_place(rng.uniform(0.5, k - 0.5), u, rand_gap(rng, delta), ymax)
    return B.sqs
