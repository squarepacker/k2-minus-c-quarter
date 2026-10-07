"""Configuration generators (exact rational squares) for the tests of Section 6. Written from scratch.
The argument paper=True of a family selects the constants of the paper (delta = 1e-5 and
alpha_F = 2 atan(7.5e-7) ~ 1.5e-6, Section 3), with a small container (k = 6 to 8);
paper=False gives scaled analogues with larger delta and alpha_F."""
import math
from fractions import Fraction as Fr
from core import Sq, disjoint, Cfg, norm_t


def rat(v, den):
    return Fr(int(round(v * den)), den)


def t_of(angle, den=10 ** 9):
    return rat(math.tan(angle / 2.0), den)


class Builder:
    def __init__(self, k):
        self.k = Fr(k)
        self.sqs = []

    def ok(self, S):
        if S.xmin < 0 or S.ymin < 0 or S.xmax > self.k or S.ymax > self.k:
            return False
        for B in self.sqs:
            if not disjoint(S, B):
                return False
        return True

    def add(self, cx, cy, t):
        S = Sq(cx, cy, t, len(self.sqs))
        if self.ok(S):
            self.sqs.append(S)
            return S
        return None

    def drop(self, cx, t, gap, floor_gap=Fr(0)):
        """drop a square of phase t at abscissa cx from above; rest on floor (+floor_gap) or a square (+gap)."""
        cx = Fr(cx)
        S0 = Sq(cx, 0, t)
        best = -S0.ymin
        on = 'floor'
        for B in self.sqs:
            lo, hi = max(S0.xmin, B.xmin), min(S0.xmax, B.xmax)
            if lo > hi:
                continue
            xs = {lo, hi}
            xs.update(v[0] for v in S0.V if lo <= v[0] <= hi)
            xs.update(v[0] for v in B.V if lo <= v[0] <= hi)
            for x in xs:
                vb, vs = B.vchord(x), S0.vchord(x)
                if vb is None or vs is None:
                    continue
                need = vb[1] - vs[0]
                if need > best:
                    best, on = need, B.idx
        cy = best + (Fr(floor_gap) if on == 'floor' else Fr(gap))
        S = Sq(cx, cy, t, len(self.sqs))
        if self.ok(S):
            self.sqs.append(S)
            return S
        return None

    def cfg(self, delta, tF, hF, name):
        return Cfg(self.k, [(S.c[0], S.c[1], S.t) for S in self.sqs], delta, tF, hF, name)


def halfwidth(t):
    t = norm_t(t)
    d = 1 + t * t
    cs, sn = (1 - t * t) / d, 2 * t / d
    return (cs + abs(sn)) / 2


# ------------------------------------------------------------------ scaled families
def fam_valley(rng, paper=False):
    if paper:
        k = rng.choice([6, 8])
        delta = Fr(1, 10 ** 5)
        tF = Fr(75, 10 ** 8)
        hF = Fr(k, 2) - 1
        den = 10 ** 15
        aF = 1.5e-6
        a_lo, a_hi = 0.3 * aF, 0.95 * aF
        w_rng = (1e-12, 4e-12)
        gam = (0.05, 0.4)
    else:
        k = rng.choice([4, 5, 6])
        delta = Fr(rng.choice([8, 15, 30, 50]), 100)
        aF = rng.uniform(0.2, 0.75)
        tF = t_of(aF)
        hF = Fr(rng.choice([Fr(k, 2) - 1, Fr(2 * k - 1, 2)]))
        den = 10 ** 9
        dl = float(delta)
        a_lo, a_hi = 0.2 * dl, min(1.3 * dl, 0.9 * aF)
        w_rng = (0.01 * dl, 0.15 * dl)
        gam = (0.03, 0.45)
    dl = float(delta)
    b = Builder(k)
    nval = rng.randint(2, 4)
    xs = sorted(rng.uniform(1.4, k - 1.4) for _ in range(nval))
    for x0 in xs:
        a1, a2 = rng.uniform(a_lo, a_hi), rng.uniform(a_lo, a_hi)
        t1, t2 = -t_of(a1, den), t_of(a2, den)
        w = rat(rng.uniform(*w_rng), den)
        fg1 = rat(rng.choice([0.0, 0.0, rng.uniform(0, 0.4) * dl]), den)
        fg2 = fg1 if rng.random() < 0.5 else rat(rng.uniform(0, 0.4) * dl, den)
        X0 = rat(x0, den)
        c1 = X0 - halfwidth(t1) - w / 2
        c2 = X0 + halfwidth(t2) + w / 2
        g1 = rat(rng.uniform(*gam) * dl, den)
        b.drop(c1, t1, g1, fg1)
        b.drop(c2, t2, g1, fg2)
        off = rat(rng.uniform(-0.3, 0.3), den)
        tz = rng.choice([Fr(0), rat(rng.uniform(-0.5, 0.5) * float(tF), den)])
        b.drop(X0 + off, tz, rat(rng.uniform(*gam) * dl, den))
        if rng.random() < 0.6:   # second-level valley on top
            a3 = rng.uniform(a_lo, a_hi)
            t3 = rng.choice([-1, 1]) * t_of(a3, den)
            b.drop(X0 + off + rat(rng.uniform(-0.6, 0.6), den), t3, rat(rng.uniform(*gam) * dl, den))
    for _ in range(rng.randint(2, 8)):
        t = rng.choice([Fr(0), t_of(rng.uniform(-1, 1) * (aF * 1.5), den)])
        b.drop(rat(rng.uniform(0.6, k - 0.6), den), t, rat(rng.uniform(*gam) * dl, den),
               rat(rng.choice([0.0, rng.uniform(0, 0.5) * dl]), den))
    return b.cfg(delta, tF, hF, 'valley' + ('_paper' if paper else ''))


def fam_zigzag(rng, paper=False):
    if paper:
        k = 8
        delta = Fr(1, 10 ** 5)
        tF = Fr(75, 10 ** 8)
        hF = Fr(k, 2) - 1
        den = 10 ** 15
        aF = 1.5e-6
        a = rng.uniform(0.3, 0.9) * aF
        w = rat(rng.uniform(1e-12, 3e-12), den)
    else:
        k = rng.choice([5, 6])
        delta = Fr(rng.choice([10, 20, 40]), 100)
        aF = rng.uniform(0.25, 0.7)
        tF = t_of(aF)
        hF = Fr(rng.choice([Fr(k, 2) - 1, Fr(2 * k - 1, 2)]))
        den = 10 ** 9
        a = min(rng.uniform(0.3, 1.2) * float(delta), 0.85 * aF)
        w = rat(rng.uniform(0.01, 0.1) * float(delta), den)
    dl = float(delta)
    b = Builder(k)
    ta = t_of(a, den)
    hw = halfwidth(ta)
    nrow = rng.randint(1, 3)
    fg = rat(rng.uniform(0, 0.3) * dl, den)
    for row in range(nrow):
        n = int((k - 0.2) / float(2 * hw + w))
        x = rat(rng.uniform(0.05, 0.2), den) + hw
        sgn = rng.choice([-1, 1])
        centers = []
        for i in range(n):
            t = sgn * ta
            sgn = -sgn
            S = b.drop(x, t, rat(rng.uniform(0.05, 0.4) * dl, den), fg)
            if S is not None:
                centers.append(x)
            x = x + 2 * hw + w
        # axis-parallel row above the junctions
        for i in range(len(centers) - 1):
            if rng.random() < 0.7:
                xm = (centers[i] + centers[i + 1]) / 2 + rat(rng.uniform(-0.1, 0.1), den)
                b.drop(xm, Fr(0), rat(rng.uniform(0.05, 0.4) * dl, den))
    return b.cfg(delta, tF, hF, 'zigzag' + ('_paper' if paper else ''))


def fam_deaths(rng, paper=False):
    if paper:
        k = rng.choice([6, 8])
        delta = Fr(1, 10 ** 5)
        tF = Fr(75, 10 ** 8)
        hF = Fr(k, 2) - 1
        den = 10 ** 15
        tilts = [lambda: Fr(0), lambda: Fr(0),
                 lambda: rat(rng.uniform(-0.9, 0.9) * 7.5e-7, den),
                 lambda: rat(rng.uniform(-3, 3) * 7.5e-7, den)]
    else:
        k = rng.choice([4, 5, 6])
        delta = Fr(rng.choice([5, 10, 20, 30]), 100)
        aF = rng.uniform(0.1, 0.6)
        tF = t_of(aF)
        hF = Fr(rng.choice([Fr(k, 2) - 1, Fr(2 * k - 1, 2), Fr(k - 1)]))
        den = 10 ** 9
        tilts = [lambda: Fr(0), lambda: Fr(0),
                 lambda: t_of(rng.uniform(-0.9, 0.9) * aF, den),
                 lambda: t_of(rng.uniform(-1.5, 1.5) * aF, den)]
    dl = delta
    fr = [Fr(1, 3), Fr(1, 2), Fr(9, 10), Fr(999999, 10 ** 6), Fr(1), Fr(1000001, 10 ** 6), Fr(1, 10)]
    b = Builder(k)
    # columns at fixed abscissae (towers) with gaps near delta (cumulative)
    ncol = rng.randint(2, int(k) - 1)
    cols = sorted(rng.uniform(0.7, k - 0.7) for _ in range(ncol))
    for cx in cols:
        h = rng.randint(2, 4)
        fgap = rng.choice(fr) * dl / 2
        for _ in range(h):
            g = rng.choice(fr) * dl / 2
            b.drop(rat(cx + rng.uniform(-0.05, 0.05), den), rng.choice(tilts)(), g, fgap)
    # a near-full horizontal row (horizontal gap of width ~ delta above it)
    if rng.random() < 0.6:
        x = Fr(1, 1000)
        while x + 1 < k:
            b.drop(x + Fr(1, 2), Fr(0), rng.choice(fr) * dl / 2, rng.choice(fr) * dl / 2)
            x += 1 + Fr(1, 997)
    for _ in range(rng.randint(2, 6)):
        b.drop(rat(rng.uniform(0.6, k - 0.6), den), rng.choice(tilts)(), rng.choice(fr) * dl / 2,
               rng.choice(fr) * dl / 2)
    return b.cfg(delta, tF, hF, 'deaths' + ('_paper' if paper else ''))


def fam_near45(rng):
    k = rng.choice([4, 5, 6])
    delta = Fr(rng.choice([10, 20, 40, 80]), 100)
    aF = rng.choice([rng.uniform(0.3, 0.7), rng.uniform(0.7, 0.785)])
    tF = t_of(aF)
    hF = Fr(rng.choice([Fr(k, 2) - 1, Fr(2 * k - 1, 2)]))
    den = 10 ** 12
    s2 = math.sqrt(2) - 1
    b = Builder(k)
    menu = [lambda: rat(s2 - rng.choice([1e-12, 1e-9, 1e-6, 1e-3]), den),      # just below 45deg
            lambda: rat(s2 + rng.choice([1e-12, 1e-9, 1e-6, 1e-3]), den),      # just above 45deg
            lambda: -rat(s2 - rng.choice([1e-12, 1e-9, 1e-6]), den),
            lambda: t_of(rng.uniform(-0.98, 0.98) * aF, den),                  # passable, maybe large
            lambda: t_of(rng.choice([-1, 1]) * aF * rng.uniform(0.97, 0.9999), den),
            lambda: Fr(0)]
    for _ in range(rng.randint(8, 20)):
        b.drop(rat(rng.uniform(0.5, k - 0.5), den), rng.choice(menu)(),
               rat(rng.uniform(0.02, 0.6) * float(delta), den),
               rat(rng.choice([0.0, rng.uniform(0, 0.6) * float(delta)]), den))
    return b.cfg(delta, tF, hF, 'near45')


def fam_wall(rng):
    """paths drifting to a wall through many squares tilted the same way (large delta allowed)."""
    k = rng.choice([4, 5, 6])
    delta = Fr(rng.choice([30, 60, 100, 200, 300]), 100)
    aF = rng.uniform(0.3, 0.75)
    tF = t_of(aF)
    hF = Fr(rng.choice([Fr(2 * k - 1, 2), Fr(k - 1), Fr(k, 2) - 1]))
    den = 10 ** 9
    b = Builder(k)
    a = aF * rng.uniform(0.7, 0.99)
    for side in (0, 1):
        if rng.random() < 0.15:
            continue
        t = t_of(a, den) * (1 if side == 0 else -1)
        hw = halfwidth(t)
        eps = rat(rng.choice([0.0, 0.0, rng.uniform(0, 0.05)]), den)
        cx = hw + eps if side == 0 else Fr(k) - hw - eps
        nst = rng.randint(1, 4)
        for i in range(nst):
            S = b.drop(cx, t, rat(rng.uniform(0.01, 0.3) * min(float(delta), 1.0), den),
                       rat(rng.choice([0.0, rng.uniform(0, 0.2)]), den))
            if S is None:
                break
            # next one shifted toward the wall a bit less (staircase) or same column
            if rng.random() < 0.5:
                cx = cx + (rat(rng.uniform(0.0, 0.4), den) if side == 0 else -rat(rng.uniform(0.0, 0.4), den))
        # a second tilted column further from the wall feeding drift paths
        if rng.random() < 0.7:
            cx2 = cx + (1 if side == 0 else -1) * rat(rng.uniform(1.0, 1.6), den)
            for i in range(rng.randint(1, 3)):
                if b.drop(cx2, t, rat(rng.uniform(0.01, 0.3) * min(float(delta), 1.0), den),
                          rat(rng.choice([0.0, rng.uniform(0, 0.2)]), den)) is None:
                    break
    for _ in range(rng.randint(0, 4)):
        b.drop(rat(rng.uniform(1.5, k - 1.5), den), rng.choice([Fr(0), t_of(rng.uniform(-a, a), den)]),
               rat(rng.uniform(0.05, 0.5) * min(float(delta), 1.0), den))
    return b.cfg(delta, tF, hF, 'wall'), a


def fam_floor(rng, paper=False):
    if paper:
        k = rng.choice([6, 8])
        delta = Fr(1, 10 ** 5)
        tF = Fr(75, 10 ** 8)
        hF = Fr(k, 2) - 1
        den = 10 ** 15
        tl = lambda: rng.choice([Fr(0), Fr(0), rat(rng.uniform(-2, 2) * 7.5e-7, den)])
    else:
        k = rng.choice([4, 5, 6])
        delta = Fr(rng.choice([10, 20, 40]), 100)
        aF = rng.uniform(0.1, 0.7)
        tF = t_of(aF)
        hF = Fr(rng.choice([Fr(k, 2) - 1, Fr(2 * k - 1, 2)]))
        den = 10 ** 9
        tl = lambda: rng.choice([Fr(0), Fr(0), t_of(rng.uniform(-1.3, 1.3) * aF, den)])
    dl = float(delta)
    b = Builder(k)
    # squares exactly on the floor (axis-parallel: whole bottom side on y=0; tilted: one vertex on y=0)
    x = rat(rng.choice([0.0, rng.uniform(0, 0.3)]), den)
    while x < k - 1:
        t = tl()
        hw = halfwidth(t)
        S = b.drop(x + hw, t, Fr(0), Fr(0))
        x = x + 2 * hw + rat(rng.choice([1e-9, 0.01, rng.uniform(0.0, 0.5)]) * (1 if not paper else 1e-4), den) \
            + Fr(1, 10 ** 12)
    for _ in range(rng.randint(3, 10)):
        b.drop(rat(rng.uniform(0.5, k - 0.5), den), tl(), rat(rng.uniform(0.01, 0.9) * dl, den),
               Fr(0))
    return b.cfg(delta, tF, hF, 'floor' + ('_paper' if paper else ''))


def fam_jam(rng, paper=False):
    if paper:
        k = rng.choice([6, 8])
        delta = Fr(1, 10 ** 5)
        tF = Fr(75, 10 ** 8)
        hF = Fr(k, 2) - 1
        den = 10 ** 15
        menu = [lambda: Fr(0), lambda: rat(rng.uniform(-1, 1) * 7.5e-7, den),
                lambda: rat(rng.uniform(-4, 4) * 7.5e-7, den), lambda: t_of(rng.uniform(-0.8, 0.8), den)]
        gaps = [lambda: rat(rng.uniform(0.01, 0.5) * 1e-5, den), lambda: rat(rng.uniform(0.9, 1.1) * 1e-5, den)]
        n = rng.randint(15, 35)
    else:
        k = rng.choice([4, 5, 6, 7])
        delta = Fr(rng.choice([5, 10, 20, 40, 70]), 100)
        aF = rng.uniform(0.1, 0.75)
        tF = t_of(aF)
        hF = Fr(rng.choice([Fr(k, 2) - 1, Fr(2 * k - 1, 2), Fr(k - 1)]))
        den = 10 ** 9
        dl = float(delta)
        menu = [lambda: Fr(0), lambda: t_of(rng.uniform(-1, 1) * aF, den),
                lambda: t_of(rng.uniform(-0.785, 0.785), den)]
        gaps = [lambda: rat(rng.uniform(0.01, 0.5) * dl, den), lambda: rat(rng.uniform(0.8, 1.2) * dl, den)]
        n = rng.randint(10, 30)
    b = Builder(k)
    for _ in range(n):
        b.drop(rat(rng.uniform(0.5, k - 0.5), den), rng.choice(menu)(), rng.choice(gaps)(),
               rng.choice([Fr(0), rng.choice(gaps)()]))
    return b.cfg(delta, tF, hF, 'jam' + ('_paper' if paper else ''))


def fam_rows(rng, paper=False):
    """dense rows spanning the width: long-lived paths, merges at converging junctions, E, T, deaths."""
    if paper:
        k = rng.choice([6, 7, 8])
        delta = Fr(1, 10 ** 5)
        tF = Fr(75, 10 ** 8)
        hF = Fr(k, 2) - 1
        den = 10 ** 16
        tl = [lambda: Fr(0), lambda: rat(rng.uniform(-1, 1) * 7.4e-7, den),
              lambda: rat(rng.choice([-1, 1]) * rng.uniform(0.5, 0.99) * 7.5e-7, den),
              lambda: rat(rng.choice([-1, 1]) * rng.uniform(1.0, 3.0) * 7.5e-7, den)]
        wl = lambda: rat(rng.choice([rng.uniform(1e-12, 5e-12), rng.uniform(1e-11, 1e-9), rng.uniform(1e-7, 1e-5)]), den)
        gl = lambda: rat(rng.uniform(0.03, 0.45) * 1e-5, den)
    else:
        k = rng.choice([4, 5, 6, 7])
        delta = Fr(rng.choice([10, 20, 35, 60]), 100)
        aF = rng.uniform(0.15, 0.75)
        tF = t_of(aF)
        hF = Fr(rng.choice([Fr(k, 2) - 1, Fr(2 * k - 1, 2), Fr(k - 1)]))
        den = 10 ** 10
        dl = float(delta)
        amax = min(0.95 * aF, 1.2 * dl)
        tl = [lambda: Fr(0), lambda: t_of(rng.uniform(-amax, amax), den),
              lambda: t_of(rng.choice([-1, 1]) * rng.uniform(0.5, 0.99) * amax, den),
              lambda: t_of(rng.choice([-1, 1]) * rng.uniform(1.0, 1.3) * aF, den)]
        wl = lambda: rat(rng.choice([rng.uniform(0.001, 0.05), rng.uniform(0.05, 0.3)]) * dl, den)
        gl = lambda: rat(rng.uniform(0.03, 0.45) * dl, den)
    b = Builder(k)
    fg = rat(rng.uniform(0, 0.4) * float(delta), den)
    for row in range(int(k) + 1):
        x = rat(rng.uniform(0.0, 0.3), den)
        while True:
            t = rng.choice(tl)()
            hw = halfwidth(t)
            if x + 2 * hw > k:
                break
            b.drop(x + hw, t, gl(), fg)
            x = x + 2 * hw + wl()
    return b.cfg(delta, tF, hF, 'rows' + ('_paper' if paper else ''))


def fam_wallmax(rng):
    """explicit drift-to-wall configuration: one tilted square per wall (vertex on floor, touching the wall),
    huge delta; returns (cfg, t_drift).  Designed to push L_W toward 2 y tan(alpha_s)."""
    k = rng.choice([6, 8, 11])
    a = rng.choice([0.03, 0.05, 0.1, 0.2, 0.3]) * rng.uniform(0.9, 1.1)
    den = 10 ** 12
    t = t_of(a, den)
    tF = t * Fr(1000001, 1000000)          # alpha_F just above the drift angle
    delta = Fr(4 * k)
    hF = Fr(2 * k - 1, 2)
    b = Builder(k)
    hw = halfwidth(t)
    b.add(hw, None or (lambda S0: -S0.ymin)(Sq(hw, 0, t)), t)             # left wall, phi=+a
    b.add(Fr(k) - hw, (lambda S0: -S0.ymin)(Sq(Fr(k) - hw, 0, -t)), -t)   # right wall, phi=-a
    if rng.random() < 0.5:   # add a second drift square higher up on each side, feeding more drift
        b.drop(hw + rat(rng.uniform(0.6, 1.2), den), t, rat(rng.uniform(0.01, 0.5), den), Fr(0))
        b.drop(Fr(k) - hw - rat(rng.uniform(0.6, 1.2), den), -t, rat(rng.uniform(0.01, 0.5), den), Fr(0))
    return b.cfg(delta, tF, hF, 'wallmax'), t


def fam_dag(rng):
    """alpha_F just below pi/4 and passable squares tilted almost alpha_F: the centre-height order
    (DAG lemma) has the smallest margin here; steep paths, E/W/T near 45 degrees."""
    k = rng.choice([4, 5, 6])
    aF = rng.uniform(0.76, 0.7853)
    den = 10 ** 12
    tF = t_of(aF, den)
    delta = Fr(rng.choice([50, 100, 200, 400]), 100)
    hF = Fr(2 * k - 1, 2)
    menu = [lambda: t_of(rng.choice([-1, 1]) * aF * rng.uniform(0.95, 0.99999), den),
            lambda: t_of(rng.choice([-1, 1]) * aF * rng.uniform(0.95, 0.99999), den),
            lambda: t_of(rng.choice([-1, 1]) * rng.uniform(aF, 0.78539816), den),
            lambda: Fr(0)]
    b = Builder(k)
    for _ in range(rng.randint(8, 24)):
        b.drop(rat(rng.uniform(0.5, k - 0.5), den), rng.choice(menu)(),
               rat(rng.uniform(0.001, 0.3) * min(float(delta), 1.0), den),
               rat(rng.choice([0.0, rng.uniform(0, 0.3)]), den))
    return b.cfg(delta, tF, hF, 'dag')


def fam_tower(rng, paper=False):
    if paper:
        k = 8
        delta = Fr(1, 10 ** 5)
        tF = Fr(75, 10 ** 8)
        hF = Fr(k, 2) - 1
        den = 10 ** 15
        tl = lambda: rat(rng.uniform(-1.2, 1.2) * 7.5e-7, den)
        gp = lambda: rat(rng.uniform(0.02, 0.45) * 1e-5, den)
    else:
        k = rng.choice([5, 6])
        delta = Fr(rng.choice([10, 30, 60]), 100)
        aF = rng.uniform(0.2, 0.7)
        tF = t_of(aF)
        hF = Fr(2 * k - 1, 2)
        den = 10 ** 9
        tl = lambda: t_of(rng.uniform(-1.1, 1.1) * aF, den)
        gp = lambda: rat(rng.uniform(0.02, 0.45) * float(delta), den)
    b = Builder(k)
    for cx in sorted(rng.uniform(0.6, k - 0.6) for _ in range(rng.randint(2, 4))):
        for _ in range(rng.randint(3, 6)):
            b.drop(rat(cx + rng.uniform(-0.02, 0.02), den), tl(), gp(), rng.choice([Fr(0), gp()]))
    return b.cfg(delta, tF, hF, 'tower' + ('_paper' if paper else ''))
