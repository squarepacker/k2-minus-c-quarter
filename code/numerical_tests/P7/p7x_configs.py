"""p7x_configs.py -- adversarial configuration families (exact rational squares), written from scratch."""
from fractions import Fraction as Fr
import math, random
from p7x_core import Sq, validate, HALF


def tq(phi, den=1000):
    """rational t ~ tan(phi/2)"""
    return Fr(int(round(math.tan(phi / 2) * den)), den)


def cs(t):
    t = Fr(t)
    return (1 - t * t) / (1 + t * t), 2 * t / (1 + t * t)


def bbox(t):
    c, s = cs(t)
    return c + abs(s)


def vprofile(S, x, upper):
    ys = []
    for (P, Q, _) in S.E.values():
        if P[0] == Q[0]:
            if P[0] == x:
                ys += [P[1], Q[1]]
        elif min(P[0], Q[0]) <= x <= max(P[0], Q[0]):
            ys.append(P[1] + (x - P[0]) * (Q[1] - P[1]) / (Q[0] - P[0]))
    if not ys:
        return None
    return max(ys) if upper else min(ys)


def drop(existing, sid, xc, t, gap, floor_gap=Fr(0)):
    """lowest position (exact) of a square of phase t centred at abscissa xc lying above every square it
    overlaps horizontally, with vertical clearance >= gap (> 0), and lowest point >= floor_gap."""
    T0 = Sq(-1, xc, 0, t)
    yreq = floor_gap - T0.ymin
    for S in existing:
        if S.xmin >= T0.xmax or S.xmax <= T0.xmin:
            continue
        lo = max(S.xmin, T0.xmin); hi = min(S.xmax, T0.xmax)
        pts = {lo, hi}
        for v in S.V + T0.V:
            if lo < v[0] < hi:
                pts.add(v[0])
        for x in pts:
            u = vprofile(S, x, True); l = vprofile(T0, x, False)
            if u is not None and l is not None:
                yreq = max(yreq, u - l + gap)
    return Sq(sid, xc, yreq, t)


# ------------------------------------------------------------ (a) blocked Z-columns
def cfg_blocked(k, tT, tZ, v1, nrows, gT, gZ, grow, flip='none', g0=Fr(1, 1000), xT0=Fr(0), xZ0=Fr(0),
                live_col=False):
    """bottom: row of T-terminators (phase tT, a >= alpha_max) on the floor; above: nrows rows of nearly
    flat squares (phase +-tZ) with lowest vertices v1, v1 + (bbox+grow), ...  live_col: leave the right-most
    unit column free of terminators and put an axis-parallel integer column there (live paths)."""
    sq = []; sid = 0
    bT = bbox(tT)
    xmaxT = k - (Fr(11, 10) if live_col else 0)
    x = xT0
    while x + bT <= xmaxT:
        sq.append(Sq(sid, x + bT / 2, g0 + bT / 2, tT)); sid += 1
        x += bT + gT
    bZ = bbox(tZ)
    rows = []
    v = Fr(v1)
    xmaxZ = k - (Fr(11, 10) if live_col else 0)
    for r in range(nrows):
        x = xZ0; i = 0
        while x + bZ <= xmaxZ:
            t = tZ
            if (flip == 'row' and r % 2) or (flip == 'col' and i % 2) or (flip == 'check' and (i + r) % 2):
                t = -tZ
            sq.append(Sq(sid, x + bZ / 2, v + bZ / 2, t)); sid += 1
            x += bZ + gZ; i += 1
        rows.append(v)
        v += bZ + grow
    if live_col:
        y = Fr(0); xc = k - HALF
        top = rows[-1] + bZ + 1
        while y + 1 <= min(top, k):
            sq.append(Sq(sid, xc, y + HALF, 0)); sid += 1
            y += 1 + Fr(1, 10 ** 5)
    return sq, rows, bZ


# ------------------------------------------------------------ (b)/(d) drift column + square on top
def cfg_drift(k, t_col, ncol, gap_col, gz, tZ, x_start, g0=Fr(0), nZ=2, zdx=Fr(-2, 5), side=1,
              filler=False):
    sq = []; sid = 0
    t = Fr(t_col) * side
    c, s = cs(t)
    n = (-s, c)
    # lowest vertex at (x_start, g0)
    if t > 0:
        cx = x_start + c / 2 - s / 2
    else:
        cx = x_start - c / 2 - s / 2
    cy = g0 + abs(s) / 2 + c / 2
    for i in range(ncol):
        sq.append(Sq(sid, cx + i * (1 + gap_col) * n[0], cy + i * (1 + gap_col) * n[1], t)); sid += 1
    top = sq[-1]
    TL = top.V[3] if t > 0 else top.V[2]
    xz = TL[0] + (zdx if t > 0 else -zdx)
    for j in range(nZ):
        Z = drop(sq, sid, xz, Fr(tZ) * (1 if j % 2 == 0 else -1), gz if j == 0 else Fr(1, 10 ** 4))
        sq.append(Z); sid += 1
    if filler:
        # axis-parallel filler columns (integer heights), skipping any square that would overlap
        from p7x_core import sat_disjoint
        base = list(sq)
        for xc in [Fr(1, 2) + i * Fr(10001, 10000) for i in range(int(k) - 1)]:
            y = Fr(0)
            while y + 1 <= Fr(k) / 2:
                T = Sq(sid, xc, y + HALF, 0)
                if all(sat_disjoint(T, S) for S in base if not (S.xmin > T.xmax or S.xmax < T.xmin or
                                                                 S.ymin > T.ymax or S.ymax < T.ymin)):
                    sq.append(T); sid += 1
                y += 1 + Fr(1, 10 ** 4)
    return sq


def cfg_grid(k, gap=Fr(1, 10 ** 4), ntilt=0, rng=None):
    """negative control: axis-parallel grid (everything reachable)."""
    sq = []; sid = 0
    n = int(k) - 1
    for i in range(n):
        for j in range(n):
            sq.append(Sq(sid, HALF + i * (1 + gap), HALF + j * (1 + gap), 0)); sid += 1
    return sq


# ------------------------------------------------------------ random Tetris-drop packings
def cfg_tetris(rng, k, nsq, tilt_mix, gap_mix, upper_tiny_from=None, tiny_max=0.004, den=2000):
    sq = []; sid = 0
    tries = 0
    while len(sq) < nsq and tries < nsq * 6:
        tries += 1
        best = None
        for _ in range(4):
            r = rng.random(); acc = 0
            for (p, lo, hi) in tilt_mix:
                acc += p
                if r <= acc:
                    phi = rng.uniform(lo, hi); break
            else:
                phi = 0.0
            if upper_tiny_from is not None and best is None:
                pass
            phi *= rng.choice((-1, 1))
            t = tq(phi, den) if phi != 0 else Fr(0)
            b = bbox(t)
            xc = Fr(rng.randrange(0, 10 ** 6), 10 ** 6) * (k - b) + b / 2
            gr = rng.random(); acc = 0
            for (p, g) in gap_mix:
                acc += p
                if gr <= acc:
                    gap = g; break
            S = drop(sq, sid, xc, t, gap)
            if upper_tiny_from is not None and S.ymin > upper_tiny_from and abs(phi) > tiny_max:
                t = tq(rng.uniform(0.0005, tiny_max) * rng.choice((-1, 1)), den)
                S = drop(sq, sid, xc, t, gap)
            if S.ymax <= k and (best is None or S.ymax < best.ymax):
                best = S
        if best is None:
            continue
        sq.append(best); sid += 1
    return sq


# ------------------------------------------------------------ merge valleys (converging beams)
def cfg_valley(rng, k, nrow, theta, gap, tiny, den=2000, v_off=Fr(0)):
    """rows alternating -theta,+theta (beams converge in pairs) dropped with small gaps; then tiny-tilt rows."""
    sq = []; sid = 0
    for r in range(nrow):
        x = Fr(rng.randrange(0, 500), 1000)
        i = 0
        while True:
            phi = theta if (i + r) % 2 else -theta
            t = tq(phi, den); b = bbox(t)
            if x + b > k:
                break
            S = drop(sq, sid, x + b / 2, t, gap, floor_gap=(v_off if r == 0 else Fr(0)))
            if S.ymax > k:
                break
            sq.append(S); sid += 1
            x += b + Fr(rng.randrange(1, 30), 10000)
            i += 1
    for r in range(3):
        x = Fr(rng.randrange(0, 300), 1000)
        while True:
            t = tq(tiny * rng.choice((-1, 1)), den); b = bbox(t)
            if x + b > k:
                break
            S = drop(sq, sid, x + b / 2, t, Fr(rng.choice((1, 5, 20, 100)), 1000))
            if S.ymax > k:
                break
            sq.append(S); sid += 1
            x += b + Fr(rng.randrange(1, 50), 10000)
    return sq


# ------------------------------------------------------------ (c) corner configurations
def cfg_corner(k, nrows, base, tZ, dipfrac, gap=Fr(1, 10 ** 6), nslots=1, slot_shift=0):
    """rows of k-nslots axis-parallel unit squares (gap between neighbours) with nslots interior slots of width
    w_s = 1 - (k-2 nslots-1) gap/nslots; above each slot a square of phase tZ>0 whose lowest corner dips by
    dipfrac*sin a cos a into the slot.  Lines through the dip: k-nslots long chords (exactly 1, E=0) + short chords."""
    k = int(k)
    sq = []; sid = 0
    tZ = abs(Fr(tZ))
    c, s = cs(tZ)
    sca = s * c; tanp = s / c
    nax = k - nslots
    ws = 1 - Fr(k - 2 * nslots - 1) * gap / nslots
    y = Fr(base)
    for r in range(nrows):
        after = sorted({((j + 1) * (nax // (nslots + 1)) + r * slot_shift) % (nax - 1) for j in range(nslots)})
        if len(after) < nslots:
            after = list(range(1, nslots + 1))
        x = Fr(0); slots = []
        for j in range(nax):
            sq.append(Sq(sid, x + HALF, y + HALF, 0)); sid += 1
            x += 1
            if j in after:
                slots.append(x); x += ws
            elif j < nax - 1:
                x += gap
        rowtop = y + 1
        dip = dipfrac * sca
        for xl in slots:
            BLx = xl + dip * tanp + Fr(1, 10 ** 5); BLy = rowtop - dip
            sq.append(Sq(sid, BLx + c / 2 - s / 2, BLy + s / 2 + c / 2, tZ)); sid += 1
        y = rowtop - dip + c + s + Fr(1, 10 ** 3)
    return sq
