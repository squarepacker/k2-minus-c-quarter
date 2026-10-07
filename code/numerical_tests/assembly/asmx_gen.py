# -*- coding: utf-8 -*-
"""asmx_gen.py -- configuration generators of the assembly test (written from scratch, 2026-10-06).
Each generator returns a list of (cx, cy, phi) in the half box y <= Hh (Hh ~ k/2 - 0.02);
compose() puts an independent design reflected into the top half."""
import math, random
from asmx_core import Sq, sat_gap, Packing


def ext_of(a):
    a = abs(a)
    return math.cos(a) + math.sin(a)


def tilted_row(k, cy, a, n, x_left, eta_h):
    """n squares of phase angle a (same sign), centres at height cy, interlocked side by side."""
    out = []
    ca = math.cos(abs(a))
    w = ext_of(a)
    step = (1 + eta_h) / ca
    x = x_left + w / 2
    for i in range(n):
        out.append((x, cy, a))
        x += step
    return out


def fits(cfg, k, Hh):
    for (cx, cy, phi) in cfg:
        s = Sq(0, cx, cy, phi)
        if s.xlo < 1e-9 or s.xhi > k - 1e-9 or s.v < -1e-12 or s.yhi > Hh:
            return False
    return True


def trim(cfg, k, Hh):
    out = []
    for (cx, cy, phi) in cfg:
        s = Sq(0, cx, cy, phi)
        if s.xlo >= 1e-9 and s.xhi <= k - 1e-9 and s.v >= -1e-12 and s.yhi <= Hh:
            out.append((cx, cy, phi))
    return out


def top_of(cfg):
    if not cfg:
        return 0.0
    return max(Sq(0, *t).yhi for t in cfg)


# ----------------------------------------------------------------------------- rows
def gen_rows(k, Hh, rng, prm):
    """rows of axis / tilted squares stacked by bounding boxes with controlled vertical gaps
    (fractional heights), big gaps (deaths), big tilts (T), tilt in (alpha(s) range) etc."""
    cfg = []
    base = 0.0 if rng.random() < 0.6 else rng.choice([0.0, 1e-6, 0.3 * rng.random()])
    amax = prm['amax']; ay1 = prm['ay1']; beta_lo = prm['beta_lo']
    while True:
        typ = rng.choices(['axis', 'tiltS', 'tiltM', 'tiltB'], [0.45, 0.25, 0.15, 0.15])[0]
        if typ == 'axis':
            a = 0.0
        elif typ == 'tiltS':
            a = rng.uniform(0.15, 0.95) * beta_lo
        elif typ == 'tiltM':
            a = rng.uniform(ay1 * 1.02, amax * 0.98)
        else:
            a = rng.uniform(amax * 1.02, min(0.6, amax * 3 + 0.05))
        a *= rng.choice([-1, 1])
        w = ext_of(a)
        eta_h = rng.choice([1e-5, 1e-4, 1e-3, 0.01])
        nmax = int((k - w) / ((1 + eta_h) / math.cos(abs(a)))) + 1
        n = max(1, min(nmax, k - 1 - rng.choice([0, 0, 1, 2])))
        span = (n - 1) * (1 + eta_h) / math.cos(abs(a)) + w
        room = k - span
        if room < 2e-6:
            n -= 1
            span = (n - 1) * (1 + eta_h) / math.cos(abs(a)) + w
            room = k - span
        if n < 1:
            break
        xl = 1e-6 + rng.random() * max(0.0, room - 2e-6)
        cy = base + w / 2
        row = tilted_row(k, cy, a, n, xl, eta_h)
        # optionally remove one square to create a hole
        if rng.random() < 0.3 and len(row) > 2:
            row.pop(rng.randrange(len(row)))
        if cy + w / 2 > Hh:
            break
        cfg += row
        gapv = rng.choices([rng.uniform(1e-7, 1e-4), rng.uniform(1e-4, 0.01), rng.uniform(0.01, 0.6)],
                           [0.55, 0.3, 0.15])[0]
        base = cy + w / 2 + gapv
    return trim(cfg, k, Hh)


# ----------------------------------------------------------------------------- nested EG strips
def gen_nested(k, Hh, rng, prm):
    """shift layer (gap > delta or big-tilt row or axis row + gap) then nested tilted rows
    (each square stacked on the one below along n) with tilt below beta (Erdos-Graham-like strip)."""
    cfg = []
    amax = prm['amax']; beta_hi = prm['beta_hi']; delta = prm['delta']
    mode = rng.choice(['gap', 'bigtilt', 'axis+gap', 'axis+tilt', 'none'])
    base = 0.0
    if mode in ('axis+gap', 'axis+tilt'):
        n = k - 1 - rng.choice([0, 1])
        xl = 1e-6 + rng.random() * (k - n * 1.0001 - 2e-6)
        cfg += tilted_row(k, 0.5, 0.0, n, xl, 1e-4)
        base = 1.0 + 1e-4
    if mode in ('bigtilt', 'axis+tilt'):
        a = rng.uniform(amax * 1.05, min(0.5, amax * 3)) * rng.choice([-1, 1])
        w = ext_of(a)
        eta_h = 1e-3
        n = int((k - w) / ((1 + eta_h) / math.cos(abs(a)))) + 1
        n = max(1, n - rng.choice([0, 1]))
        span = (n - 1) * (1 + eta_h) / math.cos(abs(a)) + w
        xl = 1e-6 + rng.random() * max(0.0, k - span - 2e-6)
        cfg += tilted_row(k, base + w / 2, a, n, xl, eta_h)
        base = base + w
    # shift gap
    shift = rng.uniform(0.15, 0.75)
    if mode in ('gap', 'axis+gap'):
        shift = max(shift, delta * 1.3)
    base += shift
    # nested rows
    a = rng.uniform(0.3, 0.97) * beta_hi
    sg = rng.choice([-1, 1])
    phi = sg * a
    eta = rng.choice([2e-6, 1e-5, 5e-5, 2e-4])
    eta_h = rng.choice([1e-6, 1e-5, 1e-4])
    w = ext_of(a)
    nrows_max = 40
    n = k - 1 - rng.choice([0, 0, 1])
    # drift per row: centre moves by (1+eta)*n = (1+eta)*(-sin phi, cos phi)
    dxr = -(1 + eta) * math.sin(phi); dyr = (1 + eta) * math.cos(phi)
    span = (n - 1) * (1 + eta_h) / math.cos(a) + w
    # number of rows that fit vertically
    m = 0
    while m < nrows_max and base + w + m * dyr <= Hh:
        m += 1
    if m == 0:
        return trim(cfg, k, Hh)
    drift = abs(dxr) * (m - 1)
    while span + drift > k - 2e-6 and n > 1:
        n -= 1
        span = (n - 1) * (1 + eta_h) / math.cos(a) + w
    room = k - span - drift
    if room < 0:
        return trim(cfg, k, Hh)
    xl = 1e-6 + rng.random() * max(0.0, room - 2e-6)
    if dxr < 0:
        xl += drift
    row0 = tilted_row(k, base + w / 2, phi, n, xl, eta_h)
    for r in range(m):
        for (cx, cy, p) in row0:
            cfg.append((cx + r * dxr, cy + r * dyr, p))
    return trim(cfg, k, Hh)


# ----------------------------------------------------------------------------- columns (mechanism A)
def gen_columns(k, Hh, rng, prm):
    """columns of squares; some columns start with a tilted 'kicker' (big tilt => T, or
    tilt in the alpha(s) range) that raises the column to a fractional height (mechanism A),
    some contain a gap (deaths), some contain small-tilt squares (Z candidates)."""
    cfg = []
    amax = prm['amax']; ay1 = prm['ay1']; beta_hi = prm['beta_hi']
    x = 1e-6
    while True:
        kind = rng.choices(['plain', 'kickT', 'kickM', 'gap', 'smalltilt'], [0.3, 0.2, 0.15, 0.15, 0.2])[0]
        items = []
        if kind == 'kickT':
            items.append(rng.uniform(amax * 1.05, min(0.6, amax * 3 + 0.1)) * rng.choice([-1, 1]))
        elif kind == 'kickM':
            for _ in range(rng.choice([1, 2])):
                items.append(rng.uniform(ay1 * 1.02, amax * 0.98) * rng.choice([-1, 1]))
        elif kind == 'gap':
            items.append(('gap', rng.uniform(0.2, 0.7)))
        tiltsmall = rng.uniform(0.2, 0.95) * beta_hi * rng.choice([-1, 1])
        widths = [ext_of(t) for t in items if not isinstance(t, tuple)]
        wcol = max([1.0] + widths + ([ext_of(tiltsmall)] if kind == 'smalltilt' else []))
        if x + wcol > k - 1e-6:
            break
        cxc = x + wcol / 2
        y = 0.0 if rng.random() < 0.7 else rng.uniform(1e-6, 0.4)
        for t in items:
            if isinstance(t, tuple):
                y += t[1]
            else:
                e = ext_of(t)
                if y + e > Hh:
                    break
                cfg.append((cxc, y + e / 2, t)); y += e + rng.choice([1e-5, 1e-3, 0.02])
        while True:
            if kind == 'smalltilt' and rng.random() < 0.6:
                t = tiltsmall
            else:
                t = 0.0
            e = ext_of(t)
            if y + e > Hh:
                break
            cfg.append((cxc, y + e / 2, t))
            y += e + rng.choices([1e-6, 1e-4, 0.01, rng.uniform(0.1, 0.5)], [0.45, 0.3, 0.15, 0.1])[0]
        x += wcol + rng.choice([1e-5, 1e-4, 1e-3, 0.02])
    return trim(cfg, k, Hh)


# ----------------------------------------------------------------------------- L-shape
def gen_lshape(k, Hh, rng, prm):
    """L-shape: an axis column along the left wall and an axis row along the floor; the rest is
    a block of slightly tilted squares at a fractional height (shadowed only through the L)."""
    cfg = []
    beta_hi = prm['beta_hi']; amax = prm['amax']
    eta = rng.choice([1e-5, 1e-3])
    # left column
    y = 0.0
    while y + 1 <= Hh:
        cfg.append((0.5 + 1e-6, y + 0.5, 0.0)); y += 1 + eta
    # bottom row
    n = k - 2
    for i in range(n):
        cfg.append((1 + 1e-6 + eta + 0.5 + i * (1 + eta), 0.5, 0.0))
    # block above the row, to the right of the column
    a = rng.uniform(0.2, 0.95) * beta_hi
    if rng.random() < 0.3:
        a = rng.uniform(amax * 1.02, min(0.5, 3 * amax))
    phi = a * rng.choice([-1, 1])
    w = ext_of(a)
    base = 1.0 + rng.uniform(0.1, 0.7)
    xl = 1.0 + 2 * eta + 1e-4
    nn = int((k - xl - w) / ((1 + 1e-4) / math.cos(a))) + 1
    while base + w <= Hh:
        cfg += tilted_row(k, base + w / 2, phi, nn, xl, 1e-4)
        base += w + rng.choice([1e-5, 1e-3, 0.05, 0.3])
    return trim(cfg, k, Hh)


# ----------------------------------------------------------------------------- valley (merges)
def gen_valley(k, Hh, rng, prm):
    """two groups with opposite small tilts (< alpha_max): paths converge into the middle,
    pass through a middle square from two sources (R1 merges) and cross in gaps (Ov_gap)."""
    cfg = []
    amax = prm['amax']; delta = prm['delta']
    a = rng.uniform(0.2, 0.6) * min(amax, delta / 2.5)
    layers = rng.choice([1, 2, 3])
    # an axis row on the floor so that the V layer sits on a flat surface
    n = k - 1
    x0 = 1e-6 + rng.random() * max(0.0, k - n * (1 + 1e-5) - 2e-6)
    cfg += tilted_row(k, 0.5, 0.0, n, x0, 1e-5)
    base = 1.0 + rng.choice([1e-7, 1e-5, 1e-3])
    for L in range(layers):
        w = ext_of(a)
        nl = (k // 2) - 1
        eta_h = rng.choice([1e-5, 1e-4, 1e-3])
        span = (nl - 1) * (1 + eta_h) / math.cos(a) + w
        mid = k / 2.0 + rng.uniform(-0.3, 0.3)
        cgap = rng.choice([1e-6, 1e-4, 1e-3, 0.01])
        xl = mid - cgap / 2 - span
        left = tilted_row(k, base + w / 2, -a, nl, xl, eta_h)       # phi=-a: n points right
        right = tilted_row(k, base + w / 2, a, nl, mid + cgap / 2, eta_h)   # phi=+a: n points left
        cfg += left + right
        base += w + rng.uniform(0.05, 0.5) * delta
        # row above catching the converging beams: axis, or big tilt (merges at T_max squares)
        if rng.random() < 0.4:
            b = rng.uniform(amax * 1.05, min(0.5, amax * 3)) * rng.choice([-1, 1])
            wb = ext_of(b)
            nb = int((k - wb) / ((1 + 1e-4) / math.cos(abs(b)))) + 1
            spb = (nb - 1) * (1 + 1e-4) / math.cos(abs(b)) + wb
            xb0 = max(1e-6, mid - spb / 2 + rng.uniform(-0.5, 0.5))
            xb0 = min(xb0, k - spb - 1e-6)
            cfg += tilted_row(k, base + wb / 2, b, nb, xb0, 1e-4)
            base += wb + rng.choice([1e-7, 1e-5, 1e-3])
        else:
            n = k - 1 - rng.choice([0, 1, 2])
            sp = n * (1 + 1e-4)
            x0 = 1e-6 + rng.random() * max(0.0, k - sp - 2e-6)
            cfg += tilted_row(k, base + 0.5, 0.0, n, x0, 1e-4)
            base += 1.0 + rng.choice([1e-7, 1e-5, 1e-3])
    # then generic rows
    rest = gen_rows(k, Hh - base, rng, prm)
    cfg += [(cx, cy + base, p) for (cx, cy, p) in rest]
    return trim(cfg, k, Hh)


# ----------------------------------------------------------------------------- floor-touching squares
def gen_floortouch(k, Hh, rng, prm):
    """tilted squares resting on the floor on a vertex (and axis squares on the floor), then rows."""
    cfg = []
    amax = prm['amax']; beta_hi = prm['beta_hi']
    x = 1e-6
    while True:
        r = rng.random()
        if r < 0.35:
            a = rng.uniform(amax * 1.02, 0.6)
        elif r < 0.6:
            a = rng.uniform(0.2, 0.95) * beta_hi
        elif r < 0.8:
            a = rng.uniform(0.3, 0.98) * amax
        else:
            a = 0.0
        a *= rng.choice([-1, 1])
        w = ext_of(a)
        if x + w > k - 1e-6:
            break
        lift = rng.choice([0.0, 0.0, 1e-9, 1e-4])
        cfg.append((x + w / 2, w / 2 + lift, a))
        x += w + rng.choice([1e-5, 1e-3, 0.05])
    base = top_of(cfg) + rng.choice([1e-4, 0.05, 0.3])
    rest = gen_rows(k, Hh - base, rng, prm)
    cfg += [(cx, cy + base, p) for (cx, cy, p) in rest]
    return trim(cfg, k, Hh)


# ----------------------------------------------------------------------------- random jammed drop
def gen_random(k, Hh, rng, prm, max_sq=None):
    amax = prm['amax']; beta_hi = prm['beta_hi']; ay1 = prm['ay1']
    placed = []
    sqs = []
    fails = 0
    if max_sq is None:
        max_sq = int(k * Hh * 0.95)
    grid = {}
    def neigh(s):
        out = set()
        for gx in range(int(math.floor(s.xlo)) - 1, int(math.floor(s.xhi)) + 2):
            for gy in range(int(math.floor(s.v)) - 1, int(math.floor(s.yhi)) + 2):
                out.update(grid.get((gx, gy), ()))
        return out
    def ok(s):
        if s.xlo < 1e-9 or s.xhi > k - 1e-9 or s.v < 0:
            return False
        for j in neigh(s):
            t = sqs[j]
            if t.xlo > s.xhi + 1e-6 or t.xhi < s.xlo - 1e-6 or t.v > s.yhi + 1e-6 or t.yhi < s.v - 1e-6:
                continue
            if sat_gap(s, t) <= 2e-9:
                return False
        return True
    while len(placed) < max_sq and fails < 60:
        r = rng.random()
        if r < 0.55:
            a = 0.0
        elif r < 0.75:
            a = rng.uniform(0.1, 0.95) * beta_hi
        elif r < 0.88:
            a = rng.uniform(ay1, amax)
        else:
            a = rng.uniform(amax, 0.785)
        a *= rng.choice([-1, 1])
        w = ext_of(a)
        cx = rng.uniform(w / 2 + 1e-6, k - w / 2 - 1e-6)
        # drop: find the lowest feasible centre height by scanning then bisecting
        y = w / 2 + 1e-7
        s = Sq(0, cx, y, a)
        if not ok(s):
            lo = y; hi = None; yy = y
            while yy < Hh:
                yy += 0.07
                if ok(Sq(0, cx, yy, a)):
                    hi = yy; break
            if hi is None:
                fails += 1; continue
            lo = hi - 0.07
            for _ in range(30):
                mid = 0.5 * (lo + hi)
                if ok(Sq(0, cx, mid, a)):
                    hi = mid
                else:
                    lo = mid
            y = hi + rng.choice([1e-7, 1e-5, 1e-3])
            s = Sq(0, cx, y, a)
            if not ok(s):
                fails += 1; continue
        if s.yhi > Hh:
            fails += 1; continue
        idx = len(sqs)
        s2 = Sq(idx, cx, y, a)
        sqs.append(s2); placed.append((cx, y, a))
        for gx in range(int(math.floor(s2.xlo)), int(math.floor(s2.xhi)) + 1):
            for gy in range(int(math.floor(s2.v)), int(math.floor(s2.yhi)) + 1):
                grid.setdefault((gx, gy), []).append(idx)
        fails = 0
    return trim(placed, k, Hh)


# ----------------------------------------------------------------------------- EG nested strip, matched
def gen_nestedEG(k, Hh, rng, prm):
    """full-width nested strip of k-1 squares tilted by a < beta(top) such that the junction
    lines have both chords short and omega < omega0 < 1 (so f = 1-omega0-E > 0 there).
    A shift layer (gap > delta, a big-tilt row, or partial) puts the junctions at fractional heights."""
    cfg = []
    amax = prm['amax']; delta = prm['delta']; a_target = prm['a_eg']
    mode = rng.choice(['gap', 'bigtilt', 'axis+gap', 'axis+tilt', 'half'])
    base = 0.0
    if mode in ('axis+gap', 'axis+tilt'):
        n = k - 1
        xl = 1e-6 + rng.random() * (k - n * 1.00001 - 2e-6)
        cfg += tilted_row(k, 0.5, 0.0, n, xl, 1e-5)
        base = 1.0 + 1e-6
    if mode in ('bigtilt', 'axis+tilt'):
        a = rng.uniform(amax * 1.05, min(0.5, amax * 3)) * rng.choice([-1, 1])
        w = ext_of(a); eta_h = 1e-4
        n = int((k - w) / ((1 + eta_h) / math.cos(abs(a)))) + 1
        span = (n - 1) * (1 + eta_h) / math.cos(abs(a)) + w
        xl = 1e-6 + rng.random() * max(0.0, k - span - 2e-6)
        cfg += tilted_row(k, base + w / 2, a, n, xl, eta_h)
        base = base + w
    if mode == 'half':
        # left part: axis column stack reaching high (alive paths), right part: big gap
        pass
    target_frac = rng.uniform(0.35, 0.65)      # where the first junction should sit
    a = a_target * rng.uniform(0.85, 1.0)
    phi = a * rng.choice([-1, 1])
    eta = rng.choice([1e-7, 1e-6, 4e-6])
    eta_h = rng.choice([1e-7, 1e-6])
    w = ext_of(a)
    dyr = (1 + eta) * math.cos(a); dxr = -(1 + eta) * math.sin(phi)
    # first junction height ~ base0 + cos a + small  -> choose base0
    shift_min = delta * 1.3 if mode in ('gap', 'axis+gap') else 1e-4
    b0 = base + shift_min
    jf = (b0 + math.cos(a) + 0.5 * math.sin(a))
    add = (target_frac - (jf - math.floor(jf))) % 1.0
    base = b0 + add
    m = 0
    while m < 60 and base + w + m * dyr <= Hh:
        m += 1
    if m == 0:
        return trim(cfg, k, Hh)
    n = k - 1
    drift = abs(dxr) * (m - 1)
    span = (n - 1) * (1 + eta_h) / math.cos(a) + w
    while span + drift > k - 2e-6 and n > 1:
        n -= 1
        span = (n - 1) * (1 + eta_h) / math.cos(a) + w
    room = k - span - drift
    xl = 1e-6 + rng.random() * max(0.0, room - 2e-6)
    if dxr < 0:
        xl += drift
    row0 = tilted_row(k, base + w / 2, phi, n, xl, eta_h)
    for r in range(m):
        for (cx, cy, p) in row0:
            cfg.append((cx + r * dxr, cy + r * dyr, p))
    return trim(cfg, k, Hh)


# ----------------------------------------------------------------------------- hill (wall losses)
def gen_hill(k, Hh, rng, prm):
    """left half tilted +a (n points left), right half tilted -a (n points right), a < alpha(y1):
    every F_s passes them and drifts into both walls -> wall loss close to 2 y tan a."""
    cfg = []
    ay1 = prm['ay1']; delta = prm['delta']
    a = rng.uniform(0.6, 0.98) * min(ay1, delta / 3.0)
    w = ext_of(a)
    eta_h = 1e-6
    nl = (k // 2) - 1
    span = (nl - 1) * (1 + eta_h) / math.cos(a) + w
    # rows stacked VERTICALLY (same x every row): the paths (direction n) drift relative to the
    # squares, left group (+a, n points left) towards x=0, right group (-a) towards x=k
    dyr = (1 + 1e-7) / math.cos(a)
    base = rng.choice([0.0, 1e-7])
    xl = rng.choice([1e-7, 1e-6, 1e-4])
    xr = k - span - rng.choice([1e-7, 1e-6, 1e-4])
    if xr < xl + span:
        return []
    L0 = tilted_row(k, base + w / 2, a, nl, xl, eta_h)
    R0 = tilted_row(k, base + w / 2, -a, nl, xr, eta_h)
    r = 0
    while base + w + r * dyr <= Hh:
        for (cx, cy, p) in L0 + R0:
            cfg.append((cx, cy + r * dyr, p))
        r += 1
    return trim(cfg, k, Hh)


def gen_hill_v1(k, Hh, rng, prm):
    """first version of gen_hill (groups stacked along n; used by campaigns A and B, which started
    before gen_hill was rewritten).  Kept so that campaigns A and B can be reproduced exactly
    (asmx_run.py --hill-v1; use_hill_v1() in asmx_replay.py)."""
    cfg = []
    ay1 = prm['ay1']
    a = rng.uniform(0.6, 0.98) * ay1
    w = ext_of(a)
    eta_h = 1e-6
    nl = (k // 2) - 1
    span = (nl - 1) * (1 + eta_h) / math.cos(a) + w
    dyr = (1 + 1e-7) * math.cos(a)
    base = rng.choice([0.0, 1e-7])
    m = 0
    while base + w + m * dyr <= Hh:
        m += 1
    if m == 0:
        return []
    drift = (m - 1) * (1 + 1e-7) * math.sin(a)
    xl = drift + 1e-6 + rng.uniform(0, 0.2)
    xr = k - drift - 1e-6 - span - rng.uniform(0, 0.2)
    if xr < xl + span:
        nl -= 1
        span = (nl - 1) * (1 + eta_h) / math.cos(a) + w
        xr = k - drift - 1e-6 - span
    L0 = tilted_row(k, base + w / 2, a, nl, xl, eta_h)
    R0 = tilted_row(k, base + w / 2, -a, nl, xr, eta_h)
    for r in range(m):
        for (cx, cy, p) in L0:
            cfg.append((cx - r * (1 + 1e-7) * math.sin(a), cy + r * dyr, p))
        for (cx, cy, p) in R0:
            cfg.append((cx + r * (1 + 1e-7) * math.sin(a), cy + r * dyr, p))
    return trim(cfg, k, Hh)


GENS = dict(rows=gen_rows, nested=gen_nested, columns=gen_columns, lshape=gen_lshape,
            valley=gen_valley, floortouch=gen_floortouch, random=gen_random, hill=gen_hill)
GENS_ALL = dict(GENS, nestedEG=gen_nestedEG)


def compose(k, bottom, top):
    cfg = list(bottom) + [(cx, k - cy, -p) for (cx, cy, p) in top]
    return cfg
