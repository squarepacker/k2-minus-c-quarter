"""p7x_run.py <family> <seed> <time_limit_s> <out.jsonl>   -- campaign driver of the test of Section 9 (P7)
families: blocked   blocked Z-columns: stacks of slightly tilted squares Z of Z(s) above a row of tilt terminators
          blocked2  the same with Z tilts >= alpha(s) (violating (C3), which Section 9 does not use) or more rows
          drift     tilted columns that carry paths sideways towards squares of Z(s) (F^0 may reach them,
                    F_s must not)
          random    random piles (tetris, layered tetris), valleys, corner slots, and an axis-parallel grid
                    ('grid_negctl', where counting all squares instead of Z(s) must exceed the bound)
          round2    zigzag rows (many merges), staggered tilted columns, L-shaped terminator rows
Configurations are generated from random.Random(seed) until the time limit; each is evaluated by p7x_eval.py.
The output (one JSON line per configuration, appended) goes to the subfolder out/ of this script's folder.
Original campaigns (2026-10-06):  A = blocked 101 1500,  B = drift 202 1500,  C = random 303 1500 (stopped early),
R2 = round2 404 1200,  D = blocked2 505 1100,  E = drift 606 1100   (out files campA.jsonl ... campE.jsonl, campR2.jsonl).
Quick run:  python p7x_run.py blocked2 505 90 quick_D.jsonl
Memory guard (Windows only; inactive elsewhere): stop if the working set exceeds 700 MB or free RAM drops below 1.2 GB
between configurations."""
import os, sys, json, random, time, math, ctypes, traceback
from fractions import Fraction as Fr
import mpmath
from p7x_core import Sq, mpf_fr
from p7x_meas import Params
from p7x_eval import evaluate
from p7x_configs import (cfg_blocked, cfg_drift, cfg_tetris, cfg_valley, cfg_corner, tq, bbox, cs)


class MEMSTAT(ctypes.Structure):
    _fields_ = [('dwLength', ctypes.c_ulong), ('dwMemoryLoad', ctypes.c_ulong),
                ('ullTotalPhys', ctypes.c_ulonglong), ('ullAvailPhys', ctypes.c_ulonglong),
                ('ullTotalPageFile', ctypes.c_ulonglong), ('ullAvailPageFile', ctypes.c_ulonglong),
                ('ullTotalVirtual', ctypes.c_ulonglong), ('ullAvailVirtual', ctypes.c_ulonglong),
                ('ullAvailExtendedVirtual', ctypes.c_ulonglong)]


class PMC(ctypes.Structure):
    _fields_ = [('cb', ctypes.c_ulong), ('PageFaultCount', ctypes.c_ulong),
                ('PeakWorkingSetSize', ctypes.c_size_t), ('WorkingSetSize', ctypes.c_size_t),
                ('QuotaPeakPagedPoolUsage', ctypes.c_size_t), ('QuotaPagedPoolUsage', ctypes.c_size_t),
                ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t), ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                ('PagefileUsage', ctypes.c_size_t), ('PeakPagefileUsage', ctypes.c_size_t)]


def free_mb():
    if not hasattr(ctypes, 'windll'):          # memory guard only on Windows
        return float('inf')
    m = MEMSTAT(); m.dwLength = ctypes.sizeof(MEMSTAT)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    return m.ullAvailPhys / 2 ** 20


def ws_mb():
    if not hasattr(ctypes, 'windll'):          # memory guard only on Windows
        return 0.0
    p = PMC(); p.cb = ctypes.sizeof(PMC)
    h = ctypes.windll.kernel32.GetCurrentProcess()
    ctypes.windll.psapi.GetProcessMemoryInfo(h, ctypes.byref(p), p.cb)
    return p.WorkingSetSize / 2 ** 20


def svalid(P, ss):
    out = []
    for s in ss:
        s = Fr(s)
        if P.y0 <= s <= P.y1 and s not in out:
            out.append(s)
    return out


def rnd_fr(rng, lo, hi, den=10 ** 6):
    return Fr(lo) + (Fr(hi) - Fr(lo)) * Fr(rng.randrange(0, den + 1), den)


# ---------------------------------------------------------------- family generators
def gen_blocked2(rng):
    """variant: Z tilts in [alpha(s), beta(y_Z)) (violates constraint (C3), which Section 9 does not use;
    Remark 9.14), or large m."""
    return gen_blocked(rng, mode=rng.choice(['viol', 'viol', 'bigm']))


def gen_blocked(rng, mode='std'):
    c0 = Fr(rng.choice([5, 8, 10, 12, 15]), 100)
    y0 = rng.choice([Fr(3, 2), Fr(2), Fr(5, 2)])
    c1 = Fr(rng.choice([2, 3, 4]), 100)
    if mode == 'viol':
        c0 = Fr(rng.choice([1, 2, 3]), 100); c1 = Fr(rng.choice([4, 5, 6]), 100)
    am = float(c0) / math.sqrt(float(y0))
    thT = am * rng.uniform(1.003, 1.25)
    tT = tq(thT, 20000)
    while 2 * math.atan(float(tT)) < am * 1.0005:
        tT += Fr(1, 20000)
    sinT = float(cs(tT)[1])
    g0 = Fr(1, 1000)
    delta = Fr(int(1000 * max(sinT + 0.001 + rng.choice([0.003, 0.01, 0.03]), rng.choice([0.05, 0.08, 0.1]))) + 1, 1000)
    m = rng.choice([3, 4, 4, 5, 6]) if mode != 'bigm' else rng.choice([7, 8])
    grow = Fr(rng.choice([1, 2, 5, 10]), 10000)
    gZ = Fr(rng.choice([1, 2, 5]), 10000)
    gT = Fr(rng.choice([1, 2, 5]), 10000)
    topT = g0 + bbox(tT)
    # estimate top witness height to bound the Z tilt
    s_est = 2.2 + m * 1.01
    amax_Z = float(c1) * s_est ** -0.75
    aZ = amax_Z * rng.uniform(0.3, 0.9)
    if mode == 'viol':
        a_s = float(c0) / math.sqrt(s_est)          # alpha(s) at the top of the stack (approx.)
        if a_s * 1.05 >= amax_Z * 0.95:
            return None
        aZ = rng.uniform(a_s * 1.05, amax_Z * 0.95)   # Z tilt >= alpha(s): Z would be a T_s terminator if entered
    tZ = tq(aZ, 40000)
    if tZ == 0:
        tZ = Fr(1, 40000)
    bZ = bbox(tZ)
    dr = float(bZ) - 1 + float(grow)
    w0g = 0.5 * float(c0) ** 2 * 1.1 * (1 + 1.1 / float(y0)) + float(delta) + 0.06
    lo_f = w0g + 0.01; hi_f = 1 - w0g - m * dr - float(bZ - 1) - 0.01
    if hi_f <= lo_f:
        return None
    f = rng.uniform(lo_f, hi_f)
    base = 1 if 1 + f >= float(topT) + 0.003 else 2
    v1 = Fr(base) + Fr(int(f * 10 ** 6), 10 ** 6)
    flip = rng.choice(['none', 'none', 'row', 'col', 'check'])
    live = rng.random() < 0.25
    sq, rows, bZ = cfg_blocked(40, tT, tZ, v1, m + 1, gT, gZ, grow, flip=flip, g0=g0, live_col=False)
    c, s_ = cs(tZ); sca = abs(s_) * c
    A = rows[0] + bZ - sca * rnd_fr(rng, Fr(1, 10), Fr(9, 10))
    sv = rows[m - 1] + sca * rnd_fr(rng, Fr(1, 10), Fr(9, 10))
    if A < y0:
        return None
    eps = 1 - A / sv
    k = int(math.ceil(2 * (float(sv) ** 2 / float(A) + 3))) + rng.choice([0, 0, 1, 2])
    if k > (46 if mode != 'bigm' else 90):
        return None
    sq, rows, bZ = cfg_blocked(k, tT, tZ, v1, m + 1, gT, gZ, grow, flip=flip, g0=g0, live_col=live)
    # optional random deletions (holes)
    if rng.random() < 0.2:
        nd = rng.randint(1, 3)
        for _ in range(nd):
            sq.pop(rng.randrange(len(sq)))
        sq = [Sq(i, S.cx, S.cy, S.t) for i, S in enumerate(sq)]
    om0 = k + 1 if rng.random() < 0.8 else Fr(rng.choice([1, 5, 10]), 10)
    P = Params(k=k, delta=delta, c0=c0, c1=c1, y0=y0, eps=eps, omega0=om0)
    ss = [sv, sv - sca / 4, sv + sca / 4] + [rnd_fr(rng, P.y0, P.y1) for _ in range(1)]
    name = '%s(m=%d,flip=%s,live=%s)' % ('blocked' if mode == 'std' else 'blocked_' + mode, m, flip, live)
    return name, sq, P, svalid(P, ss)


def gen_drift(rng):
    c0 = Fr(rng.choice([30, 40, 50, 60]), 100)
    y0 = rng.choice([Fr(3, 2), Fr(2), Fr(3)])
    am = float(c0) / math.sqrt(float(y0))
    delta = Fr(rng.choice([2, 5, 10, 15]), 100)
    c1 = Fr(rng.choice([1, 2]), 100)
    eps = Fr(rng.choice([2, 3, 5]), 10)
    mode = rng.choice(['F0only', 'F0only', 'Fs_pass', 'edge'])
    ncol = rng.randint(2, 11)
    if mode == 'F0only':
        th = am * rng.uniform(0.6, 0.995)
    elif mode == 'Fs_pass':
        th = am * rng.uniform(0.97, 0.9995)
    else:
        th = am * rng.uniform(0.3, 0.999)
    tcol = tq(th, 20000)
    if tcol == 0:
        tcol = Fr(1, 20000)
    gap_col = Fr(rng.choice([0, 1, 5, 20]) + 1, 10 ** 5)
    gz = Fr(rng.choice([2, 5, 10, 20, 40]), 1000)
    tZ = Fr(rng.choice([1, 2, 3]), 4000)
    zdx = Fr(-rng.randrange(30, 49), 100)
    yZ_est = ncol + 1.5
    k = int(math.ceil(2 * (yZ_est / (1 - float(eps)) + 3))) + 1
    k = max(k, int(2 * (float(y0) + 3)) + 1)
    if k > 60:
        return None
    side = rng.choice([1, -1])
    xs = Fr(k, 2) + (Fr(ncol, 3) if side == 1 else -Fr(ncol, 3))
    sq = cfg_drift(k, tcol, ncol, gap_col, gz, tZ, xs, g0=Fr(rng.choice([0, 1, 5]), 1000), nZ=rng.choice([1, 2, 3]),
                   zdx=zdx, side=side, filler=rng.random() < 0.3)
    P = Params(k=k, delta=delta, c0=c0, c1=c1, y0=y0, eps=eps, omega0=k + 1)
    ss = [P.y0]
    Zs = sq[ncol:ncol + 3]
    for Z in Zs:
        for yz in (Z.ymax, Z.ymin):
            ss += [yz, yz / (1 - eps) - Fr(1, 1000), (yz + yz / (1 - eps)) / 2]
    thr = (float(c0) / float(Sq(-1, 0, 0, tcol).a)) ** 2
    ss += [Fr(int(thr * 10 ** 6), 10 ** 6) - Fr(1, 10 ** 6), Fr(int(thr * 10 ** 6), 10 ** 6) + Fr(2, 10 ** 6)]
    ss = svalid(P, ss)
    rng.shuffle(ss)
    return 'drift(%s,ncol=%d)' % (mode, ncol), sq, P, ss[:6]


def gen_random(rng):
    kind = rng.choice(['tetris', 'tetris', 'tetris_layer', 'valley', 'corner', 'corner', 'grid'])
    if kind == 'grid':
        from p7x_configs import cfg_grid
        k = rng.choice([10, 12])
        P = Params(k=k, delta=Fr(1, 20), c0=Fr(1, 10), c1=Fr(3, 100), y0=2, eps=Fr(1, 2), omega0=k + 1)
        return 'grid_negctl', cfg_grid(k), P, svalid(P, [P.y0, P.y1])
    k = rng.choice([12, 14, 16, 18])
    y0 = rng.choice([Fr(1), Fr(3, 2), Fr(2)])
    c0 = Fr(rng.choice([10, 20, 30, 45]), 100)
    c1 = Fr(rng.choice([1, 2, 3, 5]), 100)
    delta = Fr(rng.choice([2, 5, 10, 20, 30]), 100)
    eps = Fr(rng.choice([2, 3, 5, 7]), 10)
    om0 = rng.choice([Fr(k + 1), Fr(1, 2), Fr(9, 10), Fr(2)])
    if kind.startswith('tetris'):
        mix = [(0.35, 0.0, 0.004), (0.2, 0.004, 0.03), (0.2, 0.03, 0.2), (0.1, 0.2, 0.7), (0.15, 0.0, 0.0)]
        gaps = [(0.4, Fr(1, 10000)), (0.25, Fr(1, 1000)), (0.15, Fr(1, 100)), (0.1, Fr(5, 100)), (0.1, Fr(1, 5))]
        nsq = rng.randint(30, int(k * k / 2.2))
        sq = cfg_tetris(rng, k, nsq, mix, gaps, upper_tiny_from=(float(y0) + 1 if kind == 'tetris_layer' else None))
    elif kind == 'valley':
        th = float(c0) / math.sqrt(float(y0)) * rng.uniform(0.2, 0.95)
        sq = cfg_valley(rng, k, rng.randint(1, 4), th, Fr(rng.choice([1, 5, 20]), 10000), rng.uniform(0.0005, 0.004),
                        v_off=Fr(rng.choice([0, 0, 1, 3]), 10))
    else:
        k = rng.choice([16, 20, 24, 28])
        y0 = rng.choice([Fr(3, 2), Fr(2)])
        c0 = Fr(1, 10); c1 = Fr(5, 100)
        delta = Fr(rng.choice([1, 2, 5]), 100)
        eps = Fr(rng.choice([2, 3, 5]), 10)
        tZ = tq(rng.uniform(0.002, 0.009), 20000)
        dipf = Fr(rng.randrange(20, 95), 100)
        sq = cfg_corner(k, rng.randint(2, 4), Fr(rng.randrange(1300, 1700), 1000), tZ, dipf,
                        nslots=rng.choice([1, 1, 2, 3]), slot_shift=rng.choice([0, 1, 3]))
        om0 = min(Fr(99, 100), 1 - dipf + Fr(rng.randrange(1, 40), 100))
    P = Params(k=k, delta=delta, c0=c0, c1=c1, y0=y0, eps=eps, omega0=om0)
    if P.y1 < P.y0:
        return None
    ss = [P.y0, P.y1] + [rnd_fr(rng, P.y0, P.y1) for _ in range(3)]
    if kind == 'corner':
        tops = sorted({S.ymax for S in sq if S.s == 0})
        ss = []
        for t in tops:
            ss += [t, t / (1 - eps) - Fr(1, 10 ** 4), t + Fr(1, 10)]
        rng.shuffle(ss)
        ss = ss[:5] + [P.y1]
    return kind, sq, P, svalid(P, ss)


def gen_round2(rng):
    """zigzag (merge-heavy), stagger (staggered tilted columns), lshape (terminators on part of the floor)."""
    from p7x_configs import drop
    kind = rng.choice(['zigzag', 'stagger', 'lshape'])
    k = rng.choice([14, 16, 18, 20])
    y0 = rng.choice([Fr(3, 2), Fr(2)])
    c0 = Fr(rng.choice([20, 30, 40, 50]), 100)
    c1 = Fr(rng.choice([1, 2, 3]), 100)
    delta = Fr(rng.choice([5, 10, 15, 20]), 100)
    eps = Fr(rng.choice([2, 3, 5]), 10)
    am = float(c0) / math.sqrt(float(y0))
    sq = []; sid = 0
    if kind == 'zigzag':
        th = am * rng.uniform(0.3, 0.97)
        nrow = rng.randint(3, 7)
        g = Fr(rng.choice([1, 5, 20, 100]), 10 ** 5)
        for r in range(nrow):
            x = Fr(rng.randrange(0, 600), 1000); i = 0
            while True:
                t = tq(th if (i + r) % 2 else -th, 20000); b = bbox(t)
                if x + b > k:
                    break
                S = drop(sq, sid, x + b / 2, t, g)
                if S.ymax > k / 2:
                    break
                sq.append(S); sid += 1
                x += b + Fr(rng.randrange(1, 40), 10 ** 5); i += 1
    elif kind == 'stagger':
        x = Fr(rng.randrange(0, 300), 1000)
        while x + Fr(3, 2) < k:
            th = am * rng.uniform(0.3, 1.0) * rng.choice((-1, 1))
            t = tq(th, 20000); b = bbox(t)
            for j in range(rng.randint(2, 7)):
                S = drop(sq, sid, x + b / 2, t, Fr(rng.choice([1, 10, 100]), 10 ** 5))
                if S.ymax > k / 2:
                    break
                sq.append(S); sid += 1
            x += b + Fr(rng.randrange(1, 200), 1000)
    else:
        wl = rng.choice([1, 2, 3])
        tT = tq(am * rng.uniform(1.01, 1.3), 20000); bT = bbox(tT)
        x = Fr(0)
        while x + bT <= k - wl - Fr(1, 10):
            sq.append(Sq(sid, x + bT / 2, Fr(1, 1000) + bT / 2, tT)); sid += 1
            x += bT + Fr(1, 1000)
        for i in range(wl):
            xc = k - Fr(1, 2) - i * Fr(10001, 10000)
            y = Fr(0)
            while y + 1 <= Fr(k) / 2:
                sq.append(Sq(sid, xc, y + Fr(1, 2), 0)); sid += 1
                y += 1 + Fr(1, 10 ** 5)
    # tiny-tilt rows on top (fractional heights arise from the drops)
    for r in range(rng.randint(2, 4)):
        x = Fr(rng.randrange(0, 300), 1000)
        while True:
            t = tq(rng.uniform(0.0005, 0.004) * rng.choice((-1, 1)), 40000); b = bbox(t)
            if x + b > k:
                break
            S = drop(sq, sid, x + b / 2, t, Fr(rng.choice([1, 50, 200, 400]), 1000))
            if S.ymax > k:
                break
            sq.append(S); sid += 1
            x += b + Fr(rng.randrange(1, 50), 10000)
    P = Params(k=k, delta=delta, c0=c0, c1=c1, y0=y0, eps=eps, omega0=rng.choice([Fr(k + 1), Fr(1, 2), Fr(9, 10)]))
    if P.y1 < P.y0:
        return None
    hs = sorted({S.ymin for S in sq if abs(S.t) < Fr(1, 400)})
    ss = [P.y0, P.y1] + [rnd_fr(rng, P.y0, P.y1) for _ in range(2)] + [h for h in rng.sample(hs, min(3, len(hs)))]
    return kind, sq, P, svalid(P, ss)


GENS = {'blocked': gen_blocked, 'blocked2': gen_blocked2, 'drift': gen_drift, 'random': gen_random,
        'round2': gen_round2}


def main():
    fam, seed, tl, outp = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    os.makedirs(outdir, exist_ok=True)
    outp = os.path.join(outdir, os.path.basename(outp))
    rng = random.Random(seed)
    t0 = time.time()
    n = 0; nerr = 0
    with open(outp, 'a', encoding='utf-8') as fo:
        while time.time() - t0 < tl:
            if ws_mb() > 700 or free_mb() < 1200:
                fo.write(json.dumps({'abort': 'memory', 'ws': ws_mb(), 'free': free_mb()}) + '\n'); break
            try:
                g = GENS[fam](rng)
                if g is None:
                    continue
                name, sq, P, ss = g
                if not ss:
                    continue
                rem = tl - (time.time() - t0)
                r = evaluate(name, sq, P, ss, rng, nys=6, xcheck_n=25, ovg_grid=20, tlimit=max(5, min(240, rem)))
                r['seed'] = seed; r['idx'] = n; r['ws_mb'] = ws_mb()
                r['sq'] = [(str(S.cx), str(S.cy), str(S.t)) for S in sq] if len(sq) <= 400 else None
                fo.write(json.dumps(r, default=str) + '\n'); fo.flush()
                n += 1
            except Exception as e:
                nerr += 1
                fo.write(json.dumps({'error': repr(e), 'tb': traceback.format_exc()[-1500:]}) + '\n'); fo.flush()
                if nerr > 50:
                    break
    print('done', fam, seed, n, nerr, time.time() - t0)


if __name__ == '__main__':
    main()
