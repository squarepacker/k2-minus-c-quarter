"""Driver: adversarial test of Theorem 6.9 (P4*), Lemma 6.7 and Proposition 6.10 with the exact tracers
of core.py.

usage: python run_tests.py unit
       python run_tests.py <families> <nconf> <seed> <budget_s> <out.jsonl> [per_conf_tl] [manyvar]

  families    comma-separated names from FAMS below (used in turn); names ending in _paper use the
              constants of the paper (generators.py)
  nconf       maximal number of configurations; budget_s: no new configuration after this many seconds
  out.jsonl   one JSON line per configuration is APPENDED (the folder is created if needed)
  per_conf_tl time limit per configuration (default 60 s; checks are cut short, flagged 'truncated')
  manyvar     more scale flows F_s per configuration: h_s at up to 8 exceptional heights of F^0 and
              tilt thresholds equal to the tilts of up to 6 passable squares

Messages of violations name the statement checked: (d) and (f) = Theorem 6.9(d),(f); (iv) = Lemma 6.7(iv);
(2) = Lemma 6.7(iii); (4) = the bound L_T <= N(s) of Theorem 6.9(f); (5) = Proposition 6.10.
"""
import os
import sys
import json
import time
import random
import bisect
from fractions import Fraction as Fr
from collections import defaultdict
from core import (Tracer, Cfg, validate, check_partition, exc_set, analyse, identity_checks, truncate,
                  point_trace, piece_events, piece_passed, ev, rss_mb, ZERO)
import generators as GEN

MEM_LIMIT_MB = 700


def fs(v):
    return None if v is None else float(v)


def choose_heights(cfg, pieces, exc, hF, rng, nrand, ncap_exc, ncap_tgt, eps):
    ys = []
    for _ in range(nrand):
        ys.append(hF * Fr(rng.randrange(1, 10 ** 9), 10 ** 9))
    el = sorted(exc)
    rng.shuffle(el)
    for e in el[:ncap_exc]:
        ys += [e, e - eps, e + eps]
    tgt = []
    for S in cfg.sq:
        for v in S.V:
            tgt += [v[1], v[1] - eps, v[1] + eps]
    smp = rng.sample(pieces, min(len(pieces), 40))
    rare = [p for p in pieces if p.kind in ('M', 'T', 'E', 'W')]
    rng.shuffle(rare)
    smp += rare[:40]
    for pc in smp:
        for xx in (pc.lo, pc.hi, (pc.lo + pc.hi) / 2):
            for (pt, lab) in piece_events(pc):
                tgt.append(ev(pt[1], xx))
    tgt = [y for y in set(tgt) if 0 < y < hF]
    rng.shuffle(tgt)
    ys += tgt[:ncap_tgt]
    out = []
    seen = set()
    for y in ys:
        if 0 < y < hF and y not in seen:
            seen.add(y)
            out.append(y)
    return out


def run_config(cfg, rng, opts, fam):
    T0 = time.time()
    tl = opts.get('tl', 90)
    ok, msg = validate(cfg.k, cfg.sq)
    if not ok:
        return {'fam': fam, 'skip': msg}
    out = {'fam': fam, 'name': cfg.name, 'k': fs(cfg.k), 'delta': fs(cfg.delta), 'tF': fs(cfg.tF),
           'hF': fs(cfg.hF), 'nsq': len(cfg.sq),
           'n_axis': sum(1 for S in cfg.sq if S.t == 0),
           'n_floor': sum(1 for S in cfg.sq if S.ymin == 0),
           'n_pass': sum(1 for S in cfg.sq if abs(S.t) < cfg.tF),
           'n_near45': sum(1 for S in cfg.sq if abs(S.t) > Fr(39, 100))}
    V = []

    def viol(kind, y, info):
        if len(V) < 40:
            V.append({'kind': kind, 'y': None if y is None else str(y), 'yf': fs(y), 'info': info})

    tr = Tracer(cfg, r1=True).run()
    P0 = tr.final
    out['tr_anom'] = tr.anom[:5]
    if tr.anom:
        viol('tracer anomaly', None, tr.anom[:5])
    pok, _ = check_partition(P0, cfg.k)
    if not pok:
        viol('partition', None, '')
    meas = defaultdict(lambda: ZERO)
    for p in P0:
        meas[p.kind] += p.hi - p.lo
    out['npieces'] = len(P0)
    out['meas'] = {kk: fs(v) for kk, v in meas.items()}
    out['minlam'] = fs(tr.minlam)
    out['minmu'] = fs(tr.minmu)
    out['merge_cells'] = tr.nentry_multi
    exc0 = exc_set(P0, cfg.hF)
    out['n_exc'] = len(exc0)
    labs = defaultdict(int)
    for e, ls in exc0.items():
        for lab in ls:
            okl = lab[0] == 'D' or (lab[0] in ('bot', 'top') and cfg.sq[lab[1]].t == 0)
            labs[lab[0] if okl else 'BAD_' + lab[0]] += 1
            if not okl:
                viol('(iv) exceptional height of unexpected type', e, str(lab))
    out['exc_labels'] = dict(labs)
    eps = opts.get('eps', Fr(1, 10 ** 12))
    ys = choose_heights(cfg, P0, exc0, cfg.hF, rng, opts.get('nrand', 25), opts.get('nexc', 20),
                        opts.get('ntgt', 35), eps)
    cache0 = {}
    st = dict(n=0, n_exc=0, n_exc_fail=0, n_exc_Bpos=0, min_margin_nonexc=None, min_margin_exc=None,
              maxrho=0.0, maxOvg=0.0, n_zero_margin=0, min_pslope=None, min_gslope=None)
    truncated = False
    for y in ys:
        if time.time() - T0 > tl:
            truncated = True
            break
        r = analyse(cfg, P0, y)
        cache0[y] = r
        st['n'] += 1
        bad = identity_checks(r, cfg.k)
        if bad:
            viol('F0 identity', y, bad)
        if r['an']:
            viol('F0 analysis anomaly', y, r['an'][:3])
        if r['Ov'] != 0:
            viol('F0 (d) Ov != 0', y, fs(r['Ov']))
        if r['overlap_pass'] != 0:
            viol('F0 pass-image overlap (R1 injectivity)', y, fs(r['overlap_pass']))
        for key in ('min_pslope', 'min_gslope'):
            if r[key] is not None:
                if r[key] < 1:
                    viol('F0 slope < 1 (' + key + ')', y, fs(r[key]))
                if st[key] is None or r[key] < st[key]:
                    st[key] = r[key]
        st['maxrho'] = max(st['maxrho'], float(r['maxrho']))
        st['maxOvg'] = max(st['maxOvg'], float(r['Ovg']))
        if y in exc0:
            st['n_exc'] += 1
            if r['B'] > 0:
                st['n_exc_Bpos'] += 1
            if r['margin'] < 0:
                st['n_exc_fail'] += 1
            if st['min_margin_exc'] is None or r['margin'] < st['min_margin_exc']:
                st['min_margin_exc'] = r['margin']
        else:
            if r['B'] != 0:
                viol('F0 B>0 at non-exceptional y', y, {'B': fs(r['B']), 'Btau': fs(r['Btau'])})
            if r['margin'] < 0:
                viol('F0 FINAL INEQUALITY FAILS at non-exceptional y', y,
                     {'margin': str(r['margin']), 'Sh': str(r['Sh']), 'LT': str(r['LT']), 'Lam': str(r['Lam'])})
            if r['margin'] != r['U']:
                viol('F0 slack != U at non-exceptional y', y, {'margin': fs(r['margin']), 'U': fs(r['U'])})
            if r['margin'] == 0:
                st['n_zero_margin'] += 1
            if st['min_margin_nonexc'] is None or r['margin'] < st['min_margin_nonexc']:
                st['min_margin_nonexc'] = r['margin']
    for kk in ('min_margin_nonexc', 'min_margin_exc', 'min_pslope', 'min_gslope'):
        st[kk] = fs(st[kk])
    out['F0'] = st
    # ---------------- scale flows F_s
    passable = sorted({abs(S.t) for S in cfg.sq if 0 < abs(S.t) < cfg.tF})
    axis_tops = sorted({S.ymax for S in cfg.sq if S.t == 0 and 0 < S.ymax <= cfg.hF})
    variants = [(cfg.tF, cfg.hF), (cfg.tF / 2, cfg.hF * 3 / 4), (cfg.tF / 4, cfg.hF / 2)]
    if passable:
        variants.append((rng.choice(passable), rng.choice(axis_tops) if axis_tops else cfg.hF * 2 / 3))
    for _ in range(opts.get('nvar_rand', 2)):
        variants.append((cfg.tF * Fr(rng.randrange(1, 1000), 1000),
                         cfg.hF * Fr(rng.randrange(300, 1001), 1000)))
    if exc0 and opts.get('var_at_exc', True):
        e = rng.choice(sorted(exc0))      # h_s placed exactly at an exceptional height of F0
        variants.append((cfg.tF * Fr(rng.randrange(1, 1000), 1000), e))
    if 'extra_variants' in opts:
        variants += opts['extra_variants']
    if opts.get('many_var'):
        el = sorted(exc0)
        rng.shuffle(el)
        for e in el[:8]:
            variants.append((cfg.tF * Fr(rng.randrange(1, 1000), 1000), e))
        for tp in passable[:6]:
            variants.append((tp, cfg.hF * Fr(rng.randrange(200, 1001), 1000)))
    fsum = []
    Dtot, Mtot, Etot = meas['D'], meas['M'], meas['E']
    for (ts, hs) in variants:
        if time.time() - T0 > tl * 1.6:
            truncated = True
            break
        Ps = truncate(cfg, P0, ts, hs)
        sv = {'ts': fs(ts), 'hs': fs(hs), 'n': 0, 'n_exc': 0, 'n_exc_fail': 0, 'min_margin_nonexc': None,
              'max_LW_ratio': 0.0, 'n_LW_pos': 0, 'min_N_minus_LT': None, 'max_Lam_ratio': 0.0}
        pok, _ = check_partition(Ps, cfg.k)
        if not pok:
            viol('Fs partition', None, str((fs(ts), fs(hs))))
        excs = exc_set(Ps, hs)
        notsub = [e for e in excs if e not in exc0 or not (0 < e < hs)]
        sv['n_exc_s'] = len(excs)
        if notsub:
            viol('(2) Exc(F_s) not subset of Exc(F0)', notsub[0], str(excs[notsub[0]]))
        N = sum((p.hi - p.lo for p in Ps if p.kind == 'T'), ZERO)
        sv['N'] = fs(N)
        tanas = 2 * ts / (1 - ts * ts)
        ys_s = [y for y in ys if y < hs]
        rng.shuffle(ys_s)
        ys_s = ys_s[:opts.get('ns', 30)]
        ys_s += [hs * Fr(rng.randrange(1, 10 ** 9), 10 ** 9) for _ in range(6)] + [hs - eps]
        own = []
        for pc in rng.sample(Ps, min(len(Ps), 15)):
            for xx in (pc.lo, pc.hi, (pc.lo + pc.hi) / 2):
                for (pt, lab) in piece_events(pc):
                    hy = ev(pt[1], xx)
                    own += [hy, hy - eps, hy + eps]
        rng.shuffle(own)
        ys_s += own[:20]
        for y in ys_s:
            if not (0 < y < hs):
                continue
            r = analyse(cfg, Ps, y)
            r0 = cache0.get(y)
            if r0 is None:
                r0 = analyse(cfg, P0, y)
                cache0[y] = r0
            sv['n'] += 1
            bad = identity_checks(r, cfg.k)
            if bad:
                viol('Fs identity', y, bad)
            if r['an']:
                viol('Fs analysis anomaly', y, r['an'][:3])
            if r['Ov'] != 0:
                viol('Fs (d) Ov != 0', y, fs(r['Ov']))
            if y in exc0:
                sv['n_exc'] += 1
                if r['margin'] < 0:
                    sv['n_exc_fail'] += 1
            else:
                if r['B'] != 0:
                    viol('Fs B>0 at y not in Exc(F0)', y, fs(r['B']))
                if r['margin'] < 0:
                    viol('Fs FINAL INEQUALITY FAILS at non-exceptional y', y, {'margin': str(r['margin'])})
                if r['margin'] != r['U']:
                    viol('Fs slack != U', y, '')
                if sv['min_margin_nonexc'] is None or r['margin'] < sv['min_margin_nonexc']:
                    sv['min_margin_nonexc'] = r['margin']
            L = r['L']
            if L['T'] > N:
                viol('(4) L_T^{F_s}(y) > N(s)', y, {'LT': str(L['T']), 'N': str(N)})
            if sv['min_N_minus_LT'] is None or N - L['T'] < sv['min_N_minus_LT']:
                sv['min_N_minus_LT'] = N - L['T']
            wb = 2 * y * tanas
            if L['W'] > wb:
                viol('(5) L_W^{F_s}(y) > 2 y tan alpha_s', y, {'LW': str(L['W']), 'bound': str(wb)})
            if L['W'] > 0:
                sv['n_LW_pos'] += 1
                sv['max_LW_ratio'] = max(sv['max_LW_ratio'], float(L['W'] / wb))
            if L['D'] > Dtot or L['M'] > Mtot or L['E'] > Etot:
                viol('Fs D/M/E loss exceeds F0 totals', y, '')
            if r['Ovg'] > r0['Ovg']:
                viol('Fs Ov_gap > F0 Ov_gap', y, {'s': str(r['Ovg']), '0': str(r0['Ovg'])})
            lam_b = Dtot + Mtot + Etot + r0['Ovg'] + wb
            if r['Lam'] > lam_b:
                viol('(5) Lambda^{F_s}(y) > Lambda_0(y-form) + 2y tan', y, {'Lam': str(r['Lam']), 'b': str(lam_b)})
            if lam_b > 0:
                sv['max_Lam_ratio'] = max(sv['max_Lam_ratio'], float(r['Lam'] / lam_b))
        for kk in ('min_margin_nonexc', 'min_N_minus_LT'):
            sv[kk] = fs(sv[kk])
        fsum.append(sv)
    out['Fs'] = fsum
    # ---------------- independent point checks
    Pso = sorted(P0, key=lambda p: p.lo)
    los = [p.lo for p in Pso]
    memo = {}
    npt = opts.get('npt', 40)
    xs = [cfg.k * Fr(rng.randrange(1, 10 ** 9), 10 ** 9) for _ in range(npt)]
    rare = [p for p in Pso if p.kind in ('M', 'T', 'E', 'W')]
    rng.shuffle(rare)
    xs += [(p.lo + p.hi) / 2 for p in rare[:opts.get('nrare', 25)]]
    pm = 0
    pc_n = 0
    pkinds = defaultdict(int)
    for x in xs:
        if time.time() - T0 > tl * 2.2:
            truncated = True
            break
        i = bisect.bisect_right(los, x) - 1
        if i < 0:
            continue
        pc = Pso[i]
        if not (pc.lo < x < pc.hi):
            continue
        try:
            kind, tau, seq = point_trace(cfg, x, memo)
        except RuntimeError as exc:
            viol('point tracer error', None, str(exc))
            continue
        pc_n += 1
        pkinds[kind] += 1
        if kind != pc.kind or tau != ev(pc.tau, x) or seq != piece_passed(pc):
            pm += 1
            viol('point/beam mismatch', None, {'x': str(x), 'point': [kind, fs(tau), list(seq)],
                                                'beam': [pc.kind, fs(ev(pc.tau, x)), list(piece_passed(pc))]})
    out['point'] = {'n': pc_n, 'mismatch': pm, 'kinds': dict(pkinds)}
    # ---------------- negative control: R1 off
    if meas['M'] > 0 and opts.get('negctl', True):
        tr2 = Tracer(cfg, r1=False).run()
        P2 = tr2.final
        exc2 = exc_set(P2, cfg.hF)
        ys2 = []
        Mp = [p for p in P0 if p.kind == 'M']
        rng.shuffle(Mp)
        for p in Mp[:12]:
            qy = ev(p.tau, (p.lo + p.hi) / 2)
            for off in (Fr(1, 10 ** 9), Fr(1, 10 ** 4), Fr(1, 10), Fr(1, 2), Fr(9, 10)):
                ys2.append(qy + off)
            ys2.append(qy + cfg.delta * Fr(rng.randrange(1, 1000), 1000))
        ys2 += [cfg.hF * Fr(rng.randrange(1, 10 ** 6), 10 ** 6) for _ in range(5)]
        nc = {'n': 0, 'n_Ov_pos': 0, 'max_Ov': 0.0, 'n_ineq_fail_nonexc': 0, 'min_margin_nonexc': None,
              'id_bad': 0}
        for y in ys2:
            if not (0 < y < cfg.hF):
                continue
            r = analyse(cfg, P2, y)
            nc['n'] += 1
            if identity_checks(r, cfg.k):
                nc['id_bad'] += 1
            if r['Ov'] > 0:
                nc['n_Ov_pos'] += 1
                nc['max_Ov'] = max(nc['max_Ov'], float(r['Ov']))
            if y not in exc2:
                if r['margin'] < 0:
                    nc['n_ineq_fail_nonexc'] += 1
                    if 'example' not in nc:
                        nc['example'] = {'y': str(y), 'margin': str(r['margin']), 'Ov': str(r['Ov']),
                                         'U': str(r['U'])}
                if nc['min_margin_nonexc'] is None or r['margin'] < nc['min_margin_nonexc']:
                    nc['min_margin_nonexc'] = r['margin']
        nc['min_margin_nonexc'] = fs(nc['min_margin_nonexc'])
        out['negctl'] = nc
    out['truncated'] = truncated
    out['viol'] = V
    out['nviol'] = len(V)
    out['time'] = round(time.time() - T0, 2)
    out['rss_mb'] = round(rss_mb(), 1)
    return out


# ------------------------------------------------------------------ unit tests (hand-computed)
def unit():
    res = []
    F = Fr
    # T1: one axis square, floor gap 1/20
    cfg = Cfg(4, [(2, F(1, 2) + F(1, 20), 0)], F(1, 10), F(1, 10), F(3, 2), 'T1')
    P = Tracer(cfg).run().final
    exc = exc_set(P, cfg.hF)
    r = analyse(cfg, P, F(1, 2))
    res.append(('T1 y=1/2', r['L']['D'] == 3 and r['summu'] == 1 and r['Fgap'] == 0 and r['omega'] == 3
                and r['Sh'] == 0 and r['margin'] == 3 and r['U'] == 3 and r['B'] == 0))
    r = analyse(cfg, P, F(108, 100))
    res.append(('T1 y=1.08', r['Fgap'] == 1 and r['omega'] == 4 and r['L']['D'] == 3 and r['U'] == 3
                and r['margin'] == 3))
    r = analyse(cfg, P, F(1, 20))
    res.append(('T1 y=1/20 (exc)', r['B'] == 1 and r['Sh'] == 1 and r['margin'] == -1 and F(1, 20) in exc))
    r = analyse(cfg, P, F(1, 10))
    res.append(('T1 y=1/10 (exc, D)', r['B'] == 3 and r['margin'] == 0 and r['U'] == 3 and F(1, 10) in exc))
    res.append(('T1 exc set', set(exc) == {F(1, 20), F(21, 20), F(1, 10), F(11, 10)}))
    # T2: the configuration of Proposition 6.17 (exceptional heights are needed)
    b = F(1, 100)
    sq = [(i + F(i + 1, 1000) + F(1, 2), b + F(1, 2), 0) for i in range(3)]
    cfg = Cfg(4, sq, F(1, 10), F(1, 10), F(7, 2), 'ex5.4')
    P = Tracer(cfg).run().final
    r = analyse(cfg, P, b)
    res.append(('T2 y=b', r['B'] == 3 and r['margin'] == -3 and r['Sh'] == 3))
    r = analyse(cfg, P, b + F(1, 10 ** 6))
    res.append(('T2 y=b+1e-6', r['B'] == 0 and r['margin'] == 0))
    # T3: symmetric valley -> merges; point tracer agreement; R1 off -> Ov>0
    a = F(1, 10)
    t1, t2 = -a / 2, a / 2   # tan(phi/2)
    hw = GEN.halfwidth(t1)
    w = F(1, 1000)
    bb = GEN.Builder(4)
    bb.drop(2 - hw - w / 2, t1, 0, 0)
    bb.drop(2 + hw + w / 2, t2, 0, 0)
    bb.drop(2, 0, F(1, 20))
    cfg = bb.cfg(F(3, 10), F(2, 10), F(3, 2), 'T3')
    tr = Tracer(cfg).run()
    P = tr.final
    M = sum((p.hi - p.lo for p in P if p.kind == 'M'), ZERO)
    res.append(('T3 merges exist', M > 0))
    rng = random.Random(1)
    o = run_config(cfg, rng, {'tl': 60, 'npt': 60}, 'unit')
    res.append(('T3 no violations', o['nviol'] == 0))
    res.append(('T3 point checks agree', o['point']['mismatch'] == 0 and o['point']['n'] > 30))
    res.append(('T3 negative control sees Ov>0', o.get('negctl', {}).get('n_Ov_pos', 0) > 0))
    for name, okk in res:
        print(('PASS ' if okk else 'FAIL ') + name)
    print(json.dumps({kk: o[kk] for kk in ('meas', 'F0', 'point', 'negctl', 'nviol')}, default=str)[:1500])
    return all(okk for _, okk in res)


FAMS = {
    'valley': lambda rng: GEN.fam_valley(rng, False),
    'valley_paper': lambda rng: GEN.fam_valley(rng, True),
    'zigzag': lambda rng: GEN.fam_zigzag(rng, False),
    'zigzag_paper': lambda rng: GEN.fam_zigzag(rng, True),
    'deaths': lambda rng: GEN.fam_deaths(rng, False),
    'deaths_paper': lambda rng: GEN.fam_deaths(rng, True),
    'near45': GEN.fam_near45,
    'wall': lambda rng: GEN.fam_wall(rng)[0],
    'floor': lambda rng: GEN.fam_floor(rng, False),
    'floor_paper': lambda rng: GEN.fam_floor(rng, True),
    'jam': lambda rng: GEN.fam_jam(rng, False),
    'jam_paper': lambda rng: GEN.fam_jam(rng, True),
    'tower': lambda rng: GEN.fam_tower(rng, False),
    'tower_paper': lambda rng: GEN.fam_tower(rng, True),
    'rows': lambda rng: GEN.fam_rows(rng, False),
    'rows_paper': lambda rng: GEN.fam_rows(rng, True),
    'wallmax': lambda rng: GEN.fam_wallmax(rng)[0],
    'dag': GEN.fam_dag,
}


def main():
    fam = sys.argv[1]
    if fam == 'unit':
        ok = unit()
        sys.exit(0 if ok else 1)
    nconf = int(sys.argv[2])
    seed = int(sys.argv[3])
    budget = float(sys.argv[4])
    outp = sys.argv[5]
    tl = float(sys.argv[6]) if len(sys.argv) > 6 else 60
    many_var = len(sys.argv) > 7 and sys.argv[7] == 'manyvar'
    fams = fam.split(',')
    rng = random.Random(seed)
    T = time.time()
    done = 0
    os.makedirs(os.path.dirname(os.path.abspath(outp)), exist_ok=True)
    with open(outp, 'a', encoding='utf-8') as fo:
        for i in range(nconf):
            if time.time() - T > budget:
                break
            if rss_mb() > MEM_LIMIT_MB:
                fo.write(json.dumps({'abort': 'memory', 'rss': rss_mb()}) + '\n')
                break
            f = fams[i % len(fams)]
            cfg = FAMS[f](rng)
            opts = {'tl': tl, 'nrand': 40, 'nexc': 40, 'ntgt': 80, 'ns': 40, 'npt': 60, 'nrare': 40,
                    'many_var': many_var}
            if many_var:
                opts['ns'] = 20
            if f.endswith('_paper'):
                opts['eps'] = Fr(1, 10 ** 13)
            try:
                o = run_config(cfg, rng, opts, f)
            except Exception as exc:
                import traceback
                o = {'fam': f, 'error': repr(exc), 'tb': traceback.format_exc()[-1500:]}
            o['seed'] = seed
            o['i'] = i
            fo.write(json.dumps(o, default=str) + '\n')
            fo.flush()
            done += 1
    print('done', done, 'configs in', round(time.time() - T, 1), 's; rss', round(rss_mb(), 1), 'MB')


if __name__ == '__main__':
    main()
