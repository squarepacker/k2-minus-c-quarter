# -*- coding: utf-8 -*-
"""asmx_run.py -- campaign driver of the assembly test.
usage: python asmx_run.py <campaign> <seed> <n_configs> <time_limit_s> [families [k_values]] [--hill-v1]
  families  comma-separated names of asmx_gen.GENS_ALL (default: all of asmx_gen.GENS); configuration i uses
            families[i % len(families)] for the bottom half (the top half is the same or a random family)
  k_values  comma-separated (default 16,20,24)
  --hill-v1 use the first version of the generator 'hill' (gen_hill_v1), as campaigns A and B did
Output (appended, one JSON line per configuration): out/<campaign>.jsonl; memory log: out/run_log.txt.
Original campaigns (2026-10-06; stopped by the time limit):
  python asmx_run.py campA 11 100000 1200 --hill-v1
  python asmx_run.py campB 12 100000 1200 nestedEG 24,28,32 --hill-v1
  python asmx_run.py campC 13 100000 1000 valley,hill,columns,floortouch,rows,lshape 20,24,28
Quick run:  python asmx_run.py quickC 13 3 150 valley,hill,columns,floortouch,rows,lshape 20,24,28
Memory guard (Windows only; skipped elsewhere): waits until 1.5 GB are free (2-min retries; after 40 min only a
tiny run)."""
import sys, os, json, time, math, random, subprocess, traceback
import numpy as np
from asmx_core import Packing
from asmx_chain import Params, LineStruct, FamilyFlows, family_links, assembly
from asmx_gen import GENS, compose

HERE = os.path.dirname(os.path.abspath(__file__))


def free_mb():
    if os.name != 'nt':                         # memory guard only on Windows
        return float('inf')
    try:
        out = subprocess.run(['powershell', '-NoProfile', '-Command',
                              '(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory'],
                             capture_output=True, text=True, timeout=60).stdout.strip()
        return int(out) / 1024.0
    except Exception:
        return -1.0


def mem_guard(logf, tag):
    t0 = time.time()
    while True:
        f = free_mb()
        with open(logf, 'a', encoding='utf-8') as fh:
            fh.write('%s %s free=%.0f MB\n' % (time.strftime('%H:%M:%S'), tag, f))
        if f >= 1536:
            return 'ok'
        if time.time() - t0 > 40 * 60:
            return 'tiny'
        time.sleep(120)


def rss_mb():
    try:
        import ctypes
        from ctypes import wintypes
        class PMC(ctypes.Structure):
            _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD),
                        ('PeakWorkingSetSize', ctypes.c_size_t), ('WorkingSetSize', ctypes.c_size_t),
                        ('QuotaPeakPagedPoolUsage', ctypes.c_size_t), ('QuotaPagedPoolUsage', ctypes.c_size_t),
                        ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t), ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                        ('PagefileUsage', ctypes.c_size_t), ('PeakPagefileUsage', ctypes.c_size_t)]
        pmc = PMC(); pmc.cb = ctypes.sizeof(PMC)
        h = ctypes.windll.kernel32.GetCurrentProcess()
        ctypes.windll.psapi.GetProcessMemoryInfo(h, ctypes.byref(pmc), pmc.cb)
        return pmc.WorkingSetSize / 1e6
    except Exception:
        return -1.0


def pick_flow_params(rng, k):
    while True:
        delta = rng.choice([0.05, 0.08, 0.1, 0.15, 0.2])
        c0 = rng.choice([0.1, 0.15, 0.2, 0.25])
        y0 = rng.choice([2, 3, 4])
        y1 = k / 2 - 3
        if y0 + 1 <= 0.5 * y1 + 0.5 or y0 == 2:
            if y0 + 1 <= 0.8 * y1:
                return delta, c0, y0


def line_variants(rng, k, delta, c0, y0):
    out = []
    tries = 0
    while len(out) < 4 and tries < 200:
        tries += 1
        c1 = rng.choice([0.02, 0.04, 0.06, 0.1, 0.15, 0.25, 0.35])
        eps = rng.choice([0.2, 0.3, 0.4, 0.5])
        om = rng.choice([0.5, 0.8, 0.9, 0.95, 0.98, 1.5, float(k)])
        prm = Params(k, delta, c0, y0, c1, eps, om)
        if prm.ok():
            continue
        out.append(prm)
    return out


def summarize_links(fl):
    per = fl['per_s']
    keys = ['r_paper', 'r_meas', 'r_sharp', 'r_nobW', 'r_shadow', 'lb1', 'lb2', 'lbZ']
    agg = {}
    for kk in keys:
        vals = [d[kk] for d in per if kk in d]
        agg[kk] = max(vals) if vals else None
    agg['nZmax'] = max(d['nZ'] for d in per)
    agg['p2v'] = sum(d['p2v'] for d in per)
    agg['outV'] = sum(d['outV'] for d in per)
    agg['tiltv'] = sum(d['tiltv'] for d in per)
    agg['step3min'] = min([d['step3min'] for d in per if 'step3min' in d] or [None]) if any('step3min' in d for d in per) else None
    agg['p4min'] = min([d['p4min'] for d in per if 'p4min' in d]) if any('p4min' in d for d in per) else None
    agg['bdmax'] = max([d['bdmax'] for d in per if 'bdmax' in d]) if any('bdmax' in d for d in per) else None
    agg['ovmax'] = max([d['ovmax'] for d in per if 'ovmax' in d]) if any('ovmax' in d for d in per) else None
    agg['LTminusN'] = max([d['LTminusN'] for d in per if 'LTminusN' in d]) if any('LTminusN' in d for d in per) else None
    agg['lam_excess'] = max([d['lam_sup_s'] - d['lam_bound'] for d in per if 'lam_sup_s' in d]) if any('lam_sup_s' in d for d in per) else None
    agg['wall_ratio'] = max([d['LW'] / d['wallb'] for d in per])
    agg['cases_Zpos'] = sum(1 for d in per if d['nZ'] > 0)
    agg['cases'] = len(per)
    for kk in ['symd_max', 'N_vs_n_max', 'N_monotone', 'T1_ratio_exact', 'T1_ratio_grid', 'T2_ratio',
               'T2_swap_rel', 'T2s_ratio', 'LBint_ratio', 'LBint_ratio_star', 'P7int_ratio_paper',
               'P7int_ratio_meas', 'chain_ratio_paper', 'chain_ratio_meas', 'chain_ratio_star_meas',
               'LHS3', 'LHS3star', 'Mb', 'Zint', 'P7int_meas', 'P7int_paper', 'RHS_paper', 'Lam0',
               'supOv0', 'intA', 'dom', 'OmegaW', 'Wall_meas', 'meas0']:
        agg[kk] = fl[kk]
    # worst P7 case detail
    best = None
    for d in per:
        v = d.get('r_sharp', d['r_meas'])
        if d['nZ'] > 0 and (best is None or v > best[0]):
            best = (v, d['s'], d['nZ'], d['N'], d.get('sharp'), d['r_paper'], d['r_meas'])
    agg['worstP7'] = best
    return agg


def run_one(cid, fam, k, rng, tlim, smoke=False):
    from asmx_gen import GENS_ALL
    eg = (fam == 'nestedEG')
    if eg:
        delta = rng.choice([0.02, 0.03, 0.05]); c0 = rng.choice([0.08, 0.1, 0.15]); y0 = rng.choice([3, 4])
        c1d = rng.choice([0.25, 0.3, 0.35]); epsd = rng.choice([0.3, 0.35])
    else:
        delta, c0, y0 = pick_flow_params(rng, k)
        c1d = rng.choice([0.04, 0.06, 0.1, 0.15, 0.25, 0.35]); epsd = 0.35
    y1 = k / 2 - 3
    top = (1 - epsd) * y1
    gprm = dict(amax=c0 / math.sqrt(y0), ay1=c0 / math.sqrt(y1), beta_lo=c1d * top ** -0.75,
                beta_hi=c1d * y0 ** -0.75, delta=delta, a_eg=0.97 * c1d * top ** -0.75)
    Hh = k / 2.0 - 0.02
    if eg:
        fam2 = 'nestedEG'
    else:
        fam2 = rng.choice(list(GENS.keys())) if rng.random() < 0.5 else fam
        if fam2 == 'random' and k > 20:
            fam2 = 'rows'
    bot = GENS_ALL[fam](k, Hh, rng, gprm)
    topd = GENS_ALL[fam2](k, Hh, rng, gprm)
    cfg = compose(k, bot, topd)
    P = Packing(k, cfg)
    ok, why = P.validate()
    rec = dict(cid=cid, fam=fam, fam_top=fam2, k=k, N=P.N, W=P.W, delta=delta, c0=c0, y0=y0, c1d=c1d)
    if not ok:
        rec['invalid'] = str(why)
        return rec, None
    s_grid = sorted(set([float(y0), float(y1)] + list(np.linspace(y0, y1, 15)[1:-1] + 1e-7 * math.pi) +
                        [rng.uniform(y0, y1) for _ in range(5)]))
    t0 = time.time()
    FFb = FamilyFlows(P, k, delta, c0, y0, s_grid, tlimit=tlim)
    Pr = P.reflect()
    FFt = FamilyFlows(Pr, k, delta, c0, y0, s_grid, tlimit=tlim)
    rec['t_flow'] = time.time() - t0
    rec['pieces'] = [FFb.npieces, FFt.npieces]
    rec['anom'] = [len(FFb.anom), len(FFt.anom)]
    rec['anom_detail'] = [str(a) for a in (FFb.anom + FFt.anom)[:5]]
    rec['dropped'] = [FFb.dropped, FFt.dropped]
    rec['total_measure_err'] = [abs(FFb.total_measure - k), abs(FFt.total_measure - k)]
    rec['F0_bd'] = [FFb.F0_bd, FFt.F0_bd]; rec['F0_ov'] = [FFb.F0_ov, FFt.F0_ov]
    rec['F0_p4min'] = [FFb.F0_p4min, FFt.F0_p4min]
    rec['R1multi'] = [FFb.nR1multi, FFt.nR1multi]
    rec['meas0_b'] = [float(x) for x in FFb.meas0]; rec['meas0_t'] = [float(x) for x in FFt.meas0]
    if eg:
        variants = []
        for om in (0.9, 0.95, 0.97, 0.99, float(k)):
            p = Params(k, delta, c0, y0, c1d, epsd, om)
            if not p.ok():
                variants.append(p)
    else:
        variants = line_variants(rng, k, delta, c0, y0)
        # always include one variant with the design c1 and a chain-friendly omega0
        for om in (0.97, float(k)):
            for eps in (0.35,):
                p = Params(k, delta, c0, y0, c1d, eps, om)
                if not p.ok():
                    variants.append(p)
    rec['nvar_ok'] = len(variants)
    rec['variants'] = []
    for prm in variants:
        tv = time.time()
        LSb = LineStruct(P, prm)
        LSt = LineStruct(Pr, prm)
        flb = family_links(FFb, LSb, prm)
        flt = family_links(FFt, LSt, prm)
        asm = assembly(LSb, LSt, flb, flt, FFb, FFt, prm, P.W)
        v = dict(prm=prm.d(), floor=summarize_links(flb), ceil=summarize_links(flt), asm=asm,
                 l35_min=[LSb.l35_min, LSt.l35_min], t=time.time() - tv)
        rec['variants'].append(v)
    # ---- negative control 1: drop the fractional constraint of H (w0 -> 0) ----
    if variants:
        import copy
        vbig = max(variants, key=lambda p: p.omega0)
        pn = copy.copy(vbig); pn.w0 = 0.0; pn.h0 = 1.0
        LSn = LineStruct(P, pn)
        rmax = 0.0; p2 = 0
        for rr in FFb.S:
            s = rr['s']; Z = LSn.Zset(s)
            den = (pn.eps * s + pn.bW) * (rr['N'] + FFb.Lam0 + rr['LW'])
            r_ = len(Z) / den if den > 0 else (0.0 if not Z else 1e9)
            rmax = max(rmax, r_)
            p2 += sum(1 for j in Z if j in rr['contacted'])
        rec['neg_w0'] = dict(rmax=rmax, p2v=p2)
    # ---- negative control 2: order "T before M" at T_max squares ----
    if FFb.nR1multi > 0:
        try:
            FFB = FamilyFlows(P, k, delta, c0, y0, s_grid, tlimit=tlim, orderB=True, nov=50)
            rec['neg_orderB'] = dict(symd=max(d['symd'] for d in FFB.S), M=float(FFB.meas0[3]),
                                     T=float(FFB.meas0[4]), M_A=float(FFb.meas0[3]))
        except Exception as ex:
            rec['neg_orderB'] = dict(error=repr(ex))
    rec['t_total'] = time.time() - t0
    rec['rss_mb'] = rss_mb()
    return rec, (P, FFb, FFt)


def main():
    argv = list(sys.argv)
    if '--hill-v1' in argv:                     # campaigns A and B were generated with gen_hill_v1
        argv.remove('--hill-v1')
        from asmx_replay import use_hill_v1
        use_hill_v1()
    camp = argv[1]; seed = int(argv[2]); ncfg = int(argv[3]); tl = float(argv[4])
    fams = argv[5].split(',') if len(argv) > 5 else list(GENS.keys())
    ks = [int(x) for x in argv[6].split(',')] if len(argv) > 6 else [16, 20, 24]
    outdir = os.path.join(HERE, 'out')
    os.makedirs(outdir, exist_ok=True)
    logf = os.path.join(outdir, 'run_log.txt')
    st = mem_guard(logf, camp + ' start')
    if st == 'tiny':
        ncfg = min(ncfg, 2); tl = min(tl, 50)
        with open(logf, 'a') as fh:
            fh.write('%s %s low memory for 40 min -> tiny run\n' % (time.strftime('%H:%M:%S'), camp))
    out = os.path.join(outdir, camp + '.jsonl')
    rng = random.Random(seed)
    T0 = time.time()
    for i in range(ncfg):
        if time.time() - T0 > tl:
            break
        fam = fams[i % len(fams)]
        k = rng.choice(ks)
        if fam == 'random' and k > 20:
            k = 20
        try:
            rec, _ = run_one('%s-%d' % (camp, i), fam, k, rng, tlim=60.0)
        except Exception as ex:
            rec = dict(cid='%s-%d' % (camp, i), fam=fam, k=k, error=repr(ex), tb=traceback.format_exc()[-1500:])
        with open(out, 'a', encoding='utf-8') as fh:
            fh.write(json.dumps(rec, default=lambda o: float(o) if isinstance(o, (np.floating,)) else str(o)) + '\n')
        if rss_mb() > 900:
            with open(logf, 'a') as fh:
                fh.write('%s %s rss>900MB stop\n' % (time.strftime('%H:%M:%S'), camp))
            break
    with open(logf, 'a', encoding='utf-8') as fh:
        fh.write('%s %s done %d configs in %.0f s, rss=%.0f MB\n' % (time.strftime('%H:%M:%S'), camp, i + 1,
                                                                   time.time() - T0, rss_mb()))


if __name__ == '__main__':
    main()
