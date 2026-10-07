"""Runner for one configuration family.

usage: python run_family.py <family> [time_limit_sec] [start]

  family          smoke, valley, bigvalley, merge, fan, three, near45, mechA, strip, jam, rowjam, regime, tri, mrand
  time_limit_sec  default 1200; when it is used up the runner prints "NEXT <index>" and stops
  start           index of the first configuration to run (default 0); a long family is run in chunks,
                  each chunk starting at the NEXT index printed by the previous one
  writes out/res_<family>_<start>.jsonl (one JSON record per configuration; overwritten).
  The family rowjam keeps its configurations in out/cache_rowjam.pkl, so that later chunks use the same
  configurations."""
import sys, json, time, math, random, traceback
from fractions import Fraction as Fr
from tracer import trace, check_config, TimeUp, out_path
from analysis import analyze
import configs as C

FAM = sys.argv[1]; TL = float(sys.argv[2]) if len(sys.argv) > 2 else 1200
START = int(sys.argv[3]) if len(sys.argv) > 3 else 0
CNT = [0]
T0 = time.time()
OUT = open(out_path(f'res_{FAM}_{START}.jsonl'), 'w')
DELTAS = [Fr(1, 10**5), Fr(1, 10**4), Fr(1, 10**3), Fr(1, 100), Fr(1, 20), Fr(1, 10)]
BIG = Fr(10)  # sF meaning: no tilt termination (|s|<1 always)

def sF_of(alpha):
    return Fr(math.tan(alpha / 2)).limit_denominator(10**12)

def run(name, sqs, k, delta, sF, extra=None):
    CNT[0] += 1
    if CNT[0] - 1 < START: return None
    if time.time() - T0 > TL:
        print('NEXT', CNT[0] - 1, flush=True); raise SystemExit(0)
    err = check_config(sqs, k)
    rec_out = dict(family=FAM, name=name, k=k, delta=str(delta), N=len(sqs), sF=str(sF), extra=extra or {})
    if err:
        rec_out['skipped'] = err; OUT.write(json.dumps(rec_out) + '\n'); OUT.flush(); return None
    h = Fr(k, 2) - 1
    t1 = time.time()
    try:
        rec = trace(sqs, k, delta, h, sF, tlimit=min(300, max(5, TL - (time.time() - T0))))
        res = analyze(rec, sqs, k, delta, h, sF, tlimit=min(300, max(5, TL - (time.time() - T0))))
    except TimeUp:
        rec_out['timeout'] = True; OUT.write(json.dumps(rec_out) + '\n'); OUT.flush(); return None
    except Exception as e:
        rec_out['error'] = traceback.format_exc()[-800:]; OUT.write(json.dumps(rec_out) + '\n'); OUT.flush()
        print('ERROR', name, rec_out['error'], flush=True); return None
    rec_out.update(res); rec_out['sec'] = time.time() - t1
    rec_out['squares'] = [(str(S.c[0]), str(S.c[1]), str(S.s)) for S in sqs] if len(sqs) <= 12 else None
    OUT.write(json.dumps(rec_out, default=str) + '\n'); OUT.flush()
    flags = []
    for key in ('n_rect_violations', 'same_source_overlaps', 'floor_square_overlaps', 'triple_regions', 'floor_merges', 'Ov_support_violations'):
        if res.get(key): flags.append(f'{key}={res[key]}')
    for key in ('merge_violations', 'E_violations', 'Ov_violations', 'dist_violations'):
        if res.get(key): flags.append(f'{key}={len(res[key])}')
    if not res['conservation_ok']: flags.append('CONSERVATION')
    if not res['death_ok']: flags.append('DEATH')
    if not res['intMg_ok']: flags.append('INTMG')
    if not res['jacobian_ok']: flags.append('JACOBIAN')
    if res['n_anom']: flags.append(f"anom={res['n_anom']}:{res['anom'][:2]}")
    print(f"{name} d={float(delta):g} N={len(sqs)} Mg={res['Mg_max']} dbl={res['n_double_regions']} pairs={res['n_pairs']} "
          f"rectA={res['rect_qa_over_we_max']:.3g} rectB={res['rect_tR_over_delta_max']:.3g} star={res['star_ratio_max']:.4g} "
          f"dist={res['dist_ratio_max']:.3g} mrg={res['merge_ratio_max']:.3g} Ypp={res['max_Y_per_pair']} E={res['E_ratio_max']:.6g} "
          f"Ov={res['Ov_ratio_max']:.3g} death={res['death_ratio_vs_waste_plus_rects']:.3g} M={res['M']:.3g} D={res['D']:.3g} "
          f"t={rec_out['sec']:.1f}s {' '.join(flags)}", flush=True)
    return res

def fam_smoke():
    sq = C.valley(0.1, Fr(1, 100))
    run('valley_smoke', sq, 6, Fr(1, 100), BIG)

def fam_valley():
    for theta in [1e-8, 1e-6, 1e-4, 1e-2, 0.05, 0.2, 0.5]:
        for delta in DELTAS:
            for gf in [Fr(0), Fr(1, 2), Fr(9, 10)]:
                if gf == 0: gfx = Fr(1, 10**6)
                else: gfx = gf
                sq = C.valley(theta, delta, gap_frac=gfx)
                run(f'valley_th{theta:g}_g{float(gfx):g}', sq, 6, delta, BIG, dict(theta=theta))
            sq = C.valley(theta, delta, base=True)
            run(f'valleybase_th{theta:g}', sq, 6, delta, BIG, dict(theta=theta))

def fam_bigvalley():
    # outside the P5 regime: large theta, alpha_F large enough that sources pass; delta below/above 0.5 cos(theta)
    for theta in [0.3, 0.6, 0.9, 1.2, 1.4, 1.5]:
        hc = 0.5 * math.cos(theta)
        for dl in sorted(set([Fr(1, 100), Fr(hc * 0.9).limit_denominator(10**6), Fr(hc * 1.1).limit_denominator(10**6),
                              Fr(min(0.45, 2 * hc)).limit_denominator(10**6), Fr(1, 10)])):
            for gf in [Fr(1, 10**6), Fr(1, 2), Fr(9, 10)]:
                sq = C.valley(theta, dl, gap_frac=gf)
                run(f'bigvalley_th{theta:g}_g{float(gf):g}', sq, 6, dl, BIG, dict(theta=theta, half_cos=hc))

def fam_merge():
    for theta in [1e-6, 1e-4, 1e-2, 0.1, 0.4]:
        for delta in DELTAS:
            for phiY in [0.0, theta / 2, 0.3, -0.3, 0.7, -0.7, 0.785]:
                for tau0 in [Fr(-9, 20), Fr(-1, 4), Fr(0), Fr(9, 20)]:
                    for ab in [(Fr(1, 2), Fr(1, 2)), (Fr(1, 10), Fr(1, 10)), (Fr(9, 10), Fr(9, 10))]:
                        sq = C.merge_cfg(theta, delta, phiY, tau0, a_frac=ab[0], b_frac=ab[1])
                        for sF in (BIG, sF_of(0.5)):
                            run(f'merge_th{theta:g}_pY{phiY:g}_t{float(tau0):g}_ab{float(ab[0]):g}_sF{float(sF):.3g}', sq, 6, delta, sF,
                                dict(theta=theta, phiY=phiY))

def fam_fan():
    for delta in DELTAS:
        for phi in [1e-6, 1e-3, 0.05, 0.3]:
            for seed in range(2):
                sq = C.rows_cfg(6, delta, [[0.0], [phi, -phi], [-phi, phi, 0.0]], seed)
                for sF in (BIG, sF_of(1.5e-6), sF_of(0.5)):
                    run(f'fan_p{phi:g}_s{seed}_sF{float(sF):.3g}', sq, 6, delta, sF, dict(phi=phi))

def fam_three():
    # three sources under one bottom: base row tilts (-p, 0, +p) tight, then Y dropped over
    for delta in DELTAS:
        for phi in [1e-6, 1e-3, 0.05, 0.2]:
            for seed in range(3):
                sq = C.rows_cfg(6, delta, [[-phi, 0.0, phi], [0.0, phi / 2, -phi / 2]], seed, spacing=1.0)
                run(f'three_p{phi:g}_s{seed}', sq, 6, delta, BIG, dict(phi=phi))

def fam_near45():
    for delta in DELTAS:
        for seed in range(3):
            def samp(r):
                return r.choice([0.0, 1e-6, -1e-6, 0.785, -0.78, 0.76, r.uniform(-0.01, 0.01)])
            sq = C.jam(22, 6, delta, samp, seed)
            for a in (1e-3, 0.7, math.pi / 4 - 0.02, None):
                sF = BIG if a is None else sF_of(a)
                run(f'near45_s{seed}_a{a}', sq, 6, delta, sF, dict(alphaF=a))

def fam_mechA():
    for delta in DELTAS:
        for phi in [1e-6, 1e-3, 0.02, 0.1]:
            sq = C.mechA(6, delta, phi, 3, 7)
            for sF in (BIG, sF_of(phi / 2)):
                run(f'mechA_p{phi:g}_sF{float(sF):.3g}', sq, 6, delta, sF, dict(phi=phi))

def fam_strip():
    for delta in DELTAS:
        for phi in [0.01, 0.1, 0.3]:
            for b in (4, 6):
                sq = C.strip_cfg(8, delta, phi, b, 11)
                for sF in (BIG, sF_of(0.005)):
                    run(f'strip_p{phi:g}_b{b}_sF{float(sF):.3g}', sq, 8, delta, sF, dict(phi=phi, b=b))

def fam_jam():
    seed = 0
    while time.time() - T0 < TL:
        for delta in DELTAS:
            seed += 1
            mode = seed % 3
            if mode == 0:
                samp = lambda r: r.choice([r.uniform(-1.5e-6, 1.5e-6), 0.0, r.uniform(-0.01, 0.01), r.uniform(-0.3, 0.3)])
                sF = sF_of(1.5e-6)
            elif mode == 1:
                samp = lambda r: r.uniform(-0.2, 0.2)
                sF = BIG
            else:
                samp = lambda r: r.uniform(-0.6, 0.6)
                sF = sF_of(0.6)
            sq = C.jam(26, 6, delta, samp, 1000 + seed, gapmax_factor=1.5)
            run(f'jam_m{mode}_s{seed}', sq, 6, delta, sF, dict(mode=mode, seed=seed))

def fam_rowjam():
    import pickle, os
    cf = out_path('cache_rowjam.pkl')
    cache = pickle.load(open(cf, 'rb')) if os.path.exists(cf) else {}
    seed = 0
    while True:
        for delta in DELTAS:
            seed += 1
            mode = seed % 5
            key = (seed, str(delta))
            if mode == 0:
                samp = lambda r: r.uniform(-1.5e-6, 1.5e-6) if r.random() < 0.85 else r.uniform(-0.3, 0.3)
                sF = sF_of(1.5e-6)
            elif mode == 1:
                samp = lambda r: r.uniform(-0.05, 0.05); sF = BIG
            elif mode == 2:
                samp = lambda r: r.uniform(-0.3, 0.3); sF = BIG
            elif mode == 3:
                samp = lambda r: r.choice([1, -1]) * r.uniform(0.6, 0.785); sF = BIG
            else:
                cnt = [0]; mag = [None]
                def samp(r, cnt=cnt, mag=mag):
                    cnt[0] += 1
                    if mag[0] is None: mag[0] = 10 ** r.uniform(-6, math.log10(0.3))
                    return mag[0] * (1 if cnt[0] % 2 else -1) * r.uniform(0.5, 1.0)
                sF = BIG
            if key in cache:
                sq = [C.Sq(i, Fr(a), Fr(b), Fr(c)) for i, (a, b, c) in enumerate(cache[key])]
            else:
                CNT[0] += 0
                if CNT[0] < START and len(cache) > 0 and key not in cache:
                    pass
                sq = C.rowjam(6, delta, 3, samp, 5000 + seed)
                cache[key] = [(str(S.c[0]), str(S.c[1]), str(S.s)) for S in sq]
                pickle.dump(cache, open(cf, 'wb'))
            run(f'rowjam_m{mode}_s{seed}', sq, 6, delta, sF, dict(mode=mode, seed=seed))

def fam_regime():
    """exact in-regime clusters: alpha_max = 1.5e-6, delta in [1e-5 .. 1e-2] (and 0.1), several valleys + merging Y's."""
    rnd = random.Random(777)
    i = 0
    while True:
        i += 1
        delta = rnd.choice([Fr(1, 10**5), Fr(1, 10**5), Fr(1, 10**4), Fr(1, 10**3), Fr(1, 100), Fr(1, 10)])
        amax = rnd.choice([1.5e-6, 1.5e-6, 1e-4, 1e-2])
        sq = C.regime_cluster(6, delta, amax, rnd, nrow=rnd.randint(2, 5))
        sF = rnd.choice([sF_of(amax), BIG])
        run(f'regime{i}_a{amax:g}', sq, 6, delta, sF, dict(amax=amax))

def fam_tri():
    """outside the P5 regime: tilts near +-45 deg, no tilt termination, delta in {0.05, 0.1}: look for M_g = 3."""
    import pickle, os
    seed = 0
    while True:
        seed += 1
        delta = Fr(1, 10) if seed % 2 else Fr(1, 20)
        samp = lambda r: r.choice([1, -1]) * r.uniform(0.70, 0.785)
        sq = C.rowjam(6, delta, 3, samp, 9000 + seed, gapmax=1.0)
        run(f'tri_s{seed}', sq, 6, delta, BIG, dict(seed=seed))

def fam_mrand():
    """adversarial random merges: valley + Y with random phase/offset; eps (floor gap) tiny."""
    rnd = random.Random(20261006)
    i = 0
    while True:
        i += 1
        theta = 10 ** rnd.uniform(-8, math.log10(0.5))
        delta = rnd.choice(DELTAS)
        phiY = rnd.choice([rnd.uniform(-0.785, 0.785), theta / 2, -theta / 2, 0.0])
        tau0 = Fr(rnd.uniform(-0.499, 0.499)).limit_denominator(1000)
        a = Fr(rnd.random()).limit_denominator(1000); b = Fr(rnd.uniform(0.01, 0.999)).limit_denominator(1000)
        gf = Fr(rnd.choice([1e-6, rnd.random()])).limit_denominator(10**7)
        eps = delta * Fr(rnd.choice([1, 10, 100, 1000]), 10000)
        sq = C.merge_cfg(theta, delta, phiY, tau0, a_frac=a, b_frac=b, gap_frac=gf, eps=eps)
        if check_config(sq, 6):
            CNT[0] += 1
            if CNT[0] - 1 >= START and time.time() - T0 > TL:
                print('NEXT', CNT[0] - 1, flush=True); raise SystemExit(0)
            continue
        sF = rnd.choice([BIG, sF_of(abs(phiY) * 0.9 + 1e-9)])
        run(f'mrand{i}_th{theta:.3g}_pY{phiY:.3g}', sq, 6, delta, sF, dict(theta=theta, phiY=phiY, tau0=str(tau0), a=str(a), b=str(b), gf=str(gf), eps=str(eps)))

{'regime': fam_regime, 'tri': fam_tri, 'rowjam': fam_rowjam, 'mrand': fam_mrand, 'smoke': fam_smoke, 'valley': fam_valley, 'bigvalley': fam_bigvalley, 'merge': fam_merge, 'fan': fam_fan,
 'three': fam_three, 'near45': fam_near45, 'mechA': fam_mechA, 'strip': fam_strip, 'jam': fam_jam}[FAM]()
print('DONE', time.time() - T0, flush=True)

