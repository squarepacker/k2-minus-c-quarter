"""p7x_eval.py -- evaluate one configuration: flow F^0, F_s for several s, all P7 checks (a)-(d):
  (a) Theorem 9.1(a): |Z(s)| <= (eps s + b_W)(N(s) + Lambda*(s))   ('a_ratio' = left side / right side), with the
      steps of its proof checked pointwise at sampled heights: Step 1 ('Z_inside_V'), Step 3 ('step3'),
      Step 4 ('P4star' = the shadow inequality of Theorem 6.9(f), 'LT_le_N', and 'DEM_le_L0', 'W_le_wall',
      'Ovg_mono' of Proposition 6.10, 'chain'); also the bookkeeping identity and the shadow identity of
      Theorem 6.9(b),(c) ('cons', 'ident');
  (b) Step 2 of that proof: no path of F_s makes an alive contact with relint bot(Z), Z in Z(s) ('b_live_on_Z'),
      and the width test of Theorem 4.11 with w(s, Z) of Lemma 9.3(e) at every alive contact ('p2_worst', 'p2_viol');
  (c) Theorem 9.2(a): I <= J <= |Z(s)| beta-hat((1-eps)s) and the resulting lower bound ('c_ok', 'c_ratio_*');
      its Step 1 uses [R, Lemma 3.5] ('lem35');
  (d) Lemma 9.6(c): {T* <= s} equals the set of x whose path is T_s-terminated in F_s ('d_Tstar_symdiff'),
      which equals the set characterised in Lemma 9.6(b),(c) ('d_Tchar_symdiff'), and T_s is nested in s
      ('d_Ts_nested'); Lemma 9.7: L_T <= N(s) ('LT_le_N').
Negative control: 'negctl_ratio_Zall' counts all squares inside V_s instead of Z(s); it may exceed 1."""
from fractions import Fraction as Fr
import random, time
import mpmath
from p7x_core import Flow, validate, trace_point, pev, mpf_fr
from p7x_meas import (Params, trunc_map, fs_pieces, recs_as_pieces, tstar_sets, height_stats, Lines,
                      p2_test, merge, measure, symdiff_measure)


def F(x):
    return float(x)


def crosscheck(fl, squares, P, n, rng):
    """compare beam records with the independent point tracer (no merging) at sampled x."""
    recs = list(fl.recs)
    rng.shuffle(recs)
    recs = [r for r in recs if r.xhi - r.xlo > 0][:n]
    bad = []; nm = 0; nr1 = 0
    for r in recs:
        x = (r.xlo + r.xhi) / 2
        typ, passed, endp, g, contacts = trace_point(x, squares, P.k, P.delta, P.alpha_max, P.hmax)
        rp = [sid for sid in r.seg if sid != -1]
        q = pev(r.way[-1], x)
        if r.typ == 'M':
            nm += 1
            j = len(rp)
            ok = passed[:j] == rp and len(contacts) > j and contacts[j][0] == r.tsq and contacts[j][1] == q
            if ok:
                gx = contacts[j][2]
                # R1: a recorded winner reaches the same point with lexicographically smaller (g, x)
                found = False
                for w in fl.recs:
                    for (sid, idx) in w.ent:
                        if sid != r.tsq:
                            continue
                        P_ = w.way[idx]
                        if P_[1] != 0:
                            xw = (q[0] - P_[0]) / P_[1]
                        elif P_[3] != 0:
                            xw = (q[1] - P_[2]) / P_[3]
                        else:
                            continue
                        if not (w.xlo < xw < w.xhi):
                            continue
                        if pev(P_, xw) != q:
                            continue
                        t2 = trace_point(xw, squares, P.k, P.delta, P.alpha_max, P.hmax)
                        gw = None
                        for (cs, cq, cg) in t2[4]:
                            if cs == r.tsq and cq == q:
                                gw = cg
                        if gw is not None and (gw, xw) < (gx, x):
                            found = True
                        break
                    if found:
                        break
                nr1 += 1
                if not found:
                    bad.append(('R1_no_winner', F(x), r.tsq))
            else:
                bad.append(('M_path_mismatch', F(x), r.tsq, passed[:6], rp[:6]))
        else:
            ok = (typ == r.typ) and passed == rp and endp == q
            if not ok:
                bad.append(('mismatch', F(x), r.typ, typ, rp[:6], passed[:6]))
    return dict(n=len(recs), nM=nm, nR1=nr1, bad=bad[:8], nbad=len(bad))


def sample_heights(lo, hi, n, rng, extra=()):
    ys = []
    for j in range(n):
        u = Fr(rng.randrange(1, 10 ** 6), 10 ** 6 + 3)
        ys.append(lo + (hi - lo) * (Fr(j) + u) / n)
    for (a, b) in extra:
        u = Fr(rng.randrange(1, 10 ** 6), 10 ** 6 + 7)
        ys.append(a + (b - a) * u)
    return ys


def evaluate(name, squares, P, s_list, rng, nys=8, xcheck_n=40, ovg_grid=30, tlimit=300):
    t0 = time.time()
    out = dict(name=name, params=P.desc(), nsq=len(squares))
    prob = validate(squares, P.k)
    if prob:
        out['invalid'] = [str(p) for p in prob[:5]]
        return out
    fl = Flow(squares, P.k, P.delta, P.alpha_max, P.hmax).run()
    out['nrec'] = len(fl.recs); out['ncand'] = len(fl.cands); out['nray'] = fl.nray
    out['anom'] = [str(a) for a in fl.anom[:10]]; out['nanom'] = len(fl.anom)
    out['near_alpha_max'] = fl.near
    tot = {t: sum((r.xhi - r.xlo for r in fl.recs if r.typ == t), Fr(0)) for t in 'DWEMTH'}
    out['F0_tot'] = {t: F(v) for t, v in tot.items()}
    L0m = tot['D'] + tot['M'] + tot['E']
    p0 = recs_as_pieces(fl)
    # sup Ov_gap of F^0 on a grid (lower estimate of the sup) + conservation on F^0
    supOvg = Fr(0); cons0 = []
    for y in sample_heights(Fr(0), P.hmax, ovg_grid, rng):
        hs = height_stats(p0, y, fl, squares, want_gap_check=False)
        supOvg = max(supOvg, hs['Ovgap'])
        if hs['cons'] != 0 or hs['ident'] != 0 or hs['Ov'] != 0:
            cons0.append((F(y), F(hs['cons']), F(hs['ident']), F(hs['Ov'])))
    out['supOvgap_grid'] = F(supOvg); out['F0_cons_fail'] = cons0[:5]
    out['xcheck'] = crosscheck(fl, squares, P, xcheck_n, rng)
    lines = Lines(squares, P)
    res = []
    for s in s_list:
        if time.time() - t0 > tlimit:
            out['timeout'] = True; break
        s = Fr(s)
        r = dict(s=F(s))
        tm, near = trunc_map(fl, P, s)
        r['near_alpha_s'] = near
        fsp = fs_pieces(fl, P, s, tm)
        bad = [p.typ for p in fsp if p.typ.startswith('BAD')]
        r['bad_pieces'] = bad[:3]
        ftot = {}
        for p in fsp:
            ftot[p.typ] = ftot.get(p.typ, Fr(0)) + (p.xhi - p.xlo)
        r['Fs_tot'] = {t: F(v) for t, v in ftot.items()}
        N = ftot.get('T', Fr(0))
        r['N'] = F(N)
        # (d) T* set equality
        Tw = merge([(p.xlo, p.xhi) for p in fsp if p.typ == 'T'])
        A, B = tstar_sets(fl, P, s, tm)
        r['d_Tstar_symdiff'] = F(symdiff_measure(Tw, A)); r['d_Tchar_symdiff'] = F(symdiff_measure(Tw, B))
        r['d_Tstar_equal_exact'] = (merge(Tw) == merge(A))
        # lines / Z(s)
        ivs, pts, bp = lines.Yb_W(s)
        Z = lines.Zset(s, ivs, pts, bp)
        r['nZ'] = len(Z)
        al = P.alpha(s)
        r['zt_over_alpha'] = max(float(fl.sq[z].a / al) for z in Z) if Z else None
        r['extra_constraint_holds'] = bool(P.beta_hat((1 - P.eps) * s) < al)
        r['Yb_W_measure'] = F(sum(((b - a) for a, b in ivs), Fr(0)))
        # (b) P2 and no live contact on Z(s)
        p2 = p2_test(fl, P, s, tm)
        r['p2_worst'] = [F(v) if not isinstance(v, int) else v for v in p2['worst']] if p2['worst'] else None
        r['p2_viol'] = p2['viol'][:5]; r['p2_ncheck'] = p2['ncheck']
        r['b_live_on_Z'] = [sid for sid in Z if sid in p2['live']]
        r['b_pass_on_Z'] = sorted({sg for p in fsp for sg in p.seg if sg != -1 and sg in Z})
        # F^0 reaches Z(s)?  (shows the test has teeth: F^0 may, F_s must not)
        f0live = {cd.sid for cd in fl.cands}
        r['F0_contact_on_Z'] = sorted(z for z in Z if z in f0live)
        # (a) main inequality
        Vlo, Vhi = P.Vs(s); Vlen = P.Vlen(s); wall = P.wall(s)
        den = Vlen * (mpf_fr(N) + mpf_fr(L0m) + wall)
        denf = Vlen * (mpf_fr(N) + mpf_fr(L0m) + mpf_fr(supOvg) + wall)
        r['a_ratio'] = float(len(Z) / den) if den > 0 else (0.0 if not Z else float('inf'))
        r['a_ratio_withOvg'] = float(len(Z) / denf) if denf > 0 else None
        inV = [S for S in squares if mpf_fr(S.ymin) >= Vlo and mpf_fr(S.ymax) <= Vhi]
        Zp = [S for S in inV if S.id not in p2['live']]
        r['nZprime'] = len(Zp); r['nZall'] = len(inV)
        r['a_ratio_Zprime'] = float(len(Zp) / den) if den > 0 else None
        r['negctl_ratio_Zall'] = float(len(inV) / den) if den > 0 else None
        r['Z_inside_V'] = all(mpf_fr(S.ymin) >= Vlo and mpf_fr(S.ymax) <= Vhi for S in squares if S.id in Z)
        r['N'] = F(N); r['L0m'] = F(L0m); r['wall'] = float(wall); r['Vlen'] = float(Vlen)
        # (c) lower bound
        I = lines.integral_I(ivs)
        J = lines.J_Z(Z, ivs)
        bh = P.beta_hat((1 - P.eps) * s)
        bnd = len(Z) * bh
        LB = mpf_fr((1 - P.eps) * s) ** mpmath.mpf(0.75) / mpf_fr(P.c1) * mpf_fr(I)
        r['c_I'] = F(I); r['c_J'] = F(J); r['c_bound'] = float(bnd)
        r['c_LB'] = float(LB)
        r['c_ok'] = bool(mpf_fr(I) <= mpf_fr(J) + mpmath.mpf('1e-30') and mpf_fr(J) <= bnd + mpmath.mpf('1e-30') and LB <= len(Z) + mpmath.mpf('1e-25'))
        r['c_ratio_IJ'] = F(I / J) if J > 0 else None
        r['c_ratio_LB_Z'] = float(LB / len(Z)) if Z else None
        # pointwise chain on sample heights in V_s
        Vlo_r = Fr(int(mpmath.ceil(Vlo * 10 ** 9)), 10 ** 9); Vhi_r = Fr(int(mpmath.floor(Vhi * 10 ** 9)), 10 ** 9)
        Vlo_r = max(Vlo_r, Fr(1, 10 ** 6))
        ys = sample_heights(Vlo_r, Vhi_r, nys, rng, extra=ivs[:4])
        worst = dict(step3=None, p4=None, chain=None, LT=None)
        fails = []
        for y in ys:
            hs = height_stats(fsp, y, fl, squares)
            hz = height_stats(p0, y, fl, squares, want_gap_check=False)
            cZ = sum((S.chord(y) for S in squares if S.id in Z), Fr(0))
            L = hs['loss']
            lam_y = L['D'] + L['E'] + L['M'] + L['W'] + hs['Ovgap']
            chk = {
                'cons': hs['cons'] == 0, 'ident': hs['ident'] == 0, 'Ov0': hs['Ov'] == 0,
                'gap_in_waste': hs.get('gap_in_square', 0) == 0,
                'muZ0': all(hs['mu'].get(z, 0) == 0 for z in Z),
                'step3': cZ <= hs['Sh'],
                'P4star': hs['Sh'] <= L['T'] + lam_y,
                'LT_le_N': L['T'] <= N,
                'DEM_le_L0': L['D'] + L['E'] + L['M'] <= L0m,
                'W_le_wall': mpf_fr(L['W']) <= wall,
                'Ovg_mono': hs['Ovgap'] <= hz['Ovgap'],
                'chain': mpf_fr(cZ) <= mpf_fr(N) + mpf_fr(L0m) + mpf_fr(max(supOvg, hz['Ovgap'])) + wall,
            }
            if 0 < y < P.k:
                om = hs['omega']; Ey = lines.E(y); ss = lines.short_sum(y)
                chk['lem35'] = ss >= 1 - om - Ey      # [R, Lemma 3.5]
            for kk, v in chk.items():
                if not v:
                    fails.append((kk, F(y)))
            if hs['degen']:
                fails.append(('degenerate_height', F(y), hs['degen']))
            rs = F(cZ / hs['Sh']) if hs['Sh'] > 0 else None
            rp = F(hs['Sh'] / (L['T'] + lam_y)) if (L['T'] + lam_y) > 0 else None
            for key, val in (('step3', rs), ('p4', rp)):
                if val is not None and (worst[key] is None or val > worst[key]):
                    worst[key] = val
        r['pointwise_fail'] = fails[:8]; r['npoint'] = len(ys); r['pointwise_worst'] = worst
        r['_Tw'] = Tw
        res.append(r)
    # Lemma 9.6(c): T_s nested in s  (N non-decreasing)
    srt = sorted(res, key=lambda z: z['s'])
    nest = True
    for a_, b_ in zip(srt, srt[1:]):
        if symdiff_measure(b_['_Tw'], merge(a_['_Tw'] + b_['_Tw'])) != 0 or a_['N'] > b_['N']:
            nest = False
    out['d_Ts_nested'] = nest
    for z in res:
        z.pop('_Tw', None)
    out['per_s'] = res
    out['time'] = time.time() - t0
    return out
