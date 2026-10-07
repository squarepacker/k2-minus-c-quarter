# -*- coding: utf-8 -*-
"""asmx_summary.py -- aggregate the campaign outputs of asmx_run.py: per link, worst LHS/RHS ratio, violation
counts, slack (the labels are explained in README.md).
usage: python asmx_summary.py [file.jsonl ...]   (default: out/camp*.jsonl; relative names are looked up in the
current directory, then in this script's folder, then in its subfolder out/)"""
import json, sys, glob, os, math
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-9


def resolve(fn):
    for cand in (fn, os.path.join(HERE, fn), os.path.join(HERE, 'out', fn)):
        if os.path.exists(cand):
            return cand
    return fn


files = [resolve(f) for f in sys.argv[1:]] if len(sys.argv) > 1 else sorted(glob.glob(os.path.join(HERE, 'out', 'camp*.jsonl')))
rows = []
for f in files:
    with open(f, encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append((os.path.basename(f), json.loads(line)))

S = defaultdict(lambda: [None, None])      # key -> [max value, where]
C = defaultdict(int)                        # counters


def upd(key, val, where, mode='max'):
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return
    cur = S[key]
    if cur[0] is None or (mode == 'max' and val > cur[0]) or (mode == 'min' and val < cur[0]):
        S[key] = [val, where]


fams = defaultdict(int)
for fn, r in rows:
    C['configs'] += 1
    if 'error' in r:
        C['errors'] += 1; print('ERROR', r.get('cid'), r['error'][:200]); continue
    if 'invalid' in r:
        C['invalid'] += 1; continue
    fams[r['fam']] += 1
    cid = r['cid']
    C['valid'] += 1
    C['anom'] += sum(r['anom'])
    upd('dropped', max(r['dropped']), cid)
    upd('total_measure_err', max(r['total_measure_err']), cid)
    upd('F0_bd', max(r['F0_bd']), cid)
    upd('F0_ov', max(r['F0_ov']), cid)
    upd('F0_p4min', min(r['F0_p4min']), cid, 'min')
    C['merges_cfg'] += 1 if (r['meas0_b'][3] + r['meas0_t'][3]) > 0 else 0
    C['T_cfg'] += 1 if (r['meas0_b'][4] + r['meas0_t'][4]) > 0 else 0
    C['W_cfg'] += 1 if (r['meas0_b'][1] + r['meas0_t'][1]) > 0 else 0
    C['E_cfg'] += 1 if (r['meas0_b'][2] + r['meas0_t'][2]) > 0 else 0
    if 'neg_w0' in r:
        C['neg_w0_runs'] += 1
        C['neg_w0_p2v'] += r['neg_w0']['p2v']
        C['neg_w0_r_gt1'] += 1 if r['neg_w0']['rmax'] > 1 + TOL else 0
        upd('neg_w0_rmax', r['neg_w0']['rmax'], cid)
    if 'neg_orderB' in r and 'symd' in r['neg_orderB']:
        C['neg_orderB_runs'] += 1
        C['neg_orderB_symd_pos'] += 1 if r['neg_orderB']['symd'] > 1e-9 else 0
        upd('neg_orderB_symd', r['neg_orderB']['symd'], cid)
    for vi, v in enumerate(r['variants']):
        where = '%s/v%d(om=%.3g,c1=%.3g,eps=%.2g)' % (cid, vi, v['prm']['omega0'], v['prm']['c1'], v['prm']['eps'])
        C['variants'] += 1
        for side in ('floor', 'ceil'):
            f = v[side]
            w2 = where + '/' + side
            C['cases_s'] += f['cases']; C['cases_Zpos'] += f['cases_Zpos']
            for kk in ('r_paper', 'r_meas', 'r_sharp', 'r_shadow', 'r_nobW', 'lbZ', 'lb1', 'lb2'):
                if f.get(kk) is not None:
                    upd('L1_' + kk, f[kk], w2)
                    if kk not in ('r_nobW',) and f[kk] > 1 + TOL:
                        C['VIOL_L1_' + kk] += 1
                    if kk == 'r_nobW' and f[kk] > 1 + TOL:
                        C['diag_nobW_gt1'] += 1
            C['VIOL_p2v'] += f['p2v']; C['VIOL_outV'] += f['outV']; C['VIOL_tiltv'] += f['tiltv']
            if f['p4min'] is not None:
                upd('L1_p4min', f['p4min'], w2, 'min'); C['VIOL_p4'] += 1 if f['p4min'] < -TOL else 0
            if f['step3min'] is not None:
                upd('L1_step3min', f['step3min'], w2, 'min'); C['VIOL_step3'] += 1 if f['step3min'] < -TOL else 0
            if f['bdmax'] is not None:
                upd('L1_bdmax', f['bdmax'], w2)
            if f['ovmax'] is not None:
                upd('L1_ovmax', f['ovmax'], w2)
            if f['LTminusN'] is not None:
                upd('L1_LTminusN', f['LTminusN'], w2); C['VIOL_LTleN'] += 1 if f['LTminusN'] > TOL else 0
            if f['lam_excess'] is not None:
                upd('L1_lam_excess', f['lam_excess'], w2); C['VIOL_lam'] += 1 if f['lam_excess'] > TOL else 0
            upd('L1_wall_ratio', f['wall_ratio'], w2); C['VIOL_wall'] += 1 if f['wall_ratio'] > 1 + TOL else 0
            upd('L1_nZmax', f['nZmax'], w2)
            # link 2
            upd('L2_symd', f['symd_max'], w2); C['VIOL_L2_symd'] += 1 if f['symd_max'] > 1e-9 else 0
            upd('L2_N_vs_n', f['N_vs_n_max'], w2)
            C['VIOL_L2_mono'] += 0 if f['N_monotone'] else 1
            upd('L2_T1_exact', f['T1_ratio_exact'], w2); C['VIOL_L2_T1'] += 1 if f['T1_ratio_exact'] > 1 + TOL else 0
            upd('L2_T1_grid', f['T1_ratio_grid'], w2)
            if f['intA'] > 0:
                C['T1_nontrivial'] += 1
            # link 3
            for kk in ('T2_ratio', 'T2s_ratio', 'LBint_ratio', 'LBint_ratio_star', 'P7int_ratio_paper',
                       'P7int_ratio_meas', 'chain_ratio_paper', 'chain_ratio_meas', 'chain_ratio_star_meas'):
                upd('L3_' + kk, f[kk], w2)
                if f[kk] > 1 + TOL:
                    C['VIOL_L3_' + kk] += 1
            upd('L3_T2_swap_rel', f['T2_swap_rel'], w2)
            if f['LHS3'] > 0:
                C['L3_f_pos'] += 1
            if f['LHS3star'] > 0:
                C['L3_fstar_pos'] += 1
        a = v['asm']
        for kk in ('ratio4_paper', 'ratio4_sharp', 'kind_i_ratio', 'iii_split_ratio', 'floor_ratio_meas',
                   'ceil_ratio_meas', 'floor_ratio_paper', 'ceil_ratio_paper'):
            upd('L4_' + kk, a[kk], where)
            if a[kk] > 1 + TOL:
                C['VIOL_L4_' + kk] += 1
        # Step 2 of the proof of Theorem 10.3, bound (10.1):
        #   l(H) >= 8 h0 ((j1+1)^{1/4} - (j0+1)^{1/4}) >= 8 h0 (((1-eps)y1)^{1/4} - (y0+1)^{1/4})
        p = v['prm']
        j0 = math.ceil(p['y0']); j1 = math.floor((1 - p['eps']) * p['y1'])
        bnd_a = 8 * p['h0'] * ((j1 + 1) ** 0.25 - (j0 + 1) ** 0.25)
        bnd_b = 8 * p['h0'] * (((1 - p['eps']) * p['y1']) ** 0.25 - (p['y0'] + 1) ** 0.25)
        if a['lH'] > 0:
            upd('L4_eq10_1_ratio', max(bnd_a, bnd_b) / a['lH'], where)
            if max(bnd_a, bnd_b) > a['lH'] * (1 + 1e-12) + 1e-12:
                C['VIOL_L4_eq10_1'] += 1
        upd('L4_Yt_reflect_diff', a['Yt_reflect_diff'], where)
        upd('L4_lH_reflect_diff', a['lH_reflect_diff'], where)
        if a['lYiii'] > 0 and (1 - v['prm']['omega0']) > 0:
            C['L4_iii_nontrivial'] += 1
        if a['fb'] + a['ft'] > 0:
            C['L4_f_pos'] += 1
        for kk in ('A_ii', 'Aprime_meas', 'CLam_meas', 'Gamma_meas'):
            upd('K_' + kk + '_max', a[kk], where)
        upd('l35_min', min(v['l35_min']), where, 'min')
        C['VIOL_l35'] += 1 if min(v['l35_min']) < -1e-9 else 0

print('files', [os.path.basename(f) for f in files])
print('families', dict(fams))
print()
for k in sorted(C):
    print('%-32s %d' % (k, C[k]))
print()
for k in sorted(S):
    v, w = S[k]
    print('%-32s %-14s %s' % (k, ('%.6g' % v) if isinstance(v, (int, float)) else str(v), w))
