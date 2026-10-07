"""p7x_summary.py file1.jsonl [file2 ...] -- aggregate the outputs of p7x_run.py and of the targeted scripts:
one line per configuration family (the part of the name before '('), then a line TOTAL.
Relative paths are looked up in the current directory, then in this script's folder, then in its subfolder out/.
The records in results/original/ store 'd_Tchar_symdiff' under an older key name (see README.md); for them the
column 'dB' stays 0 (all recorded values of that quantity are 0)."""
import os, sys, json, collections

HERE = os.path.dirname(os.path.abspath(__file__))


def resolve(fn):
    for cand in (fn, os.path.join(HERE, fn), os.path.join(HERE, 'out', fn)):
        if os.path.exists(cand):
            return cand
    return fn


def main(files):
    rows = []
    errs = 0
    for fn in [resolve(f) for f in files]:
        with open(fn, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                if 'error' in r or 'abort' in r:
                    errs += 1
                    print('ERR', str(r)[:400])
                    continue
                if 'name' not in r:
                    print('INFO', fn, str(r)[:300])
                    continue
                rows.append(r)
    fam = collections.defaultdict(lambda: dict(cfg=0, inval=0, cases=0, nZpos=0, maxZ=0, a=0.0, aZp=0.0, neg_over1=0,
                                                negmax=0.0, p2max=0.0, p2v=0, live=0, passZ=0, F0Z=0, dT=0, dB=0,
                                                cfail=0, cpos=0, cIJ=0.0, cLB=0.0, pf=collections.Counter(), anom=0,
                                                xbad=0, xn=0, xM=0, xR1=0, cons0=0, bad=0, step3=0.0, p4=0.0,
                                                best=None))
    for r in rows:
        key = r['name'].split('(')[0]
        F = fam[key]
        F['cfg'] += 1
        if 'invalid' in r:
            F['inval'] += 1; continue
        F['anom'] += r.get('nanom', 0)
        xc = r.get('xcheck', {})
        F['xbad'] += xc.get('nbad', 0); F['xn'] += xc.get('n', 0); F['xM'] += xc.get('nM', 0); F['xR1'] += xc.get('nR1', 0)
        F['cons0'] += len(r.get('F0_cons_fail', []))
        if r.get('d_Ts_nested') is False:
            F['pf']['d_Ts_not_nested'] += 1
        elif r.get('d_Ts_nested') is True:
            F['pf']['_info_nested_checked'] += 1
        m0 = r.get('F0_tot', {}).get('M', 0)
        if m0 > 0:
            F['pf']['_info_cfg_with_merges'] += 1
            F['pf']['_info_merge_measure_x1000'] += int(1000 * m0)
        for c in r.get('per_s', []):
            F['cases'] += 1
            if c['nZ'] > 0:
                F['nZpos'] += 1
            F['maxZ'] = max(F['maxZ'], c['nZ'])
            if c['a_ratio'] > F['a']:
                F['a'] = c['a_ratio']; F['best'] = (r.get('idx'), r.get('seed'), c['s'], c['nZ'], c['N'], c['L0m'], c['wall'], c['Vlen'], r['name'])
            if c.get('a_ratio_Zprime') is not None:
                F['aZp'] = max(F['aZp'], c['a_ratio_Zprime'])
            ng = c.get('negctl_ratio_Zall')
            if ng is not None:
                F['negmax'] = max(F['negmax'], ng)
                if ng > 1:
                    F['neg_over1'] += 1
            if c.get('p2_worst'):
                F['p2max'] = max(F['p2max'], c['p2_worst'][0])
            F['p2v'] += len(c.get('p2_viol', []))
            F['live'] += len(c.get('b_live_on_Z', []))
            F['passZ'] += len(c.get('b_pass_on_Z', []))
            if c.get('F0_contact_on_Z'):
                F['F0Z'] += 1
            if c.get('d_Tstar_symdiff', 0) != 0 or not c.get('d_Tstar_equal_exact', True):
                F['dT'] += 1
            if c.get('d_Tchar_symdiff', 0) != 0:
                F['dB'] += 1
            if not c.get('c_ok', True):
                F['cfail'] += 1
            if c.get('c_I', 0) > 0:
                F['cpos'] += 1
                if c.get('c_ratio_IJ'):
                    F['cIJ'] = max(F['cIJ'], c['c_ratio_IJ'])
                if c.get('c_ratio_LB_Z'):
                    F['cLB'] = max(F['cLB'], c['c_ratio_LB_Z'])
            if c['nZ'] > 0 and c.get('extra_constraint_holds') is False:
                F['pf']['_info_extra_constraint_violated_with_Z'] += 1
            if c.get('zt_over_alpha') is not None and c['zt_over_alpha'] >= 1:
                F['pf']['_info_Z_tilt_ge_alpha_s'] += 1
            for pf in c.get('pointwise_fail', []):
                F['pf'][pf[0]] += 1
            if c.get('bad_pieces'):
                F['bad'] += 1
            pw = c.get('pointwise_worst', {})
            if pw.get('step3') is not None:
                F['step3'] = max(F['step3'], pw['step3'])
            if pw.get('p4') is not None:
                F['p4'] = max(F['p4'], pw['p4'])
    print('rows', len(rows), 'errors', errs)
    for key, F in sorted(fam.items()):
        F = dict(F); F['pf'] = dict(F['pf'])
        print(key, json.dumps(F))
    # overall totals
    T = collections.Counter()
    nondecr = 0
    for r in rows:
        if 'invalid' in r:
            T['invalid'] += 1; continue
        T['configs'] += 1
        T['squares_total'] += r.get('nsq', 0)
        T['F0_pieces_total'] += r.get('nrec', 0)
        T['near_tie_alpha_max'] += len(r.get('near_alpha_max', []))
        for c in r.get('per_s', []):
            T['near_tie_alpha_s'] += len(c.get('near_alpha_s', []))
        T['F0_cons_fail'] += len(r.get('F0_cons_fail', []))
        if not r['params'].get('drift_decreasing', True):
            nondecr += 1
        T['xcheck_samples'] += r.get('xcheck', {}).get('n', 0)
        T['xcheck_bad'] += r.get('xcheck', {}).get('nbad', 0)
        T['R1_checks'] += r.get('xcheck', {}).get('nR1', 0)
        T['anomalies'] += r.get('nanom', 0)
        for c in r.get('per_s', []):
            T['cases'] += 1
            T['heights'] += c.get('npoint', 0)
            T['p2_checks'] += c.get('p2_ncheck', 0)
            T['Z_total'] += c['nZ']
            if c['nZ'] > 0:
                T['cases_Z'] += 1
            if c.get('N', 0) > 0:
                T['cases_N_pos'] += 1
            if c.get('F0_contact_on_Z'):
                T['cases_F0_reaches_Z'] += 1
            if c.get('c_I', 0) > 0:
                T['cases_I_pos'] += 1
            if c['nZ'] > 0 and c.get('extra_constraint_holds') is False:
                T['cases_Z_extra_constraint_violated'] += 1
            if c.get('zt_over_alpha') is not None and c['zt_over_alpha'] >= 1:
                T['cases_Z_tilt_ge_alpha'] += 1
            T['viol_step1_Z_not_in_V'] += int(not c.get('Z_inside_V', True))
            if r['params'].get('omega0', 99) < 1:
                T['cases_omega0_lt_1'] += 1
                if c['nZ'] > 0:
                    T['cases_omega0_lt_1_with_Z'] += 1
            T['viol_a'] += int(c['a_ratio'] > 1)
            T['viol_aZp'] += int((c.get('a_ratio_Zprime') or 0) > 1)
            T['negctl_gt1'] += int((c.get('negctl_ratio_Zall') or 0) > 1)
            T['viol_p2'] += len(c.get('p2_viol', []))
            T['viol_live_on_Z'] += len(c.get('b_live_on_Z', []))
            T['viol_pass_on_Z'] += len(c.get('b_pass_on_Z', []))
            T['viol_c'] += int(not c.get('c_ok', True))
            T['viol_d_Tstar'] += int(c.get('d_Tstar_symdiff', 0) != 0 or not c.get('d_Tstar_equal_exact', True))
            T['pointwise_fail'] += len([p for p in c.get('pointwise_fail', []) if p[0] != 'degenerate_height'])
            T['degenerate_heights'] += len([p for p in c.get('pointwise_fail', []) if p[0] == 'degenerate_height'])
            T['bad_pieces'] += int(bool(c.get('bad_pieces')))
        if r.get('d_Ts_nested') is False:
            T['viol_nested'] += 1
    print('TOTAL', json.dumps(dict(T)), 'w0_sup_nonmonotone_configs', nondecr)


if __name__ == '__main__':
    main(sys.argv[1:])
