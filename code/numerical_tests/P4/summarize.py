"""Aggregate the JSONL results of run_tests.py and ties.py (printed per family and in total).

Usage:  python summarize.py <file.jsonl> [<file.jsonl> ...]"""
import json
import sys
from collections import defaultdict


def main(files):
    agg = defaultdict(lambda: defaultdict(float))
    viols = []
    errors = []
    gl = defaultdict(float)
    gl['min_margin_nonexc'] = 1e9
    labs = defaultdict(int)
    for fn in files:
        for line in open(fn, encoding='utf-8'):
            o = json.loads(line)
            if 'abort' in o:
                errors.append((fn, 'abort', o))
                continue
            if 'error' in o:
                errors.append((fn, o['fam'], o['error'], o.get('tb', '')[-600:]))
                continue
            if 'skip' in o:
                agg[o['fam']]['skipped'] += 1
                continue
            a = agg[o['fam']]
            a['configs'] += 1
            a['squares'] += o['nsq']
            a['axis_sq'] += o['n_axis']
            a['floor_sq'] += o['n_floor']
            a['near45_sq'] += o['n_near45']
            a['pieces'] += o['npieces']
            m = o['meas']
            for kk in ('M', 'T', 'E', 'W', 'D'):
                if m.get(kk, 0) > 0:
                    a['cfg_' + kk + '>0'] += 1
            a['merge_cells'] += o['merge_cells']
            f0 = o['F0']
            a['F0_heights'] += f0['n']
            a['F0_exc_heights'] += f0['n_exc']
            a['F0_exc_ineq_fail'] += f0['n_exc_fail']
            a['F0_zero_margin'] += f0['n_zero_margin']
            if f0['min_margin_nonexc'] is not None:
                a['min_margin_nonexc'] = min(a.get('min_margin_nonexc', 1e9), f0['min_margin_nonexc'])
            a['maxrho'] = max(a['maxrho'], f0['maxrho'])
            a['maxOvg'] = max(a['maxOvg'], f0['maxOvg'])
            if f0['min_pslope'] is not None:
                a['min_pslope'] = min(a.get('min_pslope', 1e9), f0['min_pslope'])
            if f0['min_gslope'] is not None:
                a['min_gslope'] = min(a.get('min_gslope', 1e9), f0['min_gslope'])
            a['n_exc_total'] += o['n_exc']
            for kk, v in o['exc_labels'].items():
                labs[kk] += v
            for s in o['Fs']:
                a['Fs_flows'] += 1
                a['Fs_heights'] += s['n']
                a['Fs_exc_heights'] += s['n_exc']
                a['Fs_exc_ineq_fail'] += s['n_exc_fail']
                a['Fs_LW_pos'] += s['n_LW_pos']
                a['max_LW_ratio'] = max(a['max_LW_ratio'], s['max_LW_ratio'])
                a['max_Lam_ratio'] = max(a['max_Lam_ratio'], s['max_Lam_ratio'])
                if s['min_margin_nonexc'] is not None:
                    a['Fs_min_margin_nonexc'] = min(a.get('Fs_min_margin_nonexc', 1e9), s['min_margin_nonexc'])
                if s['min_N_minus_LT'] is not None:
                    a['min_N_minus_LT'] = min(a.get('min_N_minus_LT', 1e9), s['min_N_minus_LT'])
            p = o['point']
            a['point_n'] += p['n']
            a['point_mismatch'] += p['mismatch']
            for kk, v in p['kinds'].items():
                a['point_' + kk] += v
            nc = o.get('negctl')
            if nc:
                a['neg_cfgs'] += 1
                a['neg_heights'] += nc['n']
                a['neg_Ov_pos'] += nc['n_Ov_pos']
                a['neg_max_Ov'] = max(a['neg_max_Ov'], nc['max_Ov'])
                a['neg_ineq_fail'] += nc['n_ineq_fail_nonexc']
                a['neg_id_bad'] += nc['id_bad']
            a['viol'] += o['nviol']
            a['truncated'] += 1 if o.get('truncated') else 0
            a['time'] += o['time']
            a['max_rss'] = max(a['max_rss'], o.get('rss_mb', 0))
            for v in o['viol']:
                viols.append((o['fam'], o.get('seed'), o.get('i'), v))
    tot = defaultdict(float)
    for fam in sorted(agg):
        a = agg[fam]
        print('==', fam)
        print('   ' + ', '.join('%s=%s' % (kk, (round(v, 6) if isinstance(v, float) else v)) for kk, v in sorted(a.items())))
        for kk in ('configs', 'squares', 'pieces', 'F0_heights', 'F0_exc_heights', 'Fs_flows', 'Fs_heights',
                   'point_n', 'point_mismatch', 'viol', 'neg_heights', 'neg_Ov_pos', 'neg_ineq_fail', 'merge_cells',
                   'F0_exc_ineq_fail', 'Fs_exc_ineq_fail', 'n_exc_total', 'Fs_LW_pos', 'axis_sq', 'floor_sq',
                   'near45_sq', 'time'):
            tot[kk] += a.get(kk, 0)
        for kk in ('cfg_M>0', 'cfg_T>0', 'cfg_E>0', 'cfg_W>0', 'cfg_D>0'):
            tot[kk] += a.get(kk, 0)
    print('== TOTAL')
    print('   ' + ', '.join('%s=%s' % (kk, v) for kk, v in sorted(tot.items())))
    print('== exceptional-height labels (count over configs):', dict(labs))
    print('== violations:', len(viols))
    for v in viols[:30]:
        print('   ', v)
    print('== errors:', len(errors))
    for e in errors[:10]:
        print('   ', e)


if __name__ == '__main__':
    main(sys.argv[1:])
