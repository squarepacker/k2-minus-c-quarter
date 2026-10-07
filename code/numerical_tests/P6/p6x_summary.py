# p6x_summary.py -- aggregate the JSON-lines output of p6x_run.py.
# usage: python p6x_summary.py file1.jsonl [file2.jsonl ...]
#   relative paths are looked up first in the current directory, then in this script's folder,
#   then in its subfolder out/.
# The records in results/original/ store the minimum of the weak form of Lemma 8.7(c) under an older key
# name (see README.md); for them the line 'MIN L6_weak_min' is simply absent.
import os, sys, json, glob
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))


def resolve(fn):
    for cand in (fn, os.path.join(HERE, fn), os.path.join(HERE, 'out', fn)):
        if os.path.exists(cand):
            return cand
    return fn


files = [resolve(f) for f in sys.argv[1:]]
agg = Counter()
mins = {}
maxs = {}
viol = []
kinds = Counter()
by_delta = Counter()


def upd_min(k, v, tag):
    if v is None:
        return
    if k not in mins or v < mins[k][0]:
        mins[k] = (v, tag)


def upd_max(k, v, tag):
    if v is None:
        return
    if k not in maxs or v > maxs[k][0]:
        maxs[k] = (v, tag)


for fn in files:
    for line in open(fn, encoding='utf-8'):
        try:
            r = json.loads(line)
        except Exception:
            agg['partial_lines'] += 1
            continue
        if 'error' in r:
            agg['errors'] += 1
            viol.append((r['tag'], r['error']))
            continue
        if r.get('invalid'):
            agg['invalid'] += 1
            continue
        agg['configs'] += 1
        agg['squares'] += r['N']
        kinds[r['tag'].split('_')[0]] += 1
        by_delta[(r['delta'], round(r['amax'], 8))] += 1
        agg['L7_tests'] += r['L7_tests']
        agg['L7_polys'] += r['L7_polys']
        if r['nviol']:
            agg['configs_with_viol'] += 1
            viol.append((r['tag'], r['viol']))
        tag = r['tag'] + '|d=%g' % r['delta']
        for fam, s in r['fams'].items():
            for key in ('pairs', 'pairs_sq', 'L4_strict', 'inj_overlap', 'L5_tested', 'L5_sqtests', 'L6_tested',
                        'L6_ray_tests', 'L6_holes', 'chain_entries', 'pw_done', 'pw_M'):
                agg[key] += s[key]
            for key, val in s.items():
                if key.startswith('cat_'):
                    agg[key] += val
            agg['M_meas'] += s['meas'].get('M', 0.0)
            agg['T_meas'] += s['meas'].get('T', 0.0)
            if s.get('merge_at_T'):
                agg['merge_at_T_meas'] += s['merge_at_T']
            if s.get('merge_meas'):
                agg['merge_meas'] += s['merge_meas']
            if s['meas'].get('M', 0) > 0:
                agg['fams_with_M'] += 1
            upd_min('L4_slack_min', s['L4_slack_min'], tag)
            upd_min('L5_gap_min(/delta)', s['L5_gap_min'], tag)
            upd_min('L6_sharp_min(/delta*bbar)', s['L6_sharp_min'], tag)
            upd_min('L6_weak_min(/delta*bbar)', s.get('L6_weak_min'), tag)
            upd_max('L6_exc_ratio_max', s['L6_exc_ratio_max'], tag)
            upd_max('L6_gv*den/delta_max', s['L6_gv_den_max'], tag)
            upd_max('L6_gv/delta_max', s['L6_gv_max'], tag)
            upd_min('chain_min', s['chain_min'], tag)
print(json.dumps(dict(agg), indent=1))
print('kinds', dict(kinds))
print('params', {str(k): v for k, v in by_delta.items()})
for k, v in mins.items():
    print('MIN', k, v)
for k, v in maxs.items():
    print('MAX', k, v)
print('violations:', len(viol))
for v in viol[:30]:
    print(v)
