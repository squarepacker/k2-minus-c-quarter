"""Aggregate the run records res_*.jsonl of run_family.py per family (as written in the field 'family').

Usage:  python summarize.py [folder]      (default folder: out; writes out/summary.json)"""
import json, glob, sys, collections, os
from tracer import out_path
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
fams = collections.OrderedDict()
for fn in sorted(glob.glob(os.path.join(SRC, 'res_*.jsonl'))):
    for line in open(fn):
        try: r = json.loads(line)
        except Exception: continue
        fams.setdefault(r['family'], []).append(r)
keys_max = ['rect_qa_over_we_max', 'rect_tR_over_delta_max', 'star_ratio_max', 'dist_ratio_max', 'De_over_delta_we_max',
            'merge_ratio_max', 'merge_global_ratio', 'max_Y_per_pair', 'E_ratio_max', 'Ov_ratio_max',
            'death_ratio_vs_waste_plus_rects', 'intMg_over_waste_plus_rects', 'D_over_bound_global', 'Mg_max']
keys_cnt = ['n_rect_violations', 'same_source_overlaps', 'floor_square_overlaps', 'triple_regions', 'floor_merges',
            'Ov_support_violations', 'n_anom']
out = {}
for fam, rs in fams.items():
    ok = [r for r in rs if 'Mg_max' in r]
    s = dict(runs=len(rs), analysed=len(ok), skipped=sum('skipped' in r for r in rs), errors=sum('error' in r for r in rs),
             timeouts=sum('timeout' in r for r in rs) + sum(bool(r.get('timed_out')) for r in ok),
             with_doubles=sum(r['n_double_regions'] > 0 for r in ok), with_merges=sum(r['M'] > 0 for r in ok),
             with_E=sum(r['E'] > 0 for r in ok), with_D=sum(r['D'] > 0 for r in ok),
             conservation_fail=sum(not r['conservation_ok'] for r in ok), death_fail=sum(not r['death_ok'] for r in ok),
             intMg_fail=sum(not r.get('intMg_ok', True) for r in ok), jac_fail=sum(not r['jacobian_ok'] for r in ok),
             minlam_min=min([r['minlam'] for r in ok if r['minlam'] is not None], default=None),
             deltas=sorted(set(float(eval(r['delta'].replace('/', '/'))) if '/' not in r['delta'] else float(int(r['delta'].split('/')[0]) / int(r['delta'].split('/')[1])) for r in rs)))
    for k in keys_max:
        vals = [(r.get(k), r['name'], r['delta']) for r in ok if isinstance(r.get(k), (int, float))]
        if vals:
            m = max(vals, key=lambda v: v[0]); s[k] = [m[0], m[1], m[2]]
    for k in keys_cnt:
        s[k] = sum(r.get(k, 0) or 0 for r in ok)
    for k in ('merge_violations', 'E_violations', 'Ov_violations', 'dist_violations', 'rect_violations'):
        ex = [(r['name'], r['delta'], r[k][:2]) for r in ok if r.get(k)]
        s[k] = len(ex); s[k + '_ex'] = ex[:3]
    anoms = [(r['name'], r['delta'], r['anom'][:3]) for r in ok if r.get('n_anom')]
    s['anom_ex'] = anoms[:3]
    errs = [(r['name'], r['delta'], r['error'][-300:]) for r in rs if 'error' in r]
    s['err_ex'] = errs[:2]
    trip = [(r['name'], r['delta'], r.get('triple_example')) for r in ok if r.get('triple_regions')]
    s['triple_ex'] = trip[:3]
    out[fam] = s
json.dump(out, open(out_path('summary.json'), 'w'), indent=1, default=str)
for fam, s in out.items():
    print('==', fam, {k: v for k, v in s.items() if not k.endswith('_ex')})
    for k, v in s.items():
        if k.endswith('_ex') and v: print('   ', k, str(v)[:600])
