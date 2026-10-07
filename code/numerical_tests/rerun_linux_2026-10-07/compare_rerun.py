"""Compare the re-run of 7 October 2026 (console_logs.zip, outputs.zip in this folder) with the recorded results and with
the claims of the folder READMEs of code/numerical_tests. Records are compared after dropping run-time and memory fields
and mapping the keys that were renamed after the original runs (each folder README lists the renamings).

Usage:  python compare_rerun.py [--campD FILE] [--campC FILE] [--out FILE]
  --campD, --campC: the records campD.jsonl of the P7 test and campC.jsonl of the assembly test, which are not
  included in this repository (44 MB and 27 MB); without them those two rows are reported as not checked.
Writes out/compare_out.md in this folder, or FILE (a table with one row per claim). The zips are unpacked into a temporary
folder, which is removed at the end. Needs only the Python standard library."""
import os, re, sys, json, zipfile, shutil, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
NT = os.path.dirname(HERE)
X = tempfile.mkdtemp(prefix='rerun_')
for zname in ('console_logs.zip', 'outputs.zip'):
    with zipfile.ZipFile(os.path.join(HERE, zname)) as z:
        z.extractall(X)
RES = os.path.join(X, 'results'); OUT = os.path.join(X, 'outputs')
ARGS = dict(zip(sys.argv[1::2], sys.argv[2::2]))
TIMING = re.compile(r'elapsed|(^|_)time($|_)|_sec$|^sec$|rss|(^|_)mem($|_)|wall_s$')
DROPPED = set()
rows = []


def row(folder, claim, ok, detail):
    rows.append((folder, claim, 'yes' if ok is True else ('no' if ok is False else ok), detail))


def out(lane, folder, name):
    return os.path.join(OUT, lane, 'code', 'numerical_tests', folder, 'out', name)


def rec(folder, name, sub='original'):
    return os.path.join(NT, folder, 'results', sub, name)


def log(lane, label):
    return open(os.path.join(RES, f'{lane}_{label}.log'), encoding='utf-8', errors='replace').read()


def jl(text):
    return [json.loads(l) for l in text.splitlines() if l.strip()]


def jload(p):
    t = open(p, encoding='utf-8').read()
    return jl(t) if p.endswith('.jsonl') else json.loads(t)


def norm(o, ren=(), drop=None, vren=()):
    """drop timing keys (and keys matching drop), rename keys (ren) and substrings of string values (vren)"""
    if isinstance(o, dict):
        d = {}
        for k, v in o.items():
            if TIMING.search(str(k)) or (drop and re.search(drop, str(k))):
                DROPPED.add(str(k)); continue
            kk = str(k)
            for a, b in ren:
                kk = kk.replace(a, b)
            d[kk] = norm(v, ren, drop, vren)
        return d
    if isinstance(o, list):
        return [norm(v, ren, drop, vren) for v in o]
    if isinstance(o, str):
        for a, b in vren:
            o = o.replace(a, b)
    return o


def cmp(a, b, path=''):
    """(number of differing leaves, max relative float deviation among float leaves, first difference)"""
    if isinstance(a, dict) and isinstance(b, dict):
        n, m, f = 0, 0.0, None
        for k in set(a) | set(b):
            if k not in a or k not in b:
                n += 1; f = f or f'{path}/{k} only in one'; continue
            n1, m1, f1 = cmp(a[k], b[k], f'{path}/{k}')
            n, m, f = n + n1, max(m, m1), f or f1
        return n, m, f
    if isinstance(a, list) and isinstance(b, list):
        n, m, f = (0, 0.0, None) if len(a) == len(b) else (1, 0.0, f'{path} length {len(a)} vs {len(b)}')
        for i, (x, y) in enumerate(zip(a, b)):
            n1, m1, f1 = cmp(x, y, f'{path}[{i}]')
            n, m, f = n + n1, max(m, m1), f or f1
        return n, m, f
    if isinstance(a, float) and isinstance(b, float) and a != b:
        return 1, abs(a - b) / max(abs(a), abs(b)), f'{path}: {a!r} vs {b!r}'
    return (0, 0.0, None) if a == b else (1, 0.0, f'{path}: {str(a)[:60]} vs {str(b)[:60]}')


def fdev(a, b, tol, atol):
    """largest absolute deviation among float differences that are not within relative tolerance tol"""
    if isinstance(a, dict):
        return max([fdev(a[k], b[k], tol, atol) for k in a] or [0.0])
    if isinstance(a, list):
        return max([fdev(x, y, tol, atol) for x, y in zip(a, b)] or [0.0])
    if isinstance(a, float) and a != b and abs(a - b) > tol * max(abs(a), abs(b)):
        return abs(a - b)
    return 0.0


def same(new, old, ren=(), tol=0.0, drop=None, vren=(), atol=0.0):
    """(ok, detail): ok if equal after norm, or if every difference is a float difference of relative size <= tol
    or of absolute size <= atol"""
    A, B = norm(new, (), drop), norm(old, ren, drop, vren)
    n, m, f = cmp(A, B)
    if n == 0:
        return True, 'identical (timing fields dropped)'
    if atol > 0 and all_float_diffs(A, B) and fdev(A, B, tol, atol) <= atol:
        return True, (f'identical except {n} floating-point value(s): relative deviation at most {tol:.0e} or absolute '
                      f'deviation at most {fdev(A, B, tol, atol):.1e} (values near 0); first: {f}')
    if tol > 0 and m <= tol and all_float_diffs(A, B):
        return True, f'identical except {n} floating-point value(s), relative deviation at most {m:.1e}; first: {f}'
    return False, f'{n} differing values, max rel. float deviation {m:.1e}; first: {f}'


def all_float_diffs(new, old):
    """True if every difference (of already normalised objects) is a float difference"""
    def walk(a, b):
        if isinstance(a, dict) and isinstance(b, dict):
            return set(a) == set(b) and all(walk(a[k], b[k]) for k in a)
        if isinstance(a, list) and isinstance(b, list):
            return len(a) == len(b) and all(walk(x, y) for x, y in zip(a, b))
        return a == b or (isinstance(a, float) and isinstance(b, float))
    return walk(new, old)


def eol(p):
    return open(p, 'rb').read().replace(b'\r\n', b'\n')


# ---------------------------------------------------------------- integrity and exit codes
sha = open(os.path.join(RES, 'sha256_check.txt'), encoding='utf-8').read().splitlines()
row('all', 'SHA256SUMS of the uploaded repository', sum(l.endswith(': OK') for l in sha) == 176 and 'exit 0' in sha[-1],
    f'{sum(l.endswith(": OK") for l in sha)} of 176 files OK')
c = log('A', 'constants').splitlines()
row('all', 'code/k14_constants_check.py', c[-1] == 'ALL OK', f'{sum(l.startswith("OK") for l in c)} checks, last line {c[-1]}')
idx = [l for lane in 'AB' for l in open(os.path.join(RES, f'{lane}_index.txt'), encoding='utf-8') if ' done ' in l]
bad = [l.strip() for l in idx if not l.strip().endswith('EXIT 0')]
row('all', 'every command exits with code 0', not bad and len(idx) == 73, f'{len(idx)} commands; non-zero: {bad or "none"}')
tb = [f for f in os.listdir(RES) if f.endswith('.log') and 'Traceback' in open(os.path.join(RES, f), encoding='utf-8', errors='replace').read()]
row('all', 'no Python traceback in any log', not tb, ', '.join(tb) or 'none')

# ---------------------------------------------------------------- P1_P2
REN12 = (('prover', 'alphaT'), ('defs', '050001'))
for name in ('thresholds_out.json', 'reflection_out.json', 'merge_points_out.json', 'witnesses_out.json', 'side_contacts_out.json'):
    ok, d = same(jload(out('A', 'P1_P2', name)), jload(rec('P1_P2', name)), REN12)
    row('P1_P2', f'{name} equals the record' + (' (full run, not re-run before)' if 'side' in name else ''), ok, d)
row('P1_P2', 'aggregate_out.json equals the record', eol(out('A', 'P1_P2', 'aggregate_out.json')) == eol(rec('P1_P2', 'aggregate_out.json')),
    'byte-identical up to line endings' if eol(out('A', 'P1_P2', 'aggregate_out.json')) == eol(rec('P1_P2', 'aggregate_out.json')) else 'differs')
for label, summ in (('p12_smoke', 'smoke_summary.json'), ('p12_tvp', 'quick_tvp_summary.json'), ('p12_fuzz', 'fuzz_summary.json')):
    s = jload(out('A', 'P1_P2', summ))
    row('P1_P2', f'{label}: 0 failures, 0 anomalies', s['n_fail'] == 0 and s['n_anom'] == 0,
        f"{s['n_cfg']} configurations, {s['n_contacts']} contacts (time budget; the README's laptop counts are smaller)")
t = log('A', 'p12_conseq').splitlines()[-1]
row('P1_P2', 'consequence.py quick: 0 violations, 0 found by the backward search, sanity 87/87', t.startswith('done 6 ') and " 0 14688 0 " in t
    and "'checked': 87, 'matched': 87" in t, t)

# ---------------------------------------------------------------- P3_P5
REN35 = (('lemma43_fail', 'two_sources_fail'), ('lemma43', 'two_sources'), ('lemma44', 'three_sources'))
with zipfile.ZipFile(rec('P3_P5', 'runs.zip')) as z:
    chunks = {}
    for n in z.namelist():
        m = re.fullmatch(r'res_([A-Za-z0-9]+)_(\d+)\.jsonl', n)
        if m:
            chunks.setdefault(m.group(1), []).append((int(m.group(2)), jl(z.read(n).decode('utf-8'))))
VIOLK = ('rect_violations', 'n_rect_violations', 'dist_violations', 'merge_violations', 'E_violations', 'Ov_violations', 'Ov_support_violations')


def nviol(rs):
    return sum((len(r[k]) if isinstance(r.get(k), list) else int(r.get(k) or 0)) for r in rs for k in VIOLK if k in r)


def fam_check(fam, new, old, label):
    n = min(len(new), len(old))
    eq = [same(a, b, tol=1e-10) for a, b in zip(new[:n], old[:n])]
    dev = max([float(e[1].split('at most ')[1].split(';')[0]) for e in eq if e[0] and 'at most' in e[1]] or [0.0])
    sq = [a.get('squares') == b.get('squares') for a, b in zip(new[:n], old[:n])]
    nd = sum(1 for s in sq if not s)
    ok = all(e[0] for e, s in zip(eq, sq) if s) and nviol(new) == 0
    exact = sum(1 for e in eq if e[0] and e[1].startswith('identical (timing'))
    fl = sum(1 for e in eq if e[0] and not e[1].startswith('identical (timing'))
    row('P3_P5', label, ok, f'{len(new)} records on the server ({len(old)} recorded, compared {n}): {exact} identical, {fl} identical except '
        f'floating-point values (largest relative deviation {dev:.1e}), {nd} with other exact square coordinates (generated '
        f'from floating-point values); violations on the server {nviol(new)}')


for fam in ('valley', 'smoke', 'bigvalley', 'mechA', 'three', 'fan', 'strip', 'merge', 'near45', 'mrand', 'rowjam', 'regime', 'tri'):
    fam_check(fam, jload(out('A', 'P3_P5', f'res_{fam}_0.jsonl')), [r for _, rs in sorted(chunks[fam]) for r in rs],
              f'run_family.py {fam}: records equal the recorded ones')
fam_check('mrand9000', jload(out('A', 'P3_P5', 'res_mrand_9000.jsonl')), dict(chunks['mrand'])[9000], 'run_family.py mrand 20 9000: records equal res_mrand_9000.jsonl')
row('P3_P5', 'summarize.py: out/summary.json equals the record', eol(out('A', 'P3_P5', 'summary.json')) == eol(rec('P3_P5', 'summary.json')),
    'byte-identical up to line endings' if eol(out('A', 'P3_P5', 'summary.json')) == eol(rec('P3_P5', 'summary.json')) else same(jload(out('A', 'P3_P5', 'summary.json')), jload(rec('P3_P5', 'summary.json')))[1])
for name in ('degen_results.json', 'death_results.json', 'mutation_results.json', 'cross_check_results.json', 'boundary_exact.json',
             'geom_starmin_0.3.json', 'geom_star_0.3.json', 'geom_starmin_0.9.json'):
    ok, d = same(jload(out('A', 'P3_P5', name)), jload(rec('P3_P5', name)), REN35)
    row('P3_P5', f'{name} equals the record', ok, d)
for mode, vals in (('starmin', ('0.01', '1.4')), ('star', ('1e-6', '1e-3', '1.4'))):
    new = jload(out('A', 'P3_P5', f'geom_{mode}_' + '_'.join(vals) + '.json'))
    for i, v in enumerate(vals):
        old = jload(rec('P3_P5', f'geom_{mode}_{v}.json'))
        ok, d = same({k: [e for e in new[k] if abs(e['theta'] - float(v)) < 1e-15] for k in new}, old)
        row('P3_P5', f'geom_search.py {mode} theta={v} (not re-run before) equals geom_{mode}_{v}.json', ok, d)
for name, key in (('geom_two_all.json', 'two_sources'), ('geom_three_all.json', 'three_sources')):
    new, old = norm(jload(out('B', 'P3_P5', name))), norm(jload(rec('P3_P5', name)), REN35)
    m = len(old[key])
    ok, d = same({key: new[key][:m]}, old)
    row('P3_P5', f'geom_search.py {key.split("_")[0]} (not re-run before) equals {name}', ok,
        f'{len(new[key])} entries on the server, {m} recorded; the recorded ones: {d}')

# ---------------------------------------------------------------- P4
for name in ('res_negctl.json', 'res_ties.jsonl', 'res_wallsweep.json'):
    ok, d = same(jload(out('A', 'P4', name)), jload(rec('P4', name)))
    row('P4', f'{name} equals the record', ok, d)
for q, r_, n in (('quick_dag.jsonl', 'res_dag.jsonl', 6), ('quick_manyvar.jsonl', 'res_manyvar.jsonl', 6)):
    new = jload(out('A', 'P4', q)); old = jload(rec('P4', r_))[:n]
    ok, d = same(new, old, vren=(('_defs', '_paper'),))
    row('P4', f'{q} equals the first {n} records of {r_} (_defs read as _paper)', ok and len(new) == n, d)
for q in ('quick_scaled.jsonl', 'quick_paper.jsonl'):
    new = jload(out('A', 'P4', q))
    v = sum(len(r.get('violations', [])) if isinstance(r.get('violations'), list) else int(r.get('violations', 0) or 0) for r in new)
    row('P4', f'{q}: 0 violations', v == 0, f'{len(new)} configurations')
for lbl, f in (('p4_summ_all', 'summary_all.txt'), ('p4_summ_r2', 'summary_round2.txt')):
    a = log('A', lbl).replace('\r\n', '\n').strip(); b = open(rec('P4', f), encoding='utf-8-sig').read().replace('\r\n', '\n').strip()
    row('P4', f'summarize.py output equals {f}', a == b, 'identical (the record starts with a byte order mark)' if a == b else 'differs')
e, e0 = jload(out('A', 'P4', 'res_exclab.json')), jload(rec('P4', 'res_exclab.json'))
pat = lambda x: (x['table']['D']['n_fail'], x['table']['side']['n_fail'] > 0, x['table']['D']['min_margin'], x['table']['side']['min_margin'], x['examples'])
row('P4', 'exceptional_types.py 420 (full budget): the same pattern as the recorded run', pat(e) == pat(e0),
    f"{e['configs']} configurations (recorded {e0['configs']}): type D {e['table']['D']['n']} heights, 0 failures; type 'side' "
    f"{e['table']['side']['n']} heights, {e['table']['side']['n_fail']} failures, least margin {e['table']['side']['min_margin']}")
u = log('A', 'p4_unit')
row('P4', 'run_tests.py unit: 11 checks, all PASS', u.count('PASS') >= 11 and 'FAIL' not in u, f"{u.count('PASS')} PASS")

# ---------------------------------------------------------------- P6
with zipfile.ZipFile(rec('P6', 'campaigns_raw.zip')) as z:
    camp3 = jl(z.read('camp_scaled_3.jsonl').decode('utf-8'))
new = jload(out('A', 'P6', 'quick_scaled_3.jsonl'))
ok, d = same(new, camp3[:len(new)], (('L6_defs', 'L6_weak'),), drop=r'^L7_tests$|^cat_')
row('P6', 'p6x_run.py 3 25: the records equal the first 25 of campaign 3 (except L7_tests and the added cat_* counters)', ok and len(new) == 25, d)
ok, d = same(jload(out('A', 'P6', 'targeted_out.json')), jload(rec('P6', 'targeted_out.json')), (('defs', 'weak'),), tol=1e-15)
row('P6', 'targeted_out.json equals the record', ok, d)
ok, d = same(jload(out('A', 'P6', 'negctl_out.json')), jload(rec('P6', 'negctl_out.json', 'rerun_2026-10-07')), (('defs', 'weak'),))
row('P6', 'negctl_out.json equals results/rerun_2026-10-07/negctl_out.json', ok, d)
a = log('A', 'p6_ties').replace('\r\n', '\n').strip(); b = open(rec('P6', 'ties_console.txt', 'rerun_2026-10-07'), encoding='utf-8').read().replace('\r\n', '\n').strip()
row('P6', 'p6x_ties.py console output equals ties_console.txt', a == b, 'identical' if a == b else 'differs')
row('P6', 'p6x_summary.py: violations 0', log('A', 'p6_summary').strip().endswith('violations: 0'), '')

# ---------------------------------------------------------------- P7
new = jload(out('A', 'P7', 'quick_D.jsonl'))
if '--campD' in ARGS:
    campD = jload(ARGS['--campD'])
    ok, d = same(new, campD[:len(new)], (('d_lemB_symdiff', 'd_Tchar_symdiff'),))
    row('P7', 'p7x_run.py blocked2 505: the records equal the first records of campaign D (campD.jsonl, not included)', ok, f'{len(new)} records; {d}')
else:
    row('P7', 'p7x_run.py blocked2 505: the records equal the first records of campaign D (campD.jsonl, not included)', 'not checked', 'campD.jsonl not given')
with zipfile.ZipFile(rec('P7', 'records_targeted_merge.zip')) as z:
    for name in ('t7.jsonl', 'targeted.jsonl', 'targeted2.jsonl'):
        ok, d = same(jload(out('A', 'P7', name)), jl(z.read(name).decode('utf-8')), (('d_lemB_symdiff', 'd_Tchar_symdiff'),))
        row('P7', f'{name} equals the record', ok, d)
a = log('A', 'p7_example').replace('\r\n', '\n').strip().splitlines(); b = open(rec('P7', 'example_out.txt'), encoding='utf-8').read().replace('\r\n', '\n').strip().splitlines()
row('P7', 'p7x_example.py output equals example_out.txt except the first line (label)', a[1:] == b[1:], f'{len(a)} lines')
s = log('A', 'p7_summary')
row('P7', 'p7x_summary.py: every viol_* counter is 0', not re.search(r'"viol_[a-z_0-9]*": [1-9]', s), '')

# ---------------------------------------------------------------- assembly
new = jload(out('A', 'assembly', 'quickC.jsonl'))
lab = 'asmx_run.py quickC: the 40 records equal the first 40 of campaign C (campC.jsonl, not included; run name quickC read as campC)'
if '--campC' in ARGS:
    campC = jload(ARGS['--campC'])
    ok, d = same(new, campC[:len(new)], drop=r'^t$|^t_flow$|^t_total$', vren=(('campC', 'quickC'),), tol=1e-9, atol=1e-12)
    row('assembly', lab, ok and len(new) == 40, d)
else:
    row('assembly', lab, 'not checked', 'campC.jsonl not given')
for name in ('targeted_T2.json', 'replay_campC.json'):
    ok, d = same(jload(out('A', 'assembly', name)), jload(rec('assembly', name)))
    row('assembly', f'{name} equals the record', ok, d)
ok, d = same(jload(out('A', 'assembly', 'targeted_T2b.json')), jload(rec('assembly', 'targeted_T2b.json')), tol=1e-11)
row('assembly', 'targeted_T2b.json equals the record', ok, d)
s = log('A', 'asm_summary')
viol = [l for l in s.splitlines() if l.startswith('VIOL_') and l.split()[1] != '0']
row('assembly', 'asmx_summary.py: every VIOL_* counter is 0', not viol, '; '.join(viol) or '')

# ---------------------------------------------------------------- P8 (lane B: full modes, not re-run before)
REN8 = (('Thm4.3 pi*', 'Thm A pi*'), ('Thm4.4 hand', 'Thm A13'), ('main pi* c*=0.1', 'Thm A pi* c*=0.1'), ('hand c*=0.1', 'Thm A13 c*=0.1'))
REN8S = (('Thm4.3', 'Thm A'), ('Thm4.4', 'Thm A13'))
for name in ('p8x_main_out.json', 'p8x_scan_out.json', 'p8x_asym_out.json'):
    ok, d = same(jload(out('B', 'P8', name)), jload(rec('P8', name)), REN8)
    row('P8', f'full mode: {name} equals the record (labels mapped)', ok, d)
NUM = re.compile(r'-?\d+\.\d{4,}(?:e[-+]?\d+)?')
for name in ('p8x_main_out.txt', 'p8x_scan_out.txt', 'p8x_asym_out.txt'):
    a = [l for l in open(out('B', 'P8', name), encoding='utf-8').read().splitlines() if not re.match(r'\s*time |.*elapsed', l)]
    b = open(rec('P8', name), encoding='utf-8').read().split('NOTE (added 2026-10-06')[0].splitlines()
    b = [l for l in b if not re.match(r'\s*time |.*elapsed', l)]
    na, nb = NUM.findall('\n'.join(a)), NUM.findall('\n'.join(b))
    row('P8', f'full mode: every number with 4 or more decimals in {name} equals the record', na == nb, f'{len(na)} numbers')
m = jload(out('B', 'P8', 'p8x_main_out.json'))
row('P8', 'p8x_main.py (full): 623 checks, 0 fails, 0 discrepancies', (m['checks'], m['fails'], len(m['discrepancies'])) == (623, 0, 0),
    f"{m['checks']} checks, {m['fails']} fails, {len(m['discrepancies'])} discrepancies, {len(m['unsafe_rounding'])} unsafe roundings (draft values, see the README)")
for name in ('p8x_asym_fix_out.json', 'p8x_sens_out.json'):
    ok, d = same(jload(out('A', 'P8', name)), jload(rec('P8', name)), REN8S if 'sens' in name else REN8)
    row('P8', f'{name} equals the record', ok, d)

# ---------------------------------------------------------------- write
nyes = sum(r[2] == 'yes' for r in rows)
txt = ['# Re-run of 7 October 2026: comparison with the README claims', '',
       f'{nyes} of {len(rows)} claims reproduced. Fields not compared: {", ".join(sorted(DROPPED)) or "none"}.', '',
       '| folder | claim | reproduced | detail |', '|---|---|---|---|']
txt += [f'| {f} | {c} | {o} | {d.replace("|", "/")} |' for f, c, o, d in rows]
OUTF = ARGS.get('--out', os.path.join(HERE, 'out', 'compare_out.md'))
os.makedirs(os.path.dirname(os.path.abspath(OUTF)), exist_ok=True)
open(OUTF, 'w', encoding='utf-8', newline='\n').write('\n'.join(txt) + '\n')
shutil.rmtree(X, ignore_errors=True)
print('\n'.join(txt))
