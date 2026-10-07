# -*- coding: utf-8 -*-
"""asmx_replaycheck.py -- check that replaying the RNG of a campaign reproduces (k, N) of its first records.
usage: python asmx_replaycheck.py <campaign> <seed> <families> <k_values> <i_max>
reads <campaign>.jsonl from out/ (or from the current directory); campA and campB are replayed with gen_hill_v1."""
import os, sys, json, random
from asmx_replay import replay_prefix, use_hill_v1
camp = sys.argv[1]; seed = int(sys.argv[2]); fams = sys.argv[3].split(','); ks = [int(x) for x in sys.argv[4].split(',')]
if camp in ('campA', 'campB'):
    use_hill_v1()
imax = int(sys.argv[5])
rec = {}
fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', camp + '.jsonl')
if not os.path.exists(fn):
    fn = camp + '.jsonl'
with open(fn, encoding='utf-8') as fh:
    for line in fh:
        r = json.loads(line)
        rec[int(r['cid'].split('-')[-1])] = r
rng = random.Random(seed)
bad = 0
for i in range(imax + 1):
    fam = fams[i % len(fams)]
    k = rng.choice(ks)
    if fam == 'random' and k > 20:
        k = 20
    res = replay_prefix(fam, k, rng)
    r = rec.get(i)
    Nr = res[0].N if res else None
    if r is None or r.get('N') != Nr or r.get('k') != k:
        bad += 1
        if bad <= 5:
            print('DIVERGE at', i, fam, k, Nr, r.get('k') if r else None, r.get('N') if r else None)
print('checked', imax + 1, 'diverged', bad)
