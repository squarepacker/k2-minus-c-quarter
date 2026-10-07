# -*- coding: utf-8 -*-
"""list the (config, variant, side) with r_shadow > 1 in campaign files, grouped by campaign.
usage: python asmx_flagged.py out/campA.jsonl [out/campC.jsonl ...]  > flagged.json"""
import json, os, sys
res = {}
for fn in sys.argv[1:]:
    camp = os.path.basename(fn).split('.')[0]
    with open(fn, encoding='utf-8') as fh:
        for line in fh:
            r = json.loads(line)
            if 'variants' not in r:
                continue
            i = int(r['cid'].split('-')[-1])
            for vi, v in enumerate(r['variants']):
                for side in ('floor', 'ceil'):
                    x = v[side].get('r_shadow')
                    if x is not None and x > 1 + 1e-9:
                        res.setdefault(camp, {}).setdefault(str(i), []).append([vi, side])
print(json.dumps(res))
