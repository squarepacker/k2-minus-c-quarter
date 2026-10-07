"""Mutation tests of the checker: each mutation must be detected (otherwise the checker is too weak).

Usage:  python mutations.py        (writes out/mutation_results.json; under a minute)"""
from fractions import Fraction as Fr
import random, json
import tracer, analysis
import configs as C
from tracer import trace, check_config, out_path
from analysis import analyze

BIG = Fr(10)
rnd = random.Random(99)
cfgs = []
while len(cfgs) < 25:
    theta = 10 ** rnd.uniform(-7, -0.5); delta = rnd.choice([Fr(1, 10**5), Fr(1, 1000), Fr(1, 20)])
    sq = C.merge_cfg(theta, delta, rnd.choice([0.0, theta / 2, 0.3]), Fr(rnd.randint(-450, 450), 1000),
                     a_frac=Fr(rnd.randint(1, 999), 1000), b_frac=Fr(rnd.randint(1, 999), 1000), gap_frac=Fr(1, 10**6), eps=delta / 1000)
    if check_config(sq, 6) is None: cfgs.append((sq, delta))

def flags(res):
    return dict(rect=res['n_rect_violations'] > 0, merge=len(res['merge_violations']) > 0, E=len(res['E_violations']) > 0,
                same_src=res['same_source_overlaps'] > 0, conservation=not res['conservation_ok'], Ov=len(res['Ov_violations']) > 0)

out = {}
# baseline
base = [flags(analyze(trace(sq, 6, d, 2, BIG), sq, 6, d, 2, BIG)) for sq, d in cfgs]
out['baseline_any_flag'] = sum(any(f.values()) for f in base)
# MUT-A: rectangle width halved
orig = analysis.pair_geom
def half(sqs, a, b, delta):
    G = orig(sqs, a, b, delta); G['we'] = G['we'] / 2; return G
analysis.pair_geom = half
fa = [flags(analyze(trace(sq, 6, d, 2, BIG), sq, 6, d, 2, BIG)) for sq, d in cfgs]
analysis.pair_geom = orig
out['MUT-A_we_halved_detected'] = sum(f['rect'] or f['merge'] for f in fa)
# MUT-B: R1 disabled
tracer.MUT['no_r1'] = True
fb = [flags(analyze(trace(sq, 8, d, 3, BIG), sq, 8, d, 3, BIG)) for sq, d in cfgs]
tracer.MUT['no_r1'] = False
out['MUT-B_noR1_detected'] = sum(f['same_src'] or f['conservation'] for f in fb)
# MUT-C: contact before death (tie config D1)
A = C.Sq(0, Fr(3, 2), Fr(1, 2), 0); d = Fr(1, 10**5); B = C.Sq(1, Fr(17, 10), Fr(3, 2) + d, 0)
r0 = trace([A, B], 6, d, 2, BIG)
tracer.MUT['contact_before_death'] = True
r1 = trace([A, B], 6, d, 2, BIG)
tracer.MUT['contact_before_death'] = False
out['MUT-C_tie_H_normal'] = str(r0.H); out['MUT-C_tie_H_mutated'] = str(r1.H); out['MUT-C_detected'] = (r0.H != r1.H)
# MUT-D: bounds evaluated with delta/2
fd = [flags(analyze(trace(sq, 6, d, 2, BIG), sq, 6, d / 2, 2, BIG)) for sq, d in cfgs]
out['MUT-D_half_delta_detected'] = sum(f['E'] or f['rect'] for f in fd)
out['n_cfgs'] = len(cfgs)
out['n_cfgs_with_merge'] = sum(1 for sq, d in cfgs if analyze(trace(sq, 6, d, 2, BIG), sq, 6, d, 2, BIG)['M'] > 0)
print(json.dumps(out, indent=1))
json.dump(out, open(out_path('mutation_results.json'), 'w'), indent=1)
