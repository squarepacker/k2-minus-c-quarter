"""Targeted degenerate (measure-zero / tie) tests, exact arithmetic.

Usage:  python degenerate.py        (writes out/degen_results.json; a few seconds)

Note on the final count "BAD": rec.entries counts the beam pieces that enter a square, not the squares.
In D1 (gap exactly delta) the beam entering A is split at x = 1.2 (the left edge of B), so rec.entries = 2
although B is never entered (D = k, H = 0); the expectation 'no_entry_into_B' (entries == 1) is therefore
reported False for both values of delta.  See the README of this folder."""
from fractions import Fraction as Fr
import json, math
from tracer import Sq, trace, check_config, out_path
from analysis import analyze
import configs as C

BIG = Fr(10)
out = []

def go(name, sqs, delta, k=6, sF=BIG, expect=None):
    sqs = C.reindex(sqs)
    err = check_config(sqs, k)
    if err:
        out.append(dict(name=name, skipped=err)); print(name, 'SKIP', err); return None, None
    h = Fr(k, 2) - 1
    rec = trace(sqs, k, delta, h, sF)
    res = analyze(rec, sqs, k, delta, h, sF)
    r = dict(name=name, delta=str(delta), conservation=res['conservation_ok'], D=str(rec.D), W=str(rec.W), H=str(rec.H), T=str(rec.T),
             E={f'{a}->{b}': str(v) for (a, b), v in rec.E.items()}, M=[(y, w, l, str(m)) for (y, w, l, m, _) in rec.M],
             entries=rec.entries, anom=rec.anom[:5], Mg=res['Mg_max'], rect_viol=res['n_rect_violations'], star=res['star_ratio_max'],
             E_ratio=res['E_ratio_max'], E_viol=res['E_violations'], merge_ratio=res['merge_ratio_max'], merge_viol=res['merge_violations'],
             Ov=res['Ov_ratio_max'], death_ok=res['death_ok'], intMg_ok=res['intMg_ok'], dist=res['dist_ratio_max'],
             dist_viol=res['dist_violations'])
    if expect:
        r['expect'] = {k_: bool(f(rec, res)) for k_, f in expect.items()}
    out.append(r); print(json.dumps(r)[:900]); return rec, res

d = Fr(1, 10**5)
for delta in (Fr(1, 10**5), Fr(1, 100)):
    # D1: horizontal gap exactly delta between A top and B bottom: t_S == t_D on a positive-measure set -> D wins
    A = Sq(0, Fr(3, 2), Fr(1, 2), 0); B = Sq(1, Fr(17, 10), Fr(3, 2) + delta, 0)
    go(f'D1_tie_gap_eq_delta', [A, B], delta, expect={'no_entry_into_B': lambda rec, res: rec.entries == 1})
    B2 = Sq(1, Fr(17, 10), Fr(3, 2) + delta - delta / 1000, 0)
    go(f'D1b_gap_slightly_less', [A, B2], delta, expect={'entry_into_B': lambda rec, res: rec.entries >= 2})
    # D2/D5: tilted square with lowest vertex exactly on the floor: E from floor = delta*tan(phi) exactly (bound attained)
    for phi in (1e-6, 0.1, 0.7):
        S = C.sq_lowest_at(0, C.s_of(phi), (Fr(2), Fr(0)))
        go(f'D5_corner_on_floor_phi{phi}', [S], delta,
           expect={'E_ratio_le_1': lambda rec, res: res['E_ratio_max'] <= 1.0})
    # D3: vertex alignment: B left edge exactly above A left edge, gap delta/2
    A = Sq(0, Fr(3, 2), Fr(1, 2), 0); B = Sq(1, Fr(3, 2), Fr(3, 2) + delta / 2, 0)
    go('D3_aligned_stack', [A, B], delta)
    B = Sq(1, Fr(5, 2), Fr(3, 2) + delta / 2, 0)   # B's left vertex exactly above A's right vertex
    go('D3b_corner_over_corner', [A, B], delta)
    # D4: valley + Y whose bottom passes (almost) through v
    for tau0 in (Fr(-1, 2) + Fr(1, 10**6), Fr(0)):
        sq = C.merge_cfg(0.01, delta, 0.005, tau0, a_frac=Fr(0), b_frac=Fr(1, 10**6), gap_frac=Fr(1, 10**6), eps=delta / 10**4)
        go(f'D4_Y_through_v_tau{float(tau0):.3g}', sq, delta)
    # D6: valley with corner gap ~ delta*tan(theta) (boundary of dist claim)
    for gf in (Fr(999, 1000), Fr(1), Fr(1001, 1000)):
        sq = C.valley(0.02, delta, gap_frac=gf, eps=delta / 10**4)
        go(f'D6_valley_gap{float(gf)}', sq, delta)
    # D7: symmetric valley + axis-aligned Y centred on the crack (g ties at the symmetric point)
    sq = C.valley(float(delta) / 100, delta, gap_frac=Fr(1, 10**6), eps=delta / 10**4)
    L, R = sq[0], sq[1]
    v = C.vtx(R, 'TL'); w = C.vtx(L, 'TR')
    mid = ((v[0] + w[0]) / 2, max(v[1], w[1]) + delta / 3)
    Y = C.place_on_bottom(2, Fr(0), mid, Fr(0))
    go('D7_symmetric_merge', [L, R, Y], delta)
    # D8: square touching the wall x=0, tilted so its rays go to the wall
    s3 = C.s_of(0.3); tmp = Sq(0, 0, 0, s3)
    S = C.sq_vertex_at(0, s3, 'TL', (Fr(0), tmp.cs))   # TL on the wall, BL on the floor; rays lean left
    go('D8_wall', [S], delta)
    # D9: square with bottom exactly at height h (=2 for k=6): contact at h processed before H
    A = C.sq_lowest_at(0, C.s_of(4e-3), (Fr(3, 2), Fr(0)))   # k=4, h=1: A's top-left part is below h=1
    tl = C.vtx(A, 'TL')
    B = Sq(1, tl[0] + Fr(1, 10**4) - Fr(1, 2), Fr(3, 2), 0)    # B bottom exactly at y = 1 = h, right edge at TL.x + 1e-4
    go('D9_contact_at_h', [A, B], delta, k=4)
    # D10: tie between wall and contact / wall and death
    S = Sq(0, Fr(1, 2), Fr(1, 2), C.s_of(0.0))
    go('D10_square_at_wall', [S], delta)

json.dump(out, open(out_path('degen_results.json'), 'w'), indent=1, default=str)
bad = [r for r in out if r.get('expect') and not all(r['expect'].values())]
bad += [r for r in out if 'conservation' in r and (not r['conservation'] or r['anom'] or r['rect_viol'] or r['E_viol'] or r['merge_viol'] or not r['death_ok'] or not r['intMg_ok'] or r['dist_viol'])]
print('N', len(out), 'BAD', len(bad))
for r in bad: print('BAD', r['name'], r.get('expect'), r.get('anom'), r.get('E_viol'))
