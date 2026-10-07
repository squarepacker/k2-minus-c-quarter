"""(c) near-extremal death configurations (Theorem 7.1(D)): almost all strip waste swept by dying paths (exact).

Usage:  python deaths.py        (writes out/death_results.json; a few seconds)"""
from fractions import Fraction as Fr
import json
from tracer import Sq, trace, check_config, out_path
from analysis import analyze
import configs as C

out = []
for delta in (Fr(1, 10**5), Fr(1, 100), Fr(1, 10)):
    for extra in (Fr(0), delta / 1000):
        # k=4, h=1: row of 4 squares with bottoms at delta+extra (tie at extra=0 -> D wins); slits 1e-9
        g = Fr(1, 10**9)
        sqs = [Sq(i, Fr(1, 2) + i * (1 + g), Fr(1, 2) + delta + extra, 0) for i in range(4)]
        sqs = [S for S in sqs if check_config([S], 4) is None]
        err = check_config(sqs, 4)
        rec = trace(sqs, 4, delta, 1, Fr(10)); res = analyze(rec, sqs, 4, delta, 1, Fr(10))
        out.append(dict(cfg='row_above_floor', delta=str(delta), extra=str(extra), err=err, deltaD=res['deltaD'], waste=res['waste_strip'],
                        ratio=res['death_ratio_vs_waste_plus_rects'], intMg_ratio=res['intMg_over_waste_plus_rects'], ok=res['death_ok'] and res['intMg_ok']))
    # k=6, h=2: row 1 of alternating tiny tilts on the floor (valleys at the top), row 2 aligned with bottoms ~delta above row-1 tops
    for th in (1e-6, 1e-3, 0.05):
        sq = []
        x = Fr(3, 2)
        for i in range(4):
            s = C.s_of(th / 2 * (1 if i % 2 else -1), 10**15)
            if not sq:
                S = C.sq_lowest_at(0, s, (x, Fr(0)))
            else:
                P = sq[-1]
                if P.s < 0 < s:
                    w = C.vtx(P, 'TR'); tmp = Sq(0, 0, 0, s)
                    sth = P.u[0] * tmp.u[1] - P.u[1] * tmp.u[0]; cth = P.u[0] * tmp.u[0] + P.u[1] * tmp.u[1]
                    S = C.sq_vertex_at(0, s, 'TL', (w[0] + delta * sth / cth / 3, w[1]))
                else:
                    b = C.vtx(P, 'BR'); S = C.sq_vertex_at(0, s, 'BL', (b[0] + delta / 3, b[1]))
                    if min(v[1] for v in S.V) < 0: S = Sq(0, S.c[0], S.c[1] - min(v[1] for v in S.V), S.s)
            sq.append(Sq(len(sq), S.c[0], S.c[1], S.s))
        top = max(v[1] for S in sq for v in S.V)
        for j in range(5):
            sq.append(Sq(len(sq), Fr(3, 4) + j * Fr(1001, 1000), top + delta * Fr(9, 10) + Fr(1, 2), 0))
        err = check_config(sq, 6)
        if err: out.append(dict(cfg=f'two_rows_th{th}', delta=str(delta), err=err)); continue
        rec = trace(sq, 6, delta, 2, Fr(10)); res = analyze(rec, sq, 6, delta, 2, Fr(10))
        out.append(dict(cfg=f'two_rows_th{th}', delta=str(delta), deltaD=res['deltaD'], waste=res['waste_strip'], pairs=res['n_pairs'],
                        dbl=res['n_double_regions'], ratio=res['death_ratio_vs_waste_plus_rects'], intMg_ratio=res['intMg_over_waste_plus_rects'],
                        ok=res['death_ok'] and res['intMg_ok'], M=res['M'], merge_ratio=res['merge_ratio_max'], rect_viol=res['n_rect_violations']))
for o in out: print(o)
json.dump(out, open(out_path('death_results.json'), 'w'), indent=1)
