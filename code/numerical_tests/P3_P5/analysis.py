"""Analysis of an exact trace: the statements below checked exactly (Fractions) where possible.
  (a) at most two paths through a point, M_g <= 2 (Lemma 7.7); no point reached from the floor and a square
      (Lemma 7.4); no two arrivals from the same source (Lemma 5.3)
  (b) double points: q = v_1 + a u + b n with 0 <= a <= w_e = delta tan(theta)/cos(theta) and 0 <= b <= delta,
      dist(X_1, X_2) < delta tan(theta), and inequality (7.3) (Lemma 7.8)
  (c) deaths (Theorem 7.1(D), Section 7.3)       (d) merges (Theorem 7.1(M), Section 7.4)
  (e) gap overlap (Theorem 7.1(G), Section 7.5)  (f) end collisions per (source, target):
      |E_{X,Y}| <= delta tan|Delta| (Lemma 5.9(b), used for Proposition 5.10)"""
from fractions import Fraction as Fr
import math, time
import numpy as np
from tracer import dot, Sq

def cross(o, a, b): return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

def area(P):
    s = Fr(0)
    for i in range(len(P)):
        a = P[i]; b = P[(i+1) % len(P)]
        s += a[0]*b[1] - a[1]*b[0]
    return s / 2

def ccw(P):
    P = dedup(P)
    if len(P) >= 3 and area(P) < 0: P = P[::-1]
    return P

def dedup(P):
    out = []
    for p in P:
        if not out or out[-1] != p: out.append(p)
    if len(out) > 1 and out[0] == out[-1]: out.pop()
    return out

def clip(subj, clipp):
    """convex polygon intersection (both CCW), exact."""
    out = subj
    n = len(clipp)
    for i in range(n):
        A = clipp[i]; B = clipp[(i+1) % n]
        inp = out; out = []
        if not inp: break
        for j in range(len(inp)):
            P = inp[j]; Q = inp[(j+1) % len(inp)]
            cp = cross(A, B, P); cq = cross(A, B, Q)
            if cq >= 0:
                if cp < 0:
                    t = cp / (cp - cq)
                    out.append((P[0] + t*(Q[0]-P[0]), P[1] + t*(Q[1]-P[1])))
                out.append(Q)
            elif cp > 0:
                t = cp / (cp - cq)
                out.append((P[0] + t*(Q[0]-P[0]), P[1] + t*(Q[1]-P[1])))
        out = dedup(out)
    return out

def trap(gp):
    src, sa, sb, px, py, d, tf, lam, g = gp
    A = (px(sa), py(sa)); B = (px(sb), py(sb))
    ta, tb = tf(sa), tf(sb)
    P = [A, B, (B[0]+tb*d[0], B[1]+tb*d[1]), (A[0]+ta*d[0], A[1]+ta*d[1])]
    return ccw(P)

def seg_pt_d2(p, a, b):
    ab = (b[0]-a[0], b[1]-a[1]); ap = (p[0]-a[0], p[1]-a[1])
    L = dot(ab, ab); t = dot(ap, ab) / L
    t = min(max(t, Fr(0)), Fr(1))
    c = (a[0]+t*ab[0], a[1]+t*ab[1])
    return (p[0]-c[0])**2 + (p[1]-c[1])**2

def sq_dist2(A, B):
    best = None
    for P, Q in ((A, B), (B, A)):
        for v in P.V:
            for i in range(4):
                d2 = seg_pt_d2(v, Q.V[i], Q.V[(i+1) % 4])
                if best is None or d2 < best: best = d2
    return best

def frame_u(sqs, i):
    return (Fr(1), Fr(0)) if i < 0 else sqs[i].u

def pair_geom(sqs, a, b, delta):
    """R = larger phase. returns dict with exact w_e etc."""
    A, B = sqs[a], sqs[b]
    R, L = (A, B) if A.s > B.s else (B, A)
    cth = dot(R.u, L.u)
    sth = L.u[0]*R.u[1] - L.u[1]*R.u[0]   # sin(phiR - phiL)
    we = delta * sth / (cth*cth)
    h = Fr(1, 2)
    v = (R.c[0] - h*R.u[0] + h*R.n[0], R.c[1] - h*R.u[1] + h*R.n[1])
    w = (L.c[0] + h*L.u[0] + h*L.n[0], L.c[1] + h*L.u[1] + h*L.n[1])
    return dict(R=R, L=L, cth=cth, sth=sth, we=we, v=v, w=w, theta=R.phi - L.phi)

def analyze(rec, sqs, k, delta, h, sF, tlimit=900.0, ysamples=12):
    t0 = time.time()
    delta = Fr(delta); k = Fr(k); h = Fr(h)
    res = dict(anom=list(rec.anom)[:20], n_anom=len(rec.anom))
    res['conservation_ok'] = (rec.total == k)
    res['minlam'] = float(rec.minlam) if rec.minlam is not None else None
    res['D'] = float(rec.D); res['W'] = float(rec.W); res['T'] = float(rec.T); res['H'] = float(rec.H)
    res['E'] = float(sum(rec.E.values(), Fr(0))); res['M'] = float(sum((m[3] for m in rec.M), Fr(0)))
    gaps = [gp for gp in rec.gaps if gp[2] > gp[1]]
    polys = [trap(gp) for gp in gaps]
    keep = [i for i, P in enumerate(polys) if len(P) >= 3 and area(P) > 0]
    gaps = [gaps[i] for i in keep]; polys = [polys[i] for i in keep]
    n = len(polys); res['n_gap_pieces'] = n
    bb = np.array([[float(min(p[0] for p in P)), float(max(p[0] for p in P)),
                    float(min(p[1] for p in P)), float(max(p[1] for p in P))] for P in polys]) if n else np.zeros((0, 4))
    eps = 1e-12
    doubles = []   # (i, j, poly)
    same_src = 0; floor_sq = 0
    timed_out = False
    for i in range(n):
        if time.time() - t0 > tlimit: timed_out = True; break
        m = (bb[i+1:, 0] <= bb[i, 1] + eps) & (bb[i+1:, 1] >= bb[i, 0] - eps) & \
            (bb[i+1:, 2] <= bb[i, 3] + eps) & (bb[i+1:, 3] >= bb[i, 2] - eps)
        for jj in np.nonzero(m)[0]:
            j = i + 1 + int(jj)
            I = clip(polys[i], polys[j])
            if len(I) < 3 or area(I) <= 0: continue
            si, sj = gaps[i][0], gaps[j][0]
            if si == sj: same_src += 1; res.setdefault('same_src_examples', []).append((si, [(float(a), float(b)) for a, b in I][:4])); continue
            if si < 0 or sj < 0: floor_sq += 1
            doubles.append((i, j, I))
    res['timed_out_pairs'] = timed_out
    res['same_source_overlaps'] = same_src; res['floor_square_overlaps'] = floor_sq
    res['n_double_regions'] = len(doubles)
    # triples
    triples = 0; trip_ex = None
    for (i, j, I) in doubles:
        if time.time() - t0 > tlimit: timed_out = True; break
        bx = (float(min(p[0] for p in I)), float(max(p[0] for p in I)), float(min(p[1] for p in I)), float(max(p[1] for p in I)))
        m = (bb[:, 0] <= bx[1] + eps) & (bb[:, 1] >= bx[0] - eps) & (bb[:, 2] <= bx[3] + eps) & (bb[:, 3] >= bx[2] - eps)
        for kk in np.nonzero(m)[0]:
            kk = int(kk)
            if kk in (i, j): continue
            if gaps[kk][0] in (gaps[i][0], gaps[j][0]): continue
            J = clip(I, polys[kk])
            if len(J) >= 3 and area(J) > 0:
                triples += 1
                if trip_ex is None:
                    trip_ex = dict(srcs=[gaps[i][0], gaps[j][0], gaps[kk][0]], poly=[(float(a), float(b)) for a, b in J])
    res['triple_regions'] = triples; res['triple_example'] = trip_ex
    res['Mg_max'] = 3 if triples else (2 if doubles else 1)
    # (b) rectangle containment, dist, star
    pairs = {}
    worst_a = 0.0; worst_b = 0.0; worst_star = -1e9; rect_viol = []; worst_dist = 0.0
    for (i, j, I) in doubles:
        a, b = gaps[i][0], gaps[j][0]
        if a < 0 or b < 0: continue
        key = (min(a, b), max(a, b))
        if key not in pairs:
            G = pair_geom(sqs, a, b, delta)
            d2 = sq_dist2(G['R'], G['L'])
            bound2 = (delta * G['sth'] / G['cth'])**2
            G['dist_ratio'] = math.sqrt(float(d2 / bound2)) if bound2 > 0 else float('inf')
            G['dist_ok'] = d2 < bound2
            G['area'] = Fr(0)
            pairs[key] = G
        G = pairs[key]; R, L = G['R'], G['L']
        G['area'] += area(I)
        worst_dist = max(worst_dist, G['dist_ratio'])
        for z in I:
            qa = dot((z[0]-G['v'][0], z[1]-G['v'][1]), R.u)
            tR = dot((z[0]-G['v'][0], z[1]-G['v'][1]), R.n)
            m_ = dot((G['w'][0]-z[0], G['w'][1]-z[1]), L.u)
            tL = dot((z[0]-G['w'][0], z[1]-G['w'][1]), L.n)
            ra = float(qa / G['we']) if G['we'] > 0 else float('inf'); rb = float(tR / delta)
            worst_a = max(worst_a, ra); worst_b = max(worst_b, rb)
            okr = (0 <= qa <= G['we']) and (0 <= tR <= delta)
            if not okr: rect_viol.append(dict(pair=key, z=(float(z[0]), float(z[1])), qa_over_we=ra, tR_over_delta=rb, qa=float(qa)))
            rhs = G['sth'] * max(tR, tL / G['cth'])
            lhs = G['cth'] * qa + m_
            if rhs > 0: worst_star = max(worst_star, float(lhs / rhs))
            elif lhs > 0: worst_star = float('inf')
    res['n_pairs'] = len(pairs)
    res['rect_qa_over_we_max'] = worst_a; res['rect_tR_over_delta_max'] = worst_b
    res['rect_violations'] = rect_viol[:10]; res['n_rect_violations'] = len(rect_viol)
    res['dist_ratio_max'] = worst_dist
    res['dist_violations'] = [k_ for k_, G in pairs.items() if not G['dist_ok']]
    res['star_ratio_max'] = worst_star
    res['De_over_delta_we_max'] = max([float(G['area'] / (delta * G['we'])) for G in pairs.values()], default=0.0)
    sum_we = sum((G['we'] for G in pairs.values()), Fr(0))
    res['sum_we'] = float(sum_we)
    # (c) death
    waste = k * h
    strip = [(Fr(0), Fr(0)), (k, Fr(0)), (k, h), (Fr(0), h)]
    for S in sqs:
        I = clip(ccw(list(S.V)), strip)
        if len(I) >= 3: waste -= area(I)
    res['waste_strip'] = float(waste)
    dD = delta * rec.D
    res['deltaD'] = float(dD); res['Dcells_area'] = float(rec.Dcells_area)
    res['jacobian_ok'] = rec.Dcells_area >= dD
    rhs = waste + delta * sum_we
    res['death_ratio_vs_waste_plus_rects'] = float(rec.Dcells_area / rhs) if rhs > 0 else None
    res['death_ok'] = rec.Dcells_area <= rhs
    intMg = sum((area(P) for P in polys), Fr(0))      # = integral of M_g over the waste (all live gap segments)
    res['intMg_over_waste_plus_rects'] = float(intMg / rhs) if rhs > 0 else None
    res['intMg_ok'] = intMg <= rhs
    N = len(sqs); Wtot = k*k - N
    res['W_container'] = float(Wtot)
    res['D_over_bound_global'] = float(rec.D / ((1 + Fr(5625, 100) * delta**2) * Wtot / delta))
    # (d) merges
    mg = {}
    floor_merge = 0
    for (Y, ws, ls, meas, slen) in rec.M:
        if ws < 0 or ls < 0: floor_merge += 1; continue
        key = (min(ws, ls), max(ws, ls))
        mg.setdefault((Y, key), Fr(0)); mg[(Y, key)] += meas
    worst_m = 0.0; m_viol = []; ycount = {}
    for (Y, key), meas in mg.items():
        G = pairs.get(key)
        if G is None:
            G = pair_geom(sqs, key[0], key[1], delta)
        cyr = abs(dot(sqs[Y].u, G['R'].u))
        bnd = G['we'] / cyr
        r = float(meas / bnd) if bnd > 0 else float('inf')
        worst_m = max(worst_m, r)
        if meas > bnd: m_viol.append(dict(Y=Y, pair=key, meas=float(meas), bound=float(bnd)))
        ycount.setdefault(key, set()).add(Y)
    res['floor_merges'] = floor_merge
    res['merge_ratio_max'] = worst_m; res['merge_violations'] = m_viol[:10]
    res['max_Y_per_pair'] = max([len(v) for v in ycount.values()], default=0)
    amax = max([abs(S.phi) for S in sqs if abs(S.s) < Fr(sF)] + [0.0])
    gb = 4 / math.cos(math.pi/4 + amax) * float(sum_we) if sum_we > 0 else 0.0
    res['merge_global_ratio'] = (res['M'] / gb) if gb > 0 else (0.0 if res['M'] == 0 else float('inf'))
    # (f) E per (src, target)
    worst_e = 0.0; e_viol = []
    for (src, Y), meas in rec.E.items():
        ux = frame_u(sqs, src); uy = sqs[Y].u
        cD = dot(ux, uy); sD = ux[0]*uy[1] - ux[1]*uy[0]
        if cD <= 0: e_viol.append(dict(src=src, Y=Y, note='cos<=0')); continue
        bnd = delta * abs(sD) / cD
        r = float(meas / bnd) if bnd > 0 else float('inf')
        worst_e = max(worst_e, r)
        if meas > bnd: e_viol.append(dict(src=src, Y=Y, meas=float(meas), bound=float(bnd)))
    res['E_ratio_max'] = worst_e; res['E_violations'] = e_viol[:10]
    # (e) Ov_gap(y) on sampled heights
    ys = set()
    for (i, j, I) in doubles:
        yl = sorted(set(p[1] for p in I))
        for a_, b_ in zip(yl[:-1], yl[1:]): ys.add((a_ + b_) / 2)
    ys = sorted(ys)
    if len(ys) > ysamples:
        step = len(ys) / ysamples; ys = [ys[int(i*step)] for i in range(ysamples)]
    ys += [h * Fr(i, 7) + Fr(1, 997) for i in range(1, 7)]
    ov_worst = 0.0; ov_viol = []; supp_viol = 0
    bound_ov = sum((G['we'] / G['R'].cs for G in pairs.values()), Fr(0))
    for y in ys:
        if time.time() - t0 > tlimit: timed_out = True; break
        ivs = []
        for gp in gaps:
            src, sa, sb, px, py, d, tf, lam, g = gp
            # t(x) = (y - py(x))/d1 ; need 0 < t < tf
            tx = (py.scale(-1) + y).scale(1 / d[1])
            lo, hi = sa, sb
            ok = True
            for f in (tx, tf - tx):   # f(x) > 0
                if f.c1 == 0:
                    if f.c0 <= 0: ok = False
                else:
                    r_ = f.root()
                    if f.c1 > 0: lo = max(lo, r_)
                    else: hi = min(hi, r_)
            if not ok or lo >= hi: continue
            X = px + tx.scale(d[0])
            X1, X2 = sorted([X(lo), X(hi)])
            dens = abs(1 / X.c1)
            ivs.append((X1, X2, dens))
        pts = sorted(set([a for a, b, c in ivs] + [b for a, b, c in ivs]))
        ov = Fr(0)
        for a_, b_ in zip(pts[:-1], pts[1:]):
            rho = sum((c for (A_, B_, c) in ivs if A_ <= a_ and b_ <= B_), Fr(0))
            if rho > 1:
                ov += (rho - 1) * (b_ - a_)
                zm = ((a_ + b_) / 2, y)
                inside = False
                for G in pairs.values():
                    qa = dot((zm[0]-G['v'][0], zm[1]-G['v'][1]), G['R'].u)
                    tR = dot((zm[0]-G['v'][0], zm[1]-G['v'][1]), G['R'].n)
                    if 0 <= qa <= G['we'] and 0 <= tR <= delta: inside = True; break
                if not inside: supp_viol += 1
        r = float(ov / bound_ov) if bound_ov > 0 else (0.0 if ov == 0 else float('inf'))
        ov_worst = max(ov_worst, r)
        if ov > bound_ov: ov_viol.append(dict(y=float(y), ov=float(ov), bound=float(bound_ov)))
    res['Ov_ratio_max'] = ov_worst; res['Ov_violations'] = ov_viol[:5]; res['Ov_support_violations'] = supp_viol
    res['timed_out'] = timed_out
    res['analysis_sec'] = time.time() - t0
    return res
