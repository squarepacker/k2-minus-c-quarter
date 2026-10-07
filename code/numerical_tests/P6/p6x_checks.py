# p6x_checks.py -- checks of Section 8 of the manuscript.  The labels L4..L7 are the names used in
# the output files; they refer to the following statements:
#   L4 = Lemma 8.5 (density): w_X is injective on P_e and |I_e| >= m_e;
#   L5 = Lemma 8.6 (the hull slit is waste): (i) the hull endpoints lie in top(X) / bot(Y) and
#        0 <= t(w) <= delta; (ii) the open set O_e meets no square;
#   L6 = Lemma 8.7 (vertical comparison P6.V) and Remark 8.8 (the exceptional set);
#   L7 = Lemma 8.9 (disjoint vertical slits: (a) floor and ceiling families share no square,
#        (c) the slits S_e are pairwise disjoint);
#   chain = Step 1 of the proof of Theorem 8.3 (sum of theta'_e along a path >= min(beta-bar, a(Z_j)));
# plus a pointwise audit of rule R1 (Definition 3.14) with an independent point tracer.
from fractions import Fraction as Fr
import bisect
import random
import mpmath as mp
from p6x_core import (HALF, mpf, interiors_disjoint, sat_gap, poly_area2, dedup, bbox, bb_hit,
                      hit, trace_point, theta)


def frame(flow, X):
    if X == -1:
        return (Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(0), Fr(1))
    S = flow.byid[X]
    o = (S.c[0] + HALF * S.n[0], S.c[1] + HALF * S.n[1])
    return o, S.u, S.n


def pt(o, u, n, w, s):
    return (o[0] + w * u[0] + s * n[0], o[1] + w * u[1] + s * n[1])


def check_pairs(flow, tiny=False, polys_out=None, k=None, mut=None):
    """returns dict of stats + list of violations. polys_out: list to append (fam,(X,Y),poly_real)."""
    delta = flow.delta
    bbar = mp.atan(mpf(flow.Ta))          # beta-bar := alpha_max (in the scaled analogue beta-bar is replaced
                                          # by the tilt threshold alpha_max of the flow; manuscript: 1.5e-6)
    viol = []
    st = dict(pairs=0, pairs_sq=0, L4_slack_min=None, L4_strict=0, inj_overlap=0,
              L5_tested=0, L5_sqtests=0, L5_gap_min=None, L5_gap_min_where=None,
              L6_tested=0, L6_sharp_min=None, L6_weak_min=None, L6_exc_ratio_max=None,
              L6_gv_den_max=None, L6_gv_max=None, L6_ray_tests=0, L6_holes=0)
    sqs = flow.sqs
    for (X, Yid), pieces in flow.pairs.items():
        st['pairs'] += 1
        if X != -1:
            st['pairs_sq'] += 1
        Y = flow.byid[Yid]
        o, u, n = frame(flow, X)
        cosf, sinf = u[0], u[1]
        # ---- L4 (Lemma 8.5): pieces = exit-coordinate intervals of the paths of P_e
        m = sum(b - a for (a, b, _, _) in pieces)
        ws = sorted((wa, wb) for (_, _, wa, wb) in pieces)
        for i in range(len(ws) - 1):
            if ws[i][1] > ws[i + 1][0]:
                st['inj_overlap'] += 1
                viol.append(('L4_injectivity', X, Yid, float(ws[i][1] - ws[i + 1][0])))
        I = sum(b - a for (a, b) in ws)
        if m > I:
            viol.append(('L4_density', X, Yid, float(m - I)))
        sl = float(I - m)
        if I > m:
            st['L4_strict'] += 1
        if st['L4_slack_min'] is None or sl < st['L4_slack_min']:
            st['L4_slack_min'] = sl
        w1, w2 = ws[0][0], ws[-1][1]
        if mut == 'widen':      # mutation test: widen hull by 0.15 each side (should be caught)
            w1, w2 = w1 - Fr(15, 100), w2 + Fr(15, 100)
        # holes in I_e?
        merged_gaps = 0
        cur = ws[0][1]
        for (a, b) in ws[1:]:
            if a > cur:
                merged_gaps += 1
            cur = max(cur, b)
        if merged_gaps:
            st['L6_holes'] += 1
        # categories (coverage statistics)
        if X != -1 and abs(u[1] / u[0]) > flow.Ta / 2:
            st['cat_Xtilt_gt_half_alpha'] = st.get('cat_Xtilt_gt_half_alpha', 0) + 1
        if Y.tan_abs > Fr(4, 5):
            st['cat_Y_near45'] = st.get('cat_Y_near45', 0) + 1
        if Y.tan_abs >= flow.Ta:
            st['cat_Y_terminator'] = st.get('cat_Y_terminator', 0) + 1
        if X == -1:
            st['cat_floor_pair'] = st.get('cat_floor_pair', 0) + 1
        if len(ws) > 1:
            st['cat_multi_piece'] = st.get('cat_multi_piece', 0) + 1
        # ---- t(w) of (8.3): the line l_w meets the line L_Y of bot(Y) at o + w u + t(w) n
        nn = n[0] * Y.n[0] + n[1] * Y.n[1]
        un = u[0] * Y.n[0] + u[1] * Y.n[1]
        oc = (o[0] - Y.c[0]) * Y.n[0] + (o[1] - Y.c[1]) * Y.n[1]
        t0 = (-HALF - oc) / nn
        kap = -un / nn

        def t(w):
            return t0 + kap * w
        if not (t(w1) >= 0 and t(w2) >= 0 and t(w1) <= delta and t(w2) <= delta):
            viol.append(('L5_t_range', X, Yid, float(t(w1)), float(t(w2))))
        # hull endpoints in closed top(X) / bot(Y)
        if X != -1 and not (w1 >= -HALF and w2 <= HALF):
            viol.append(('L5_hull_outside_topX', X, Yid))
        mY = (Y.c[0] - HALF * Y.n[0], Y.c[1] - HALF * Y.n[1])
        for w in (w1, w2):
            q = pt(o, u, n, w, t(w))
            tau = (q[0] - mY[0]) * Y.u[0] + (q[1] - mY[1]) * Y.u[1]
            if abs(tau) > HALF:
                viol.append(('L5_hull_outside_botY', X, Yid, float(tau)))
        # ---- L5 (Lemma 8.6(ii)): the polygon O_e over the closed hull [w1, w2] must meet no square
        Om = dedup([pt(o, u, n, w1, 0), pt(o, u, n, w2, 0), pt(o, u, n, w2, t(w2)), pt(o, u, n, w1, t(w1))])
        if len(Om) >= 3 and poly_area2(Om) != 0:
            st['L5_tested'] += 1
            for p_ in Om:
                if not (0 <= p_[0] <= flow.k and 0 <= p_[1] <= flow.k):
                    viol.append(('L5_outside_box', X, Yid))
            bo = bbox(Om)
            for Z in sqs:
                if not bb_hit(bo, Z.bb, 0.25):
                    continue
                st['L5_sqtests'] += 1
                if not interiors_disjoint(Om, Z.V):
                    viol.append(('L5_square_in_Omega', X, Yid, Z.id, sat_gap(Om, Z.V)))
                if Z.id not in (X, Yid):
                    gp = sat_gap(Om, Z.V) / float(delta)
                    if gp < 0.01:
                        st['cat_foreign_within_0.01delta'] = st.get('cat_foreign_within_0.01delta', 0) + 1
                    if st['L5_gap_min'] is None or gp < st['L5_gap_min']:
                        st['L5_gap_min'] = gp
                        st['L5_gap_min_where'] = (X, Yid, Z.id)
        # ---- L6 (Lemma 8.7, P6.V): den = chi_e, De = Delta_e, G = G_e, Jp = |J'_e|;
        #      'sharp'  : |J'_e| >= cos(phi_X) (m_e - Delta_e)              (first inequality of 8.7(c))
        #      'weak'   : |J'_e| >= cos(beta-bar) m_e - 2.1 delta beta-bar   (a weaker form of the second
        #                 inequality of 8.7(c); the manuscript has 1.0000017 in place of 2.1)
        #      'exc'    : excluded end interval <= cos(phi_X) Delta_e         (Remark 8.8)
        #      'gv'     : g^v_e < delta / chi_e                               (Lemma 8.7(a))
        den = cosf - kap * sinf
        if den <= 0:
            viol.append(('L6_den_nonpos', X, Yid, float(den)))
            continue
        De = delta * abs(sinf) / den
        a_ = cosf / den
        b_ = sinf * t0 / den       # w_c(w0) = a_ w0 + b_
        glo = max(w1, (w1 - b_) / a_)
        ghi = min(w2, (w2 - b_) / a_)
        if mut == 'noDtrim':    # mutation test: forget the exceptional end (should be caught)
            glo, ghi = w1, w2
        G = max(Fr(0), ghi - glo)
        Jp = cosf * G
        st['L6_tested'] += 1
        sharp = Jp - cosf * (m - De)
        if sharp < 0:
            viol.append(('L6_sharp', X, Yid, float(sharp)))
        scale = mpf(delta) * bbar
        sm = mpf(sharp) / scale
        if st['L6_sharp_min'] is None or sm < st['L6_sharp_min']:
            st['L6_sharp_min'] = sm
        weak = mpf(Jp) - (mp.cos(bbar) * mpf(m) - mp.mpf('2.1') * mpf(delta) * bbar)
        if weak < 0:
            viol.append(('L6_weak_form', X, Yid, float(weak)))
        dm = weak / scale
        if st['L6_weak_min'] is None or dm < st['L6_weak_min']:
            st['L6_weak_min'] = dm
        if De > 0:
            exc = (cosf * (w2 - w1) - Jp) / (cosf * De)
            if exc > 1:
                viol.append(('L6_exc_gt_bound', X, Yid, float(exc)))
            if st['L6_exc_ratio_max'] is None or exc > st['L6_exc_ratio_max']:
                st['L6_exc_ratio_max'] = float(exc)
        if G > 0:
            lam_lo = t(glo) / den
            lam_hi = t(ghi) / den
            for lam in (lam_lo, lam_hi):
                r = lam * den / delta
                if r > 1:
                    viol.append(('L6_gv_gt_delta_over_den', X, Yid, float(r)))
                if st['L6_gv_den_max'] is None or r > st['L6_gv_den_max']:
                    st['L6_gv_den_max'] = float(r)
                r2 = lam / delta
                if st['L6_gv_max'] is None or r2 > st['L6_gv_max']:
                    st['L6_gv_max'] = float(r2)
                if tiny and r2 >= Fr(10000016, 10000000):
                    viol.append(('L6_gv_ge_1.0000016delta', X, Yid, float(r2)))
            b1, b2 = pt(o, u, n, glo, 0), pt(o, u, n, ghi, 0)
            c1, c2 = (b1[0], b1[1] + lam_lo), (b2[0], b2[1] + lam_hi)
            # endpoints of c in closed bot(Y)
            for c in (c1, c2):
                tau = (c[0] - mY[0]) * Y.u[0] + (c[1] - mY[1]) * Y.u[1]
                on = (c[0] - mY[0]) * Y.n[0] + (c[1] - mY[1]) * Y.n[1]
                if on != 0 or abs(tau) > HALF:
                    viol.append(('L6_c_not_on_botY', X, Yid, float(tau), float(on)))
            V = dedup([b1, b2, c2, c1])
            if len(V) >= 3 and poly_area2(V) != 0:
                bv = bbox(V)
                for Z in sqs:
                    if bb_hit(bv, Z.bb, 1e-6) and not interiors_disjoint(V, Z.V):
                        viol.append(('L6_V_hits_square', X, Yid, Z.id))
                if polys_out is not None:
                    if flow.fam == 'floor':
                        Vr = V
                    else:
                        Vr = [(p_[0], k - p_[1]) for p_ in V]
                    polys_out.append((flow.fam, (X, Yid), Vr))
            # independent vertical ray casts from sample points of G
            for fr_ in (Fr(1, 10 ** 6), Fr(1, 7), Fr(1, 2), Fr(5, 6), 1 - Fr(1, 10 ** 6)):
                w0 = glo + fr_ * (ghi - glo)
                b = pt(o, u, n, w0, 0)
                st['L6_ray_tests'] += 1
                best = None
                for Z in sqs:
                    if Z.id == X:
                        continue
                    if not (Z.bb[0] - 1e-6 <= float(b[0]) <= Z.bb[2] + 1e-6):
                        continue
                    if Z.bb[3] < float(b[1]) - 1e-6:
                        continue
                    try:
                        r = hit(Z, b, (Fr(0), Fr(1)))
                    except RuntimeError:
                        viol.append(('L6_raycast_start_inside', X, Yid, Z.id))
                        r = None
                    if r is not None and (best is None or r[0] < best[0]):
                        best = (r[0], Z.id, r[1], r[2])
                lam = t(w0) / den
                if best is None or best[1] != Yid or best[2] != 'bot' or best[3] or best[0] != lam:
                    viol.append(('L6_raycast', X, Yid, None if best is None else (best[1], best[2], float(best[0])),
                                 float(lam)))
                # b is upper end of X's vertical chord
                if X != -1:
                    S = flow.byid[X]
                    eps = Fr(1, 10 ** 12)
                    up = (b[0], b[1] + eps)
                    dn = (b[0], b[1] - eps)

                    def inside(S, p):
                        return all(nu[0] * p[0] + nu[1] * p[1] <= ce for (_, nu, ce) in S.E)
                    if inside(S, up) or not inside(S, dn):
                        viol.append(('L6_b_not_upper_end', X, Yid))
                else:
                    for Z in sqs:
                        if Z.id != Yid and all(nu[0] * b[0] + nu[1] * b[1] <= ce for (_, nu, ce) in Z.E):
                            viol.append(('L6_floor_b_in_square', Yid, Z.id))
    return st, viol


def check_L7(polys, sqsets):
    """Lemma 8.9: (c) the slits S_e of all pairs of both families are pairwise disjoint (exact test of
    every two slits whose x-ranges overlap); (a) no square is a member of a floor pair and of a ceiling pair.
    polys: list of (fam,(X,Y),poly real coords). sqsets: dict fam->set of square ids used in pairs."""
    viol = []
    tests = 0
    n = len(polys)
    bbs = [bbox(p[2]) for p in polys]
    order = sorted(range(n), key=lambda i: bbs[i][0])
    for ii in range(n):
        i = order[ii]
        for jj in range(ii + 1, n):
            j = order[jj]
            if bbs[j][0] > bbs[i][2] + 1e-9:
                break
            # all pairs of slits sharing vertical lines (x-ranges overlap) are tested exactly
            if polys[i][:2] == polys[j][:2]:
                continue
            tests += 1
            if not interiors_disjoint(polys[i][2], polys[j][2]):
                viol.append(('L7_overlap', polys[i][:2], polys[j][:2]))
    sh = sqsets.get('floor', set()) & sqsets.get('ceil', set())
    if sh:
        viol.append(('L7_families_share', sorted(sh)))
    return tests, viol


def check_chain(flow):
    """Step 1 of the proof of Theorem 8.3: for every R1-passed entry (Z_j, q_j) of every path,
    sum_{i<=j} min(theta(Z_{i-1}, Z_i), beta-bar) >= min(beta-bar, a(Z_j)).  Also returns sum_e theta'_e m_e."""
    bbar = mp.atan(mpf(flow.Ta))
    viol = []
    mn = None
    phi = {S.id: S.phi for S in flow.sqs}
    phi[-1] = mp.mpf(0)
    cnt = 0
    lhs = mp.mpf(0)
    for (Yid, xa, xb, hist) in flow.entries:
        cnt += 1
        S = mp.mpf(0)
        prev = -1
        for z in hist:
            S += min(theta(phi[prev], phi[z]), bbar)
            prev = z
        beta = min(bbar, abs(phi[Yid]))
        mg = S - beta
        if mg < -mp.mpf(10) ** -30:
            viol.append(('chain', Yid, float(mg)))
        if mn is None or mg < mn:
            mn = mg
    # aggregate: sum theta'_e m_e
    Sm = mp.mpf(0)
    for (X, Y), pieces in flow.pairs.items():
        m = sum(b - a for (a, b, _, _) in pieces)
        Sm += min(theta(phi[X], phi[Y]), bbar) * mpf(m)
    return cnt, mn, viol, Sm


def pointwise_audit(flow, nsamp, rng):
    """compare beam tracer with independent pointwise tracer; verify R1 winners (Definition 3.14:
    the winner at an entry candidate has the lexicographically smallest (g, x))."""
    sqs = flow.sqs
    terms = sorted(flow.terms, key=lambda t: t[1])
    starts = [t[1] for t in terms]
    viol = []
    nM = 0
    done = 0
    xs_list = [Fr(rng.randrange(1, 10 ** 9), 10 ** 9) * flow.k for _ in range(nsamp)]
    Mp = [t for t in terms if t[0] == 'M']
    rng.shuffle(Mp)
    for t in Mp[:30]:
        xs_list.append(t[1] + (t[2] - t[1]) * Fr(rng.randrange(1, 10 ** 6), 10 ** 6))
    Tp = [t for t in terms if t[0] in ('T', 'E')]
    rng.shuffle(Tp)
    for t in Tp[:15]:
        xs_list.append(t[1] + (t[2] - t[1]) * Fr(rng.randrange(1, 10 ** 6), 10 ** 6))
    for x in xs_list:
        i = bisect.bisect_right(starts, x) - 1
        if i < 0:
            continue
        kind, a, b, hist, z = terms[i]
        if not (a < x < b):
            continue          # measure-zero boundary
        done += 1
        kd, ents, zz = trace_point(sqs, x, flow.k, flow.delta, flow.Ta, flow.h)
        seq = [e[0] for e in ents]
        if kind != 'M':
            if kd != kind or tuple(seq) != tuple(hist) or (kind in ('E', 'T') and zz != z):
                viol.append(('pointwise_mismatch', float(x), kind, kd, hist, tuple(seq)))
                continue
        else:
            nM += 1
            L = len(hist)
            if tuple(seq[:L]) != tuple(hist) or len(seq) <= L or seq[L] != z:
                viol.append(('pointwise_M_mismatch', float(x), hist, tuple(seq)))
                continue
            Yid, q, gx = ents[L]
            Y = flow.byid[Yid]
            mY = (Y.c[0] - HALF * Y.n[0], Y.c[1] - HALF * Y.n[1])
            tau = (q[0] - mY[0]) * Y.u[0] + (q[1] - mY[1]) * Y.u[1]
            wpc = [w for w in flow.winners[Yid] if w[0] < tau < w[1]]
            if len(wpc) != 1:
                viol.append(('R1_no_unique_winner', float(x), Yid, len(wpc)))
                continue
            _, _, xf, gf = wpc[0]
            xs = xf[0] + xf[1] * tau
            gs = gf[0] + gf[1] * tau
            if not ((gs, xs) < (gx, x)):
                viol.append(('R1_winner_not_lexmin', float(x), float(xs)))
            kd2, ents2, _ = trace_point(sqs, xs, flow.k, flow.delta, flow.Ta, flow.h)
            if not any(e[0] == Yid and e[1] == q and e[2] == gs for e in ents2):
                viol.append(('R1_winner_not_arriving', float(x), float(xs)))
        # every entry of x (winner) must be its own winner piece
        for (Yid, q, gx) in ents[:len(hist)]:
            Y = flow.byid[Yid]
            mY = (Y.c[0] - HALF * Y.n[0], Y.c[1] - HALF * Y.n[1])
            tau = (q[0] - mY[0]) * Y.u[0] + (q[1] - mY[1]) * Y.u[1]
            wpc = [w for w in flow.winners[Yid] if w[0] < tau < w[1]]
            if len(wpc) != 1 or wpc[0][2][0] + wpc[0][2][1] * tau != x:
                viol.append(('R1_entry_not_self_winner', float(x), Yid))
    return done, nM, viol
