"""Exact (mpmath 50 digits) re-check of boundary configurations found by geom_search.py (outside the regime
of Section 7, where delta <= 1e-5).
(1) The conclusion of Lemma 7.5 fails just above delta = cos(theta)/2  (theta = 0.6); cf. Remark 7.6.
(2) Inequality (7.3) of Lemma 7.8 fails for delta >= ~cos(theta)  (theta = 1.4, delta = 0.45).
Each: two closed unit squares R (phase theta), L (phase 0), point q with source delays tR, tL; check disjointness by SAT
with a strictly positive margin, ray clearness, and the violated inequality.

Usage:  python boundary_exact.py        (writes out/boundary_exact.json; about 1 s)"""
import mpmath as mp, json, os
mp.mp.dps = 50


def out_path(name):
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)

def frame(phi):
    return (mp.cos(phi), mp.sin(phi)), (-mp.sin(phi), mp.cos(phi))

def sq(q, phi, tau, t):
    u, n = frame(phi)
    c = (q[0] - (mp.mpf(1)/2 + t)*n[0] - tau*u[0], q[1] - (mp.mpf(1)/2 + t)*n[1] - tau*u[1])
    V = [(c[0] + a*u[0]/2 + b*n[0]/2, c[1] + a*u[1]/2 + b*n[1]/2) for a, b in ((-1,-1),(1,-1),(1,1),(-1,1))]
    return V, u, n

def sep(A, B, extra_axes=()):
    best = mp.mpf(-10)
    for ax in [A[1], A[2], B[1], B[2]] + list(extra_axes):
        pa = [v[0]*ax[0] + v[1]*ax[1] for v in A[0]]; pb = [v[0]*ax[0] + v[1]*ax[1] for v in B[0]]
        best = max(best, max(min(pb) - max(pa), min(pa) - max(pb)))
    return best

def seg_sep(P, Q, B):
    d = (Q[0]-P[0], Q[1]-P[1]); L = mp.sqrt(d[0]**2 + d[1]**2); d = (d[0]/L, d[1]/L); nn = (-d[1], d[0])
    A = ([P, Q], d, nn)
    return sep(A, B)

out = {}
# (1) Lemma 7.5: theta=0.6, case A (qa >= sin(theta)/2) with tR slightly above cos(theta)/2
th = mp.mpf('0.6'); q = (mp.mpf(0), mp.mpf(0))
tR = mp.cos(th)/2 * (1 + mp.mpf('1e-6')); tL = mp.mpf('0.2')
tauR = -mp.mpf(1)/2 + mp.sin(th)/2 * (1 + mp.mpf('1e-9'))     # qa = tauR + 1/2 slightly above sin/2
best = None
for tauL in [mp.mpf(1)/2 - mp.mpf(10)**(-e) for e in range(3, 13)]:
    R = sq(q, th, tauR, tR); L = sq(q, mp.mpf(0), tauL, tL)
    s = sep(R, L)
    if best is None or s > best[0]: best = (s, tauL)
s, tauL = best
R = sq(q, th, tauR, tR); L = sq(q, mp.mpf(0), tauL, tL)
pR = (q[0] - tR*R[2][0], q[1] - tR*R[2][1]); pL = (q[0] - tL*L[2][0], q[1] - tL*L[2][1])
out['two_sources_fail'] = dict(theta=str(th), delta_needed=mp.nstr(max(tR, tL), 20), half_cos=mp.nstr(mp.cos(th)/2, 20),
                           sep=mp.nstr(s, 10), qa=mp.nstr(tauR + mp.mpf(1)/2, 15), half_sin=mp.nstr(mp.sin(th)/2, 15),
                           ray_R_clear_of_L=mp.nstr(seg_sep(pR, q, L), 8), ray_L_clear_of_R=mp.nstr(seg_sep(pL, q, R), 8),
                           R_vertices=[(mp.nstr(v[0], 20), mp.nstr(v[1], 20)) for v in R[0]],
                           L_vertices=[(mp.nstr(v[0], 20), mp.nstr(v[1], 20)) for v in L[0]])
# (2) star: theta=1.4, from geom_search.py star (delta=0.45)
th = mp.mpf('1.4'); d = mp.mpf('0.45')
x = [mp.mpf('0.49991869552904733'), mp.mpf('0.21587702672659614'), mp.mpf('0.99992509328384')*d, mp.mpf('8.719748499014246e-05')*d]
R = sq(q, th, x[0], x[2]); L = sq(q, mp.mpf(0), x[1], x[3])
k, sg = mp.cos(th), mp.sin(th)
qa = x[0] + mp.mpf(1)/2; m = mp.mpf(1)/2 - x[1]
F = k*qa + m - sg*max(x[2], x[3]/k)
pR = (q[0] - x[2]*R[2][0], q[1] - x[2]*R[2][1]); pL = (q[0] - x[3]*L[2][0], q[1] - x[3]*L[2][1])
out['star_fail'] = dict(theta='1.4', delta='0.45', cos_theta=mp.nstr(k, 15), tR=mp.nstr(x[2], 15), tL=mp.nstr(x[3], 15),
                        sep=mp.nstr(sep(R, L), 10), star_lhs_minus_rhs=mp.nstr(F, 15),
                        ray_R_clear_of_L=mp.nstr(seg_sep(pR, q, L), 8), ray_L_clear_of_R=mp.nstr(seg_sep(pL, q, R), 8))
print(json.dumps(out, indent=1))
json.dump(out, open(out_path('boundary_exact.json'), 'w'), indent=1)
