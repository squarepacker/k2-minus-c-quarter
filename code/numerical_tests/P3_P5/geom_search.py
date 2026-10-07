"""Lemma-level boundary search (written from scratch, floats; differential evolution from scipy).
Only geometry: squares X_i with phase phi_i, source point p_i in relint top(X_i), q = p_i + t_i n_i.

usage: python geom_search.py <mode> [time_limit_sec] [comma-separated parameter values]
  two      Lemma 7.5 (two sources -> near corners): smallest delay max(t_R, t_L) at which its conclusion
           fails, per phase difference theta, compared with cos(theta)/2 (Remark 7.6)
  three    Lemma 7.7 (three sources): smallest delay for three disjoint squares with phases in [-alpha, alpha]
           sharing a point, compared with cos(2 alpha)/2
  star     inequality (7.3) of Lemma 7.8: largest excess (lhs - rhs) / (delta tan theta) per theta and delta
  starmin  smallest delay at which (7.3) fails, per theta
  writes out/geom_<mode>_<values or all>.json
"""
import numpy as np, json, sys, time, math, os
from scipy.optimize import differential_evolution
import mpmath as mp

T0 = time.time(); TLIM = float(sys.argv[2]) if len(sys.argv) > 2 else 1500.0


def out_path(name):
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)

def frame(phi):
    u = np.array([math.cos(phi), math.sin(phi)]); n = np.array([-math.sin(phi), math.cos(phi)])
    return u, n

def square_from_source(q, phi, tau, t):
    u, n = frame(phi)
    c = q - (0.5 + t) * n - tau * u
    V = [c + a*0.5*u + b*0.5*n for a, b in ((-1,-1),(1,-1),(1,1),(-1,1))]
    return np.array(V), u, n

def sep(VA, axA, VB, axB):
    """SAT separation (>0 iff disjoint closed convex sets); axes given as list of unit vectors."""
    best = -1e9
    for a in list(axA) + list(axB):
        pa = VA @ a; pb = VB @ a
        g = max(pb.min() - pa.max(), pa.min() - pb.max())
        best = max(best, g)
    return best

def seg_axes(P, Q):
    d = Q - P; L = np.linalg.norm(d)
    if L < 1e-15: return [np.array([1.0, 0.0]), np.array([0.0, 1.0])]
    d = d / L; return [d, np.array([-d[1], d[0]])]

PEN = 1e3

def two_sources_obj(x, theta, case, clear):
    tauR, tauL, tR, tL = x
    q = np.zeros(2)
    VR, uR, nR = square_from_source(q, theta, tauR, tR)
    VL, uL, nL = square_from_source(q, 0.0, tauL, tL)
    s = sep(VR, [uR, nR], VL, [uL, nL])
    pen = max(0.0, 1e-9 - s)
    sig = math.sin(theta)
    qa = tauR + 0.5; m = 0.5 - tauL
    if case == 'A': pen += max(0.0, 0.5*sig - qa)
    else: pen += max(0.0, 0.5*sig - m)
    if clear:
        pR = q - tR*nR; pL = q - tL*nL
        segR = np.array([pR, q]); segL = np.array([pL, q])
        pen += max(0.0, 1e-9 - sep(segR, seg_axes(pR, q), VL, [uL, nL]))
        pen += max(0.0, 1e-9 - sep(segL, seg_axes(pL, q), VR, [uR, nR]))
    return max(tR, tL) + PEN * pen

def three_sources_obj(x, alpha, clear):
    tA, tB, tC, uA, uB, uC, b1, b2 = x
    # phases: A = alpha, B = alpha - b1*2alpha, C = alpha - (b1+b2)*2alpha*... keep within (-alpha, alpha)
    phA = alpha; phB = alpha - 2*alpha*b1; phC = phB - 2*alpha*b2
    pen = max(0.0, -alpha - phC + 1e-12) * PEN
    q = np.zeros(2)
    sq = [square_from_source(q, ph, tau, t) for ph, tau, t in ((phA, uA, tA), (phB, uB, tB), (phC, uC, tC))]
    for i in range(3):
        for j in range(i+1, 3):
            s = sep(sq[i][0], sq[i][1:], sq[j][0], sq[j][1:])
            pen += PEN * max(0.0, 1e-9 - s)
    if clear:
        ts = (tA, tB, tC)
        for i in range(3):
            p = q - ts[i]*sq[i][2]
            seg = np.array([p, q])
            for j in range(3):
                if j == i: continue
                pen += PEN * max(0.0, 1e-9 - sep(seg, seg_axes(p, q), sq[j][0], sq[j][1:]))
    return max(tA, tB, tC) + pen

def star_obj(x, theta, delta):
    """return -(star excess)/(delta tan theta) with penalty; sources R(phi=theta), L(phi=0); t in (0,delta)."""
    tauR, tauL, tR, tL = x
    tR *= delta; tL *= delta
    q = np.zeros(2)
    VR, uR, nR = square_from_source(q, theta, tauR, tR)
    VL, uL, nL = square_from_source(q, 0.0, tauL, tL)
    s = sep(VR, [uR, nR], VL, [uL, nL])
    k, sg = math.cos(theta), math.sin(theta)
    qa = tauR + 0.5; m = 0.5 - tauL
    F = k*qa + m - sg*max(tR, tL/k)
    if s <= 0.0:
        return 1e6 + 1e6 * min(1.0, -s)      # hard rejection: R, L must be disjoint
    return -F/(delta*math.tan(theta))

def starmin_obj(x, theta, clear):
    """minimize max(tR,tL) subject to (star) violated (lhs >= rhs + 1e-9) and R, L disjoint."""
    tauR, tauL, tR, tL = x
    q = np.zeros(2)
    VR, uR, nR = square_from_source(q, theta, tauR, tR)
    VL, uL, nL = square_from_source(q, 0.0, tauL, tL)
    s = sep(VR, [uR, nR], VL, [uL, nL])
    k, sg = math.cos(theta), math.sin(theta)
    qa = tauR + 0.5; m = 0.5 - tauL
    T = max(tR, tL, 1e-15)
    rhs = sg*max(tR, tL/k)
    viol_rel = (k*qa + m) / max(rhs, 1e-300) - 1.0
    pen = max(0.0, -s) / T + max(0.0, 1e-6 - viol_rel)
    if clear:
        pR = q - tR*nR; pL = q - tL*nL
        pen += max(0.0, -sep(np.array([pR, q]), seg_axes(pR, q), VL, [uL, nL])) / T
        pen += max(0.0, -sep(np.array([pL, q]), seg_axes(pL, q), VR, [uR, nR])) / T
    return T + PEN*pen

def run():
    out = {}
    mode = sys.argv[1]
    if mode == 'two':
        res = []
        for theta in [0.01, 0.1, 0.3, 0.6, 0.9, 1.2, 1.4, 1.5]:
            for case in 'AB':
                for clear in (False, True):
                    if time.time() - T0 > TLIM: break
                    r = differential_evolution(two_sources_obj, [(-0.5, 0.5), (-0.5, 0.5), (0, 1.2), (0, 1.2)],
                                               args=(theta, case, clear), seed=1, maxiter=600, popsize=40, tol=1e-12, polish=True)
                    res.append(dict(theta=theta, case=case, clear=clear, Tmin=float(r.fun), x=[float(v) for v in r.x],
                                    half_cos=0.5*math.cos(theta)))
                    print(res[-1], flush=True)
        out['two_sources'] = res
    elif mode == 'three':
        res = []
        for alpha in [0.01, 0.1, 0.3, 0.5, 0.6, 0.7, 0.75, 0.785398]:
            for clear in (False, True):
                if time.time() - T0 > TLIM: break
                r = differential_evolution(three_sources_obj, [(0, 1.5)]*3 + [(-0.5, 0.5)]*3 + [(0, 1), (0, 1)],
                                           args=(alpha, clear), seed=2, maxiter=800, popsize=40, tol=1e-12, polish=True)
                res.append(dict(alpha=alpha, clear=clear, Tmin=float(r.fun), x=[float(v) for v in r.x],
                                half_cos2a=0.5*math.cos(2*alpha), cos2a=math.cos(2*alpha)))
                print(res[-1], flush=True)
        out['three_sources'] = res
    elif mode == 'star':
        res = []
        for theta in ([float(v) for v in sys.argv[3].split(',')] if len(sys.argv) > 3 else [1e-6, 1e-3, 0.05, 0.3, 0.8, 1.2]):
            for delta in [1e-5, 1e-2, 0.1, 0.3, 0.45]:
                if time.time() - T0 > TLIM: break
                r = differential_evolution(star_obj, [(-0.5, 0.5), (-0.5, 0.5), (1e-9, 1), (1e-9, 1)],
                                           args=(theta, delta), seed=3, maxiter=200, popsize=20, tol=1e-14, polish=True)
                res.append(dict(theta=theta, delta=delta, max_excess_ratio=float(-r.fun), x=[float(v) for v in r.x],
                                half_cos=0.5*math.cos(theta)))
                print(res[-1], flush=True)
        out['star'] = res
    elif mode == 'starmin':
        res = []
        for theta in ([float(v) for v in sys.argv[3].split(',')] if len(sys.argv) > 3 else [0.01, 0.1, 0.3, 0.6, 0.9, 1.2, 1.4]):
            for clear in (False, True):
                if time.time() - T0 > TLIM: break
                r = differential_evolution(starmin_obj, [(-0.5, 0.5), (-0.5, 0.5), (0, 1.2), (0, 1.2)],
                                           args=(theta, clear), seed=4, maxiter=250, popsize=25, tol=1e-12, polish=True)
                res.append(dict(theta=theta, clear=clear, Tmin=float(r.fun), x=[float(v) for v in r.x], half_cos=0.5*math.cos(theta)))
                print(res[-1], flush=True)
        out['starmin'] = res
    json.dump(out, open(out_path(f'geom_{mode}_' + (sys.argv[3].replace(',', '_') if len(sys.argv) > 3 else 'all') + '.json'), 'w'), indent=1)

if __name__ == '__main__':
    run()
