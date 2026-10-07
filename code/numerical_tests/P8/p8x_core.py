# -*- coding: utf-8 -*-
"""
p8x_core.py - independent re-implementation of the constants of Sections 10-12 of the manuscript and of
the explicit bound Psi_pi(k) of (10.3).  Written from scratch for this test (2026-10-06); it does not
import or run any other code (in particular not the script k14_constants_check.py of the manuscript).

Formulas are transcribed from the definitions only:
  Definition 3.1 (alpha, beta-hat, b*, w0, h0), Definition 3.2 (constraints (C0)-(C7), (Q1), (Q2)),
  Definition 3.6 (A, A', C_Lambda, invF, Gamma, Omega_W, K_W, B),
  Theorem 10.3 ((10.1) lower bound for l(H), (10.2) upper bound for Omega_W, (10.3) Psi_pi(k)),
  Proposition 11.2 (D1 and the derivative of Psi_pi(k) - c* k^{1/4}).
Every function takes a context `ctx` which is either mpmath.mp (50 digits) or
mpmath.iv (50-digit interval arithmetic, outward rounding).
Variants of the fixed constants: 'exact9' / 'exact13' = A, A', C_Lambda with K = 9 / K = 13 (Table 1);
'round9' / 'round13' = the upper bounds A-bar, A'-bar, C_Lambda-bar used in Sections 11.1 and 11.5.
"""
from mpmath import mp, iv, mpf

mp.dps = 50
iv.dps = 50


def num(ctx, s):
    """decimal string -> number of the context (interval encloses the decimal)."""
    return ctx.mpf(s)


def pw(ctx, x, p, q):
    """x**(p/q) for x>0, computed as exp((p/q) log x) (q a power of 2 -> exact exponent)."""
    return ctx.exp(ctx.log(x) * ctx.mpf(p) / q)


def lo(x):
    """lower end (mp.mpf) of a number or interval"""
    if hasattr(x, '_mpi_'):
        return mp.make_mpf(x._mpi_[0])
    return x


def hi(x):
    if hasattr(x, '_mpi_'):
        return mp.make_mpf(x._mpi_[1])
    return x


def constants(ctx, variant):
    K = {}
    K['delta'] = num(ctx, '1e-5')
    K['betabar'] = num(ctx, '1.5e-6')
    K['bW'] = num(ctx, '2.0002')
    Q = num(ctx, '0.32')
    if variant in ('exact9', 'exact13'):
        Kw = 9 if variant == 'exact9' else 13
        A = 4 * ctx.sqrt(ctx.mpf(Kw) / Q) / (2 - num(ctx, '3e-6'))
        Ap = A / ctx.cos(K['betabar']) + num(ctx, '1e-8')
        CL = 1 + (num(ctx, '1185.6') if Kw == 9 else num(ctx, '1712.6')) * K['delta'] ** 2
    elif variant == 'round9':
        A, Ap, CL = num(ctx, '10.6067'), num(ctx, '10.6068'), num(ctx, '1.0000002')
    elif variant == 'round13':
        A, Ap, CL = num(ctx, '12.7476'), num(ctx, '12.7476'), num(ctx, '1.0000002')
    else:
        raise ValueError(variant)
    K['A'], K['Ap'], K['CL'] = A, Ap, CL
    K['variant'] = variant
    return K


def params(ctx, c0, c1, y0, eps, om0):
    return dict(c0=num(ctx, c0), c1=num(ctx, c1), y0=num(ctx, y0),
                eps=num(ctx, eps), om0=num(ctx, om0),
                raw=(c0, c1, y0, eps, om0))


def derived(ctx, P, K):
    c0, c1, y0, eps, om0 = P['c0'], P['c1'], P['y0'], P['eps'], P['om0']
    delta, bW, A, Ap, CL = K['delta'], K['bW'], K['A'], K['Ap'], K['CL']
    D = {}
    D['alpha0'] = c0 / ctx.sqrt(y0)                       # alpha(y0) = c0 y0^{-1/2}
    D['beta0'] = c1 * pw(ctx, y0, -3, 4)                   # beta(y0)  = c1 y0^{-3/4}
    D['beta_eps0'] = c1 * pw(ctx, (1 - eps) * y0, -3, 4)   # beta((1-eps) y0)
    # Definition 3.1: w0 = 0.50001 c0^2 (1+1.0001/y0) + delta + 1.0001 c1 ((1-eps) y0)^{-3/4} + 1e-12
    D['w0'] = (num(ctx, '0.50001') * c0 ** 2 * (1 + num(ctx, '1.0001') / y0) + delta
               + num(ctx, '1.0001') * D['beta_eps0'] + num(ctx, '1e-12'))
    D['h0'] = 1 - 2 * D['w0']
    e34 = pw(ctx, 1 - eps, 3, 4)
    # Definition 3.6: invF = 4 (1-eps)^{3/4} (1-(1-eps)^{3/4}) / (3 eps (1 + bW/(eps y0)))
    D['invF'] = 4 * e34 * (1 - e34) / (3 * eps * (1 + bW / (eps * y0)))
    D['epsy0'] = eps * y0
    # Gamma = 1.00002 (2 c1 A y0^{-1/2} + 4/5 c1^2 delta^{-1} y0^{-5/4})
    D['Gamma'] = num(ctx, '1.00002') * (2 * c1 * A / ctx.sqrt(y0)
                                        + ctx.mpf(4) / 5 * c1 ** 2 / delta * pw(ctx, y0, -5, 4))
    D['T1'] = 1 / (om0 * pw(ctx, y0, 3, 4))
    D['T2'] = A / (2 * c1)
    D['T3'] = 2 * c1 / (c0 * D['invF']) * (Ap + D['alpha0'] * CL / delta)
    D['B'] = D['T1'] + D['T2'] + D['T3'] + D['Gamma']
    D['Kcoef'] = 4 * c1 / (c0 * D['invF'])                  # C_wall, the coefficient of Omega_W in (star)
    D['KW'] = D['Kcoef'] * num(ctx, '1.000001') * c0 ** 2   # K_W of (10.3)
    D['P'] = D['h0'] * (1 - om0) * pw(ctx, 1 - eps, 1, 4) / D['B']
    D['rho'] = 8 * pw(ctx, ctx.mpf(2), -1, 4) * D['P']       # lim Psi/k^{1/4} = Psi_infty(pi) of Section 11.4
    return D


def at_k(ctx, P, K, D, k, cstar=None):
    """Psi_pi(k) of (10.3) and its pieces; k may be a string or number."""
    c0, y0, eps, om0 = P['c0'], P['y0'], P['eps'], P['om0']
    kk = num(ctx, k) if isinstance(k, str) else k
    u = kk / 2 - 3
    R = {}
    R['u'] = u
    R['lH'] = 8 * D['h0'] * (pw(ctx, (1 - eps) * u, 1, 4) - pw(ctx, y0 + 1, 1, 4))  # (10.1)
    R['OmUB'] = num(ctx, '1.000001') * c0 ** 2 * (ctx.log(u / y0) + 2 / y0)        # (10.2)
    R['wall'] = D['Kcoef'] * R['OmUB']
    R['Psi'] = ((1 - om0) * R['lH'] - R['wall']) / D['B']                            # (10.3)
    R['k14'] = pw(ctx, kk, 1, 4)
    R['ratio'] = R['Psi'] / R['k14']
    R['C7'] = (1 - eps) * u - (y0 + 1)        # >= 0 needed
    if cstar is not None:
        cs = num(ctx, cstar)
        R['target'] = cs * R['k14']
        R['D1'] = D['P'] - cs / (4 * pw(ctx, ctx.mpf(2), 3, 4))
        R['u14D1'] = pw(ctx, u, 1, 4) * R['D1']
        R['KW2B'] = D['KW'] / (2 * D['B'])
    return R


def psi_deriv(ctx, P, D, k, cstar):
    """analytic phi'(k) of the proof of Proposition 11.2, phi = Psi - c* k^{1/4}"""
    u = k / 2 - 3
    cs = num(ctx, cstar)
    return (D['P'] * pw(ctx, u, -3, 4) - D['KW'] / (2 * D['B'] * u)
            - cs / 4 * pw(ctx, 2 * u + 6, -3, 4))


def constraints(ctx, P, K, D, k=None):
    """returns list of (name, margin) where margin >= 0 (or > 0 for strict) means satisfied"""
    c0, c1, y0, eps, om0 = P['c0'], P['c1'], P['y0'], P['eps'], P['om0']
    out = []
    out.append(('C0 c0>0', c0))
    out.append(('C0 c1>0', c1))
    out.append(('C0 eps>0', eps))
    out.append(('C0 1-eps>0', 1 - eps))
    out.append(('C0 om0>0', om0))
    out.append(('C1 betabar-alpha(y0)>=0', K['betabar'] - D['alpha0']))
    out.append(('C2 1.5e-6-beta(y0)>=0', num(ctx, '1.5e-6') - D['beta0']))
    out.append(('C3 c0 y0^{1/4}-c1(1-eps)^{-3/4}>0', c0 * pw(ctx, y0, 1, 4) - c1 * pw(ctx, 1 - eps, -3, 4)))
    out.append(('C4 1/2-om0>0', ctx.mpf(1) / 2 - om0))
    out.append(('C5 1/2-w0>0', ctx.mpf(1) / 2 - D['w0']))
    out.append(('C6 cos(2alpha0)/2-delta>=0', ctx.cos(2 * D['alpha0']) / 2 - K['delta']))
    # further side conditions: X1 = (Q1) and X2 = (Q2) of Definition 3.2; X3 = the range of the inequality
    # sin a + 1 - cos a <= 1.0001 a used in Lemma 9.3(e) (Remark 9.4); X4 = (C5) restated; X5, X6 = Lemma 10.2
    out.append(('X1 1e-3-alpha(y0)>=0 (sec x-1<=0.50001x^2 range)', num(ctx, '1e-3') - D['alpha0']))
    out.append(('X2 1e-4-beta((1-eps)y0)>=0', num(ctx, '1e-4') - D['beta_eps0']))
    out.append(('X3 1.5e-6-beta((1-eps)y0)>=0 (range of sin a+1-cos a<=1.0001a)', num(ctx, '1.5e-6') - D['beta_eps0']))
    out.append(('X4 h0>0', D['h0']))
    out.append(('X5 invF>0', D['invF']))
    out.append(('X6 1-invF>=0', 1 - D['invF']))
    if k is not None:
        u = (num(ctx, k) if isinstance(k, str) else k) / 2 - 3
        out.append(('C7 (1-eps)y1-(y0+1)>=0', (1 - eps) * u - (y0 + 1)))
    return out
