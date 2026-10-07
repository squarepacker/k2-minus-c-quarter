# -*- coding: utf-8 -*-
"""Exact check of Lemma 4.7 (reflection), parts (a)-(c), on random rational squares:
the vertices of rho Z, v_{rho Z} = k - (v_Z + cos a + sin a), closure of R(rho Z) = k - closure of R(Z),
dist(k - y, Z) = dist(y, Z), and the symmetry of the window {y} in [w0, 1 - w0].

Usage:  python reflection.py          (writes out/reflection_out.json; about 1 s)
"""
import random, json, time
from fractions import Fraction as Fr
from tracer import Sq, HALF, dist_int
from common import out_path

T0 = time.time()
r = random.Random(7)
bad = []
N = 3000
for i in range(N):
    k = r.randint(4, 60)
    t = Fr(r.randint(-414, 414), 1000)
    c0 = Fr(r.randint(800, k * 1000 - 800), 1000)
    c1 = Fr(r.randint(800, k * 1000 - 800), 1000)
    Z = Sq(c0, c1, t)
    R = Z.reflect(k)
    # vertex sets: reflection of Z's vertices == R's vertices
    vz = sorted((p[0], k - p[1]) for p in Z.verts)
    vr = sorted(R.verts)
    if vz != vr:
        bad.append(("verts", i))
    if R.v != k - (Z.v + Z.ca + Z.sa):
        bad.append(("v", i))
    rz = sorted((k - hi, k - lo) for (lo, hi) in Z.ramps_closure())
    rr = sorted(R.ramps_closure())
    if t != 0 and rz != rr:
        bad.append(("ramps", i))
    for (lo, hi) in Z.ramps_closure():
        for y in (lo, hi, (lo + hi) / 2):
            if dist_int(k - y) != dist_int(y):
                bad.append(("dist", i))
    # frac in [w0, 1-w0] symmetric
    w0 = Fr(r.randint(1, 499), 1000)
    for y in (Fr(r.randint(0, 10 ** 6), 10 ** 5) for _ in range(3)):
        fy = y - (y.numerator // y.denominator)
        fky = (k - y) - ((k - y).numerator // (k - y).denominator)
        if (w0 <= fy <= 1 - w0) != (w0 <= fky <= 1 - w0):
            bad.append(("frac", i))
out = dict(n=N, bad=bad[:20], n_bad=len(bad), elapsed=time.time() - T0)
with open(out_path("reflection_out.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, default=str)
print(out)
