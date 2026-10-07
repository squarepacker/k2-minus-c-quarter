# -*- coding: utf-8 -*-
"""Shared helpers for the tests of Section 4 (P1 and P2).

* out_path(name): path of an output file in the folder out/ next to the scripts
  (created on first use).
* rss_mb(): memory self-check used by the drivers to abort above 700 MB
  (Windows: current working set; other systems: peak resident set size;
  -1 if it cannot be determined).
* Exact configuration builders: place, ok_with, place_min_gap, build_column.
"""
import os
import sys
from fractions import Fraction as Fr
from tracer import Sq, rot, disjoint, HALF

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "out")


def out_path(name):
    os.makedirs(OUT_DIR, exist_ok=True)
    return os.path.join(OUT_DIR, name)


# ------------------------------------------------------------ self memory check
_WIN = None


def _win_setup():
    import ctypes
    import ctypes.wintypes

    class PMC(ctypes.Structure):
        _fields_ = [("cb", ctypes.wintypes.DWORD), ("PageFaultCount", ctypes.wintypes.DWORD),
                    ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]

    k32 = ctypes.windll.kernel32
    k32.GetCurrentProcess.restype = ctypes.wintypes.HANDLE
    psapi = ctypes.windll.psapi
    psapi.GetProcessMemoryInfo.argtypes = [ctypes.wintypes.HANDLE, ctypes.POINTER(PMC), ctypes.wintypes.DWORD]
    psapi.GetProcessMemoryInfo.restype = ctypes.wintypes.BOOL
    return ctypes, PMC, k32, psapi


def rss_mb():
    global _WIN
    if sys.platform == "win32":
        try:
            if _WIN is None:
                _WIN = _win_setup()
            ctypes, PMC, k32, psapi = _WIN
            c = PMC()
            c.cb = ctypes.sizeof(PMC)
            psapi.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(c), c.cb)
            return c.WorkingSetSize / 2 ** 20
        except Exception:
            return -1.0
    try:
        import resource
        r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return r / 2 ** 20 if sys.platform == "darwin" else r / 1024.0
    except Exception:
        return -1.0


# ------------------------------------------------------------ builders
def place(p, d, gap, t, xi):
    q = (p[0] + gap * d[0], p[1] + gap * d[1])
    c, s = rot(t)
    cx = q[0] - s / 2 - xi * c
    cy = q[1] + c / 2 - xi * s
    return Sq(cx, cy, t)


def ok_with(S, existing, k):
    for v in S.verts:
        if not (0 <= v[0] <= k and 0 <= v[1] <= k):
            return False
    for E in existing:
        if not disjoint(S, E):
            return False
    return True


def place_min_gap(existing, p, d, t, xi, gbase, k, gcap):
    """smallest gap >= gbase (geometric ladder then bisection) making S valid"""
    gbase = Fr(gbase)
    S = place(p, d, gbase, t, xi)
    if ok_with(S, existing, k):
        return S, gbase
    lo, hi = gbase, None
    gtry = max(gbase * 2, Fr(1, 10 ** 12))
    while gtry <= gcap:
        S = place(p, d, gtry, t, xi)
        if ok_with(S, existing, k):
            hi = gtry
            break
        lo = gtry
        gtry = gtry * 2
    if hi is None:
        return None, None
    for _ in range(30):
        mid = (lo + hi) / 2
        mid = Fr(mid).limit_denominator(10 ** 15) if mid.denominator > 10 ** 15 else mid
        if not (lo < mid < hi):
            break
        S = place(p, d, mid, t, xi)
        if ok_with(S, existing, k):
            hi = mid
        else:
            lo = mid
    return place(p, d, hi, t, xi), hi


def build_column(x0, specs, gaps, k, existing=(), gcap=Fr(1, 2)):
    """specs: list of (t, xi); gaps: list of base gaps g_0..g_{m-1} (g_0 from floor).
    Returns (squares, final exit point, final direction, list of used gaps) or None."""
    sq = []
    p = (Fr(x0), Fr(0))
    d = (Fr(0), Fr(1))
    used = []
    for i, (t, xi) in enumerate(specs):
        S, g = place_min_gap(list(existing) + sq, p, d, t, xi, gaps[i], k, gcap)
        if S is None:
            return None
        sq.append(S)
        used.append(g)
        q = S.pt(Fr(xi), -HALF)
        nn = S.n()
        p = (q[0] + nn[0], q[1] + nn[1])
        d = nn
    return sq, p, d, used
