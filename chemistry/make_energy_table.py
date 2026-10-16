#!/usr/bin/env python3
"""CHEMCOOL energy -> (temperature, gamma) lookup table.

Specific energy per particle (in K) of an H/He/H2 gas:
    ekn = x * e_H2(T) + (1 - x) * 1.5 T,     x = abh2 / (1 + abhe - abh2)
x is the H2 share of all heavy particles, so a table in (ekn, x) is exact for any abhe.
gamma = 1 + 1/c_v,  c_v = d ekn / dT.

  make_energy_table.py --mode x    --dx 0.05 --out energy-table-x.txt     (abhe-independent, new)
  make_energy_table.py --mode abh2 --abhe 0.1 --out energy-table.txt     (legacy layout, abh2 = 0.1..0.5)

Input: eh2-table.txt (T, e_H2(T) per molecule in K, T = 1..50000). That table switches from a partition-function
energy (T < 20000) to analytic 3T (T > 20000) with a 45 K step; by default each side gets its own spline so the
step does not ring into gamma. --legacy-seam splines across it and reproduces the original GIZMO table byte for byte.
"""
import argparse
import numpy as np
from scipy.interpolate import CubicSpline

N_EKN = 50000       # ekn = 1..50000 per block, as calc_temp.F expects
T_SEAM = 20000.     # eh2 regime switch


class SeamSpline:
    """separate cubic splines below and above the eh2 seam"""
    def __init__(self, T, e):
        lo, hi = T < T_SEAM, T > T_SEAM
        self.lo, self.hi = CubicSpline(T[lo], e[lo]), CubicSpline(T[hi], e[hi])

    def __call__(self, T, nu=0):
        T = np.asarray(T, dtype=float)
        return np.where(T < T_SEAM, self.lo(T, nu), self.hi(T, nu))


def load_eh2(path, legacy_seam=False):
    d = np.loadtxt(path)
    bad = (d[:, 0] == 20000) & (abs(d[:, 1] - 1.5 * d[:, 0]) < 1e-6)   # known glitch row: 1.5T instead of ~3T
    if bad.any(): print(f'eh2: dropping glitch row(s) at T={d[bad, 0]}')
    d = d[~bad]
    return CubicSpline(d[:, 0], d[:, 1]) if legacy_seam else SeamSpline(d[:, 0], d[:, 1])


def block(eh2, x, ekn):
    """T(ekn) and gamma(ekn) at fixed H2 particle share x"""
    T = np.concatenate([np.linspace(0.05, 1, 200, endpoint=False), np.arange(1, 50001, 0.25)])
    T = T[T <= 50000]
    eT = np.where(T >= 1, eh2(np.clip(T, 1, 5e4)), 1.5 * T)   # below 1 K: H2 internal dof frozen
    dT = np.where(T >= 1, eh2(np.clip(T, 1, 5e4), 1), 1.5)
    f = x * eT + (1 - x) * 1.5 * T
    cv = x * dT + (1 - x) * 1.5
    Tk = np.interp(ekn, f, T)
    gk = 1 + 1 / np.interp(ekn, f, cv)
    return Tk, gk


def make(eh2, xs):
    ekn = np.arange(1, N_EKN + 1, dtype=float)
    rows = []
    for x in xs:
        T, g = block(eh2, x, ekn)
        rows.append(np.column_stack([ekn, T, g]))
    return np.vstack(rows)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--eh2', default='eh2-table.txt')
    ap.add_argument('--mode', choices=['x', 'abh2'], default='x')
    ap.add_argument('--dx', type=float, default=0.05)
    ap.add_argument('--abhe', type=float, default=0.1)
    ap.add_argument('--legacy-seam', action='store_true', help='spline across the 20000 K seam (original GIZMO table)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    eh2 = load_eh2(a.eh2, a.legacy_seam)
    if a.mode == 'x':
        n = int(round(1 / a.dx))
        xs = np.linspace(0, 1, n + 1)
        hdr = f'energy table in (ekn, x), x = abh2/(1+abhe-abh2); nblocks={n + 1} dx={1 / n:.10g} nekn={N_EKN}'
    else:
        abh2 = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        xs = abh2 / (1 + a.abhe - abh2)
        hdr = f'legacy energy table in (ekn, abh2=0.1..0.5) for abhe={a.abhe}'
    tab = make(eh2, xs)
    np.savetxt(a.out, tab, fmt=['%d', '%.14g', '%.14g'], delimiter='  ')
    print(f'{hdr}: wrote {len(tab)} rows to {a.out}')
