#!/usr/bin/env python3
"""
Invariant tests for the stellar tables under interpolation.

Each test_* function returns (passed, message). Run as a script to get a
PASS/FAIL summary. Designed to catch table bugs that only show up away from
grid points — i.e., at the runtime interpolation locations the FB code uses.
"""

import sys
import numpy as np
from interp import StellarTable, load_default
from interp import REM_NAMES, SN_TYPES, NO_EJECTA_TYPES, REM_WD

FLOAT_TOL    = 1e-12         # for grid-point exactness
INTERP_TOL   = 1e-9          # for closure under interpolation
PHYS_TOL     = 1e-6          # for physics invariants (some interp rounding ok)

N_LM_FINE = 200
N_LZ_FINE = 50
N_AGE_FINE = 100


# ── helpers ────────────────────────────────────────────────────────────

def _fine_grid(T, n_lm=N_LM_FINE, n_lz=N_LZ_FINE):
    lm = np.linspace(T.logM[0],  T.logM[-1],  n_lm)
    lz = np.linspace(T.logZ[0],  T.logZ[-1],  n_lz)
    return lm, lz


def _fine_age(T, n=N_AGE_FINE):
    return np.linspace(T.log_age[0], T.log_age[-1], n)


# ── invariant 1: closure  M_preSN - rem - Σ elem_SN ≈ 0  for SN types ──

def test_closure_SN_interpolated(T):
    lm_fine, lz_fine = _fine_grid(T)
    sum_SN_field = T.elem_ej_SN.sum(axis=2)
    worst = (0.0, None)
    n_checked = 0
    for lm in lm_fine:
        for lz in lz_fine:
            rt = T.remnant_type(lm, lz)
            if rt not in SN_TYPES:
                continue
            mp = T.m_preSN(lm, lz)
            rm = T.remnant_mass(lm, lz)
            ss = T._bilinear_rt_aware(sum_SN_field, lm, lz)
            r = abs(mp - rm - ss)
            n_checked += 1
            if r > worst[0]:
                worst = (r, (10**lm, 10**lz, REM_NAMES[rt], mp, rm, ss))
    if worst[0] < INTERP_TOL:
        return True, f"max|M_preSN-rem-ΣSN| = {worst[0]:.2e} over {n_checked} interp pts"
    return False, (f"max|M_preSN-rem-ΣSN| = {worst[0]:.4e} > {INTERP_TOL:.0e}\n"
                   f"    worst at M={worst[1][0]:.2f} Z={worst[1][1]:.5f} "
                   f"type={worst[1][2]} M_preSN={worst[1][3]:.3f} "
                   f"rem={worst[1][4]:.3f} ΣSN={worst[1][5]:.3f}")


# ── invariant 2: FSN/DBH have rem == M_preSN, Σ SN == 0 ────────────────

def test_FSN_rem_equals_M_preSN(T):
    # check grid points (interpolation may average across boundaries)
    worst = 0.0
    n = 0
    for iz in range(T.NZ):
        for im in range(T.NM):
            if int(T.rem_type[iz, im]) not in NO_EJECTA_TYPES:
                continue
            d = abs(T.M_preSN[iz, im] - T.rem_mass[iz, im])
            worst = max(worst, d)
            n += 1
    if n == 0:
        return True, "no FSN/DBH grid points to check"
    return (worst < FLOAT_TOL,
            f"max|rem-M_preSN| over {n} FSN/DBH grid pts = {worst:.2e}")


def test_FSN_yields_zero(T):
    worst = 0.0
    n = 0
    for iz in range(T.NZ):
        for im in range(T.NM):
            if int(T.rem_type[iz, im]) not in NO_EJECTA_TYPES:
                continue
            s = float(T.elem_ej_SN[iz, im].sum())
            worst = max(worst, abs(s))
            n += 1
    return worst < FLOAT_TOL, f"max ΣSN over {n} FSN/DBH grid pts = {worst:.2e}"


# ── invariant 3: closure for WD via AGB  M_preSN - rem - Σ AGB ≈ 0 ─────

def test_closure_AGB_interpolated(T):
    lm_fine, lz_fine = _fine_grid(T)
    sum_AGB_field = T.elem_ej_AGB.sum(axis=2)
    worst = 0.0
    n_checked = 0
    for lm in lm_fine:
        for lz in lz_fine:
            rt = T.remnant_type(lm, lz)
            if rt != REM_WD:
                continue
            mp = T.m_preSN(lm, lz)
            rm = T.remnant_mass(lm, lz)
            ss = T._bilinear_rt_aware(sum_AGB_field, lm, lz)
            worst = max(worst, abs(mp - rm - ss))
            n_checked += 1
    return (worst < INTERP_TOL,
            f"max|M_preSN-rem-ΣAGB|={worst:.2e} over {n_checked} WD interp pts")


# ── invariant 4: grid-point round-trip ─────────────────────────────────

def test_grid_point_exactness(T):
    worst = 0.0
    for iz in range(T.NZ):
        for im in range(T.NM):
            for name, field, expected in (
                ("rem_mass", T.rem_mass, T.rem_mass[iz, im]),
                ("M_preSN",  T.M_preSN,  T.M_preSN[iz, im]),
                ("lifetime", T.lifetime, T.lifetime[iz, im]),
            ):
                v = T._bilinear(field, T.logM[im], T.logZ[iz])
                worst = max(worst, abs(v - expected))
    return worst < FLOAT_TOL, f"max grid round-trip error = {worst:.2e}"


# ── invariant 5: Σ cum_wind monotonic in log_age ────────────────────────

def test_cum_wind_monotonic_total(T):
    if T.cum_wind is None:
        return True, "cum_wind data not present — skipped"
    sum_cw = T.cum_wind.sum(axis=3)  # (NZ, NM, NAGE)
    diffs = np.diff(sum_cw, axis=2)
    bad_count = int(np.sum(diffs < -1e-10))
    if bad_count == 0:
        return True, "Σ cum_wind is non-decreasing in age (all grid samples)"
    min_diff = float(diffs.min())
    return False, f"{bad_count} samples with Σ cum_wind decreasing; min Δ = {min_diff:.2e}"


# ── invariant 6: cum_wind at lifetime = M_init - M_preSN ──────────────

def test_cum_wind_endpoint(T):
    if T.cum_wind is None:
        return True, "cum_wind data not present — skipped"
    sum_cw_end = T.cum_wind[:, :, -1, :].sum(axis=2)  # (NZ, NM)
    expected   = T.M_init[None, :] - T.M_preSN
    diff = sum_cw_end - expected
    worst = float(np.max(np.abs(diff)))
    return (worst < PHYS_TOL,
            f"max |Σcum_wind(t_end) - (M_init-M_preSN)| = {worst:.2e}")


# ── invariant 7: per-element cum_wind monotonic ────────────────────────

def test_cum_wind_per_elem_monotonic(T):
    if T.cum_wind is None:
        return True, "cum_wind data not present — skipped"
    diffs = np.diff(T.cum_wind, axis=2)  # (NZ, NM, NAGE-1, NELEM)
    bad = int(np.sum(diffs < -1e-10))
    if bad == 0:
        return True, "per-element cum_wind non-decreasing"
    min_diff = float(diffs.min())
    return False, f"{bad} per-element samples decreasing; min Δ = {min_diff:.2e}"


# ── invariant 8: M_current + Σ cum_wind == M_init  (per (Z,M,age)) ─────

def test_M_current_plus_wind_ledger(T):
    if T.cum_wind is None:
        return True, "cum_wind data not present — skipped"
    # Restrict to ages ≤ lifetime: post-death M_current is sentinel-zero, so the
    # invariant only applies while the star is alive.
    log_lifetime = np.log10(np.maximum(T.lifetime, 1.0))   # (NZ, NM)
    age_3d = np.broadcast_to(T.log_age[None, None, :], T.M_current.shape)
    alive = (T.M_current > 0) & (age_3d <= log_lifetime[:, :, None])
    sum_cw = T.cum_wind.sum(axis=3)
    diff = (T.M_current + sum_cw - T.M_init[None, :, None])
    if not alive.any():
        return True, "no alive samples"
    worst = float(np.max(np.abs(diff[alive])))
    mean_ = float(np.mean(np.abs(diff[alive])))
    return (worst < PHYS_TOL,
            f"max |M_current + Σ cum_wind - M_init| (alive) = {worst:.2e}, mean = {mean_:.2e}")


# ── invariants 9-11: sanity bounds ─────────────────────────────────────

def test_elem_SN_nonneg(T):
    mn = float(T.elem_ej_SN.min())
    return mn >= 0.0, f"min elem_ej_SN = {mn:.2e}"


def test_rem_mass_bounds(T):
    if (T.rem_mass < -FLOAT_TOL).any():
        return False, f"rem_mass has negative entries: min={T.rem_mass.min():.2e}"
    over = T.rem_mass - T.M_preSN
    worst = float(over.max())
    return (worst < PHYS_TOL,
            f"max (rem - M_preSN) = {worst:.2e}  [should be ≤ 0]")


def test_M_current_bounds(T):
    if (T.M_current < -FLOAT_TOL).any():
        return False, f"M_current has negative entries"
    over = T.M_current - T.M_init[None, :, None]
    worst = float(over.max())
    return (worst < PHYS_TOL,
            f"max (M_current - M_init) = {worst:.2e}  [should be ≤ 0]")


# ── invariant 12: rem_type behavior at grid points ─────────────────────

def test_rem_type_integer(T):
    rts = set(int(x) for x in T.rem_type.flatten())
    bad = rts - set(REM_NAMES.keys())
    return len(bad) == 0, f"rem_type values seen: {sorted(rts)}"


# ── invariant 13: edge / extrapolation sanity ──────────────────────────

def test_edge_no_nan(T):
    for lm in (T.logM[0], T.logM[-1]):
        for lz in (T.logZ[0], T.logZ[-1]):
            for name, fn in (("rem_mass", T.remnant_mass),
                             ("M_preSN",  T.m_preSN),
                             ("lifetime", T.lifetime_yr)):
                v = fn(lm, lz)
                if not np.isfinite(v):
                    return False, f"{name}({10**lm:.2f},{10**lz:.5f}) = {v}"
    return True, "no NaN/inf at grid corners"


# ── driver ─────────────────────────────────────────────────────────────

TESTS = [
    ("1.  SN closure interp",          test_closure_SN_interpolated),
    ("2a. FSN/DBH rem == M_preSN",     test_FSN_rem_equals_M_preSN),
    ("2b. FSN/DBH Σ SN == 0",          test_FSN_yields_zero),
    ("3.  WD/AGB closure interp",      test_closure_AGB_interpolated),
    ("4.  grid-point round-trip",      test_grid_point_exactness),
    ("5.  Σ cum_wind monotonic",       test_cum_wind_monotonic_total),
    ("6.  cum_wind endpoint",          test_cum_wind_endpoint),
    ("7.  per-elem cum_wind mono",     test_cum_wind_per_elem_monotonic),
    ("8.  M_current + wind ledger",    test_M_current_plus_wind_ledger),
    ("9.  elem_ej_SN ≥ 0",             test_elem_SN_nonneg),
    ("10. 0 ≤ rem_mass ≤ M_preSN",     test_rem_mass_bounds),
    ("11. 0 ≤ M_current ≤ M_init",     test_M_current_bounds),
    ("12. rem_type integer-valued",    test_rem_type_integer),
    ("13. edge no NaN",                test_edge_no_nan),
]


def main():
    T = load_default()
    print(f"Loaded table: {T.NZ} Z x {T.NM} M x {T.NAGE} ages, {T.NELEM} elements")
    print(f"cum_wind data: {'present' if T.cum_wind is not None else 'MISSING'}\n")
    n_pass = n_fail = 0
    for label, fn in TESTS:
        try:
            ok, msg = fn(T)
        except Exception as e:
            ok, msg = False, f"EXCEPTION: {e!r}"
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {label}: {msg}")
        n_pass += ok
        n_fail += not ok
    print(f"\n{n_pass}/{len(TESTS)} passed, {n_fail} FAILED")
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
