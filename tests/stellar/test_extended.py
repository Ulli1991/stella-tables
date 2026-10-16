#!/usr/bin/env python3
"""Extended checks on the built stellar tables: catalogue stitch, scaled-solar mixture, spectral LUT,
band/bolometric consistency, ionizing-photon ordering, and dataset fingerprints against a reference build.

  test_extended.py BUILD_DIR [--write-fingerprints]     exit 0 = all PASS (expected failures do not fail the run)
"""
import os, sys, json, hashlib
import numpy as np, h5py

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
FP_FILE = os.path.join(HERE, 'reference_fingerprints.json')
TESTS = []
def test(xfail=None):
    def deco(fn): TESTS.append((fn, xfail)); return fn
    return deco


def fingerprints(table, wind):
    fp = {}
    for tag, path in (('table', table), ('wind', wind)):
        with h5py.File(path, 'r') as f:
            def v(n, o):
                if isinstance(o, h5py.Dataset): fp['%s:%s' % (tag, n)] = hashlib.sha256(np.ascontiguousarray(o[()]).tobytes()).hexdigest()
            f.visititems(v)
    return fp


def seam_excess(y, M, src, ok=None):
    """per track-source seam: |step - local slope x dlogM| (dex), local slope from the two neighbouring intervals"""
    lm = np.log10(M); worst = (0., None)
    for k in np.where(src[1:] != src[:-1])[0]:
        if k < 1 or k + 2 >= len(M) or (ok is not None and not ok[k - 1:k + 3].all()): continue
        slope = 0.5 * ((y[k] - y[k - 1]) / (lm[k] - lm[k - 1]) + (y[k + 2] - y[k + 1]) / (lm[k + 2] - lm[k + 1]))
        e = abs(y[k + 1] - y[k] - slope * (lm[k + 1] - lm[k]))
        if np.isfinite(e) and e > worst[0]: worst = (e, (M[k], M[k + 1]))
    return worst


def per_z(t, fn):
    worst = (0., None)
    for iz in range(t['Z'].shape[0]):
        e, at = fn(iz)
        if e > worst[0]: worst = (e, (t['Z'][iz],) + at)
    return worst


@test(xfail='open item: MIST->PARSEC v2 stitch at 2 Msun steps lifetime up to 0.06 dex (MIST nearest-[Fe/H], and MIST itself non-monotonic near 1.95 Msun at Z=0.006)')
def stitch_continuity_lifetime(t, w, lut):
    """log lifetime(M) continuous across track-source seams (lifetimes below a Hubble time)"""
    lt, src, M = np.log10(t['lifetime_yr'][()]), t['track_source'][()], t['M_init'][()]
    e, at = per_z(t, lambda iz: seam_excess(lt[iz], M, src[iz], lt[iz] < np.log10(1.4e10)))
    return e < 0.05, 'worst excess seam step %.3f dex at Z,M=%s' % (e, at)


def ms_logL(t):
    """log L_bol at 30 % of each star's lifetime (main sequence)"""
    L, lt, ages = t['log_L_bol'][()], t['lifetime_yr'][()], t['log_age_yr'][()]
    ia = np.clip(np.searchsorted(ages, np.log10(0.3 * lt)), 0, len(ages) - 1)
    return np.take_along_axis(L, ia[..., None], 2)[..., 0]


@test()
def stitch_continuity_ms_luminosity(t, w, lut):
    """main-sequence log L_bol(M) continuous across track-source seams"""
    L, src, M = ms_logL(t), t['track_source'][()], t['M_init'][()]
    e, at = per_z(t, lambda iz: seam_excess(L[iz], M, src[iz]))
    return e < 0.1, 'worst excess seam step %.3f dex at Z,M=%s' % (e, at)


@test(xfail='MIST and PARSEC v2 treat the pre-main sequence differently: L_bol at 1e5 yr jumps 0.1-0.3 dex at 2 Msun')
def stitch_continuity_pms_luminosity(t, w, lut):
    """pre-main-sequence (1e5 yr) log L_bol(M) continuous across track-source seams"""
    L, src, M = t['log_L_bol'][()], t['track_source'][()], t['M_init'][()]
    ia = int(np.searchsorted(t['log_age_yr'][()], 5.0))
    e, at = per_z(t, lambda iz: seam_excess(L[iz, :, ia], M, src[iz]))
    return e < 0.1, 'worst excess seam step %.3f dex at Z,M=%s' % (e, at)


@test()
def scaled_solar_metal_pattern(t, w, lut):
    """the metal mixture in X_init is the same scaled pattern at every Z"""
    X = t['X_init'][()]; met = X[:, 2:]; pat = met / met.sum(1, keepdims=True)
    dev = np.max(np.abs(pat / pat[len(pat) // 2] - 1)[:, pat[len(pat) // 2] > 1e-6])
    return dev < 1e-6, 'max relative pattern deviation across Z = %.2e' % dev


@test()
def surface_abundances_sum_to_one(t, w, lut):
    """surface mass fractions sum to 1 wherever set, so the wind mass ledger equals the track mass loss"""
    s = t['surface_abundances'][()].astype(np.float64).sum(-1); m = s > 0
    dev = np.max(np.abs(s[m] - 1)); return dev < 1e-6, 'max |sum_k surface - 1| = %.1e over %d alive entries' % (dev, m.sum())


@test()
def x_init_sums_to_one(t, w, lut):
    s = t['X_init'][()].sum(1); return np.max(np.abs(s - 1)) < 1e-6, 'max|sum X_init - 1| = %.2e' % np.max(np.abs(s - 1))


@test()
def lut_fractions_physical(t, w, lut):
    """band energy fractions finite, in [0,1], and summing to <= 1 wherever the LUT has data"""
    has = lut['has_data'][()]; tot = 0; bad = 0
    for k in lut['energy_fracs']:
        f = lut['energy_fracs'][k][()]
        bad += int(np.sum(~np.isfinite(f[has]) | (f[has] < -1e-12) | (f[has] > 1 + 1e-9))); tot = tot + np.where(has, f, 0)
    ion = sum(lut['energy_fracs'][k][()] for k in ('ion_H0', 'ion_He0', 'ion_He1', 'ion_He2'))
    over = np.max(np.where(has, tot - ion, 0))   # ion_H0..He2 are sub-bands of the ionizing band: count the non-ionizing bands + ion once
    return bad == 0 and over <= 1 + 1e-6, 'out-of-range fractions: %d, max non-ion band sum %.4f' % (bad, over)


@test()
def band_le_bolometric(t, w, lut):
    """no band luminosity exceeds L_bol"""
    Lb = t['log_L_bol'][()]; worst = -np.inf
    for k in ('log_L_FUV_total', 'log_L_NUV', 'log_L_OPT_NIR', 'log_L_ion_tot'):
        d = t[k][()] - Lb; m = np.isfinite(d); worst = max(worst, np.max(d[m]) if m.any() else -np.inf)
    return worst <= 1e-6, 'max log(L_band/L_bol) = %.3e' % worst


def q_drops(t, lo, hi):
    Q, M = t['log_Q_ion'][()], t['M_init'][()]; ia = int(np.searchsorted(t['log_age_yr'][()], 5.0)); sel = (M >= lo) & (M <= hi)
    out = []
    for iz in range(Q.shape[0]):
        q = Q[iz, sel, ia]; ms = M[sel]; d = np.diff(q)
        out += [(t['Z'][iz], ms[k], ms[k + 1], round(-d[k], 2)) for k in np.where(d < -0.05)[0]]
    return out


@test()
def ionizing_rate_rises_with_mass(t, w, lut):
    """ionizing photon rate at 1e5 yr increases with mass for 15-150 Msun at every Z"""
    b = q_drops(t, 15, 150); return not b, '%d decreasing steps (> 0.05 dex): %s' % (len(b), b[:4])


@test(xfail='PARSEC v2 very massive stars at Z >= 0.017: log Q falls 0.2-0.6 dex from 250 to 500 Msun (origin not verified)')
def ionizing_rate_rises_with_mass_vms(t, w, lut):
    """ionizing photon rate at 1e5 yr increases with mass for 150-500 Msun at every Z"""
    b = q_drops(t, 150, 500); return not b, '%d decreasing steps (> 0.05 dex): %s' % (len(b), b[:4])


@test()
def wind_grid_matches_table(t, w, lut):
    ok = np.array_equal(t['Z'][()], w['Z'][()]) and np.array_equal(t['M_init'][()], w['M_init'][()]) and np.array_equal(t['log_age_yr'][()], w['log_age_yr'][()])
    return ok, 'Z / M_init / log_age grids identical in table and wind file'


@test()
def fingerprints_match_reference(t, w, lut):
    if not os.path.exists(FP_FILE): return False, 'no reference_fingerprints.json'
    ref = json.load(open(FP_FILE)); cur = fingerprints(t.filename, w.filename)
    diff = sorted(k for k in set(ref) | set(cur) if ref.get(k) != cur.get(k))
    return not diff, ('all %d datasets identical to the reference build' % len(ref)) if not diff else 'differ: %s' % diff[:8]


def main():
    out = sys.argv[1]
    table, wind = os.path.join(out, 'stellar_tables_unified.hdf5'), os.path.join(out, 'stellar_wind_cumulative.hdf5')
    if '--write-fingerprints' in sys.argv:
        json.dump(fingerprints(table, wind), open(FP_FILE, 'w'), indent=1, sort_keys=True); print('wrote', FP_FILE); return
    lut = h5py.File(os.environ.get('STELLA_SPECTRAL_LUT', os.path.join(REPO, 'stellar', 'data', 'lut', 'spectral_band_lut.hdf5')), 'r')
    with h5py.File(table, 'r') as t, h5py.File(wind, 'r') as w:
        nf = 0
        for fn, xf in TESTS:
            try: ok, msg = fn(t, w, lut)
            except Exception as e: ok, msg = False, 'EXCEPTION %r' % e
            tag = ('PASS' if ok else 'FAIL') if not xf else ('XPASS' if ok else 'XFAIL')
            if tag == 'FAIL': nf += 1
            print('  [%s] %s: %s%s' % (tag, fn.__name__, msg, ('  (expected: %s)' % xf) if xf else ''))
    print('\n%d/%d failed' % (nf, len(TESTS))); sys.exit(1 if nf else 0)


if __name__ == '__main__':
    main()
