#!/usr/bin/env python3
"""Radiation checks of the stellar table against analytic spectra and the raw PARSEC v2 files.

  test_radiation.py BUILD_DIR      exit 0 = all PASS (raw inputs from BUILD_DIR/data_root)

- band integration: a blackbody sampled like ATLAS9 gives the analytic Planck photon rates (band edges included)
- the table names the radiation source it was built with (no "Planck" label on a model-atmosphere table)
- time integrals: post-MS and lifetime integral of L_bol = the raw PARSEC v2 track's (age sampling)
- O stars: the hotstar_stack Q(H)/L and Q(HeI)/L at the Teff and log g of Martins+05 Tables 1-3 (CMFGEN with winds,
  Galactic) within 0.15 and 0.5 dex at every point, Q(H) within 0.05 dex on average
- PARSEC v2 photon tables (Costa+25, same tracks, WM-basic/PoWR/ATLAS9 atmospheres): lifetime-integrated Q(H) within
  0.3 dex at 20-100 Msun, Q(HeI) within 0.35 dex at 40-100 Msun; Q(HeI) at 20 Msun, Q(HeII) and LW printed only.
  At late O (28-34 kK) PARSEC's WM-basic rates are 0.25 dex (Q(H)) and 0.6 dex (Q(HeI)) below Martins+05, while
  TLUSTY, CMFGEN and PoWR agree with it to 0.05 dex (Q(H)); so 20 Msun at Z >= 0.01 sits +0.2-0.3 dex above PARSEC
"""
import os, re, sys
from glob import glob
import numpy as np, h5py
trapz = getattr(np, 'trapezoid', None) or np.trapz   # numpy >= 2 renamed it

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'stellar'))
TESTS = []
def test(fn): TESTS.append(fn); return fn
h, c, k = 6.62607015e-27, 2.99792458e10, 1.380649e-16
POINTS = [(0.001, 20.), (0.001, 40.), (0.001, 100.), (0.014, 20.), (0.014, 40.), (0.014, 100.)]


@test
def band_integration_planck(t, w, root):
    import build_spectral_lut as S
    from scipy.integrate import quad
    worst = (0., None)
    for T in (20000., 39000., 60000.):
        lam = np.concatenate([np.arange(90, 1000, 2.), np.arange(1000, 30000, 50.), np.arange(30000, 2e6, 2000.)])
        B = lambda l: 2 * h * c**2 / (l * 1e-8)**5 / np.expm1(h * c / (l * 1e-8 * k * T)) * 1e-8
        F = B(lam); _, qf = S.integrate_bands(lam, F); Ft = trapz(F, lam)
        for name, Elo, Ehi in S.PHOTON_RATE_BANDS:
            if name == 'Q_He1' and T < 30000: continue     # 1e-9 of Q_H: trapezoid curvature of the Wien tail
            lo, hi = (S.HC_eV_A / Ehi if Ehi else lam[0]), S.HC_eV_A / Elo
            qa = quad(lambda l: B(l) * l * 1e-8 / (h * c), lo, hi, limit=500)[0] / Ft
            e = abs(qf[name] / qa - 1)
            if e > worst[0]: worst = (e, (T, name))
    return worst[0] < 0.01, 'max |Q_band/Planck - 1| = %.2e at %s' % worst


@test
def radiation_source_labelled(t, w, root):
    opt, src = str(t.attrs.get('radiation_source_option', '')), str(t.attrs.get('radiation_source', ''))
    ok = opt in ('atlas9', 'planck', 'atlas9_nonionizing', 'hotstar_stack') and 'TODO' not in src and \
        (opt == 'planck') == src.startswith('Planck')
    return ok, 'radiation_source_option = %r, radiation_source = %r' % (opt, src[:90])


@test
def ostar_rates_vs_martins05(t, w, root):
    import build_stellar_tables as B
    B.RADIATION['source'] = 'hotstar_stack'
    d = np.genfromtxt(os.path.join(HERE, 'data', 'martins05_ostar_theoretical.txt'), dtype=None, encoding=None)
    T = np.array([r[2] for r in d], float); lg = np.array([r[3] for r in d], float)
    lT = np.log10(T); sT4 = 5.670374e-5 * T ** 4
    rad = B.compute_radiation_from_lut(lT, np.zeros(len(T)), lg, np.full(len(T), np.log10(B.Z_SOLAR)))
    L1 = B.L_sun                                    # log L = 0: rates per L_sun
    d0 = np.log10(rad['Q_ion'] / L1) - (np.array([r[4] for r in d], float) - np.log10(sT4))
    d1 = np.log10(rad['Q_He0'] / L1) - (np.array([r[5] for r in d], float) - np.log10(sT4))
    k0, k1 = int(np.argmax(abs(d0))), int(np.argmax(abs(d1)))
    ok = abs(d0).max() < 0.15 and abs(d1).max() < 0.5 and abs(d0.mean()) < 0.05
    return ok, 'stack - Martins+05 (%d O stars): dlog Q_H/L mean %+.3f, max %+.2f (%s %s); dlog Q_HeI/L mean %+.2f, max %+.2f (%s %s)' % (
        len(T), d0.mean(), d0[k0], d[k0][1], d[k0][0], d1.mean(), d1[k1], d[k1][1], d[k1][0])


def _raw_track(Z, M):
    d = glob(os.path.join(root_g, 'parsec2_vms', 'tracks_Z%s' % ('%g' % Z), '*_M%.1f.TAB' % M))
    if not d or os.path.getsize(d[0]) == 0: return None
    L = open(d[0]).read().splitlines(); n = L[2].split()
    a = np.array([[float(x) for x in l.split()] for l in L[3:] if len(l.split()) == len(n)])
    return {k_: a[:, n.index(k_)] for k_ in ('AGE', 'LOG_L', 'XCEN')}


def _table_integral(t, key, iz, im, a0=0.0):
    """integral of 10^key over table age (linear in log value vs log age, as GIZMO interpolates), from table age a0"""
    la, v = t['log_age_yr'][()], t[key][iz, im].astype(float)
    alive = v > -90
    if not alive.any(): return 0.0
    la, v = la[alive], v[alive]
    x = np.linspace(la[0], la[-1], 20000); y = 10 ** np.interp(x, la, v); age = 10 ** x
    m = (age >= a0) & (age <= t['lifetime_yr'][iz, im])   # GIZMO ends the star at its lifetime
    return trapz(y[m], age[m])


@test
def lbol_time_integral_vs_raw_track(t, w, root):
    Zt, Mt = t['Z'][()], t['M_init'][()]
    rows, worst = [], 0.
    for Z, M in POINTS + [(0.014, 5.)]:
        iz, im = int(np.argmin(abs(Zt - Z))), int(np.argmin(abs(Mt - M)))
        tr = _raw_track(Z, M)
        if tr is None: continue
        i0 = int(np.argmax(tr['XCEN'] < tr['XCEN'][0] - 1e-4)); age = tr['AGE'] - tr['AGE'][i0]   # ZAMS as the builder
        Lr = 10 ** tr['LOG_L'] * 3.828e33
        ms_end = age[int(np.argmax(tr['XCEN'] < 1e-4))]
        raw_ms, raw_post = trapz(Lr[(age >= 0) & (age <= ms_end)], age[(age >= 0) & (age <= ms_end)]), trapz(Lr[age >= ms_end], age[age >= ms_end])
        tab_all, tab_post = _table_integral(t, 'log_L_bol', iz, im), _table_integral(t, 'log_L_bol', iz, im, ms_end)
        r_post, r_all = tab_post / raw_post, tab_all / (raw_ms + raw_post)
        worst = max(worst, abs(np.log10(r_post)), abs(np.log10(r_all)))
        rows.append('Z=%g M=%g post-MS %.2f life %.2f' % (Z, M, r_post, r_all))
    return bool(rows) and worst < np.log10(1.15), 'table/raw integral of L_bol: ' + '; '.join(rows)


@test
def photons_vs_parsec2_tables(t, w, root):
    Zt, Mt = t['Z'][()], t['M_init'][()]
    rows, worst = [], {'H': 0., 'HeI': 0.}
    for Z, M in POINTS:
        f = glob(os.path.join(root, 'parsec2_vms', 'photons', '*', 'Z%gY*_M%.1f.TAB.QH' % (Z, M)))
        if not f or os.path.getsize(f[0]) == 0: continue
        a = np.loadtxt(f[0], comments='#', skiprows=3)
        age = a[:, 0] * 1e6
        iz, im = int(np.argmin(abs(Zt - Z))), int(np.argmin(abs(Mt - M)))
        par = {'H': trapz(10 ** a[:, 5], age), 'HeI': trapz(10 ** a[:, 6], age), 'HeII': trapz(10 ** a[:, 7], age),
               'LW': trapz(np.maximum(10 ** a[:, 9] - 10 ** a[:, 5], 0), age)}
        tab = {'H': _table_integral(t, 'log_Q_ion', iz, im), 'HeI': _table_integral(t, 'log_Q_He0', iz, im),
               'HeII': _table_integral(t, 'log_Q_He1', iz, im),
               'LW': _table_integral(t, 'log_L_LW', iz, im) / (12.4 * 1.602e-12)}   # LW photons at the band mean energy
        d = {q: np.log10(tab[q] / par[q]) if tab[q] > 0 and par[q] > 0 else np.nan for q in par}
        worst['H'] = max(worst['H'], abs(d['H']))
        if M >= 40: worst['HeI'] = max(worst['HeI'], abs(d['HeI']))   # late-O Q(HeI): atmosphere families differ by ~1 dex
        rows.append('Z=%g M=%g dlog Q_H %+.2f Q_HeI %+.2f Q_HeII %+.2f LW %+.2f' % (Z, M, d['H'], d['HeI'], d['HeII'], d['LW']))
    ok = bool(rows) and worst['H'] < 0.3 and worst['HeI'] < 0.35
    return ok, 'table - PARSEC v2 (lifetime integrals): ' + '; '.join(rows)


def main():
    global root_g
    b = sys.argv[1]; root = root_g = os.path.join(b, 'data_root')
    with h5py.File(os.path.join(b, 'stellar_tables_unified.hdf5'), 'r') as t, \
         h5py.File(os.path.join(b, 'stellar_wind_cumulative.hdf5'), 'r') as w:
        nf = 0
        for fn in TESTS:
            try: ok, msg = fn(t, w, root)
            except Exception as e: ok, msg = False, 'EXCEPTION %r' % e
            nf += not ok
            print('  [%s] %s: %s' % ('PASS' if ok else 'FAIL', fn.__name__, msg))
    print('\n%d/%d failed' % (nf, len(TESTS))); sys.exit(1 if nf else 0)


if __name__ == '__main__':
    main()
