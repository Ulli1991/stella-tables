#!/usr/bin/env python3
"""Physics, schema and provenance tests for stellar_tables_unified.hdf5.

  test_table_physics.py [--table T.hdf5]                     run all checks, exit 0 = all PASS
  test_table_physics.py --compare A.hdf5 B.hdf5              dataset-by-dataset reproducibility diff

Complements test_invariants.py (mass/wind closure); this file covers what the
consumer (GIZMO resolvedism_stellar_tables.cc) and the physics rely on.
"""
import argparse, os, sys
import numpy as np, h5py

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'build')
ELEMS = ['H', 'He', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
         'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
IA_EJECTA_MASS = 1.378          # GIZMO resolvedism_fb_shared.h
# datasets GIZMO's loader opens by name (grep of resolvedism_stellar_tables.cc)
LOADER_DATASETS = ['M_CO_core', 'M_He_core', 'M_current', 'M_init', 'M_preSN', 'X_init', 'Z',
                   'elem_ej_AGB_mass', 'elem_ej_SN_mass', 'lifetime_yr', 'logR_cm', 'log_L_FUV',
                   'log_L_FUV_total', 'log_L_LW', 'log_L_NUV', 'log_L_NUV_lo', 'log_L_OPT_NIR',
                   'log_L_bol', 'log_L_ion_tot', 'log_M_init', 'log_Q_ion', 'log_Z', 'log_age_yr',
                   'phase_transitions/t_AGB_end', 'phase_transitions/t_AGB_start',
                   'phase_transitions/t_CHeB_start', 'phase_transitions/t_MS_end',
                   'phase_transitions/t_MS_start', 'phase_transitions/t_postAGB_start',
                   'phase_transitions/t_RGB_start', 'phase_transitions/t_WR_start',
                   'remnant_mass', 'remnant_type', 'surface_abundances', 't_PMS',
                   'type_ia_yields', 'v_wind']


class T:
    def __init__(self, path):
        self.f = h5py.File(path, 'r'); f = self.f
        self.Z, self.M = f['Z'][()], f['M_init'][()]
        self.nZ, self.nM, self.nA = len(self.Z), len(self.M), len(f['log_age_yr'])
        self.iZsun = int(np.argmin(abs(self.Z - 0.014)))
    def __getitem__(self, k): return self.f[k][()]
    def im(self, m): return int(np.argmin(abs(self.M - m)))


TESTS = []
def test(fn): TESTS.append(fn); return fn
XFAIL = {   # known open table issues: reported, not counted as failures
    'lifetime_monotonic_in_mass': 'MIST v1.2 tracks themselves: lifetime rises ~4% at one mass step near 1.9 Msun (1.88->1.90 at [Fe/H]=-0.5, 1.94->1.96 at -0.25)',
}


@test
def schema_loader_datasets(t):
    """every dataset the GIZMO loader reads exists with the expected shape"""
    missing = [d for d in LOADER_DATASETS if d not in t.f]
    want = {1: None, 2: (t.nZ, t.nM), 3: (t.nZ, t.nM, t.nA)}
    bad = []
    for d in LOADER_DATASETS:
        if d in t.f and t.f[d].ndim == 3 and d not in ('elem_ej_AGB_mass', 'elem_ej_SN_mass') \
                and t.f[d].shape != want[3]: bad.append((d, t.f[d].shape))
    return not missing and not bad, f'missing={missing} badshape={bad}'

@test
def schema_element_order(t):
    """27 elements in the H..Zn order GIZMO's METALS_INDIVIDUAL indices assume"""
    e = [x.decode().strip() for x in t['elements']]
    return e == ELEMS, f'{e}'

@test
def finite_everywhere(t):
    """no NaN/inf in any float dataset (log fields may be -inf only where documented: none)"""
    bad = []
    def v(n, o):
        if isinstance(o, h5py.Dataset) and o.dtype.kind == 'f' and not np.all(np.isfinite(o[()])):
            bad.append((n, int((~np.isfinite(o[()])).sum())))
    t.f.visititems(v)
    return not bad, f'nonfinite={bad}'

@test
def ia_yields_w7(t):
    """Type Ia: H=He=0, sum = IA_EJECTA_MASS, Fe=0.744 and Ni~0.14 (Iwamoto+99 W7)"""
    y = t['type_ia_yields']; d = dict(zip(ELEMS, y))
    ok = (d['H'] == 0 and d['He'] == 0 and abs(y.sum() - IA_EJECTA_MASS) < 1e-6
          and abs(d['Fe'] - 0.744) < 1e-6 and 0.12 < d['Ni'] < 0.16 and np.all(y >= 0))
    return ok, f"sum={y.sum():.6f} Fe={d['Fe']} Ni={d['Ni']:.4f} Si={d['Si']} S={d['S']}"

@test
def lifetime_monotonic_in_mass(t):
    """lifetime (t_PMS+lifetime) decreases with mass at every Z (above 0.9 Msun, alive tracks)"""
    life = t['lifetime_yr'] + t['t_PMS']; viol = []
    sel = t.M >= 0.9
    for iz in range(t.nZ):
        l = life[iz, sel]; ok = l > 0
        d = np.diff(l[ok])
        if np.any(d > 0.02 * l[ok][1:]): viol.append((float(t.Z[iz]), int((d > 0.02 * l[ok][1:]).sum())))
    return not viol, f'Z with >2% upticks: {viol}'

@test
def lifetime_anchors_solar(t):
    """solar-Z lifetimes: 1 Msun 8-13 Gyr, 8 Msun 30-50 Myr, 40 Msun 3.5-6 Myr"""
    life = t['lifetime_yr'][t.iZsun] + t['t_PMS'][t.iZsun]
    v = {m: life[t.im(m)] for m in (1, 8, 40)}
    ok = 8e9 < v[1] < 1.3e10 and 3e7 < v[8] < 5e7 and 3.5e6 < v[40] < 6e6
    return ok, ' '.join(f'{m}Msun={x:.3g}yr(M={t.M[t.im(m)]:.3g})' for m, x in v.items())

@test
def pms_positive_and_short(t):
    """0 < t_PMS < lifetime for every star"""
    p, l = t['t_PMS'], t['lifetime_yr']
    bad = int(((p <= 0) | (p >= l)).sum())
    return bad == 0, f'violations={bad} min t_PMS={p.min():.3g} max t_PMS/lifetime={np.max(p / l):.3g}'

@test
def remnant_types_by_mass(t):
    """WD below 7 Msun, no WD above 12 Msun; WD masses 0.45-1.4 Msun"""
    rt, rm = t['remnant_type'], t['remnant_mass']
    lo = t.M < 7; hi = t.M > 12
    wd = rt == 0
    bad_lo = int((~wd[:, lo]).sum()); bad_hi = int((wd[:, hi]).sum())
    wdm = rm[wd & (t.M[None, :] > 0.9)]
    ok = bad_lo == 0 and bad_hi == 0 and wdm.min() > 0.45 and wdm.max() < 1.4
    return ok, f'nonWD<7={bad_lo} WD>12={bad_hi} WDmass=[{wdm.min():.3f},{wdm.max():.3f}]'

@test
def ifmr_anchor_solar(t):
    """solar-Z initial-final mass: 1 Msun->0.5-0.6, 3 Msun->0.65-0.85 (Cummings+18 band)"""
    rm = t['remnant_mass'][t.iZsun]
    a, b = rm[t.im(1)], rm[t.im(3)]
    return 0.5 < a < 0.62 and 0.65 < b < 0.85, f'M_WD(1)={a:.3f} M_WD(3)={b:.3f}'

def _zams(t, key, m):
    """value at the first alive age after t_PMS (ZAMS) at solar Z"""
    iz, i = t.iZsun, t.im(m)
    age = 10 ** t['log_age_yr']; k = int(np.searchsorted(age, t['t_PMS'][iz, i] * 1.05))
    return t.f[key][iz, i, k]

@test
def zams_luminosity_anchors(t):
    """solar-Z ZAMS log L_bol: 1 Msun 33.2-33.8 (L_sun=33.58), 10 Msun 36.9-37.6"""
    a, b = _zams(t, 'log_L_bol', 1), _zams(t, 'log_L_bol', 10)
    return 33.2 < a < 33.8 and 36.9 < b < 37.6, f'logL(1)={a:.2f} logL(10)={b:.2f}'

@test
def ionizing_photon_anchor(t):
    """solar-Z ZAMS log Q_ion: 40 Msun 48.7-49.6 (Sternberg+03/Martins+05 band)"""
    q = _zams(t, 'log_Q_ion', 40)
    return 48.7 < q < 49.6, f'logQ(40)={q:.2f}'

@test
def x_init_consistent_with_Z(t):
    """X_init metals (C..Zn) sum to Z at every grid metallicity (1e-6 rel)"""
    x = t['X_init']; zsum = x[:, 2:].sum(1)
    err = np.max(abs(zsum / t.Z - 1))
    return err < 1e-6 and np.allclose(x.sum(1), 1, atol=1e-8), f'max|ΣZ/Z-1|={err:.2e} max|ΣX-1|={np.max(abs(x.sum(1)-1)):.2e}'

@test
def imf_integrated_yields(t):
    """Kroupa IMF, 0.08-100 Msun, solar Z: CCSN+AGB metal yield 0.005-0.05 per Msun formed; Fe 3e-4-3e-3"""
    m = t.M; iz = t.iZsun; s = (m >= 0.08) & (m <= 100)
    xi = np.where(m < 0.5, m ** -1.3, 0.5 * m ** -2.3)
    w = np.gradient(m) * xi; w = w * s; w /= np.sum(w * m)
    ej = t['elem_ej_SN_mass'][iz] + t['elem_ej_AGB_mass'][iz]
    xin = t['X_init'][iz]
    net = ej - np.outer(ej.sum(1), xin)          # net = ejected minus initial composition
    yZ = np.sum(w[:, None] * net[:, 2:]); yFe = np.sum(w * net[:, ELEMS.index('Fe')])
    return 0.005 < yZ < 0.05 and 3e-4 < yFe < 3e-3, f'y_Z={yZ:.4f} y_Fe={yFe:.2e} per Msun formed'

@test
def provenance_attrs(t):
    """file records yield/evolution sources and radiation source"""
    need = ['yield_sources', 'evolution_sources', 'radiation_source']
    miss = [a for a in need if a not in t.f.attrs]
    return not miss, f"missing={miss} radiation_source={t.f.attrs.get('radiation_source')}"


def recompute_radiation(t, source, lookup):
    """band outputs recomputed by the builder from the table's own Teff, L_bol, log g and Z"""
    import build_stellar_tables as B
    B.RADIATION.update(source=source, lookup=lookup)
    lT, lL, lR, Mc = t['logTeff'], t['log_L_bol'], t['logR_cm'], t['M_current']
    alive = (lL > -50) & (Mc > 0)
    lZ = np.broadcast_to(np.log10(t.Z)[:, None, None], lT.shape)
    track = None
    if source == 'hotstar_stack':
        sa = t['surface_abundances']; age = 10 ** t['log_age_yr']
        mdot = -np.gradient(Mc, age, axis=2)                 # Msun/yr from the resampled M_current
        track = {'X': sa[..., 0][alive], 'C': sa[..., 2][alive], 'N': sa[..., 3][alive], 'R_cm': 10 ** lR[alive],
                 'log_Mdot': np.log10(np.maximum(mdot, 1e-30))[alive]}
    r = B.compute_radiation_from_lut(lT[alive], lL[alive] - np.log10(B.L_sun), B.log_g(Mc[alive], 10 ** lR[alive]), lZ[alive], track=track)
    return alive, r

def radiation_residuals(t, source, lookup):
    alive, r = recompute_radiation(t, source, lookup); out = {}
    for key, ds in (('Q_ion_H0', 'log_Q_ion_H0'), ('L_LW', 'log_L_LW'), ('L_FUV', 'log_L_FUV')):
        tab = t[ds][alive]; new = np.log10(np.maximum(r[key], 1e-300))
        m = (tab > 20) & (new > 20)
        out[key] = np.percentile(np.abs(new[m] - tab[m]), [50, 90])
    return out

@test
def photon_rates_nested(t):
    """recorded source on a Teff (3.5-200 kK) x log g (0-6) x Z lattice: Q_ion >= Q_He0 >= Q_He1 and Q_ion >= Q_ion_H0"""
    import build_stellar_tables as B
    a = t.f.attrs; B.RADIATION.update(source=a['radiation_source_option'], lookup=a.get('atlas9_lookup_option'))
    lT, lg, lZ = np.meshgrid(np.linspace(np.log10(3500), np.log10(2e5), 60), np.linspace(0, 6, 13), np.log10(t.Z), indexing='ij')
    lT, lg, lZ = lT.ravel(), lg.ravel(), lZ.ravel()
    r = B.compute_radiation_from_lut(lT, np.full(lT.size, 5.0), lg, lZ)
    worst = max(float(np.max(r[x] - r[y] * (1 + 1e-4))) for y, x in (('Q_ion', 'Q_He0'), ('Q_He0', 'Q_He1'), ('Q_ion', 'Q_ion_H0')))
    return worst <= 0, 'max (inner - outer) photon rate on the lattice: %.3g /s (L = 1e5 Lsun)' % worst

def _hotstar(fn):
    """evaluate the hotstar_stack source directly (independent of the table's own source)"""
    import build_stellar_tables as B
    saved = dict(B.RADIATION); B.RADIATION.update(source='hotstar_stack', lookup=None)
    try: return fn(B)
    finally: B.RADIATION.update(saved)

def _ladder(B, T, lg=4.0, z=0.014, track=None):
    lT = np.log10(T); r = B.compute_radiation_from_lut(lT, np.full(lT.size, 5.0), np.full(lT.size, lg), np.full(lT.size, np.log10(z)), track=track)
    L = 10 ** 5.0 * B.L_sun
    return {k: r[k] / L for k in ('Q_ion_H0', 'Q_He0', 'Q_He1', 'L_FUV', 'L_LW')}

@test
def hotstar_seams_continuous(t):
    """hotstar_stack has no jump at its seams: |log v(T(1+1e-5)) - log v(T(1-1e-5))| < 0.01 dex for Q_H0, Q_He0, Q_He1,
    L_FUV, L_LW at the ATLAS9/TLUSTY blend edges (15, 20 kK), the BSTAR/OSTAR blend (27, 30 kK), the TLUSTY hot edge
    and its blackbody blend (55, 65 kK) and the PoWR OB He+ edge and blend (62, 73 kK); log g 4, three Z"""
    def run(B):
        worst = (0., None)
        for z in (0.014, 0.0014, 0.00014):
            for Ts in (15000, 20000, 27000, 30000, 55000, 65000, 62000, 62000 * 65 / 55.):
                r = _ladder(B, np.array([Ts * (1 - 1e-5), Ts * (1 + 1e-5)]), 4.0, z)
                for k, v in r.items():
                    if min(v) <= 0: continue
                    d = abs(np.log10(v[1] / v[0]))
                    if d > worst[0]: worst = (d, '%s at %.0f K, Z=%g' % (k, Ts, z))
        return worst
    w = _hotstar(run)
    return w[0] < 0.01, 'largest jump %.2e dex (%s)' % w

@test
def hotstar_q_h0_matches_tlusty_nodes(t):
    """at TLUSTY model nodes >= 20 kK (cell class 0) the source returns the LUT Q_H0 per erg to 1e-9; surface q0 at
    30 kK / log g 4 / Zsun (BSTAR -> OSTAR blend) lies between Lanz & Hubeny 2007 Table 4 (BSTAR, log 22.90) and the
    OSTAR2002 model (log 23.01)"""
    def run(B):
        H = B.get_hotstar_lut(); tl = H.tl; i0 = H.QB.index('Q_ion_H0'); worst = 0.; n = 0
        for i, lT in enumerate(tl['log_Teff']):
            if 10 ** lT < 20000: continue
            for j, lg in enumerate(tl['logg']):
                for k, zr in enumerate(tl['Z_rel']):
                    if tl['cell_class'][i, j, k] != 0: continue
                    r = B.compute_radiation_from_lut(np.array([lT]), np.array([5.0]), np.array([lg]), np.array([np.log10(zr * B.Z_SOLAR)]))
                    q = r['Q_ion_H0'][0] / (1e5 * B.L_sun); ref = tl['photon_rates_per_L'][i, j, k, i0]
                    worst = max(worst, abs(q / ref - 1)); n += 1
        T = 30000.; r = _ladder(B, np.array([T])); q0 = np.log10(r['Q_ion_H0'][0] * 5.670374419e-5 * T ** 4)
        return worst, n, q0
    worst, n, q0 = _hotstar(run)
    return worst < 1e-9 and 22.89 <= q0 <= 23.02, '%d nodes, max rel dev %.1e; log q0(30 kK, 4.0, Zsun) = %.3f (BSTAR L&H07 22.90, OSTAR 23.01)' % (n, worst, q0)

@test
def hotstar_he_plus_flagged(t):
    """hotstar_stack documents Q_He1 as order-of-magnitude (wind-set) in the module and in the table attributes"""
    import hotstar_radiation as HS
    ok = 'order-of-magnitude' in HS.__doc__
    if t.f.attrs.get('radiation_source_option') == 'hotstar_stack':
        ok &= 'order of magnitude' in str(t.f.attrs.get('hotstar_rules')) and 'order of magnitude' in str(t.f['log_Q_He1'].attrs.get('note', ''))
    return ok, 'Q_He1 flagged (module%s)' % (', table attrs' if t.f.attrs.get('radiation_source_option') == 'hotstar_stack' else '')

@test
def hotstar_nested_and_wr(t):
    """hotstar_stack: Q_ion >= Q_He0 >= Q_He1 on a lattice incl. WR points (X = 0, Mdot 1e-5) at 4 Z; WR He+ is never blackbody"""
    def run(B):
        import hotstar_radiation as HS
        lT, lg, lZ = np.meshgrid(np.linspace(np.log10(3500), np.log10(2e5), 50), np.linspace(1, 6, 11), np.log10([0.0002, 0.0014, 0.006, 0.02]), indexing='ij')
        lT, lg, lZ = lT.ravel(), lg.ravel(), lZ.ravel(); n = lT.size
        out = []
        for X, C in ((0.7, 0.002), (0.0, 0.002), (0.0, 0.3), (0.2, 0.001)):
            tr = {'X': np.full(n, X), 'C': np.full(n, C), 'N': np.full(n, 0.01), 'R_cm': np.full(n, 2 * 6.957e10), 'log_Mdot': np.full(n, -5.0)}
            r = B.compute_radiation_from_lut(lT, np.full(n, 5.0), lg, lZ, track=tr)
            out.append(max(float(np.max(r[x] - r[y] * (1 + 1e-4))) for y, x in (('Q_ion', 'Q_He0'), ('Q_He0', 'Q_He1'), ('Q_ion', 'Q_ion_H0'))))
            E, Q, Cc = HS.hotstar_band_fractions(B.get_hotstar_lut(), lT, lg, lZ, B.Z_SOLAR, B.spectral_band_fractions, B.planck_band_fractions, tr)
            wr = (tr['X'] < HS.WR_X_MAX) & (lT > HS.WR_LOGTEFF_MIN)
            out.append(int(np.sum(Cc[wr, HS.GROUPS.index('Q_He1')] == 3)))
        return out
    o = _hotstar(run)
    return max(o[0::2]) <= 0 and sum(o[1::2]) == 0, 'max inner-outer %.3g /s; WR points with blackbody He+: %d' % (max(o[0::2]), sum(o[1::2]))

@test
def hotstar_cmfgen_continuous(t):
    """hotstar_stack with the CMFGEN layer has no jump: dense scans in log Teff (20-100 kK, steps of 5e-5 dex) at log g 3.0,
    3.5, 4.0, 4.3, 4.5 and in log g (2.0-4.8, steps of 5e-4) at 35, 45, 52, 58, 66 kK, Z = 0.0002, 0.002, 0.014; at the 20
    largest steps of each line and band (log Q_H0, Q_He0, Q_He1, L_FUV, L_LW per erg) the change over +-1e-7 dex is < 0.01 dex"""
    def run(B):
        import hotstar_radiation as HS
        H = B.get_hotstar_lut(); worst = (0., None); big = (0., None)
        def bands(lT, lg, z):
            E, Q, C = HS.hotstar_band_fractions(H, lT, lg, np.full(lT.size, np.log10(z)), B.Z_SOLAR, B.spectral_band_fractions, B.planck_band_fractions)
            return {'Q_H0': Q[:, H.QB.index('Q_ion_H0')], 'Q_He0': Q[:, H.QB.index('Q_He0')], 'Q_He1': Q[:, H.QB.index('Q_He1')],
                    'L_FUV': E[:, H.EB.index('NUV_hi')] + E[:, H.EB.index('FUV_mid')], 'L_LW': E[:, H.EB.index('LW')]}
        lines = [(np.arange(np.log10(2e4), 5.0, 5e-5), np.full(1, g), 'log g %.1f' % g, 0) for g in (3.0, 3.5, 4.0, 4.3, 4.5)]
        lines += [(np.full(1, np.log10(T)), np.arange(2.0, 4.8, 5e-4), 'T %d kK' % (T / 1e3), 1) for T in (35e3, 45e3, 52e3, 58e3, 66e3)]
        for z in (0.0002, 0.002, 0.014):
            for lT, lg, name, ax in lines:
                lT, lg = np.broadcast_arrays(lT, lg); lT, lg = lT.ravel(), lg.ravel()
                for k, a in bands(lT, lg, z).items():
                    d = np.abs(np.diff(np.log10(np.maximum(a, 1e-300)))); d[(a[1:] < 1e-200) | (a[:-1] < 1e-200)] = 0
                    if d.max() > big[0]: big = (float(d.max()), '%s, %s, Z=%g' % (k, name, z))
                    for i in np.argsort(d)[-20:]:
                        x = np.array([lT[i], lT[i + 1]]) if ax == 0 else np.array([lg[i], lg[i + 1]])
                        x0 = x[0] + (x[1] - x[0]) * np.linspace(0, 1, 9)        # refine the step, take its steepest piece
                        p = (x0, np.full(9, lg[i])) if ax == 0 else (np.full(9, lT[i]), x0)
                        v = np.log10(np.maximum(bands(p[0], p[1], z)[k], 1e-300)); j = int(np.argmax(np.abs(np.diff(v))))
                        c = 0.5 * (x0[j] + x0[j + 1]); q = (np.array([c - 1e-7, c + 1e-7]), np.full(2, lg[i])) if ax == 0 else (np.full(2, lT[i]), np.array([c - 1e-7, c + 1e-7]))
                        w = bands(q[0], q[1], z)[k]
                        if min(w) > 1e-200 and abs(np.log10(w[1] / w[0])) > worst[0]:
                            worst = (abs(float(np.log10(w[1] / w[0]))), '%s, %s, Z=%g at %s %.6f' % (k, name, z, 'log T' if ax == 0 else 'log g', c))
        return worst, big
    w, big = _hotstar(run)
    return w[0] < 0.01, 'largest change over 2e-7 dex at the steepest scan steps %.2e (%s); largest scan step %.2f dex (%s)' % (w + big)

@test
def hotstar_cmfgen_nodes(t):
    """every CMFGEN model is reproduced by its mesh (rel. 1e-9, all bands); where the stack takes CMFGEN fully (W = 1, a node
    of a set at its own Z, inside the hull) hotstar_stack returns the node's Q_H0, Q_He0, Q_He1, L_FUV_mid, L_LW per erg to 1e-9"""
    def run(B):
        import hotstar_radiation as HS
        H = B.get_hotstar_lut(); w1 = 0.; w2 = 0.; n2 = 0
        for c in H.cm:
            v, a, lq, out = HS.cmf_set_eval(c, c['log_Teff'], c['logg'])
            pos = c['vals'] > 1e-200
            w1 = max(w1, float(np.max(np.abs(v[pos] / c['vals'][pos] - 1))))
            dg = {}; lZ = np.full(len(c['logg']), np.log10(c['Z_rel'] * B.Z_SOLAR))
            E, Q, C = HS.hotstar_band_fractions(H, c['log_Teff'], c['logg'], lZ, B.Z_SOLAR, B.spectral_band_fractions, B.planck_band_fractions, diag=dg)
            m = dg['cmfgen_W'] >= 1 - 1e-12; n2 += int(m.sum())
            for col, ref in ((Q[:, H.QB.index('Q_ion_H0')], c['photon_rates_per_L'][:, H.QB.index('Q_ion_H0')]),
                             (Q[:, H.QB.index('Q_He0')], c['photon_rates_per_L'][:, H.QB.index('Q_He0')]),
                             (Q[:, H.QB.index('Q_He1')], c['photon_rates_per_L'][:, H.QB.index('Q_He1')]),
                             (E[:, H.EB.index('FUV_mid')], c['energy_fracs'][:, H.EB.index('FUV_mid')]), (E[:, H.EB.index('LW')], c['energy_fracs'][:, H.EB.index('LW')])):
                k = m & (ref > 1e-200)
                if k.any(): w2 = max(w2, float(np.max(np.abs(col[k] / ref[k] - 1))))
        return w1, w2, n2, sum(len(c['logg']) for c in H.cm)
    w1, w2, n2, n = _hotstar(run)
    return w1 < 1e-9 and w2 < 1e-9 and n2 > 0, 'mesh: %d nodes, max rel dev %.1e; %d nodes taken fully by the stack, max rel dev %.1e' % (n, w1, n2, w2)

@test
def hotstar_cmfgen_only_fills(t):
    """the CMFGEN layer changes no point whose TLUSTY interpolation uses only non-filled cells at or below the 55 kK edge (lattice
    20-55 kK x log g 2-4.75 x 4 Z): those points are bit-identical with and without the CMFGEN sets"""
    def run(B):
        import hotstar_radiation as HS
        H = B.get_hotstar_lut(); tl = H.tl
        lT, lg, lZ = np.meshgrid(np.linspace(np.log10(2e4), np.log10(5.5e4), 60), np.linspace(2, 4.75, 45), np.log10([0.0002, 0.002, 0.014, 0.03]), indexing='ij')
        lT, lg, lZ = lT.ravel(), lg.ravel(), lZ.ravel()
        args = (lT, lg, lZ, B.Z_SOLAR, B.spectral_band_fractions, B.planck_band_fractions)
        E1, Q1, C1 = HS.hotstar_band_fractions(H, *args)
        cm = H.cm; H.cm = []
        try: E0, Q0, C0 = HS.hotstar_band_fractions(H, *args)
        finally: H.cm = cm
        (f,), _, _ = HS._interp((tl['log_Teff'], tl['logg'], np.log10(tl['Z_rel'])), ((tl['cell_class'] == 2).astype(float)[..., None],),
                                np.zeros_like(tl['cell_class']), np.column_stack([lT, lg, lZ - np.log10(B.Z_SOLAR)]))
        clean = f[:, 0] == 0; changed = ~(np.all(E1 == E0, 1) & np.all(Q1 == Q0, 1) & np.all(C1 == C0, 1))
        return int(clean.sum()), int((clean & changed).sum()), int((~clean & changed).sum())
    n, bad, used = _hotstar(run)
    return bad == 0 and used > 0, '%d points without filled cells, %d changed; %d hole points take CMFGEN' % (n, bad, used)

@test
def radiation_source_recorded(t):
    """table names its radiation source and options; atlas9 builds record the LUT checksum and coverage"""
    a = t.f.attrs; src = a.get('radiation_source_option')
    if src not in ('atlas9', 'planck', 'atlas9_nonionizing', 'hotstar_stack'): return False, f'radiation_source_option={src!r}'
    ok = 'TODO' not in str(a.get('radiation_source', 'TODO'))
    if src != 'planck':
        import hashlib, build_stellar_tables as B
        sha = hashlib.sha256(open(B.SPECTRAL_LUT_PATH, 'rb').read()).hexdigest()
        ok &= a.get('spectral_lut_sha256') == sha and len(a.get('radiation_coverage_track_points', [])) == (4 if src == 'hotstar_stack' else 3)
        ok &= (a.get('atlas9_lookup_option') in ('legacy', 'clamp', 'masked_planck')) == (src == 'atlas9')
        keys = ('Q_H0', 'Q_He0', 'Q_He1', 'L_FUV', 'L_LW') if src == 'hotstar_stack' else ('Q_H0_planck', 'L_FUV', 'L_LW')
        ok &= all(abs(sum(a.get('radiation_coverage_%s_frac' % k, [0])) - 1) < 1e-9 for k in keys)
        if src == 'hotstar_stack':
            ok &= a.get('hotstar_lut_sha256') == hashlib.sha256(open(B.HOTSTAR_LUT_PATH, 'rb').read()).hexdigest()
            ok &= all(len(a['radiation_coverage_%s_frac' % k]) == 4 for k in keys)
    return bool(ok), f"{src} / {a.get('atlas9_lookup_option')}: {a.get('radiation_source')}"

@test
def radiation_no_silent_fallback(t):
    """builder refuses to compute radiation without a source, and atlas9 without the LUT file"""
    import build_stellar_tables as B
    x = np.array([4.0]); errs = []
    saved = dict(B.RADIATION), B.SPECTRAL_LUT, B.SPECTRAL_LUT_PATH
    try:
        for src, look, path in ((None, None, saved[2]), ('atlas9', None, saved[2]), ('atlas9', 'clamp', saved[2] + '.missing')):
            B.RADIATION.update(source=src, lookup=look); B.SPECTRAL_LUT, B.SPECTRAL_LUT_PATH = None, path
            try: B.compute_radiation_from_lut(x, x, x, x); errs.append('no error')
            except (ValueError, FileNotFoundError) as e: errs.append(type(e).__name__)
    finally:
        B.RADIATION.update(saved[0]); B.SPECTRAL_LUT, B.SPECTRAL_LUT_PATH = saved[1], saved[2]
    return 'no error' not in errs, 'no source / atlas9 without lookup / missing LUT -> %s' % errs

@test
def radiation_matches_recorded_source(t):
    """Q_H0, L_LW, L_FUV recomputed from the table's (age-resampled) Teff/L/g/Z with the recorded source match the table
    (median < 0.02 dex), and the other source is >= 10x further off"""
    a = t.f.attrs; src, look = a['radiation_source_option'], a.get('atlas9_lookup_option')
    same = radiation_residuals(t, src, look)
    other = radiation_residuals(t, *{'atlas9': ('planck', None), 'planck': ('atlas9', 'clamp'),
                                     'atlas9_nonionizing': ('atlas9', 'legacy'), 'hotstar_stack': ('planck', None)}[src])
    ok = max(v[0] for v in same.values()) < 0.02 and max(v[0] for v in other.values()) >= 10 * max(v[0] for v in same.values())
    return ok, 'recorded: ' + ' '.join(f'{k} {v[0]:.4f}/{v[1]:.3f}' for k, v in same.items()) + \
               ' | other source: ' + ' '.join(f'{k} {v[0]:.3f}' for k, v in other.items()) + ' (median/p90 |dex|)'


def compare(a, b):
    fa, fb = h5py.File(a, 'r'), h5py.File(b, 'r')
    da, db = {}, {}
    fa.visititems(lambda n, o: da.__setitem__(n, o) if isinstance(o, h5py.Dataset) else None)
    fb.visititems(lambda n, o: db.__setitem__(n, o) if isinstance(o, h5py.Dataset) else None)
    nd = 0
    for k in sorted(set(da) ^ set(db)): print(f'  ONLY IN {"A" if k in da else "B"}: {k}'); nd += 1
    for k in sorted(set(da) & set(db)):
        x, y = da[k][()], db[k][()]
        if x.shape != y.shape: print(f'  SHAPE {k} {x.shape} vs {y.shape}'); nd += 1
        elif not np.array_equal(x, y):
            if x.dtype.kind in 'fi':
                fin = np.isfinite(x) & np.isfinite(y)
                d = np.max(abs(x[fin] - y[fin])) if fin.any() else np.nan
                print(f'  DIFF {k}: {int((x != y).sum())} entries, max|Δ|={d:.3g}')
            else: print(f'  DIFF {k}')
            nd += 1
    ka, kb = dict(fa.attrs), dict(fb.attrs)
    for k in sorted(set(ka) | set(kb)):
        if str(ka.get(k)) != str(kb.get(k)): print(f'  ATTR {k}: {ka.get(k)!r} vs {kb.get(k)!r}'); nd += 1
    print(f'{"IDENTICAL" if nd == 0 else f"{nd} differences"}: {a} vs {b}')
    return nd == 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--table', default=os.path.join(ROOT, 'stellar_tables_unified.hdf5'))
    ap.add_argument('--compare', nargs=2)
    a = ap.parse_args()
    if a.compare: sys.exit(0 if compare(*a.compare) else 1)
    t = T(a.table); print(f'table {a.table}: {t.nZ} Z x {t.nM} M x {t.nA} ages')
    nf = 0
    for fn in TESTS:
        try: ok, msg = fn(t)
        except Exception as e: ok, msg = False, f'EXCEPTION {e!r}'
        xf = XFAIL.get(fn.__name__)
        tag = ('PASS' if ok else 'FAIL') if not xf else ('XPASS' if ok else 'XFAIL')
        nf += tag == 'FAIL'
        print(f"  [{tag}] {fn.__name__}: {msg}" + (f'  (expected: {xf})' if xf else ''))
    print(f'\n{nf}/{len(TESTS)} failed'); sys.exit(1 if nf else 0)
