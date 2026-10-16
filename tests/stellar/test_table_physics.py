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
