#!/usr/bin/env python3
"""Per-element checks of the stellar table against the raw track and yield files.

  test_yields_raw.py BUILD_DIR      exit 0 = all PASS (raw inputs from BUILD_DIR/data_root)

- pre-SN mass = last MASS of the raw PARSEC v2 track (all wind mass reaches the wind table) for M >= 8
- surface composition at the first alive age = birth composition X_init, every element, every (Z, M)
- AGB: Karakas 2010 Z = 1e-4, 6 Msun is the VW93 model alone; KL16 elements that AGB stars do not process leave at
  their birth fraction
- 8-14 Msun at Z = 0.014: wind + SN = X_init (M - rem) + Limongi+2025 net yield (total - X Mej) per Limongi element
"""
import os, re, sys
from glob import glob
import numpy as np, h5py

ELEMS = ['H', 'He', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
         'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
TESTS = []
def test(fn): TESTS.append(fn); return fn


@test
def preSN_mass_is_raw_track_final_mass(t, w, root):
    Zt, Mt, Mp = t['Z'][()], t['M_init'][()], t['M_preSN'][()]
    worst, n = (0., None), 0
    for d in glob(os.path.join(root, 'parsec2_vms', 'tracks_Z*')):
        Zp = float(re.search(r'tracks_Z([\d.Ee-]+)', d).group(1))
        iz = np.where(np.isclose(Zt, Zp, rtol=1e-9))[0]
        if not iz.size: continue
        for f in glob(os.path.join(d, '*.TAB')):
            M = float(re.search(r'_M([\d.]+)\.TAB$', f).group(1))
            im = np.where(Mt == M)[0]
            if M < 8 or not im.size or os.path.getsize(f) == 0: continue
            lines = open(f).read().splitlines(); names = lines[2].split()
            rows = [l.split() for l in lines[3:] if len(l.split()) == len(names)]
            if len(rows) < 50: continue   # the builder skips truncated tracks
            m_end = float(rows[-1][names.index('MASS')])
            e = abs(Mp[iz[0], im[0]] - m_end) / M; n += 1
            if e > worst[0]: worst = (e, (Zp, M, Mp[iz[0], im[0]], m_end))
    return n > 0 and worst[0] < 1e-4, '%d PARSEC v2 tracks; worst |M_preSN - M_end|/M = %.2e at (Z, M, table, raw) = %s' % (n, worst[0], worst[1])


@test
def zams_surface_is_birth_composition(t, w, root):
    S, X, Mc = t['surface_abundances'], t['X_init'][()], t['M_current'][()]
    worst, nbad = (0., None), 0
    for iz in range(S.shape[0]):
        for im in range(S.shape[1]):
            ia = np.where(Mc[iz, im] > 0)[0]
            if not ia.size: continue
            s0 = S[iz, im, ia[0], :].astype(np.float64)
            e = np.max(np.abs(s0 / X[iz] - 1))
            nbad += e > 1e-3
            if e > worst[0]: worst = (e, (t['Z'][iz], t['M_init'][im], ELEMS[int(np.argmax(np.abs(s0 / X[iz] - 1)))]))
    return nbad == 0, 'max |X_surf(ZAMS)/X_init - 1| = %.2e at (Z, M, elem) = %s; %d (Z, M) above 1e-3' % (worst[0], worst[1], nbad)


def _agb_net(t, iz, im):
    """AGB lifetime net yield as GIZMO ejects it (no winds below 8 Msun: elem_ej_AGB scaled to M - rem)"""
    M, rm = t['M_init'][im], t['remnant_mass'][iz, im]
    ej = t['elem_ej_AGB_mass'][iz, im].astype(np.float64)
    return ej / ej.sum() * (M - rm) - t['X_init'][iz] * (M - rm)


@test
def karakas2010_z1e4_6msun_is_vw93(t, w, root):
    f = os.path.join(root, 'karakas2010', 'tablea5.dat')
    nN = {}
    for ln in open(f):
        if not ln.strip(): continue
        M0, Z0, Mf, iso = float(ln[0:4]), float(ln[5:11]), float(ln[12:17]), ln[21:25].strip().lower()
        if (M0, Z0) == (6.0, 0.0001) and iso in ('n13', 'n14', 'n15'):
            nN[Mf] = nN.get(Mf, 0.0) + float(ln[34:48])
    iz, im = int(np.argmin(abs(t['Z'][()] - 1e-4))), int(np.argmin(abs(t['M_init'][()] - 6.0)))
    have = _agb_net(t, iz, im)[ELEMS.index('N')]
    return abs(have / nN[1.006] - 1) < 0.02, 'net N: table %.4e, VW93 %.4e, Reimers %.4e' % (have, nN[1.006], nN[1.008])


def _kl16_models(path):
    """{M: (M_mix, Mexp, {elem: [X/H]})} of a KL16 yield file, the standard M_mix per mass (KL16 Table 2)"""
    std = lambda M: 2e-3 if M <= 3 else 1e-3 if M <= 4 else 1e-4 if M < 5 else 0.0
    out = {}
    for blk in open(path).read().split('# Initial mass')[1:]:
        L = blk.splitlines()
        M = float(L[0].split('=')[1].split(',')[0]); mmix = float(re.search(r'M_mix = ([\d.Ee+-]+)', L[0]).group(1))
        mexp = float(re.search(r'Mass expelled =\s*([\d.]+)', L[1]).group(1))
        xh = {l.split()[0]: float(l.split()[3]) for l in L[3:] if len(l.split()) >= 7}
        if M not in out or abs(mmix - std(M)) < abs(out[M][0] - std(M)): out[M] = (mmix, mexp, xh)
    return out


@test
def kl16_unprocessed_elements_follow_raw(t, w, root):
    """KL16 net yields against the model's own initial composition: S, Ca, Ti, Cr, Fe leave at birth fraction times
    the raw production factor 10^([X/H] - [Ca/H]) (Ca is not processed by AGB stars)"""
    Zt, Mt, rt, rm = t['Z'][()], t['M_init'][()], t['remnant_type'][()], t['remnant_mass'][()]
    worst = (0., None)
    for zs, Z in (('007', 0.007), ('014', 0.014), ('03', 0.03)):
        iz = np.where(np.isclose(Zt, Z))[0]
        if not iz.size: continue
        iz = iz[0]
        for M, (_, mexp, xh) in _kl16_models(os.path.join(root, 'karakas_lugaro2016', 'yields', 'yield_z%s.dat' % zs)).items():
            im = np.where(Mt == M)[0]
            if not im.size or rt[iz, im[0]] != 0: continue
            im = im[0]; ej = t['elem_ej_AGB_mass'][iz, im].astype(np.float64); fr = ej / ej.sum()
            for e in ('S', 'Ca', 'Ti', 'Cr', 'Fe'):
                k = ELEMS.index(e)
                pf = 10 ** (xh[e.lower()] - xh['ca'])
                want = 1 + (pf - 1) * mexp / (M - rm[iz, im])
                r = abs(fr[k] / t['X_init'][iz, k] - want)
                if r > worst[0]: worst = (r, (Z, M, e, fr[k] / t['X_init'][iz, k], want))
    return worst[0] < 0.01, 'max |X_AGB/X_init - raw production factor| = %.4f at (Z, M, elem, table, raw) = %s' % worst


@test
def limongi_lifetime_budget_solar(t, w, root):
    import build_stellar_tables as B
    ms, rems, tab = B.get_limongi25_yields()
    x0 = B._parse_parsec2_init_comp(os.path.join(root, 'parsec2_vms'))[0.014]
    Zt, Mt = t['Z'][()], t['M_init'][()]
    iz = int(np.argmin(abs(Zt - 0.014)))
    cw, SN, X, rm = w['elem_ej_wind_cumulative'], t['elem_ej_SN_mass'][()], t['X_init'][()], t['remnant_mass'][()]
    worst = (0., None)
    for M in (11., 12.):
        im = int(np.where(Mt == M)[0][0]); mej = M - np.interp(M, ms, rems)
        for e, ys in tab.items():
            if e == 'H': continue        # H takes the sum of all net yields (mass closure, tested in test_invariants)
            k = ELEMS.index(e)
            want = np.interp(M, ms, ys) - x0[e] * mej
            have = float(cw[iz, im, -1, k]) + SN[iz, im, k] - X[iz, k] * (M - rm[iz, im])
            err = abs(have - want) / (1e-3 * X[iz, k] * M + 1e-3 * abs(want) + 1e-9)
            if err > worst[0]: worst = (err, (M, e, have, want))
    return worst[0] < 1, 'Limongi elements except H at 11, 12 Msun: worst |dnet|/tol = %.2g at (M, elem, table, raw) = %s' % worst


@test
def remnant_type_is_raw_parsec2_fate(t, w, root):
    """remnant type at PARSEC v2 grid points (Z, M >= 14) = the raw total-ejecta SNT column (e.g. Z = 0.03: failed SNe
    only at 26-30 Msun, core collapse from 35 Msun; Z = 0.04, 0.06 take the Z = 0.03 fates)"""
    code = {'CCSN': 2, 'FSN': 3, 'PPISN': 4, 'PISN': 5, 'DBH': 6, 'NOT': 0}
    Zt, Mt, rt = t['Z'][()], t['M_init'][()], t['remnant_type'][()]
    bad, n = [], 0
    for f in glob(os.path.join(root, 'parsec2_vms', 'ejecta_Z*', '*total_ejecta.dat')):
        Zp = float(re.search(r'ejecta_Z([\d.Ee-]+)', f).group(1))
        izs = [i for i, z in enumerate(Zt) if np.isclose(z, Zp) or (Zp == 0.03 and z > 0.03)]
        rows = {}
        for ln in open(f).read().splitlines()[4:]:
            p = ln.split()
            if len(p) > 7: rows.setdefault(float(p[0]), p[6])
        for iz in izs:
            for im, M in enumerate(Mt):
                if M >= 14 and M in rows:
                    n += 1
                    if rt[iz, im] != code[rows[M]]: bad.append((Zt[iz], M, int(rt[iz, im]), rows[M]))
    return n > 0 and not bad, '%d grid points; mismatches (Z, M, table, raw): %s' % (n, bad[:6])


@test
def yield_source_mass_ranges_stated(t, w, root):
    src = str(t.attrs.get('yield_sources', ''))
    return 'Limongi+(2025) 8-14 Msun' in src and '8-13' not in src, 'yield_sources = %r' % src


def main():
    b = sys.argv[1]; root = os.path.join(b, 'data_root')
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'stellar'))
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
