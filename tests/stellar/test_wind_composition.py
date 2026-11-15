#!/usr/bin/env python3
"""Wind composition and per-star element budget against the raw catalogues.

  test_wind_composition.py BUILD_DIR      exit 0 = all PASS (raw inputs read from BUILD_DIR/data_root, as build.sh lays them out)

- BoOST surface columns: the builder's isotope -> element map equals the names in fdat-columns.txt.
- Massive-star winds (M >= 14, PARSEC v2 metallicity): the end-of-wind element masses carry the PARSEC v2 wind
  ejecta enrichment for all 27 elements, on the table's birth composition:
      cum_wind_k(end) / W - X_init_k  ==  ejecta_k / (M - Mfin) - X_PARSEC,init_k      (W = M_init - M_preSN)
- Per-star budget (PARSEC v2 CCSN/PPISN/PISN points): wind + SN ejecta = X_init (M - rem) + PARSEC v2 total net yield,
  per element, i.e. the wind-phase production is not counted again at the SN.
Only grid points where the table Z and M are PARSEC v2 grid values are compared (no interpolation involved).
"""
import os, re, sys
from glob import glob
import numpy as np, h5py

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'stellar'))
ELEMS = ['H', 'He', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
         'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
# float32 cumulative wind: element mass known to ~1e-7 relative; compare enrichments to 1e-5 of the element's
# birth fraction plus enrichment, plus 1e-9 absolute (trace elements, e.g. Sc ~ 1e-8 by mass)
REL_TOL, ABS_TOL = 1e-5, 1e-9
TESTS = []
def test(fn): TESTS.append(fn); return fn


def elem_of(iso):
    if iso.upper() == 'AL26': return 'Mg'   # 26Al decays to 26Mg (0.72 Myr)
    e = re.sub(r'\d+$', '', iso)
    e = e[0] + e[1:].lower()
    return e if e in ELEMS else None


def read_parsec2(path, n_lead):
    """rows {Min: (row dict, {isotope: mass}, unfilled isotopes)} and the initial composition {isotope: X} of a PARSEC
    v2 ejecta file; an isotope column that is exactly 0 while present initially is not filled by PARSEC v2 and is
    set to initial fraction x ejected mass (no net yield)"""
    lines = open(path).read().splitlines()
    names = lines[3].split()
    tok = lines[1].lstrip('#').split()
    init_tok = tok[tok.index('SN') + 1:] if 'SN' in tok else tok[n_lead:]
    iso_names = names[n_lead + (1 if 'SNT' in names else 0):]
    x0i = dict(zip(iso_names, (float(v) for v in init_tok)))
    rows = {}
    for ln in lines[4:]:
        p = ln.split()
        if len(p) != len(names): continue
        r = dict(zip(names, p))
        mej = float(r['Min']) - float(r['Mbar'] if 'Mbar' in r else r['Mfin'])
        iso, unf = {}, set()
        for n in iso_names:
            v = float(r[n])
            if v == 0 and x0i[n] > 0: v = x0i[n] * mej; unf.add(n)
            iso[n] = v
        rows.setdefault(float(r['Min']), (r, iso, unf))   # first row wins on a duplicated mass (identical rows)
    return rows, x0i


def to_elems(iso):
    el = dict.fromkeys(ELEMS, 0.0)
    for n, v in iso.items():
        if elem_of(n): el[elem_of(n)] += v
    return el


def parsec2_points(t, root):
    """(iz, im, Zp, (wind row, wind elems), (total row, total elems), X_P elems) at table grid points that are PARSEC v2
    grid points with M >= 14. Where the wind file leaves an isotope unfilled, a filled total holds the explosion part
    only: the wind part (initial fraction x (M - Mfin)) is added."""
    Zt, Mt = t['Z'][()], t['M_init'][()]
    out = []
    for d in sorted(glob(os.path.join(root, 'parsec2_vms', 'ejecta_Z*'))):
        Zp = float(re.search(r'ejecta_Z([\d.Ee-]+)', d).group(1))
        iz = np.where(np.isclose(Zt, Zp, rtol=1e-9))[0]
        if not iz.size: continue
        wf, tf = glob(os.path.join(d, '*winds_ejecta.dat')), glob(os.path.join(d, '*total_ejecta.dat'))
        if not wf or not tf: continue
        wind, x0 = read_parsec2(wf[0], 4)
        tot, _ = read_parsec2(tf[0], 6)
        for im, M in enumerate(Mt):
            if M >= 14 and M in wind and M in tot:
                (wr, wiso, wunf), (tr, tiso, tunf) = wind[M], tot[M]
                tiso = dict(tiso)
                for n in wunf - tunf:
                    tiso[n] += x0[n] * (M - float(wr['Mfin']))
                out.append((int(iz[0]), im, Zp, (wr, to_elems(wiso)), (tr, to_elems(tiso)), to_elems(x0)))
    return out


@test
def boost_isotope_map_matches_fdat_columns(t, w, root):
    import build_stellar_tables as B
    hdr = glob(os.path.join(root, 'boost', 'BoOST-stellarmodels-*', 'fdat-columns.txt'))
    hdr = [h for h in hdr if os.path.getsize(h) > 0]
    if not hdr: return False, 'no fdat-columns.txt in data_root/boost'
    cols = {}
    for tok in open(hdr[0]).read().split():
        m = re.match(r'(\d+):s([A-Z][a-z]?)(\d+)$', tok)   # surface isotope columns, 1-indexed
        if m and m.group(2) in ELEMS and m.group(2) + m.group(3) not in ('C11', 'N12'):   # half-lives 20 min / 11 ms
            el = 'Mg' if m.group(2) + m.group(3) == 'Al26' else m.group(2)                # 26Al decays to 26Mg
            cols.setdefault(el, []).append(int(m.group(1)) - 1)
    bad = [(e, B.BOOST_SURFACE_ELEM.get(e), cols.get(e, [])) for e in sorted(set(cols) | set(B.BOOST_SURFACE_ELEM))
           if sorted(B.BOOST_SURFACE_ELEM.get(e, [])) != sorted(cols.get(e, []))]
    return not bad, 'builder map == fdat-columns.txt for %d elements' % len(cols) if not bad else 'differ (elem, builder, file): %s' % bad


@test
def boost_track_surface_na_mg(t, w, root):
    import build_stellar_tables as B
    fs = [f for f in sorted(glob(os.path.join(root, 'boost', 'BoOST-stellarmodels-MW', 'f0*-100.mw.dat'))) if os.path.getsize(f) > 0]
    if not fs: return False, 'no BoOST MW track with content in data_root'
    B.RADIATION.update(source='planck', lookup=None)   # radiation is not under test here
    raw = np.loadtxt(fs[-1]); tr = B.boost_track_to_common(fs[-1])
    iNa, iMg = B.TRACKED_ELEMENTS.index('Na'), B.TRACKED_ELEMENTS.index('Mg')
    n = len(tr['surface'])     # trim_to_ZAMS drops leading rows only
    na_raw, mg_raw = raw[-n:, 50], raw[-n:, 51] + raw[-n:, 52] + raw[-n:, 53] + raw[-n:, 54]
    e = max(np.max(abs(tr['surface'][:, iNa] - na_raw)), np.max(abs(tr['surface'][:, iMg] - mg_raw)))
    return e < 1e-12, '%s: max |builder - sNa23|, |builder - sMg24..26 + sAl26| = %.2e' % (os.path.basename(fs[-1]), e)


@test
def wind_end_matches_parsec2_wind_ejecta(t, w, root):
    cw, X, Mp = w['elem_ej_wind_cumulative'], t['X_init'][()], t['M_preSN'][()]
    Mt, Zt = t['M_init'][()], t['Z'][()]
    pts = parsec2_points(t, root)
    if not pts: return False, 'no PARSEC v2 grid points in the table grid'
    worst, nbad, zs = (0., None), 0, set()
    for iz, im, Zp, (wr, wel), _, xw in pts:
        W = Mt[im] - Mp[iz, im]
        if W <= 0: continue
        Wp = Mt[im] - float(wr['Mfin'])
        end = cw[iz, im, -1, :].astype(np.float64)
        # PARSEC v2 enrichment; data rules: unfilled (exactly 0) isotope columns = net 0, Li/Be and rounding to H,
        # depletion capped at our birth fraction with the excess to H
        dx = np.array([wel[e] / Wp - xw[e] for e in ELEMS])
        dx[0] -= dx.sum(); dx = np.maximum(dx, -X[iz]); dx[0] -= dx.sum()
        for k, e in enumerate(ELEMS):
            want, have = dx[k], end[k] / W - X[iz, k]
            err = abs(have - want) / (REL_TOL * (X[iz, k] + abs(want)) + ABS_TOL)
            if err > 1: nbad += 1
            if err > worst[0]: worst = (err, (Zt[iz], Mt[im], e, have, want))
        zs.add(Zp)
    msg = '%d (Z,M) points at Z=%s x 27 elements; worst |dx|/tol = %.2g at %s' % (len(pts), sorted(zs), worst[0], worst[1])
    return nbad == 0, msg + ('' if not nbad else '; %d element values off' % nbad)


@test
def wind_na_al_follow_parsec2(t, w, root):
    """wind Na/Mg, Al/Mg at the end of the wind = birth + PARSEC v2 wind enrichment (2%), for Z = 0.001, 0.01 and 20/40/100 Msun"""
    cw, X, Mt, Zt, Mp = w['elem_ej_wind_cumulative'], t['X_init'][()], t['M_init'][()], t['Z'][()], t['M_preSN'][()]
    iNa, iAl, iMg = ELEMS.index('Na'), ELEMS.index('Al'), ELEMS.index('Mg')
    rows, ok = [], True
    for iz, im, Zp, (wr, wel), _, xw in parsec2_points(t, root):
        if Zp not in (0.001, 0.01) or Mt[im] not in (20., 40., 100.): continue
        Wp = Mt[im] - float(wr['Mfin'])
        want = {k: X[iz, k] + wel[ELEMS[k]] / Wp - xw[ELEMS[k]] for k in (iNa, iAl, iMg)}
        e = cw[iz, im, -1, :].astype(np.float64)
        for a_, nm in ((iNa, 'Na'), (iAl, 'Al')):
            have, exp = e[a_] / e[iMg], want[a_] / want[iMg]
            ok &= abs(have / exp - 1) < 0.02
            rows.append('Z=%g M=%g %s/Mg %.3f (raw %.3f, birth %.3f)' % (Zp, Mt[im], nm, have, exp, X[iz, a_] / X[iz, iMg]))
    return ok and len(rows) == 12, '; '.join(rows)


@test
def lifetime_budget_wind_plus_sn(t, w, root):
    cw, X, SN = w['elem_ej_wind_cumulative'], t['X_init'][()], t['elem_ej_SN_mass'][()]
    rt, rm, Mt, Zt = t['remnant_type'][()], t['remnant_mass'][()], t['M_init'][()], t['Z'][()]
    worst, nbad, n = (0., None), 0, 0
    for iz, im, Zp, _, (tr, tel), xw in parsec2_points(t, root):
        if int(rt[iz, im]) not in (2, 4) or tr['SNT'] not in ('CCSN', 'PPISN'): continue   # explosions with a remnant
        M = Mt[im]; Mej_p = M - float(tr['Mbar'])
        end = cw[iz, im, -1, :].astype(np.float64)
        n += 1
        want = np.array([tel[e] - xw[e] * Mej_p for e in ELEMS])   # PARSEC v2 lifetime net yield
        want[0] -= want.sum()            # the source's sum of net yields (Li/Be, rounding) is taken up by H
        for k, e in enumerate(ELEMS):
            have = end[k] + SN[iz, im, k] - X[iz, k] * (M - rm[iz, im])
            want_k = want[k]
            err = abs(have - want_k) / (1e-4 * (abs(want_k) + X[iz, k] * M) + 1e-9)
            if err > 1: nbad += 1
            if err > worst[0]: worst = (err, (Zt[iz], M, e, have, want_k))
    if not n: return False, 'no CCSN/PPISN PARSEC v2 grid points'
    return nbad == 0, '%d stars x 27 elements; worst |dnet|/tol = %.2g at %s%s' % (n, worst[0], worst[1], '' if not nbad else '; %d off' % nbad)


def main():
    b = sys.argv[1]; root = os.path.join(b, 'data_root')
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
