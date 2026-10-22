#!/usr/bin/env python3
"""Solar composition: A09 transcription, normalization, table X_init, and GIZMO's init.cc array all agree.

  test_solar.py BUILD_DIR [--gizmo GIZMO_SRC]      exit 0 = all PASS
"""
import os, re, sys, argparse
import numpy as np, h5py

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'stellar'))
import solar_abundances as sa

TESTS = []
def test(fn): TESTS.append(fn); return fn


@test
def a09_transcription(a):
    """log eps table reproduces A09's quoted present-day Z = 0.0134 and Z/X = 0.0181 (Sect. 3.12)"""
    X, Y, Z = sa.a09_xyz()
    return abs(Z - 0.0134) < 5e-5 and abs(Z / X - 0.0181) < 5e-5, 'X=%.4f Y=%.4f Z=%.5f Z/X=%.5f' % (X, Y, Z, Z / X)


@test
def pattern_normalization(a):
    p = sa.solar_pattern(); zs = sum(p[e] for e in sa.ELEMENTS[2:])
    ok = abs(zs - sa.Z_SUN) < 1e-15 and abs(sum(p.values()) - 1) < 1e-15 and p['He'] == sa.Y_SUN
    return ok, 'metals sum %.15f, all sum %.15f, He %.4f' % (zs, sum(p.values()), p['He'])


@test
def table_x_init_is_pattern(a):
    """X_init[Z] = pattern x Z/Z_sun for metals (sums to Z), solar He, H = remainder"""
    p = sa.solar_pattern()
    with h5py.File(os.path.join(a.build, 'stellar_tables_unified.hdf5'), 'r') as f:
        X, Z = f['X_init'][()], f['Z'][()]
    ref = np.array([[p[e] * z / sa.Z_SUN for e in sa.ELEMENTS] for z in Z])
    ref[:, 1] = sa.Y_SUN; ref[:, 0] = 1 - sa.Y_SUN - Z
    dev = np.max(np.abs(X - ref) / np.maximum(ref, 1e-30)); zdev = np.max(np.abs(X[:, 2:].sum(1) / Z - 1))
    return dev < 1e-12 and zdev < 1e-12, 'max rel dev from pattern %.1e, max |sum metals/Z - 1| %.1e' % (dev, zdev)


@test
def gizmo_init_array(a):
    """GIZMO core/init.cc solar_elem[27] equals the generator (needs --gizmo)"""
    if not a.gizmo: return True, 'skipped (no --gizmo)'
    src = open(os.path.join(a.gizmo, 'core', 'init.cc')).read()
    m = re.search(r'static const double solar_elem\[27\]\s*=\s*\{(.*?)\};', src, re.S)
    if not m: return False, 'solar_elem[27] not found in core/init.cc'
    vals = [float(v) for v in re.findall(r'([0-9.]+e[-+][0-9]+)', re.sub(r'/\*.*?\*/', '', m.group(1), flags=re.S))]
    p = sa.solar_pattern(); ref = [p[e] for e in sa.ELEMENTS]
    if len(vals) != 27: return False, 'found %d values' % len(vals)
    dev = max(abs(v / r - 1) for v, r in zip(vals, ref))
    return dev < 1e-9, 'max rel dev init.cc vs generator %.1e' % dev


@test
def no_stray_solar_constants(a):
    """builder scripts take Z_sun from solar_abundances.py: no literal 0.0134 / 0.0142 or division by a literal 0.014"""
    pat = re.compile(r'0\.0134\b|0\.0142\b|\s/\s*0\.014\b|log10\(0\.014\)')
    sdir = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'stellar'); hits = []
    for fn in sorted(os.listdir(sdir)):
        if fn.endswith('.py') and fn != 'solar_abundances.py':
            for i, line in enumerate(open(os.path.join(sdir, fn)), 1):
                if pat.search(line.split('#')[0]): hits.append('%s:%d' % (fn, i))
    return not hits, 'stray solar constants: %s' % (hits or 'none')


@test
def lut_metallicity_on_z_sun(a):
    """spectral LUT [M/H] maps to Z with the tables' Z_sun (builder lookup and LUT generator)"""
    import h5py as _h, build_stellar_tables as B, build_spectral_lut as S
    with _h.File(B.SPECTRAL_LUT_PATH, 'r') as f: mh = f['metallicity_MH'][()]
    lz = B.load_spectral_lut()[2]
    d1, d2 = np.max(np.abs(lz - np.log10(sa.Z_SUN) - mh)), abs(S.mh_to_Z(0.0) / sa.Z_SUN - 1)
    return max(d1, d2) < 1e-12, 'builder max |log Z - log Z_sun - [M/H]| = %.1e dex, generator |Z([M/H]=0)/Z_sun - 1| = %.1e' % (d1, d2)


@test
def mist_grid_z_matches_headers(a):
    """MIST grid Z used for log-Z interpolation equals the Zinit in every (non-empty) MIST track header"""
    import glob, build_stellar_tables as B
    root = os.path.join(a.build, 'data_root', 'mist'); worst = (0., None); n = 0
    for key, z in B.MIST_Z_GRID.items():
        for fn in glob.glob(os.path.join(root, 'MIST_v1.2_feh_%s_*' % key, '*.track.eep')):
            if os.path.getsize(fn) == 0: continue
            with open(fn) as fh: zh = float([next(fh) for _ in range(5)][4].split()[2])
            n += 1; d = abs(z / zh - 1)
            if d > worst[0]: worst = (d, os.path.relpath(fn, root))
    return n > 0 and worst[0] < 1e-5, '%d tracks, max rel dev %.1e %s' % (n, worst[0], worst[1] or '')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('build'); ap.add_argument('--gizmo'); a = ap.parse_args()
    nf = 0
    for fn in TESTS:
        try: ok, msg = fn(a)
        except Exception as e: ok, msg = False, 'EXCEPTION %r' % e
        nf += not ok; print('  [%s] %s: %s' % ('PASS' if ok else 'FAIL', fn.__name__, msg))
    print('\n%d/%d failed' % (nf, len(TESTS))); sys.exit(1 if nf else 0)


if __name__ == '__main__':
    main()
