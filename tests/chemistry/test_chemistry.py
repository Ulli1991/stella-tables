#!/usr/bin/env python3
"""Checks on the built SGCHEM tables (no GIZMO tree needed).

  test_chemistry.py BUILD_DIR      exit 0 = all PASS (expected failures do not fail the run)

eh2-table.txt and H2-cooling-ratios.dat are inherited tables of unknown provenance; wss_z_collis.txt is the
Wiersma, Schaye & Smith (2009) CIE table (official release). All three are checksum-locked here so any change is
deliberate, and sanity-checked. energy-table-x.txt is regenerated from eh2.
"""
import os, sys, hashlib, importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
spec = importlib.util.spec_from_file_location('m', os.path.join(REPO, 'chemistry', 'make_energy_table.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
SHA = {   # sha256 of the built tables (energy-table-x.txt = seam-split generator output)
    'eh2-table.txt': 'd423c682d6e9036edc9100e9e2ee990084ee191d500289dc485ce0a486f4450b',
    'H2-cooling-ratios.dat': '13cf3c756a6b9e4e993ba9fb5fa93883e1dedd788d6415149497ef1fee893a0d',
    'wss_z_collis.txt': 'eaa290e0895be4846147e99880f18e5b981d9293cde632e4b3c4be71f10f8372',
    'energy-table-x.txt': 'e1b03445ee95b020e611f73d3d84ca395378c607517e3367ceb18d95a712b96d',
}
ORIG = os.path.join(REPO, 'chemistry', 'reference', 'energy-table-x.gizmo-original.txt.xz')   # table GIZMO shipped before the seam fix
ORIG_SHA = 'bde7be8337a1ddc6f92b2da6ee43b051c9221112a7ef3d250fa48734e7e72983'
NBLK, NEKN = 21, 50000   # gamma.h: netab_x = 21, 50000 rows per block
SEAM = (19995, 20005)    # eh2 regime switch, see eh2_seam_20000K
TESTS = []
def test(xfail=None):
    def deco(fn): TESTS.append((fn, xfail)); return fn
    return deco


def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()


@test()
def checksums(d):
    bad = [f for f, h in SHA.items() if sha(os.path.join(d, f)) != h]
    return not bad, 'all %d tables match the locked sha256' % len(SHA) if not bad else 'changed: %s' % bad


@test()
def eh2_grid_and_finite(d):
    e = np.loadtxt(os.path.join(d, 'eh2-table.txt'))
    ok = e.shape == (50000, 2) and np.array_equal(e[:, 0], np.arange(1, 50001)) and np.all(np.isfinite(e[:, 1]))
    return ok, 'shape %s, T = 1..50000 step 1, all finite' % (e.shape,)


@test()
def eh2_low_T_translational(d):
    """below ~10 K the H2 rotational levels are frozen: e_H2 = 1.5 T"""
    e = np.loadtxt(os.path.join(d, 'eh2-table.txt')); s = e[:, 0] <= 10
    dev = np.max(np.abs(e[s, 1] / (1.5 * e[s, 0]) - 1)); return dev < 1e-4, 'max |e/1.5T - 1| (T<=10 K) = %.1e' % dev


@test()
def eh2_monotonic_heat_capacity_bounds(d):
    """e_H2 rises with T and its heat capacity stays in [1.5, 4.5] (away from the 20000 K seam)"""
    e = np.loadtxt(os.path.join(d, 'eh2-table.txt')); s = (e[:, 0] < SEAM[0]) | (e[:, 0] > SEAM[1])
    cv = np.diff(e[s, 1]) / np.diff(e[s, 0]); cv = cv[np.diff(e[s, 0]) == 1]
    return np.all(cv > 0) and cv.min() > 1.5 - 1e-3 and cv.max() < 4.5, 'c_v in [%.3f, %.3f] outside %s K' % (cv.min(), cv.max(), SEAM)


@test()
def eh2_seam_20000K(d):
    """inherited eh2 jumps ~45 K at 20000 K (partition function -> analytic 3T); the built table keeps gamma in [1+1/4.5, 5/3] there"""
    e = np.loadtxt(os.path.join(d, 'eh2-table.txt')); t = etab(d); g, T = t[:, :, 2], t[:, :, 1]
    step = e[e[:, 0] == 20001, 1][0] - e[e[:, 0] == 19999, 1][0]
    w = (T > SEAM[0]) & (T < SEAM[1]); lo, hi = g[w].min(), g[w].max()
    return lo >= 1 + 1 / 4.5 - 1e-9 and hi <= 5 / 3. + 1e-9, 'eh2 step e(20001)-e(19999) = %.1f K; built gamma in [%.4f, %.4f] in %s K' % (step, lo, hi, SEAM)


def orig_tab():
    import lzma
    return np.loadtxt(lzma.open(ORIG, 'rt')).reshape(NBLK, NEKN, 3)


@test()
def original_table_preserved(d):
    import lzma
    h = hashlib.sha256(lzma.open(ORIG).read()).hexdigest()
    return h == ORIG_SHA, 'stored original GIZMO energy table sha256 %s' % ('intact' if h == ORIG_SHA else 'CHANGED ' + h)


@test()
def legacy_seam_reproduces_original(d):
    tmp = os.path.join(d, '.regen_legacy.txt')
    np.savetxt(tmp, m.make(m.load_eh2(os.path.join(d, 'eh2-table.txt'), legacy_seam=True), np.linspace(0, 1, NBLK)), fmt=['%d', '%.14g', '%.14g'], delimiter='  ')
    ok = sha(tmp) == ORIG_SHA; os.remove(tmp)
    return ok, '--legacy-seam output byte-identical to the original GIZMO table'


@test()
def new_vs_original_only_seam_changes(d):
    """the seam fix changes rows only near 20000 K; everywhere else the built table equals the original exactly"""
    o, n = orig_tab(), etab(d); ch = np.any(o != n, axis=2); T = o[:, :, 1]
    far = (T < 19900) | (T > 20100)
    ok = np.array_equal(o[:, :, 0], n[:, :, 0]) and not np.any(ch & far)
    return ok, '%d of %d rows changed, all at T = %.0f..%.0f K; rows outside 19900-20100 K identical: %s' % (
        ch.sum(), ch.size, T[ch].min() if ch.any() else 0, T[ch].max() if ch.any() else 0, not np.any(ch & far))


@test()
def original_gamma_defect_documented(d):
    """the original table has gamma outside [1+1/4.5, 5/3] only inside the seam window (what the fix removes)"""
    o = orig_tab(); g, T = o[:, :, 2], o[:, :, 1]; bad = (g < 1 + 1 / 4.5 - 1e-9) | (g > 5 / 3. + 1e-9)
    return bad.any() and T[bad].min() > SEAM[0] and T[bad].max() < SEAM[1], 'original: %d out-of-bound gamma rows, gamma %.3f..%.3f, all in %s K' % (bad.sum(), g[bad].min(), g[bad].max(), SEAM)


@test()
def eh2_known_glitch_present_and_dropped(d):
    """the inherited table carries one bad row (T=20000, e=1.5T); the generator must drop exactly that row"""
    e = np.loadtxt(os.path.join(d, 'eh2-table.txt')); r = e[e[:, 0] == 20000, 1][0]
    nb = 0.5 * (e[e[:, 0] == 19999, 1][0] + e[e[:, 0] == 20001, 1][0])
    return abs(r - 30000) < 1e-6 and abs(nb / r - 1) > 0.5, 'row T=20000 has %.0f vs neighbours %.0f' % (r, nb)


@test()
def energy_table_regenerates(d):
    tmp = os.path.join(d, '.regen_energy-table-x.txt')
    e = m.load_eh2(os.path.join(d, 'eh2-table.txt'))
    np.savetxt(tmp, m.make(e, np.linspace(0, 1, NBLK)), fmt=['%d', '%.14g', '%.14g'], delimiter='  ')
    ok = sha(tmp) == sha(os.path.join(d, 'energy-table-x.txt')); os.remove(tmp)
    return ok, 'generator output byte-identical to the built energy-table-x.txt'


def etab(d):
    t = np.loadtxt(os.path.join(d, 'energy-table-x.txt')); return t.reshape(NBLK, NEKN, 3)


@test()
def energy_table_layout(d):
    t = etab(d); ok = all(np.array_equal(t[b, :, 0], np.arange(1, NEKN + 1)) for b in range(NBLK)) and np.all(np.isfinite(t))
    return ok, '%d blocks x %d rows, ekn = 1..%d in each, all finite' % (NBLK, NEKN, NEKN)


@test()
def energy_table_atomic_block(d):
    """x = 0 (no H2): T = ekn/1.5 and gamma = 5/3 exactly"""
    t = etab(d)[0]; dT = np.max(np.abs(t[:, 1] / (t[:, 0] / 1.5) - 1)); dg = np.max(np.abs(t[:, 2] - 5 / 3.))
    return dT < 1e-10 and dg < 1e-10, 'max rel T dev %.1e, max |gamma - 5/3| %.1e' % (dT, dg)


@test()
def energy_table_monotonic_and_gamma_bounds(d):
    """T rises with ekn in every block; T falls with x at fixed ekn (H2 soaks up energy); 1+1/4.5 <= gamma <= 5/3 off the seam"""
    t = etab(d); T, g = t[:, :, 1], t[:, :, 2]
    a = np.all(np.diff(T, axis=1) > 0); b = np.all(np.diff(T, axis=0) <= 1e-9 * T[1:])
    off = (T < SEAM[0]) | (T > SEAM[1]); g = g[off]
    return a and b and g.min() >= 1 + 1 / 4.5 - 1e-9 and g.max() <= 5 / 3. + 1e-9, 'gamma in [%.4f, %.4f] outside %s K; T(ekn) increasing: %s; T(x) non-increasing: %s' % (g.min(), g.max(), SEAM, a, b)


@test()
def energy_table_pure_h2_consistent(d):
    """x = 1 block inverts e_H2: e_H2(T(ekn)) = ekn"""
    t = etab(d)[-1]; e = m.load_eh2(os.path.join(d, 'eh2-table.txt'))
    s = (t[:, 1] >= 1) & (t[:, 1] <= 5e4); dev = np.max(np.abs(e(t[s, 1]) / t[s, 0] - 1))
    return dev < 1e-3, 'max |e_H2(T)/ekn - 1| = %.1e' % dev


@test()
def h2_cooling_ratios_grid(d):
    """H2-cooling-ratios.dat: full rectangular (col1, col2) grid of finite values"""
    r = np.loadtxt(os.path.join(d, 'H2-cooling-ratios.dat')); u1, u2 = np.unique(r[:, 0]), np.unique(r[:, 1])
    ok = len(r) == len(u1) * len(u2) and np.all(np.isfinite(r))
    return ok, '%d rows = %d x %d grid, col1 %.2f..%.2f, col2 %.2f..%.2f, col3 %.3g..%.3g' % (len(r), len(u1), len(u2), u1[0], u1[-1], u2[0], u2[-1], r[:, 2].min(), r[:, 2].max())


@test()
def wss_table_grid(d):
    """wss_z_collis.txt: 25 columns, finite, temperature column strictly increasing and log-spaced"""
    w = np.loadtxt(os.path.join(d, 'wss_z_collis.txt')); T = w[:, 0]; dl = np.diff(np.log10(T))
    ok = w.shape[1] == 25 and np.all(np.isfinite(w)) and np.all(dl > 0) and np.ptp(dl) < 1e-3
    return ok, 'shape %s, T %.0f..%.3g K, dlog10T = %.5f +- %.1e' % (w.shape, T[0], T[-1], dl.mean(), np.ptp(dl))


def main():
    d = os.path.join(os.path.abspath(sys.argv[1]), 'chemistry'); nf = 0
    for fn, xf in TESTS:
        try: ok, msg = fn(d)
        except Exception as e: ok, msg = False, 'EXCEPTION %r' % e
        tag = ('PASS' if ok else 'FAIL') if not xf else ('XPASS' if ok else 'XFAIL')
        if tag == 'FAIL': nf += 1
        print('  [%s] %s: %s%s' % (tag, fn.__name__, msg, ('  (expected: %s)' % xf) if xf else ''))
    print('\n%d/%d failed' % (nf, len(TESTS))); sys.exit(1 if nf else 0)


if __name__ == '__main__':
    main()
