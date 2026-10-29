#!/usr/bin/env python3
"""Round trip of the built energy-table-x.txt through GIZMO's own lookups (Fortran calc_temp.F and the C++ port).

  run_lookup_test.py GIZMO_SRC BUILD_DIR [WORK_DIR]

GIZMO_SRC  GIZMO source tree (only cooling/sgchem/{calc_temp.F,*.h,cxx/} are read; cooling/chemcool before the rename)
BUILD_DIR  output of build.sh (chemistry/energy-table-x.txt, chemistry/eh2-table.txt)
Exact T -> ekn via ekn = x e_H2(T) + (1-x) 1.5 T, x = abh2/(1+abhe-abh2); each backend must return T to <2e-3
for every abhe (0, 0.079, 0.1, 0.12), and the two backends must agree. The same points go through the original GIZMO
table (chemistry/reference/): results must be identical off the 20000 K seam, and the built gamma bounded on it.
Nothing is written to GIZMO_SRC. Exit 0 = PASS.
"""
import os, lzma, shutil, subprocess, sys, tempfile, importlib.util
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
spec = importlib.util.spec_from_file_location('m', os.path.join(REPO, 'chemistry', 'make_energy_table.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
TOL = 2e-3
MODS = 'module purge && module load gcc openmpi gsl hwloc fftw-mpi hdf5-mpi'
CONFIG = '#define %s\n#define CHEMISTRYNETWORK 5\n#define WSS_CIE_COOL\n'   # production chemistry (option name per tree)


def sh(cmd, cwd):
    r = subprocess.run(['bash', '-lc', '%s && %s' % (MODS, cmd)], cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    if r.returncode: print(r.stdout); sys.exit('FAILED: ' + cmd)


def backends(src, table, eh2_file, work, pts):
    """T, gamma from GIZMO's Fortran and C++ lookups with TABLE staged where calc_temp.F opens it"""
    name = 'sgchem' if os.path.isdir(os.path.join(src, 'cooling', 'sgchem')) else 'chemcool'   # CHEMCOOL -> SGCHEM rename, 2026-10
    cc_src, cc = os.path.join(src, 'cooling', name), os.path.join(work, 'gizmo', 'cooling', name)
    os.makedirs(cc, exist_ok=True)
    open(os.path.join(work, 'gizmo', 'GIZMO_config.h'), 'w').write(CONFIG % name.upper())
    for f in ('calc_temp.F', 'cool.h', 'gamma.h', name + '_consts.h', 'f2c.h'): shutil.copy(os.path.join(cc_src, f), cc)
    shutil.copy(table, os.path.join(cc, 'energy-table-x.txt')); shutil.copy(eh2_file, os.path.join(cc, 'eh2-table.txt'))
    shutil.copy(os.path.join(HERE, 'test_lookup.F'), work)
    sh('mpifort -cpp -O2 -I%s -I%s test_lookup.F %s/calc_temp.F -o test_lookup.x' % (os.path.join(work, 'gizmo'), cc, cc), work)
    sh('g++ -std=c++17 -O2 -I%s/cxx %s %s/cxx/chem_calc_temp.cc -o test_lookup_cxx.x' % (cc_src, os.path.join(HERE, 'test_lookup_cxx.cc'), cc_src), work)
    np.savetxt(os.path.join(work, 'points.txt'), pts[:, :3], fmt='%.16e')
    sh('./test_lookup.x', work)
    sh('./test_lookup_cxx.x gizmo/cooling/%s/energy-table-x.txt points.txt > cxx_out.txt' % name, work)
    fo, co = np.loadtxt(os.path.join(work, 'lookup_out.txt')), np.loadtxt(os.path.join(work, 'cxx_out.txt'))
    return {'fortran': (fo[:, 3], fo[:, 4]), 'c++': (co[:, 0], co[:, 1])}


def main():
    src, build = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    work = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else tempfile.mkdtemp(prefix='lookup_')
    eh2_file = os.path.join(build, 'chemistry', 'eh2-table.txt')
    eh2 = m.load_eh2(eh2_file)
    rng = np.random.RandomState(7); pts = []
    for abhe in (0.0, 0.079, 0.1, 0.12):
        abh2 = rng.uniform(0.005, 0.5, 4000); T = 10 ** rng.uniform(1, np.log10(15000), 4000)
        x = abh2 / (1 + abhe - abh2)
        ekn = x * eh2(T) + (1 - x) * 1.5 * T
        keep = ekn < 4.99e4
        pts += [(a, abhe, e, t) for a, e, t in zip(abh2[keep], ekn[keep], T[keep])]
    for x in np.linspace(0.02, 0.6, 30):          # seam probes: T(ekn) across 19900..20100 K at abhe = 0.1
        for T in np.linspace(19900, 20100, 41):
            pts.append((x * 1.1 / (1 + x), 0.1, x * eh2(T) + (1 - x) * 1.5 * T, T))
    pts = np.array(pts); seam = (pts[:, 3] > 19800) & (pts[:, 3] < 20200)

    os.makedirs(work, exist_ok=True); orig = os.path.join(work, 'energy-table-x.gizmo-original.txt')
    with lzma.open(os.path.join(REPO, 'chemistry', 'reference', 'energy-table-x.gizmo-original.txt.xz')) as fi, open(orig, 'wb') as fo:
        shutil.copyfileobj(fi, fo)
    new = backends(src, os.path.join(build, 'chemistry', 'energy-table-x.txt'), eh2_file, os.path.join(work, 'built'), pts)
    old = backends(src, orig, eh2_file, os.path.join(work, 'original'), pts)
    ok = True
    for tab, res in (('built', new), ('original', old)):
        for name in ('fortran', 'c++'):
            err = abs(res[name][0] / pts[:, 3] - 1)
            for abhe in (0.0, 0.079, 0.1, 0.12):
                s = (pts[:, 1] == abhe) & ~seam; e = err[s]; p = e.max() < TOL; ok &= p
                print('  %-8s %-7s abhe=%-5s n=%d median reldT=%.1e p99=%.1e max=%.1e  %s' % (tab, name, abhe, s.sum(), np.median(e), np.percentile(e, 99), e.max(), 'PASS' if p else 'FAIL'))
        par = max(np.max(abs(res['c++'][0] / res['fortran'][0] - 1)), np.max(abs(res['c++'][1] / res['fortran'][1] - 1)))
        ok &= par < 1e-5
        print('  %-8s fortran vs c++ max rel diff (T, gamma) = %.1e  %s' % (tab, par, 'PASS' if par < 1e-5 else 'FAIL'))
    for name in ('fortran', 'c++'):   # built vs original through GIZMO: identical off the seam, bounded gamma on it
        dT = abs(new[name][0] / old[name][0] - 1); dg = abs(new[name][1] - old[name][1])
        same = dT[~seam].max() == 0 and dg[~seam].max() == 0
        gn, go = new[name][1][seam], old[name][1][seam]
        bnd = gn.min() >= 1 + 1 / 4.5 - 1e-9 and gn.max() <= 5 / 3. + 1e-9
        ok &= same and bnd
        print('  %-7s built vs original: off-seam identical %s (n=%d); seam (n=%d) max reldT %.1e, gamma built %.3f..%.3f vs original %.3f..%.3f  %s' % (
            name, same, (~seam).sum(), seam.sum(), dT[seam].max(), gn.min(), gn.max(), go.min(), go.max(), 'PASS' if same and bnd else 'FAIL'))
    print('energy-table lookup round trip:', 'PASS' if ok else 'FAIL'); sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
