#!/usr/bin/env python3
"""Run every table test against one build directory.

  run_all.py BUILD_DIR [--gizmo GIZMO_SRC]

BUILD_DIR   output of build.sh (stellar_tables_unified.hdf5, stellar_wind_cumulative.hdf5, chemistry/)
--gizmo     also round-trip the energy table through GIZMO's Fortran and C++ lookups (needs the MPCDF modules)
Exit 0 = every suite passed.
"""
import os, sys, subprocess, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('build'); ap.add_argument('--gizmo')
    a = ap.parse_args(); b = os.path.abspath(a.build)
    table, wind = os.path.join(b, 'stellar_tables_unified.hdf5'), os.path.join(b, 'stellar_wind_cumulative.hdf5')
    env = dict(os.environ, STELLA_TABLE=table, STELLA_WIND=wind,
               PYTHONPATH=os.pathsep.join(p for p in (os.path.join(REPO, 'stellar'), os.environ.get('PYTHONPATH')) if p))
    py = sys.executable
    suites = [
        ('stellar: invariants', [py, os.path.join(HERE, 'stellar', 'test_invariants.py')]),
        ('stellar: physics/schema', [py, os.path.join(HERE, 'stellar', 'test_table_physics.py'), '--table', table]),
        ('stellar: extended', [py, os.path.join(HERE, 'stellar', 'test_extended.py'), b]),
        ('chemistry: tables', [py, os.path.join(HERE, 'chemistry', 'test_chemistry.py'), b]),
    ]
    if a.gizmo:
        suites.append(('chemistry: GIZMO lookups', [py, os.path.join(HERE, 'chemistry', 'run_lookup_test.py'), a.gizmo, b,
                                                    os.path.join(b, 'lookup_test')]))
    missing = [p for p in (table, wind) if not os.path.exists(p)]
    if missing: sys.exit('missing build outputs: %s' % missing)
    res = []
    for name, cmd in suites:
        print('\n==== %s' % name); sys.stdout.flush()
        res.append((name, subprocess.call(cmd, env=env, cwd=HERE) == 0))
    print('\n==== summary')
    for name, ok in res: print('  %-28s %s' % (name, 'PASS' if ok else 'FAIL'))
    if not a.gizmo: print('  %-28s skipped (no --gizmo)' % 'chemistry: GIZMO lookups')
    sys.exit(0 if all(ok for _, ok in res) else 1)


if __name__ == '__main__':
    main()
