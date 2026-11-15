#!/usr/bin/env python3
"""Compare the running Python environment with requirements-build.txt; print the mismatches, exit 1 if any."""
import os, sys, platform
want = {}
for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'requirements-build.txt')):
    line = line.split('#')[0].strip()
    if '==' in line:
        k, v = line.split('=='); want[k.strip()] = v.strip()
have = {'python': platform.python_version()}
for mod in ('numpy', 'scipy', 'h5py'):
    try: have[mod] = __import__(mod).__version__
    except ImportError: have[mod] = 'missing'
bad = ['%s %s (need %s)' % (k, have.get(k), v) for k, v in want.items() if have.get(k) != v]
print('build environment: ' + ', '.join('%s %s' % kv for kv in have.items()) + ('' if not bad else '; differs: ' + ', '.join(bad)))
sys.exit(1 if bad else 0)
