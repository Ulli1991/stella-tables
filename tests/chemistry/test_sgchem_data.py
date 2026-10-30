#!/usr/bin/env python3
"""SGChem data: net-18 and chem-tables binaries are byte-identical to the pinned hashes (the ones GIZMO checks),
the chem-tables binary round-trips to the legacy header, and the frozen reference data match data.sha256.
  test_sgchem_data.py BUILD_DIR [--gizmo GIZMO_SRC]
--gizmo: also compare with GIZMO's own entry in galaxy_sf/stella_tables_check.h and its chem_tables_data.h."""
import hashlib, os, re, struct, sys, argparse

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SG = os.path.join(REPO, 'chemistry', 'sgchem')
fails = []


def check(ok, msg):
    print(('PASS ' if ok else 'FAIL ') + msg)
    if not ok: fails.append(msg)


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()


def read_bin(p):
    d = open(p, 'rb').read(); assert d[:8] == b'SGCT0001'
    n, = struct.unpack_from('<Q', d, 8); o = 16; out = {}
    for _ in range(n):
        name = d[o:o + 48].rstrip(b'\0').decode(); m, = struct.unpack_from('<Q', d, o + 48); o += 56
        out[name] = struct.unpack_from('<%dd' % m, d, o); o += 8 * m
    assert o == len(d)
    return out


def header_arrays(p):
    t = open(p).read(); out = {}
    for m in re.finditer(r'static const double (\w+)\[(\d+)\] = \{(.*?)\};', t, re.S):
        v = [float(x) for x in m.group(3).replace('\n', ' ').split(',') if x.strip()]
        assert len(v) == int(m.group(2)); out[m.group(1)] = tuple(v)
    return out


ap = argparse.ArgumentParser(); ap.add_argument('build'); ap.add_argument('--gizmo')
a = ap.parse_args(); ch = os.path.join(os.path.abspath(a.build), 'chemistry')
pinned = {l.split()[1]: l.split()[0] for l in open(os.path.join(SG, 'expected_sha256.txt')) if l.strip()}
for f, want in pinned.items():
    check(sha(os.path.join(ch, f)) == want, '%s equals pinned sha256' % f)
n18 = os.path.getsize(os.path.join(ch, 'sgchem_net18_tables.bin'))
check(n18 == 8 * 75321, 'net18 bin size = 75321 float64 (%d bytes)' % n18)
b = read_bin(os.path.join(ch, 'sgchem_chem_tables.bin'))
check(b == header_arrays(os.path.join(ch, 'chem_tables_data.h')), 'chem tables bin == legacy header, %d arrays' % len(b))
# frozen reference data vs data.sha256
rd = os.path.join(SG, 'reference_data'); bad = 0; n = 0
for l in open(os.path.join(rd, 'data.sha256')):
    if not l.strip(): continue
    h, p = l.split(None, 1); p = p.strip(); n += 1
    if not os.path.exists(os.path.join(rd, p)) or sha(os.path.join(rd, p)) != h: bad += 1; print('  mismatch', p)
check(bad == 0 and n > 0, 'reference_data matches data.sha256 (%d files)' % n)
listed = {l.split(None, 1)[1].strip() for l in open(os.path.join(rd, 'data.sha256')) if l.strip()}
found = {os.path.relpath(os.path.join(r, f), rd) for r, _, fs in os.walk(rd) for f in fs if f != 'data.sha256'}
check({os.path.normpath(x) for x in listed} == found, 'reference_data has no unlisted or missing files')
if a.gizmo:
    hdr = open(os.path.join(a.gizmo, 'galaxy_sf', 'stella_tables_check.h')).read()
    m = re.search(r'\{"sgchem_net18_tables.bin",\s*"([0-9a-f]{64})"\}', hdr)
    check(bool(m) and m.group(1) == pinned['sgchem_net18_tables.bin'], 'net18 hash equals GIZMO stella_tables_check.h')
    g = os.path.join(a.gizmo, 'cooling', 'sgchem', 'cxx', 'chem_tables_data.h')
    if os.path.exists(g): check(sha(g) == pinned['chem_tables_data.h'], 'chem_tables_data.h equals GIZMO source copy')
    gd = os.path.join(a.gizmo, 'cooling', 'sgchem', 'cxx', 'fortran_reference', 'data.sha256')
    if os.path.exists(gd): check(sha(gd) == sha(os.path.join(rd, 'data.sha256')), 'data.sha256 equals GIZMO copy')
sys.exit(1 if fails else 0)
