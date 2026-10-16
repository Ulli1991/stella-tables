#!/usr/bin/env python3
"""Reduce the full stellar-track catalogues to exactly what build_stellar_tables.py reads.

  reduce_tracks.py RAW_ROOT OUT_ROOT [--manifest stellar/data/reduced/MANIFEST.txt]

RAW_ROOT holds mist/, parsec2_vms/, parsec1.2/, boost/ as downloaded (see fetch/). Every file in the
manifest is copied to OUT_ROOT under the same relative path. MIST *.track.eep and PARSEC v2 *.TAB tracks
keep their header lines but only the columns whose names appear as string literals in
build_stellar_tables.py, as the original text tokens, so the builder parses identical numbers; rows the
builder would skip (wrong column count) are dropped.
All other files are copied verbatim. sha256 of each original is checked against the manifest.
"""
import os, re, sys, hashlib, shutil, argparse

# every quoted identifier in the builder: a header column is kept iff the builder can name it
BUILDER_NAMES = set(re.findall(r"""['"]([A-Za-z0-9_]+)['"]""",
                               open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'build_stellar_tables.py')).read()))


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''): h.update(chunk)
    return h.hexdigest()


def reduce_mist(src, dst):
    names, idx, out = None, None, []
    for line in open(src):
        if line.startswith('#'):
            if 'star_age' in line:
                names = line.strip().lstrip('#').split()
                idx = [i for i, c in enumerate(names) if c in BUILDER_NAMES]
                line = '#' + ' '.join(names[i] for i in idx) + '\n'
            out.append(line); continue
        vals = line.split()
        if names and len(vals) == len(names): out.append(' '.join(vals[i] for i in idx) + '\n')
    open(dst, 'w').writelines(out)


def reduce_parsec2(src, dst):
    lines = open(src).readlines()
    names = lines[2].split(); idx = [i for i, c in enumerate(names) if c in BUILDER_NAMES]
    out = lines[:2] + [' '.join(names[i] for i in idx) + '\n']
    for line in lines[3:]:
        vals = line.split()
        if len(vals) == len(names): out.append(' '.join(vals[i] for i in idx) + '\n')
    open(dst, 'w').writelines(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('raw_root'); ap.add_argument('out_root')
    ap.add_argument('--manifest', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'reduced', 'MANIFEST.txt'))
    a = ap.parse_args()
    bad = 0
    for line in open(a.manifest):
        if not line.strip() or line.startswith('#'): continue
        digest, rel = line.split(None, 1); rel = rel.strip()
        src, dst = os.path.join(a.raw_root, rel), os.path.join(a.out_root, rel)
        if digest == 'placeholder':  # discovered by filename only: keep the name so track discovery is unchanged
            os.makedirs(os.path.dirname(dst), exist_ok=True); open(dst, 'w').close(); continue
        if sha256(src) != digest: print('CHECKSUM MISMATCH', rel); bad += 1; continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if rel.startswith('mist/') and rel.endswith('.track.eep'): reduce_mist(src, dst)
        elif rel.startswith('parsec2_vms/') and rel.endswith('.TAB'): reduce_parsec2(src, dst)
        else: shutil.copy2(src, dst)
    if bad: sys.exit('%d input files differ from the manifest' % bad)


if __name__ == '__main__':
    main()
