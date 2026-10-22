#!/usr/bin/env python3
"""Download the hot-star atmosphere grids behind stellar/data/lut/hotstar_lut.hdf5 into RAW_ROOT/hotstar.

  fetch_hotstar.py RAW_ROOT              fetch what is missing, then write RAW_ROOT/hotstar/MANIFEST.txt
  fetch_hotstar.py RAW_ROOT --manifest   only (re)write the manifest

svo/        TLUSTY OSTAR2002 (Lanz & Hubeny 2003, ApJS 146, 417) and BSTAR2006 (Lanz & Hubeny 2007, ApJS 169, 83),
            every spectrum listed by the SVO theory-server SSAP service (tlusty_ostarbin, tlusty_bstarbin)
powr/       PoWR (https://www.astro.physik.uni-potsdam.de/PoWR/): OB Vd3 grids GAL/LMC/SMC/subSMC (Hainich+2019,
            Pauli+2026; recomputed after the 2025-12-15 bug fix) parameters and log Q tables; WR grids MW/LMC/SMC/0.07 Zsun
            WNE, WNL-Hxx, WC (Todt+2015, Hamann+2006, Sander+2012, Hainich+2014/2015) parameters, log Q tables where
            served, model lists, every SED (text) and model header
gotberg18/  Gotberg+2018 (A&A 615, A78; CDS J/A+A/615/A78) SED.txt of every model, ReadMe, and the arXiv:1802.03018
            source (parameter tables)
The servers offer only their current versions; the manifest pins what was used (date, sha256).
"""
import os, re, sys, glob, time, hashlib, subprocess, datetime

SVO = 'http://svo2.cab.inta-csic.es/theory/newov2/ssap.php?model=%s&'
POWR = 'https://www.astro.physik.uni-potsdam.de/PoWR/'
CDS = 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/615/A78'
OB = ['gal-ob-vd3', 'lmc-ob-vd3', 'smc-ob-vd3', 'subsmc-ob-vd3']
WR = ['wne', 'wnl', 'wnl-h50', 'wc', 'lmc-wne', 'lmc-wnl-h20', 'lmc-wnl-h40', 'lmc-wc', 'smc-wne', 'smc-wnl-h20',
      'smc-wnl-h40', 'smc-wnl-h60', 'smc-wc-2021', '007-wne-2015', '007-wnl-h20-2015', '007-wnl-h40-2015',
      '007-wnl-h60-2015', '007-wc-2021']
CITE = {'svo/tlusty_ostarbin': 'Lanz & Hubeny 2003, ApJS 146, 417 (via SVO theory server)',
        'svo/tlusty_bstarbin': 'Lanz & Hubeny 2007, ApJS 169, 83 (via SVO theory server)',
        'powr/gal-ob-vd3': 'Pauli et al. 2026 (arXiv:2608.25003); PoWR', 'powr/lmc-ob-vd3': 'Pauli et al. 2026; PoWR',
        'powr/smc-ob-vd3': 'Pauli et al. 2026; PoWR', 'powr/subsmc-ob-vd3': 'Pauli et al. 2026; PoWR',
        'powr/': 'Hamann & Graefener 2004; Todt et al. 2015; Sander et al. 2012; Hainich et al. 2014, 2015 (PoWR WR grids)',
        'gotberg18/': 'Gotberg et al. 2018, A&A 615, A78'}


def curl(out, url, data=None):
    if os.path.exists(out) and os.path.getsize(out) > 0: return
    os.makedirs(os.path.dirname(out), exist_ok=True)
    for _ in range(3):
        subprocess.call(['curl', '-s', '-m', '600', '-o', out] + (['-X', 'POST', '-d', data] if data else []) + [url])
        if os.path.exists(out) and os.path.getsize(out) > 0: return
        time.sleep(5)


def fetch(root):
    for m, key in (('tlusty_ostarbin', 'ostar'), ('tlusty_bstarbin', 'bstar')):
        lst = os.path.join(root, 'svo', key + '_ssap.xml'); curl(lst, SVO % m)
        for u in re.findall(r'CDATA\[(.*?)\]', open(lst).read()):
            curl(os.path.join(root, 'svo', m, 'fid%s.xml' % u.split('fid=')[1]), u)
    for g in OB + WR:
        for f in ('parameters.txt', 'logQ_HII.txt', 'logQ_HeI.txt', 'logQ_HeII.txt'):
            out = os.path.join(root, 'powr', g, f); curl(out, POWR + 'isoplots/%s/%s' % (g, f))
            if os.path.exists(out) and '<html' in open(out, errors='ignore').read(500).lower(): os.rename(out, out + '.notfound.html')
    for g in WR:
        lst = os.path.join(root, 'powr', g, 'modellist.html'); curl(lst, POWR + 'powrgrid2.php', 'type=wr&grid=%s&Submit=Next' % g)
        mods = re.findall(r'id="modlink([0-9-]+)"[^>]*\s*title="T= ([0-9.]+) kK&#013;log Rt= ([-0-9.]+)"', open(lst).read())
        with open(os.path.join(root, 'powr', g, 'models.txt'), 'w') as fh:
            for m, T, R in mods: fh.write('%s %s %s\n' % (m, T, R))
        for m, _, _ in mods:
            curl(os.path.join(root, 'powr', g, 'sed', m + '.txt'), POWR + 'download.php',
                 'content=specdl&grid=%s&model=%s&datatype=sed&dlformat=txt&specdownload=Download' % (g, m))
            curl(os.path.join(root, 'powr', g, 'sed', m + '.modinfo'), POWR + 'download.php?content=modinfo&grid=%s&model=%s' % (g, m))
    gb = os.path.join(root, 'gotberg18'); curl(os.path.join(gb, 'ReadMe'), CDS + '/ReadMe')
    for g in ('grid_014', 'grid_006', 'grid_002', 'grid_0002'):
        idx = os.path.join(gb, g, 'index.html'); curl(idx, CDS + '/%s/' % g)
        for m in sorted(set(re.findall(r'href="(M1_[^"/]*)', open(idx).read()))):
            curl(os.path.join(gb, g, m, 'SED.txt'), CDS + '/%s/%s/SED.txt' % (g, m))
    src = os.path.join(gb, 'paper', 'src.tgz'); curl(src, 'https://arxiv.org/e-print/1802.03018')
    os.makedirs(os.path.join(gb, 'paper', 'x'), exist_ok=True)
    subprocess.call(['tar', '-xzf', src, '-C', os.path.join(gb, 'paper', 'x')])


def manifest(root):
    lines = ['# hot-star grids: path  sha256  bytes  mtime  source', '# written %s' % datetime.date.today().isoformat()]
    for f in sorted(glob.glob(os.path.join(root, '**', '*'), recursive=True)):
        if not os.path.isfile(f) or f.endswith('MANIFEST.txt'): continue
        rel = os.path.relpath(f, root); h = hashlib.sha256(open(f, 'rb').read()).hexdigest()
        cite = next((v for k, v in sorted(CITE.items(), key=lambda kv: -len(kv[0])) if rel.startswith(k)), '')
        lines.append('%s  %s  %d  %s  %s' % (rel, h, os.path.getsize(f),
                                             datetime.datetime.fromtimestamp(os.path.getmtime(f)).strftime('%Y-%m-%d'), cite))
    open(os.path.join(root, 'MANIFEST.txt'), 'w').write('\n'.join(lines) + '\n')
    print('MANIFEST: %d files' % (len(lines) - 2))


if __name__ == '__main__':
    root = os.path.join(sys.argv[1], 'hotstar')
    if '--manifest' not in sys.argv: fetch(root)
    manifest(root)
