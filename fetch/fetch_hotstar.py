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
pollux/     POLLUX (https://pollux.oreme.org, Palacios+2010) CMFGEN SED collections CMFGEN-SED (Martins & Palacios 2017, 2021),
            VMS-ZLMC-SED (Martins & Palacios 2022), VMS-ZSMC/Z0p1Zsun/Z0p01Zsun-SED (Martins+2025), each zip as served and
            unpacked; the arXiv:2505.02993 source (Martins+2025 parameter tables)
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
        'gotberg18/': 'Gotberg et al. 2018, A&A 615, A78',
        'pollux/VMS-ZLMC-SED': 'Martins & Palacios 2022, A&A 659, A163 (CMFGEN, via POLLUX: Palacios et al. 2010, A&A 516, A13)',
        'pollux/VMS-Z': 'Martins, Palacios, Schaerer & Marques-Chaves 2025, A&A 698, A262 (CMFGEN, via POLLUX: Palacios et al. 2010)',
        'pollux/CMFGEN-SED': 'Martins & Palacios 2017, 2021 (CMFGEN, via POLLUX: Palacios et al. 2010)',
        'pollux/paper_martins2025': 'Martins et al. 2025, A&A 698, A262 (arXiv:2505.02993 source)'}
POLLUX_URL = 'https://pollux.oreme.org'
POLLUX = ['VMS-ZLMC-SED', 'VMS-ZSMC-SED', 'VMS-Z0p1Zsun-SED', 'VMS-Z0p01Zsun-SED', 'CMFGEN-SED']


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
    for m in ('m014.dat', 'm006.dat', 'm002.dat', 'm0002.dat'): curl(os.path.join(gb, m), CDS + '/' + m)
    for g in ('grid_014', 'grid_006', 'grid_002', 'grid_0002'):
        idx = os.path.join(gb, g, 'index.html'); curl(idx, CDS + '/%s/' % g)
        for m in sorted(set(re.findall(r'href="(M1_[^"/]*)', open(idx).read()))):
            curl(os.path.join(gb, g, m, 'SED.txt'), CDS + '/%s/%s/SED.txt' % (g, m))
    src = os.path.join(gb, 'paper', 'src.tgz'); curl(src, 'https://arxiv.org/e-print/1802.03018')
    os.makedirs(os.path.join(gb, 'paper', 'x'), exist_ok=True)
    subprocess.call(['tar', '-xzf', src, '-C', os.path.join(gb, 'paper', 'x')])
    pollux(os.path.join(root, 'pollux'))


def pollux(d):
    """POLLUX collection zips (POST with the page's CSRF token), each unpacked to <collection>/ (one zip per SED inside)"""
    os.makedirs(d, exist_ok=True); jar = os.path.join(d, '.cookies')
    for c in POLLUX:
        z = os.path.join(d, c + '.zip')
        if not (os.path.exists(z) and os.path.getsize(z) > 0):
            page = subprocess.check_output(['curl', '-s', '-m', '120', '-c', jar, POLLUX_URL + '/collections/']).decode()
            tok = re.search(r'csrfmiddlewaretoken" value="([^"]+)', page).group(1)
            subprocess.call(['curl', '-s', '-m', '3600', '-b', jar, '-c', jar, '-e', POLLUX_URL + '/collections/', '-X', 'POST',
                             '-d', 'csrfmiddlewaretoken=' + tok, '-o', z, POLLUX_URL + '/collections/%s/download/' % c])
        x = os.path.join(d, c)
        if not os.path.isdir(x):
            os.makedirs(x); subprocess.call(['unzip', '-q', '-o', z, '-d', x])
            for inner in glob.glob(os.path.join(x, '*.zip')):
                subprocess.call(['unzip', '-q', '-o', inner, '-d', x]); os.remove(inner)
    if os.path.exists(jar): os.remove(jar)
    src = os.path.join(d, 'paper_martins2025', 'src.tgz'); curl(src, 'https://arxiv.org/e-print/2505.02993')   # tables A.1-A.6
    os.makedirs(os.path.join(d, 'paper_martins2025', 'x'), exist_ok=True)
    subprocess.call(['tar', '-xzf', src, '-C', os.path.join(d, 'paper_martins2025', 'x')])


def url_of(rel):
    """where a raw file came from (POST requests: the form data in brackets)"""
    p = rel.split('/')
    if p[0] == 'svo':
        if len(p) == 2: return SVO % ('tlusty_%sbin' % p[1].split('_')[0])
        return (SVO % p[1]) + 'fid=' + p[2][3:-4]
    if p[0] == 'powr':
        g = p[1]
        if p[2] == 'modellist.html': return POWR + 'powrgrid2.php [POST type=wr&grid=%s&Submit=Next]' % g
        if p[2] == 'models.txt': return 'parsed from powr/%s/modellist.html' % g
        if p[2] == 'sed' and p[3].endswith('.modinfo'): return POWR + 'download.php?content=modinfo&grid=%s&model=%s' % (g, p[3][:-8])
        if p[2] == 'sed': return POWR + 'download.php [POST content=specdl&grid=%s&model=%s&datatype=sed&dlformat=txt]' % (g, p[3][:-4])
        return POWR + 'isoplots/%s/%s' % (g, p[2].replace('.notfound.html', ''))
    if p[0] == 'gotberg18':
        if p[1] == 'paper': return 'https://arxiv.org/e-print/1802.03018' if p[2] == 'src.tgz' else 'unpacked from gotberg18/paper/src.tgz'
        if len(p) == 3 and p[2] == 'index.html': return CDS + '/%s/' % p[1]
        return CDS + '/' + '/'.join(p[1:])
    if p[0] == 'pollux':
        if p[1] == 'paper_martins2025': return 'https://arxiv.org/e-print/2505.02993' if p[2] == 'src.tgz' else 'unpacked from pollux/paper_martins2025/src.tgz'
        if len(p) == 2: return POLLUX_URL + '/collections/%s/download/ [POST csrfmiddlewaretoken]' % p[1][:-4]
        return 'unpacked from pollux/%s.zip' % p[1]
    return 'unknown'


def manifest(root):
    lines = ['# hot-star grids: path  sha256  bytes  mtime (download date; archive date of unpacked files)  url  source', '# written %s' % datetime.date.today().isoformat()]
    for f in sorted(glob.glob(os.path.join(root, '**', '*'), recursive=True)):
        if not os.path.isfile(f) or f.endswith('MANIFEST.txt'): continue
        rel = os.path.relpath(f, root); h = hashlib.sha256(open(f, 'rb').read()).hexdigest()
        cite = next((v for k, v in sorted(CITE.items(), key=lambda kv: -len(kv[0])) if rel.startswith(k)), '')
        lines.append('%s  %s  %d  %s  %s  %s' % (rel, h, os.path.getsize(f),
                                                 datetime.datetime.fromtimestamp(os.path.getmtime(f)).strftime('%Y-%m-%d'), url_of(rel), cite))
    open(os.path.join(root, 'MANIFEST.txt'), 'w').write('\n'.join(lines) + '\n')
    print('MANIFEST: %d files' % (len(lines) - 2))


if __name__ == '__main__':
    root = os.path.join(sys.argv[1], 'hotstar')
    if '--manifest' not in sys.argv: fetch(root)
    manifest(root)
