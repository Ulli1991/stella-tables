#!/usr/bin/env python3
"""Hot-star band LUT for the hotstar_stack radiation source, from the raw grids (fetch/fetch_hotstar.sh).

  build_hotstar_lut.py RAW_ROOT [OUT.hdf5]      (default OUT: data/lut/hotstar_lut.hdf5)

RAW_ROOT/svo/tlusty_{ostar,bstar}bin/fid*.xml   TLUSTY OSTAR2002 (Lanz & Hubeny 2003) and BSTAR2006 (2007) from SVO
RAW_ROOT/powr/<grid>/                            PoWR OB Vd3 log Q tables, WR SEDs (sed/*.txt, *.modinfo), models.txt
RAW_ROOT/gotberg18/grid_*/<model>/SED.txt        Gotberg+2018 stripped stars (CDS J/A+A/615/A78), paper tables in paper/x
RAW_ROOT/pollux/<collection>/*.sed(.txt)        CMFGEN O stars and very massive stars (Martins & Palacios 2017-2022, Martins+2025)
                                                 from POLLUX; Martins+2025 tables in pollux/paper_martins2025/x

Every spectrum goes through build_spectral_lut.integrate_bands: band energy fractions L_band/L_bol and photon rates
per erg/s of L_bol, the same convention as the ATLAS9 LUT. Cell class per grid cell: 0 model, 1 interpolated in Z
towards the grid's Z = 0 set (BSTAR below 0.1 Zsun), 2 filled from the nearest model (Eddington-limit holes and
missing models: nearest log g at the same Teff; PoWR: nearest valid node).
"""
import os, re, sys, glob
import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from build_stellar_tables import ENERGY_BANDS, PHOTON_RATE_BANDS, L_sun
from build_spectral_lut import integrate_bands

EB = [n for n, _, _ in ENERGY_BANDS]
QB = [n for n, _, _ in PHOTON_RATE_BANDS]


def bands(w, F):
    o = np.argsort(w); ef, qf = integrate_bands(w[o], F[o])
    return np.array([ef[n] for n in EB]), np.array([qf[n] for n in QB])


def read_svo(path):
    s = open(path).read()
    p = {k: float(v) for k, v in re.findall(r'<PARAM name="(teff|logg|z)"[^>]*value="([^"]*)"', s)}
    d = np.array(re.findall(r'<TR>\s*<TD>([^<]+)</TD>\s*<TD>([^<]+)</TD>\s*</TR>', s), float)
    return p, d[:, 0], d[:, 1]


def fill_logg(e, q, cls, have):
    """holes at fixed (Teff, Z): copy the nearest log g that has a model (class 2)"""
    nT, ng, nZ = have.shape
    for i in range(nT):
        for k in range(nZ):
            ok = np.where(have[i, :, k])[0]
            for j in range(ng):
                if not have[i, j, k] and len(ok):
                    s = ok[np.argmin(abs(ok - j))]
                    e[i, j, k], q[i, j, k], cls[i, j, k] = e[i, s, k], q[i, s, k], 2


BLEND_BO = (27000.0, 30000.0)   # BSTAR -> OSTAR blend window [K]


# ── TLUSTY: BSTAR2006 up to 30 kK, OSTAR2002 above ─────────────────────────────
def tlusty(root):
    mods = {}
    for grid in ('bstar', 'ostar'):
        for f in sorted(glob.glob(os.path.join(root, 'svo', 'tlusty_%sbin' % grid, 'fid*.xml'))):
            p, w, F = read_svo(f)
            mods[(grid, p['teff'], p['logg'], p['z'])] = bands(w, F)
    TB = sorted({t for g, t, _, _ in mods if g == 'bstar'}); TO = sorted({t for g, t, _, _ in mods if g == 'ostar' and t > TB[-1]})
    T = np.array(TB + TO); G = np.array(sorted({lg for _, _, lg, _ in mods}))
    ZO = sorted({z for g, _, _, z in mods if g == 'ostar' and z > 0}); ZB = sorted({z for g, _, _, z in mods if g == 'bstar'})
    Z = np.array(ZO)                                       # Z/Zsun axis (nonzero OSTAR compositions)
    nT, ng, nZ = len(T), len(G), len(Z)
    e = np.zeros((nT, ng, nZ, len(EB))); q = np.zeros((nT, ng, nZ, len(QB)))
    cls = np.zeros((nT, ng, nZ), np.int8); have = np.zeros((nT, ng, nZ), bool)
    for i, t in enumerate(T):
        grid = 'bstar' if t <= TB[-1] else 'ostar'
        for j, lg in enumerate(G):
            for k, z in enumerate(Z):
                m = mods.get((grid, t, lg, z))
                if m is None and grid == 'bstar' and z < min(x for x in ZB if x > 0):
                    lo, hi = mods.get((grid, t, lg, 0.0)), mods.get((grid, t, lg, min(x for x in ZB if x > 0)))
                    if lo is not None and hi is not None:
                        a = z / min(x for x in ZB if x > 0)
                        m = ((1 - a) * lo[0] + a * hi[0], (1 - a) * lo[1] + a * hi[1]); cls[i, j, k] = 1
                if m is not None:
                    e[i, j, k], q[i, j, k], have[i, j, k] = m[0], m[1], True
    # BSTAR -> OSTAR: blend linearly in Teff over BLEND_BO (OSTAR linear between its 27.5 and 30 kK nodes)
    for i, t in enumerate(T):
        if not (BLEND_BO[0] < t <= BLEND_BO[1]): continue
        wO = (t - BLEND_BO[0]) / (BLEND_BO[1] - BLEND_BO[0])
        for j, lg in enumerate(G):
            for k, z in enumerate(Z):
                o1, o2 = mods.get(('ostar', 27500.0, lg, z)), mods.get(('ostar', 30000.0, lg, z))
                if not have[i, j, k] or o1 is None or o2 is None: continue
                a = np.clip((t - 27500.0) / 2500.0, 0, 1)
                e[i, j, k] = (1 - wO) * e[i, j, k] + wO * ((1 - a) * o1[0] + a * o2[0])
                q[i, j, k] = (1 - wO) * q[i, j, k] + wO * ((1 - a) * o1[1] + a * o2[1]); cls[i, j, k] = 1
    fill_logg(e, q, cls, have)
    return dict(log_Teff=np.log10(T), logg=G, Z_rel=Z, energy_fracs=e, photon_rates_per_L=q, cell_class=cls,
                n_models=len(mods))


# ── PoWR ───────────────────────────────────────────────────────────────────────
def powr_params(path):
    out = {}
    for blk in open(path).read().split('MODEL:')[1:]:
        lines = blk.strip().splitlines(); mid = lines[0].strip(); d = {}
        for ln in lines[1:]:
            m = re.match(r'\s*(.+?)\s*=\s*([-0-9.Ee+]+)', ln)
            if m: d[m.group(1).strip()] = float(m.group(2))
        out[mid] = d
    return out


def powr_q(path):
    if not os.path.exists(path): return {}
    return {a: float(b) for a, b in (ln.split() for ln in open(path) if ln.strip())}


OB_GRIDS = {'gal-ob-vd3': 1.0, 'lmc-ob-vd3': 0.5, 'smc-ob-vd3': 1 / 7., 'subsmc-ob-vd3': 1 / 31.}  # Z/Zsun from the grid pages


def powr_ob(root):
    out = {}
    for g, zrel in OB_GRIDS.items():
        P = powr_params(os.path.join(root, 'powr', g, 'parameters.txt'))
        Q = {k: powr_q(os.path.join(root, 'powr', g, 'logQ_%s.txt' % k)) for k in ('HII', 'HeI', 'HeII')}
        T = np.array(sorted({p['Teff  [kK]'] for p in P.values()})); G = np.array(sorted({round(p['log g [cgs]'], 2) for p in P.values()}))
        arr = np.full((3, len(T), len(G)), np.nan); cls = np.zeros((len(T), len(G)), np.int8)
        for mid, p in P.items():
            i, j = int(np.searchsorted(T, p['Teff  [kK]'])), int(np.searchsorted(G, round(p['log g [cgs]'], 2)))
            lL = np.log10(L_sun) + p['log L [L_sun]']
            for a, k in enumerate(('HII', 'HeI', 'HeII')):
                if mid in Q[k]: arr[a, i, j] = (Q[k][mid] if Q[k][mid] > 0 else -60.0) - lL   # log photons/s per erg/s; 0.00 = none tabulated
        have = np.isfinite(arr[2])
        for i in range(len(T)):
            ok = np.where(have[i])[0]
            for j in range(len(G)):
                if not have[i, j]:
                    if len(ok): s = (i, ok[np.argmin(abs(ok - j))])
                    else:
                        ii, jj = np.where(have); s = (ii[np.argmin(abs(ii - i) * 100 + abs(jj - j))], jj[np.argmin(abs(ii - i) * 100 + abs(jj - j))])
                    arr[:, i, j] = arr[:, s[0], s[1]]; cls[i, j] = 2
        out[g] = dict(Teff_kK=T, logg=G, log_q_per_L=arr, cell_class=cls, Z_rel=zrel, n_models=len(P))
    return out


WR_GRIDS = ['wne', 'wnl', 'wnl-h50', 'wc', 'lmc-wne', 'lmc-wnl-h20', 'lmc-wnl-h40', 'lmc-wc', 'smc-wne', 'smc-wnl-h20',
            'smc-wnl-h40', 'smc-wnl-h60', 'smc-wc-2021', '007-wne-2015', '007-wnl-h20-2015', '007-wnl-h40-2015',
            '007-wnl-h60-2015', '007-wc-2021']


# iron-group mass fraction in the model headers (*.modinfo): MW 1.4E-3, LMC 7E-4, SMC 3E-4, 0.07 Zsun 9.2E-5
WR_FE = {'': 1.4e-3, 'lmc': 7e-4, 'smc': 3e-4, '007': 9.2e-5}


NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+)(?:[EeDd][-+]?\d+)?')


def read_powr_sed(path):
    """log lambda [A], log F_lambda [erg/cm2/s/A at 10 pc]; fixed-width columns can touch (0.26E-01-100.0)"""
    return np.array([[float(x.replace('D', 'E')) for x in NUM.findall(ln)] for ln in open(path) if len(NUM.findall(ln)) == 2])


def powr_wr(root):
    out = {}
    for g in WR_GRIDS:
        d = os.path.join(root, 'powr', g)
        mlist = [ln.split() for ln in open(os.path.join(d, 'models.txt'))]
        P = powr_params(os.path.join(d, 'parameters.txt'))
        T = np.array(sorted({float(t) for _, t, _ in mlist})); R = np.array(sorted({float(r) for _, _, r in mlist}))
        e = np.zeros((len(T), len(R), len(EB))); q = np.zeros((len(T), len(R), len(QB)))
        have = np.zeros((len(T), len(R)), bool); qcheck = []
        Qt = {k: powr_q(os.path.join(d, 'logQ_%s.txt' % k)) for k in ('HII', 'HeI')}
        for mid, t, r in mlist:
            f = os.path.join(d, 'sed', mid + '.txt')
            if not os.path.exists(f) or os.path.getsize(f) == 0: continue
            a = read_powr_sed(f); w = 10 ** a[:, 0]; F = np.where(a[:, 1] > -99, 10 ** a[:, 1], 0.0)
            i, j = int(np.searchsorted(T, float(t))), int(np.searchsorted(R, float(r)))
            e[i, j], q[i, j] = bands(w, F); have[i, j] = True
            if mid in P and mid in Qt['HII'] and Qt['HII'][mid] > 0:
                lL = np.log10(L_sun) + P[mid]['log L [L_sun]']
                qcheck.append(np.log10(max(q[i, j, QB.index('Q_ion')], 1e-300)) + lL - Qt['HII'][mid])
        cls = np.where(have, 0, 2).astype(np.int8)
        for i in range(len(T)):          # nearest model in (T, Rt) index space
            for j in range(len(R)):
                if not have[i, j]:
                    ii, jj = np.where(have); k = np.argmin((ii - i) ** 2 + (jj - j) ** 2)
                    e[i, j], q[i, j] = e[ii[k], jj[k]], q[ii[k], jj[k]]
        m = re.search(r'-h(\d\d)', g)
        xh = 0.0 if ('wne' in g or 'wc' in g) else (float(m.group(1)) / 100 if m else 0.2)   # MW 'wnl' is WNL-H20
        fe = WR_FE[g.split('-')[0] if g.split('-')[0] in WR_FE else '']
        out[g] = dict(T_kK=T, log_Rt=R, energy_fracs=e, photon_rates_per_L=q, cell_class=cls, X_H=xh,
                      family='WC' if 'wc' in g else ('WNE' if 'wne' in g else 'WNL'), Fe=fe, Z_rel=fe / WR_FE[''],
                      vinf=float(np.median([p['v_inf [km/s]'] for p in P.values()])), D=float(np.median([p['D_inf'] for p in P.values()])),
                      n_models=int(have.sum()), qcheck_max_dex=float(np.max(np.abs(qcheck))) if qcheck else np.nan)
    return out


# ── Gotberg+2018 ───────────────────────────────────────────────────────────────
GB = {'0.014': ('grid_014', 'table_param.tex'), '0.006': ('grid_006', 'table_param_006.tex'),
      '0.002': ('grid_002', 'table_param_002.tex'), '0.0002': ('grid_0002', 'table_param_0002.tex')}


def gb_table(fn):
    rows = []
    for line in open(fn):
        line = line.replace('\\midrule', '').replace('\\bottomrule', '')
        if '&' not in line or 'Group' in line or '[$' in line: continue
        c = [x.strip() for x in line.replace('\\\\', '').split('&')]
        m = re.match(r'([0-9.]+)', c[1])
        if m: rows.append([float(m.group(1))] + [float(re.match(r'[-0-9.]+', x).group(0)) for x in c[2:]])
    return np.array(rows)   # M_init P M_strip M_H logL logL_H T* Teff logg R X_s Y_s logMdot Gamma vinf logQ0 logQ1 logQ2


def gotberg(root):
    out = {}
    for z, (gdir, tab) in GB.items():
        t = gb_table(os.path.join(root, 'gotberg18', 'paper', 'x', tab)); rec = []
        for f in glob.glob(os.path.join(root, 'gotberg18', gdir, '*', 'SED.txt')):
            Mi = float(re.search(r'M1_([0-9.]+)q', f).group(1)); row = t[np.argmin(abs(t[:, 0] - Mi))]
            a = np.loadtxt(f, comments='#'); ef, qf = bands(a[:, 0], a[:, 1])
            L = 4 * np.pi * (3.0857e21) ** 2 * np.trapz(a[:, 1], a[:, 0])
            rec.append((row[7], row[8], row[10], Mi, ef, qf, np.log10(L / L_sun) - row[4], np.log10(qf[QB.index('Q_ion')] * L) - row[15]))
            if abs(rec[-1][6]) > 0.1: print('  Gotberg Z=%s M_init=%g: SED log L %.2f vs table %.2f (Teff from the table, fractions per L)' % (z, Mi, np.log10(L / L_sun), row[4]))
        rec.sort(key=lambda r: r[0])
        out[z] = dict(Teff_kK=np.array([r[0] for r in rec]), logg=np.array([r[1] for r in rec]), X_s=np.array([r[2] for r in rec]),
                      M_init=np.array([r[3] for r in rec]), energy_fracs=np.array([r[4] for r in rec]),
                      photon_rates_per_L=np.array([r[5] for r in rec]), Z=float(z),
                      logL_check_max_dex=float(np.max(np.abs([r[6] for r in rec]))), logQ0_check_max_dex=float(np.max(np.abs([r[7] for r in rec]))))
    return out


# ── CMFGEN (Martins & Palacios) via POLLUX ────────────────────────────────────
# Z/Zsun of each set (relative to the grids' own solar), POLLUX collection and file prefix
CMF_SETS = [(1.0, 'CMFGEN-SED', 'MW_'), (0.4, 'VMS-ZLMC-SED', 'Zlmc_'), (0.2, 'CMFGEN-SED', 'Zsmc_'), (0.2, 'VMS-ZSMC-SED', 'Zsmc_'),
            (0.1, 'VMS-Z0p1Zsun-SED', 'Z0p1_'), (1 / 30., 'CMFGEN-SED', 'Z1on30_'), (0.01, 'VMS-Z0p01Zsun-SED', 'Z0p01_')]
CMF_TABLE = {'VMS-ZSMC-SED': 'zsmc', 'VMS-Z0p1Zsun-SED': 'z0p1', 'VMS-Z0p01Zsun-SED': 'z0p01'}
CMF_X_MIN = 0.4            # surface X below this is a WR phase in hotstar_radiation (PoWR WR / Gotberg)


def pollux_meta(path):
    return {k: v.strip() for k, v in re.findall(r"^\s*(\w+)\s*=\s*'([^']*)'", open(path, errors='ignore').read(), re.M)}


def martins25_rows(tex):
    """Martins+2025 Tables A.1-A.6 per model: M, Teff, log g, log L, log Mdot, X; set x0 (VMS Mdot not scaled with Z) or x1 (scaled with Z/Z_LMC)"""
    rows = []
    f = lambda x: float(x.rstrip('.')) if x not in ('-', '') else np.nan
    for m in re.finditer(r'\\label\{tab_(grz?)_(zsmc|z0p1|z0p01)\}(.*?)\\end\{tabular\}', open(tex).read(), re.S):
        s, z, body = m.groups()
        for ln in body.splitlines():
            c = [x.strip() for x in ln.replace('\\\\', '').split('&')]
            if len(c) < 10 or not re.match(r'^[0-9.]+$', c[0]): continue
            rows.append(dict(set='x1' if s == 'grz' else 'x0', z=z, M=f(c[0]), T=f(c[1]), g=f(c[2]), logL=f(c[3]), log_Mdot=f(c[5]),
                             X=f(c[8]) if c[8] != '-' else 1 - f(c[9])))
    return rows


def _closest(rows, T, lg, key, fn):
    """parameter row of the model named fn (Teff T, log g lg from the name): Teff within 300 K, log g within 0.25 (typos in names)"""
    c = [r for r in rows if abs(key(r)[0] - T) < 300 and abs(key(r)[1] - lg) < 0.25]
    if not c: raise ValueError('no parameter row for %s' % fn)
    return min(c, key=lambda r: (abs(key(r)[0] - T) / 100) ** 2 + ((key(r)[1] - lg) / 0.01) ** 2)


def cmfgen(root):
    """per Z set: scattered (log Teff, log g) CMFGEN models with band fractions. Coordinates are the evolutionary-model Teff and
    log g each spectrum was computed for: the Martins+2025 tables for the low-Z VMS, else the file name (several POLLUX
    headers carry a neighbouring model's parameters; X and Mdot come from the header row whose Teff and log g match).
    Low-Z VMS: the Z-scaled mass-loss sequences (x1) only. Models with X < CMF_X_MIN are dropped."""
    d = os.path.join(root, 'pollux'); rows = martins25_rows(os.path.join(d, 'paper_martins2025', 'x', 'vms_lowZ.tex'))
    out = {}
    for zrel, coll, pre in CMF_SETS:
        files = sorted(glob.glob(os.path.join(d, coll, '*.sed')))
        metas = [pollux_meta(f + '.txt') for f in files]
        files, metas = zip(*[(f, m) for f, m in zip(files, metas) if m['filename_atm'].startswith(pre)])
        for f, me in zip(files, metas):
            fn = me['filename_atm']
            t = re.search(r'_T(\d+)p(\d+)', fn); g = re.search(r'_hy?d?ro(\d)p(\d+)', fn)
            T, lg = float(t.group(1) + '.' + t.group(2)) * 1e3, float(g.group(1) + '.' + g.group(2))
            if coll in CMF_TABLE:
                r = _closest([r for r in rows if r['z'] == CMF_TABLE[coll]], T, lg, lambda r: (r['T'], r['g']), fn)
                if r['set'] != 'x1': continue
                T, lg, X, lMd, src = r['T'], r['g'], r['X'], r['log_Mdot'], 'Martins+2025'
            else:
                hm = _closest(metas, T, lg, lambda m: (float(m['Teff']), float(m['logg'])), fn)
                T, lg = float(hm['Teff']), float(hm['logg'])
                gchk = np.log10(4 * np.pi * 5.670374e-5 * 6.674e-8 * float(hm['mass']) * 1.989e33 * T ** 4 / (10 ** float(hm['lum']) * L_sun))
                if abs(gchk - lg) > 0.03: raise ValueError('%s: header log g %.2f vs GM/R^2 %.2f' % (fn, lg, gchk))
                X = 1 / (1 + 4 * 10 ** (float(hm['Abu_02']) - 12)); lMd = float(hm['Mdot'])
                src = 'Martins & Palacios 2022' if coll == 'VMS-ZLMC-SED' else ('Martins & Palacios 2017' if pre == 'MW_' else 'Martins & Palacios 2021')
            if X < CMF_X_MIN: continue
            a = np.loadtxt(f); ef, qf = bands(a[:, 0], a[:, 1])
            tsed = (np.trapz(a[:, 1], a[:, 0]) / 5.670374e-5) ** 0.25      # flux is given below the photosphere: T_sed >= Teff
            if not 0.97 < tsed / T < 1.6: raise ValueError('%s: SED temperature %.0f K vs Teff %.0f K' % (fn, tsed, T))
            o = out.setdefault('%.4f' % zrel, dict(Z_rel=zrel, rec=[]))
            o['rec'].append((np.log10(T), lg, X, lMd, ef, qf, os.path.basename(f), src))
    res = {}
    for k, o in out.items():
        rec = sorted(o['rec'], key=lambda r: (r[0], r[1])); keep = []
        for r in rec:                  # one model per (Teff, log g) position
            if not any(abs(r[0] - q[0]) < 1e-4 and abs(r[1] - q[1]) < 0.006 for q in keep): keep.append(r)
        res[k] = dict(log_Teff=np.array([r[0] for r in keep]), logg=np.array([r[1] for r in keep]), X_surf=np.array([r[2] for r in keep]),
                      log_Mdot=np.array([r[3] for r in keep]), energy_fracs=np.array([r[4] for r in keep]),
                      photon_rates_per_L=np.array([r[5] for r in keep]), files=np.array([r[6] for r in keep], dtype='S'),
                      Z_rel=o['Z_rel'], sources=', '.join(sorted({r[7] for r in keep})), n_models=len(keep),
                      n_duplicates_dropped=len(rec) - len(keep))
    return res


def write(path, tl, ob, wr, gb, cm):
    def put(grp, d):
        for k, v in d.items():
            if isinstance(v, np.ndarray): grp.create_dataset(k, data=v, compression='gzip', compression_opts=4)
            else: grp.attrs[k] = v
    with h5py.File(path, 'w') as f:
        f.attrs['energy_bands'] = EB; f.attrs['photon_rate_bands'] = QB
        f.attrs['convention'] = 'energy_fracs = L_band/L_bol; photon_rates_per_L = photons/s per erg/s of L_bol; cell_class 0 model, 1 Z-interpolated towards Z=0 or BSTAR/OSTAR blend (27-30 kK), 2 filled from nearest model'
        put(f.create_group('tlusty'), tl)
        f['tlusty'].attrs['source'] = 'TLUSTY BSTAR2006 (<= 30 kK; Lanz & Hubeny 2007) + OSTAR2002 (> 30 kK; Lanz & Hubeny 2003), SVO theory server tlusty_bstarbin / tlusty_ostarbin'
        for g, d in ob.items(): put(f.create_group('powr_ob/' + g), d)
        f['powr_ob'].attrs['source'] = 'PoWR OB Vd3 grids (Hainich+2019, Pauli+2026), log Q tables per model divided by the model L; grids recomputed 2025-12-15'
        for g, d in wr.items(): put(f.create_group('powr_wr/' + g), d)
        f['powr_wr'].attrs['source'] = 'PoWR WR grids (Todt+2015, Hamann+2006, Sander+2012, Hainich+2014/2015), SEDs per model'
        for z, d in gb.items(): put(f.create_group('gotberg/' + z), d)
        f['gotberg'].attrs['source'] = 'Gotberg+2018 A&A 615 A78 stripped-star CMFGEN SEDs; Teff, log g, X_s from the paper tables'
        for z, d in cm.items(): put(f.create_group('cmfgen/' + z), d)
        f['cmfgen'].attrs['source'] = ('CMFGEN O-star and very-massive-star SEDs along evolutionary tracks (Martins & Palacios 2017, 2021, 2022; '
                                       'Martins+2025, Z-scaled VMS mass loss), POLLUX (Palacios+2010); scattered in (log Teff, log g) per Z')


if __name__ == '__main__':
    root = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'data', 'lut', 'hotstar_lut.hdf5')
    tl = tlusty(root); print('TLUSTY: %d models, Teff %d-%d K, cells by class %s' % (tl['n_models'], 10 ** tl['log_Teff'][0], 10 ** tl['log_Teff'][-1], np.bincount(tl['cell_class'].ravel())), flush=True)
    ob = powr_ob(root); print('PoWR OB:', {g: (d['n_models'], int((d['cell_class'] == 2).sum())) for g, d in ob.items()}, flush=True)
    wr = powr_wr(root); print('PoWR WR:', {g: (d['family'], round(d['Z_rel'], 3), d['X_H'], d['n_models'], round(d['qcheck_max_dex'], 3)) for g, d in wr.items()}, flush=True)
    gb = gotberg(root); print('Gotberg:', {z: (len(d['Teff_kK']), round(d['logL_check_max_dex'], 3), round(d['logQ0_check_max_dex'], 3)) for z, d in gb.items()}, flush=True)
    cm = cmfgen(root); print('CMFGEN:', {z: (d['n_models'], d['n_duplicates_dropped'], d['sources']) for z, d in cm.items()}, flush=True)
    write(out, tl, ob, wr, gb, cm); print('wrote', out)
