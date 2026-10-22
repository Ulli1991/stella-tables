"""Band fractions for the hotstar_stack radiation source (LUT from build_hotstar_lut.py).

Per track point, in this order:
  WR phase (surface X < WR_X_MAX and log Teff > WR_LOGTEFF_MIN):
      Z >= GOTBERG_BELOW_ZREL Zsun: PoWR WR grid (WNL if X > WNL_X_MIN, else WC if surface C > N, else WNE), entered at
          T* = Teff and the transformed radius Rt = R [(v_inf / 2500 km/s) / (Mdot sqrt(D) / 1e-4 Msun/yr)]^(2/3) with the
          grid's own v_inf and D; linear in surface X between WNL hydrogen grids, linear in log Z between metallicities
          (lowest 0.07 Zsun); every band from the SED.
      Z <  GOTBERG_BELOW_ZREL Zsun: Gotberg+2018 stripped-star SEDs, linear in log Teff along each Z sequence,
          linear in log Z between sequences (lowest Z = 0.0002).
  other stars, Teff >= TLUSTY_SEAM[0]: TLUSTY BSTAR2006/OSTAR2002, trilinear in (log Teff, log g, log Z/Zsun),
      blended with ATLAS9 across TLUSTY_SEAM (weight linear in log Teff; energy fractions linearly, photon rates
      geometrically), the 55 kK edge models blended into a blackbody over
      55-65 kK, blackbody above;
      He+-ionizing photons (Q_He1) from the PoWR OB Vd3 log Q(HeII) tables (per erg/s, bilinear in Teff, log g, linear
      in log Z between 1/31 and 1 Zsun), blended in log q with the ATLAS9 side over TLUSTY_SEAM and into the blackbody
      over 62-73 kK; the He+ energy bands follow that rate.
  cooler stars: ATLAS9 where every interpolation corner is a model, blackbody elsewhere (as masked_planck).
Coverage class per point and band group: 0 model, 1 blended (ATLAS9/TLUSTY seam or a Z-interpolated BSTAR cell),
2 clamped (point outside a grid's range, or a filled grid cell, or a WR / O-star He+ metallicity outside the grids),
3 blackbody.  Q_He1 is order-of-magnitude: it is set by the wind, not the photosphere.
"""
import numpy as np
import h5py

WR_X_MAX = 0.4
WR_LOGTEFF_MIN = 4.45
WNL_X_MIN = 0.05
TLUSTY_SEAM = (15000.0, 20000.0)
GOTBERG_BELOW_ZREL = 0.1
POWR_OB_TMAX_KK = 62.0
EDGE_BLEND = 65000.0 / 55000.0       # TLUSTY hot edge -> blackbody over 55-65 kK, PoWR OB He+ edge over 62-73 kK (factor in T)
CLASSES = ('model', 'blended', 'clamped', 'blackbody')
GROUPS = ('Q_H0', 'Q_He0', 'Q_He1', 'L_FUV', 'L_LW')
E_GROUP = {'OPT_NIR': 'L_FUV', 'NUV_lo': 'L_FUV', 'NUV_hi': 'L_FUV', 'FUV_mid': 'L_FUV', 'LW': 'L_LW',
           'ion_H0': 'Q_H0', 'ion_He0': 'Q_He0', 'ion_He1': 'Q_He1', 'ion_He2': 'Q_He1'}
Q_GROUP = {'Q_ion': 'Q_H0', 'Q_ion_H0': 'Q_H0', 'Q_He0': 'Q_He0', 'Q_He1': 'Q_He1'}
E_HE1 = 54.4 * 1.602176634e-12      # erg, He+ edge


class HotLUT:
    def __init__(self, path):
        with h5py.File(path, 'r') as f:
            self.EB = [x.decode() if isinstance(x, bytes) else x for x in f.attrs['energy_bands']]
            self.QB = [x.decode() if isinstance(x, bytes) else x for x in f.attrs['photon_rate_bands']]
            t = f['tlusty']
            self.tl = {k: t[k][()] for k in ('log_Teff', 'logg', 'Z_rel', 'energy_fracs', 'photon_rates_per_L', 'cell_class')}
            self.ob = {}
            for g in f['powr_ob']:
                self.ob[g] = {k: f['powr_ob'][g][k][()] for k in ('Teff_kK', 'logg', 'log_q_per_L', 'cell_class')}
                self.ob[g]['Z_rel'] = float(f['powr_ob'][g].attrs['Z_rel'])
            self.wr = {}
            for g in f['powr_wr']:
                d = f['powr_wr'][g]
                self.wr[g] = {k: d[k][()] for k in ('T_kK', 'log_Rt', 'energy_fracs', 'photon_rates_per_L', 'cell_class')}
                self.wr[g].update({k: d.attrs[k] for k in ('X_H', 'family', 'Z_rel', 'vinf', 'D')})
                self.wr[g]['family'] = self.wr[g]['family'].decode() if isinstance(self.wr[g]['family'], bytes) else self.wr[g]['family']
            self.gb = {float(z): {k: f['gotberg'][z][k][()] for k in ('Teff_kK', 'energy_fracs', 'photon_rates_per_L')} for z in f['gotberg']}


def _corners(axes, pts):
    """per-axis lower index, weight and an out-of-range flag; points are clamped to the axes"""
    idx, wts, out = [], [], np.zeros(len(pts), bool)
    for d, g in enumerate(axes):
        x = pts[:, d]; out |= (x < g[0] - 1e-9) | (x > g[-1] + 1e-9); x = np.clip(x, g[0], g[-1])
        if len(g) == 1:
            idx.append(np.zeros(len(x), int)); wts.append(np.zeros(len(x))); continue
        i = np.clip(np.searchsorted(g, x) - 1, 0, len(g) - 2); idx.append(i); wts.append((x - g[i]) / (g[i + 1] - g[i]))
    return idx, wts, out


def _interp(axes, tables, cell_class, pts):
    """multilinear interpolation of each table (shape axes + (n,)); returns values, worst corner class, out flag"""
    idx, wts, out = _corners(axes, pts); nd = len(axes)
    vals = [np.zeros((len(pts), t.shape[-1])) for t in tables]; worst = np.zeros(len(pts), np.int8)
    for corner in range(2 ** nd):
        bits = [(corner >> d) & 1 for d in range(nd)]
        w = np.ones(len(pts)); ii = []
        for d in range(nd):
            w = w * (wts[d] if bits[d] else 1 - wts[d]); ii.append(np.minimum(idx[d] + bits[d], len(axes[d]) - 1))
        ii = tuple(ii)
        for v, t in zip(vals, tables): v += w[:, None] * t[ii]
        worst = np.maximum(worst, np.where(w > 0, cell_class[ii], 0))
    return vals, worst, out


def hotstar_band_fractions(H, lT, lg, lZ, z_sun, atlas9_fn, planck_fn, track=None):
    """energy fractions E (N, nE), photon rates per erg Q (N, nQ), coverage class per group C (N, 5)"""
    EB, QB = H.EB, H.QB; N = len(lT)
    pe, pq = planck_fn(lT)
    PE = np.column_stack([pe[n] for n in EB]); PQ = np.column_stack([pq[n] for n in QB])
    E, Q = PE.copy(), PQ.copy(); C = np.full((N, len(GROUPS)), 3, np.int8)
    zrel = 10 ** lZ / z_sun
    tr = track or {}
    X = np.asarray(tr.get('X', np.full(N, 0.7))); Xc = np.asarray(tr.get('C', np.zeros(N))); Xn = np.asarray(tr.get('N', np.zeros(N)))
    wr = (X < WR_X_MAX) & (lT > WR_LOGTEFF_MIN) & ('log_Mdot' in tr)
    T = 10 ** lT

    # ATLAS9 (masked): cool stars and the cool side of the seam
    lo = ~wr & (T < TLUSTY_SEAM[1])
    if lo.any():
        ef, qf, cls = atlas9_fn(lT[lo], lg[lo], lZ[lo])
        use = cls == 0
        E[lo] = np.where(use[:, None], np.column_stack([ef[n] for n in EB]), PE[lo])
        Q[lo] = np.where(use[:, None], np.column_stack([qf[n] for n in QB]), PQ[lo])
        C[lo] = np.where(use, 0, 3)[:, None]

    # TLUSTY
    hot = ~wr & (T >= TLUSTY_SEAM[0])
    tl = H.tl
    tmax = tl['log_Teff'][-1]; ledge = np.log10(EDGE_BLEND)
    inside = hot & (lT <= tmax + ledge)
    if inside.any():
        pts = np.column_stack([np.minimum(lT[inside], tmax), lg[inside], np.log10(np.maximum(zrel[inside], 1e-30))])
        (e, q), worst, out = _interp((tl['log_Teff'], tl['logg'], np.log10(tl['Z_rel'])), (tl['energy_fracs'], tl['photon_rates_per_L']), tl['cell_class'], pts)
        c = np.where(out | (worst == 2), 2, np.where(worst == 1, 1, 0))
        Tin = T[inside]; w = np.clip((np.log10(Tin) - np.log10(TLUSTY_SEAM[0])) / np.log10(TLUSTY_SEAM[1] / TLUSTY_SEAM[0]), 0, 1)
        Ein, Qin = E[inside], Q[inside]
        E[inside] = np.where((w < 1)[:, None], (1 - w)[:, None] * Ein + w[:, None] * e, e)
        lgb = (1 - w)[:, None] * np.log(np.maximum(Qin, 1e-300)) + w[:, None] * np.log(np.maximum(q, 1e-300))   # photon rates: geometric blend
        Q[inside] = np.where((w < 1)[:, None], np.exp(lgb), q)
        C[inside] = np.where(w < 1, np.maximum(c, 1), c)[:, None]
        b = np.clip((lT[inside] - tmax) / ledge, 0, 1)       # beyond the hot edge: edge model -> blackbody
        E[inside] = (1 - b)[:, None] * E[inside] + b[:, None] * PE[inside]; Q[inside] = (1 - b)[:, None] * Q[inside] + b[:, None] * PQ[inside]
        C[inside] = np.where((b > 0)[:, None], np.maximum(C[inside], 1), C[inside])
    above = hot & ~inside            # blackbody beyond TLUSTY and the blend: E, Q already Planck, class 3
    C[above] = 3

    # O/B-star He+ from PoWR OB Vd3
    heb = hot & (T <= POWR_OB_TMAX_KK * 1e3 * EDGE_BLEND)
    if heb.any():
        lzr = np.log10(np.maximum(zrel[heb], 1e-30)); grids = sorted(H.ob.values(), key=lambda g: g['Z_rel'])
        lzg = np.log10([g['Z_rel'] for g in grids]); vals, wst = [], []
        for g in grids:
            pts = np.column_stack([np.minimum(T[heb] / 1e3, POWR_OB_TMAX_KK), lg[heb]])
            (v,), worst, out = _interp((g['Teff_kK'], g['logg']), (g['log_q_per_L'][2][..., None],), g['cell_class'], pts)
            vals.append(v[:, 0]); wst.append(np.where(out | (worst == 2), 2, 0))
        vals, wst = np.array(vals), np.array(wst)
        x = np.clip(lzr, lzg[0], lzg[-1]); i = np.clip(np.searchsorted(lzg, x) - 1, 0, len(lzg) - 2); a = (x - lzg[i]) / (lzg[i + 1] - lzg[i])
        k = np.arange(heb.sum())
        lq = (1 - a) * vals[i, k] + a * vals[i + 1, k]
        ch = np.maximum(wst[i, k], wst[i + 1, k]); ch = np.where((lzr < lzg[0] - 1e-9) | (lzr > lzg[-1] + 1e-9), 2, ch)
        bb = np.clip(np.log10(T[heb] / (POWR_OB_TMAX_KK * 1e3)) / ledge, 0, 1)     # He+ edge -> blackbody, in log q
        lq = (1 - bb) * lq + bb * np.log10(np.maximum(PQ[heb, QB.index('Q_He1')], 1e-300)); ch = np.where(bb > 0, np.maximum(ch, 1), ch)
        ws = np.clip(np.log10(T[heb] / TLUSTY_SEAM[0]) / np.log10(TLUSTY_SEAM[1] / TLUSTY_SEAM[0]), 0, 1)   # ATLAS9 side of the seam
        lq = ws * lq + (1 - ws) * np.log10(np.maximum(Q[heb, QB.index('Q_He1')], 1e-300)); ch = np.where(ws < 1, np.maximum(ch, 1), ch)
        qn = 10 ** lq; _set_he1(E, Q, heb, qn, EB, QB); C[heb, GROUPS.index('Q_He1')] = ch

    # WR phases
    if wr.any():
        hi = wr & (zrel >= GOTBERG_BELOW_ZREL)
        if hi.any(): _wr_powr(H, E, Q, C, hi, T, X, Xc, Xn, zrel, tr)
        lo = wr & (zrel < GOTBERG_BELOW_ZREL)
        if lo.any(): _wr_gotberg(H, E, Q, C, lo, T, 10 ** lZ)
    return E, Q, C


def _set_he1(E, Q, m, qn, EB, QB):
    """replace the He+-ionizing photon rate by qn (per erg) and keep the band ledger consistent"""
    i1, i2, iq1, iq0, iqi = EB.index('ion_He1'), EB.index('ion_He2'), QB.index('Q_He1'), QB.index('Q_He0'), QB.index('Q_ion')
    qo = Q[m, iq1]; eo = E[m, i1] + E[m, i2]
    r = np.where(qo > 0, qn / np.where(qo > 0, qo, 1), 0)
    E[m, i1] = np.where(qo > 0, E[m, i1] * r, qn * 1.1 * E_HE1); E[m, i2] = np.where(qo > 0, E[m, i2] * r, 0)
    Q[m, iq0] += qn - qo; Q[m, iqi] += qn - qo; Q[m, iq1] = qn


def _wr_rt(R_cm, log_Mdot, vinf, D):
    """transformed radius [Rsun]"""
    Mdot = 10 ** log_Mdot
    return (R_cm / 6.957e10) * ((vinf / 2500.0) / (np.maximum(Mdot, 1e-30) * np.sqrt(D) / 1e-4)) ** (2.0 / 3.0)


def _wr_powr(H, E, Q, C, m, T, X, Xc, Xn, zrel, tr):
    idx = np.where(m)[0]
    fam = np.where(X[idx] > WNL_X_MIN, 'WNL', np.where(Xc[idx] > Xn[idx], 'WC', 'WNE'))
    for f in ('WNL', 'WNE', 'WC'):
        sel = idx[fam == f]
        if not len(sel): continue
        grids = [g for g in H.wr.values() if g['family'] == f]
        zs = np.array(sorted({round(float(g['Z_rel']), 3) for g in grids}))
        per_z = []
        for z in zs:
            gz = sorted([g for g in grids if round(float(g['Z_rel']), 3) == z], key=lambda g: g['X_H'])
            res = [_wr_grid(g, T[sel], tr, sel) for g in gz]
            if len(gz) == 1 or f != 'WNL': per_z.append(res[0]); continue
            xs = np.array([g['X_H'] for g in gz]); x = np.clip(X[sel], xs[0], xs[-1])
            i = np.clip(np.searchsorted(xs, x) - 1, 0, len(xs) - 2); a = (x - xs[i]) / (xs[i + 1] - xs[i]); k = np.arange(len(sel))
            Es = np.array([r[0] for r in res]); Qs = np.array([r[1] for r in res]); Cs = np.array([r[2] for r in res])
            cl = np.maximum(Cs[i, k], Cs[i + 1, k]); cl = np.where((X[sel] < xs[0]) | (X[sel] > xs[-1]), 2, cl)
            per_z.append(((1 - a)[:, None] * Es[i, k] + a[:, None] * Es[i + 1, k], (1 - a)[:, None] * Qs[i, k] + a[:, None] * Qs[i + 1, k], cl))
        lzg = np.log10(zs); lz = np.log10(np.maximum(zrel[sel], 1e-30)); x = np.clip(lz, lzg[0], lzg[-1])
        if len(zs) > 1:
            i = np.clip(np.searchsorted(lzg, x) - 1, 0, len(lzg) - 2); a = (x - lzg[i]) / (lzg[i + 1] - lzg[i])
        else:
            i = np.zeros(len(sel), int); a = np.zeros(len(sel))
        j = np.minimum(i + 1, len(zs) - 1); k = np.arange(len(sel))
        Es = np.array([p[0] for p in per_z]); Qs = np.array([p[1] for p in per_z]); Cs = np.array([p[2] for p in per_z])
        E[sel] = (1 - a)[:, None] * Es[i, k] + a[:, None] * Es[j, k]; Q[sel] = (1 - a)[:, None] * Qs[i, k] + a[:, None] * Qs[j, k]
        cl = np.maximum(Cs[i, k], Cs[j, k]); cl = np.where((lz < lzg[0] - 1e-9) | (lz > lzg[-1] + 1e-9), 2, cl)
        C[sel] = cl[:, None]


def _wr_grid(g, T, tr, sel):
    rt = _wr_rt(np.asarray(tr['R_cm'])[sel], np.asarray(tr['log_Mdot'])[sel], g['vinf'], g['D'])
    pts = np.column_stack([T / 1e3, np.log10(np.maximum(rt, 1e-30))])
    # PoWR T axis is logarithmic: interpolate in log T*
    (e, q), worst, out = _interp((np.log10(g['T_kK']), g['log_Rt']), (g['energy_fracs'], g['photon_rates_per_L']), g['cell_class'],
                                 np.column_stack([np.log10(pts[:, 0]), pts[:, 1]]))
    return e, q, np.where(out | (worst == 2), 2, 0).astype(np.int8)


def _wr_gotberg(H, E, Q, C, m, T, Z):
    zs = np.array(sorted(H.gb)); lzg = np.log10(zs); sel = np.where(m)[0]
    per_z = []
    for z in zs:
        g = H.gb[z]; o = np.argsort(g['Teff_kK']); lt = np.log10(g['Teff_kK'][o]); x = np.log10(T[sel] / 1e3)
        xc = np.clip(x, lt[0], lt[-1]); e = np.column_stack([np.interp(xc, lt, g['energy_fracs'][o][:, b]) for b in range(g['energy_fracs'].shape[1])])
        q = np.column_stack([np.interp(xc, lt, g['photon_rates_per_L'][o][:, b]) for b in range(g['photon_rates_per_L'].shape[1])])
        per_z.append((e, q, np.where((x < lt[0]) | (x > lt[-1]), 2, 0)))
    lz = np.log10(Z[sel]); x = np.clip(lz, lzg[0], lzg[-1])
    i = np.clip(np.searchsorted(lzg, x) - 1, 0, len(lzg) - 2); a = (x - lzg[i]) / (lzg[i + 1] - lzg[i]); k = np.arange(len(sel))
    Es = np.array([p[0] for p in per_z]); Qs = np.array([p[1] for p in per_z]); Cs = np.array([p[2] for p in per_z])
    E[sel] = (1 - a)[:, None] * Es[i, k] + a[:, None] * Es[i + 1, k]; Q[sel] = (1 - a)[:, None] * Qs[i, k] + a[:, None] * Qs[i + 1, k]
    cl = np.maximum(Cs[i, k], Cs[i + 1, k]); cl = np.where(lz < lzg[0] - 1e-9, 2, cl)
    C[sel] = cl[:, None]
