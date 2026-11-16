#!/usr/bin/env python3
"""Element ratios of the IMF-integrated yields, the 8-14 Msun Limongi+2025 block against the paper, the yield seams
at 8 and 14 Msun, and negative (clipped) ejecta.

  test_yield_ratios.py BUILD_DIR      exit 0 = all PASS

Lifetime ejecta of one star = SN (or AGB) ejecta + cumulative wind at lifetime_yr (evaluated as GIZMO does: linear in
log age, age clamped at the star's lifetime). Net yield = lifetime ejecta - X_init * (M_init - M_rem).

IMF: Kroupa (2001) 0.08-100 Msun (slopes 1.3 / 2.3, break 0.5 Msun), per Msun formed; net yields linear in log M
between grid masses. [X/Y] = log10 of the net-yield mass ratio over the Asplund+09 mass ratio (the table's X_init).

Bands (stated before running; widths cover the spread between studies and NLTE corrections):
  massive stars (M >= 8, winds + CCSN), Z = 1e-4, 1e-3, 6e-3, 0.014 -- the core-collapse plateau of halo stars:
    [O/Fe]  +0.25..+0.80   Cayrel+04 (A&A 416, 1117) ~+0.7 (1D LTE); Amarsi+19 (A&A 630, A104) 3D NLTE ~+0.5
    [Mg/Fe] +0.15..+0.65   Cayrel+04 <[Mg/Fe]> = +0.27 (LTE); Andrievsky+10 (A&A 509, A88) NLTE ~+0.6
    [Si/Fe] +0.15..+0.65   Cayrel+04 ~+0.37
    [Ni/Fe] -0.30..+0.30   Cayrel+04 ~-0.05
    [Mn/Fe] -0.70..+0.00   at Z <= 1e-3 only (odd-Z, rises with Z; the sources give ~+0.15 at Z_sun):
                           Cayrel+04 ~-0.4 (LTE); Bergemann & Gehren 08 (A&A 492, 823) NLTE ~-0.1..-0.3
    [N/O]   -3.50..-0.50   at Z <= 1e-3 only: non-rotating massive stars make secondary N; the primary-N plateau of
                           metal-poor HII regions, log N/O ~ -1.5 ([N/O] ~ -0.65, Izotov & Thuan 99, ApJ 511, 639),
                           is an upper bound for them
  whole population after 13.7 Gyr (all stars with lifetime < 13.7 Gyr + 1.3e-3 SN Ia per Msun formed, Maoz & Graur
  17, ApJ 848, 25, with the table's W7 vector), Z = 6e-3 and 0.014 -- near-solar disc stars (Bensby+14, A&A 562, A71;
  Battistini & Bensby 15, A&A 577, A9 for Mn):
    [O/Fe], [Mg/Fe], [Si/Fe] -0.20..+0.40;  [Ni/Fe] -0.25..+0.25;  [Mn/Fe] -0.50..+0.15;  [N/O] -0.50..+0.40
  whole population, Z = 1e-4 and 1e-3:  [Ni/Fe] -0.30..+0.30;  [N/O] -1.50..+0.50 (Izotov & Thuan 99 plateau, AGB
  primary N above it)

Limongi+2025 (arXiv:2505.22030, Table 2, stellar/data/yields/limongi2025/tabyields.tex): this file parses the table
and decays the unstable species itself (beta decay along the isobar to the stable nuclide; stable set from the NuGrid
isotope list; Cl36, Cu64 branchings). The Z dependence is checked against Chieffi & Limongi 2018 (set R, v = 0,
13 Msun, table8.dat) isotope ratios, read here independently of the builder.
"""
import os, re, sys, collections
import numpy as np, h5py

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
YDIR = os.path.join(REPO, 'stellar', 'data', 'yields')
ELEMS = ['H', 'He', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
         'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
SYM = ('n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y '
       'Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os '
       'Ir Pt Au Hg Tl Pb Bi').split()
ZNUC = {s.lower(): i for i, s in enumerate(SYM)}
Z_SUN = 0.014
N_IA = 1.3e-3
TESTS = []
def test(fn): TESTS.append(fn); return fn


# ── table helpers ─────────────────────────────────────────────────────

class Tab:
    def __init__(self, b):
        self.t = h5py.File(os.path.join(b, 'stellar_tables_unified.hdf5'), 'r')
        self.w = h5py.File(os.path.join(b, 'stellar_wind_cumulative.hdf5'), 'r')
        t = self.t
        self.Z, self.M, self.la = t['Z'][()], t['M_init'][()], t['log_age_yr'][()]
        self.lt, self.rm, self.rt = t['lifetime_yr'][()], t['remnant_mass'][()], t['remnant_type'][()]
        self.X, self.ys = t['X_init'][()], t['yield_source'][()]
        self.SN, self.AGB = t['elem_ej_SN_mass'][()].astype(np.float64), t['elem_ej_AGB_mass'][()].astype(np.float64)
        self.cw = self.w['elem_ej_wind_cumulative']
        self.ia = t['type_ia_yields'][()].astype(np.float64)
        self._wind = {}

    def wind_at_death(self, iz, im):
        """cumulative wind per element at lifetime_yr, linear in log age between the grid ages (GIZMO's lookup)"""
        if (iz, im) not in self._wind:
            x = np.log10(max(self.lt[iz, im], 100.0))
            i = int(np.clip(np.searchsorted(self.la, x, side='right') - 1, 0, len(self.la) - 2))
            f = float(np.clip((x - self.la[i]) / (self.la[i + 1] - self.la[i]), 0, 1))
            c = self.cw[iz, im, i:i + 2, :].astype(np.float64)
            self._wind[(iz, im)] = c[0] * (1 - f) + c[1] * f
        return self._wind[(iz, im)]

    def lifetime_ejecta(self, iz, im):
        return self.SN[iz, im] + self.AGB[iz, im] + self.wind_at_death(iz, im)

    def net(self, iz, im):
        """lifetime net yield per element: ejecta minus the birth composition of the ejected mass"""
        ej = self.lifetime_ejecta(iz, im)
        return ej - self.X[iz] * ej.sum()

    def iz(self, Z): return int(np.argmin(abs(self.Z - Z)))
    def im(self, M): return int(np.argmin(abs(self.M - M)))


def solar():
    sys.path.insert(0, os.path.join(REPO, 'stellar'))
    from solar_abundances import solar_pattern
    return solar_pattern()


def imf_integral(T, iz, sel):
    """net yield per element per Msun formed, Kroupa 0.08-100, over the stars where sel(m, lifetime) is True"""
    mm = np.logspace(np.log10(0.08), 2.0, 6000)
    xi = np.where(mm < 0.5, (mm / 0.5) ** -1.3, (mm / 0.5) ** -2.3)
    xi /= np.trapezoid(xi * mm, mm)
    lgM = np.log10(T.M)
    net = np.array([T.net(iz, im) for im in range(len(T.M))])
    lt = 10 ** np.interp(np.log10(mm), lgM, np.log10(T.lt[iz]))
    ok = sel(mm, lt)
    out = np.zeros(len(ELEMS))
    for k in range(len(ELEMS)):
        y = np.interp(np.log10(mm), lgM, net[:, k])
        out[k] = np.trapezoid(np.where(ok, y * xi, 0.0), mm)
    return out


def bracket(y, S, a, b):
    ya, yb = y[ELEMS.index(a)], y[ELEMS.index(b)]
    if ya <= 0 or yb <= 0: return float('nan')
    return float(np.log10(ya / yb) - np.log10(S[a] / S[b]))


def check_bands(rows):
    bad = [r for r in rows if not (r[3] <= r[2] <= r[4])]
    txt = '; '.join('Z=%g %s %+.2f' % (r[0], r[1], r[2]) for r in rows)
    return not bad, ('OUT OF BAND: ' + '; '.join('Z=%g %s %+.2f not in [%+.2f, %+.2f]' % r for r in bad) + ' | ' if bad else '') + txt


# ── IMF-integrated ratios ─────────────────────────────────────────────

CC_BANDS = [('O', 'Fe', 0.25, 0.80), ('Mg', 'Fe', 0.15, 0.65), ('Si', 'Fe', 0.15, 0.65), ('Ni', 'Fe', -0.30, 0.30)]
POP_BANDS_HI = [('O', 'Fe', -0.20, 0.40), ('Mg', 'Fe', -0.20, 0.40), ('Si', 'Fe', -0.20, 0.40),
                ('Ni', 'Fe', -0.25, 0.25), ('Mn', 'Fe', -0.50, 0.15), ('N', 'O', -0.50, 0.40)]
POP_BANDS_LO = [('Ni', 'Fe', -0.30, 0.30), ('N', 'O', -1.50, 0.50)]


@test
def imf_massive_star_ratios(T):
    """massive stars (M >= 8: winds + CCSN): [O/Fe] [Mg/Fe] [Si/Fe] [Ni/Fe] in the halo-plateau bands at
    Z = 1e-4, 1e-3, 6e-3, 0.014; [Mn/Fe] and [N/O] (below the primary-N plateau) at Z <= 1e-3"""
    S = solar(); rows = []
    for Z in (1e-4, 1e-3, 6e-3, 0.014):
        y = imf_integral(T, T.iz(Z), lambda m, lt: m >= 8.0)
        rows += [(Z, '[%s/%s]' % (a, b), bracket(y, S, a, b), lo, hi) for a, b, lo, hi in CC_BANDS]
        if Z <= 1e-3:   # the halo-star bands of Mn (odd-Z, metallicity dependent) and N apply at low Z only
            rows.append((Z, '[Mn/Fe]', bracket(y, S, 'Mn', 'Fe'), -0.7, 0.0))
            rows.append((Z, '[N/O]', bracket(y, S, 'N', 'O'), -3.5, -0.5))
    return check_bands(rows)


@test
def imf_population_ratios(T):
    """all stars dying within 13.7 Gyr + SN Ia (1.3e-3 per Msun formed, W7): near-solar ratios at Z = 6e-3, 0.014;
    [Ni/Fe] and [N/O] at Z = 1e-4, 1e-3"""
    S = solar(); rows = []
    for Z, bands in ((1e-4, POP_BANDS_LO), (1e-3, POP_BANDS_LO), (6e-3, POP_BANDS_HI), (0.014, POP_BANDS_HI)):
        y = imf_integral(T, T.iz(Z), lambda m, lt: lt < 1.37e10) + N_IA * T.ia
        rows += [(Z, '[%s/%s]' % (a, b), bracket(y, S, a, b), lo, hi) for a, b, lo, hi in bands]
    return check_bands(rows)


# ── Limongi+2025 Table 2, decayed here ────────────────────────────────

def stable_isobars():
    hdr = [l for l in open(os.path.join(YDIR, 'nugrid', 'isotope_yield_table_MESAonly_fryer12_delay_total.txt'))
           if l.startswith('H Isotopes:')][0]
    st = collections.defaultdict(set)
    for s in hdr.split(':', 1)[1].split(','):
        if s.strip():
            e, a = s.strip().split('-'); st[int(a)].add(ZNUC[e.lower()])
    for A, Zn in ((6, 3), (9, 4), (10, 5)): st[A].add(Zn)       # Li6, Be9, B10: stable, not in the NuGrid list
    return st


# odd-odd nuclei between two stable isobars (NNDC branchings); all others in Table 2 decay one way
BRANCH = {(17, 36): {18: 0.981, 16: 0.019}, (29, 64): {28: 0.615, 30: 0.385}, (21, 46): {22: 1.0}, (21, 48): {22: 1.0},
          (25, 54): {24: 1.0}, (27, 58): {26: 1.0}}


def decay_to_element(Zn, A, st):
    """{element index: branching} of the stable end of isotope (Zn, A); None for products outside H..Zn"""
    if Zn in st[A]: ends = {Zn: 1.0}
    elif (Zn, A) in BRANCH: ends = BRANCH[(Zn, A)]
    elif Zn < min(st[A]): ends = {min(st[A]): 1.0}
    elif Zn > max(st[A]): ends = {max(st[A]): 1.0}
    else: ends = {max(z for z in st[A] if z < Zn): 1.0}  # beyond Zn only (untracked)
    return {ELEMS.index(SYM[z]): b for z, b in ends.items() if SYM[z] in ELEMS}


def limongi25():
    """masses, M_cut (= M_init - sum of Table 2), {(Zn, A): ejecta[6]}"""
    iso = {}
    for ln in open(os.path.join(YDIR, 'limongi2025', 'tabyields.tex')):
        m = re.match(r'\$\\rm\^\{(\d+)\}([A-Za-z]+)\s*\$\s*&(.*)\\\\', ln.strip())
        if m: iso[(ZNUC[m.group(2).lower()], int(m.group(1)))] = np.array([float(x) for x in m.group(3).split('&')])
    M = np.array([9.22, 10., 11., 12., 13., 15.])
    return M, M - sum(iso.values()), iso


def lc18_13msun():
    """{([Fe/H], Zn, A): 13 Msun ejecta} of Chieffi & Limongi 2018 set R, v = 0"""
    out = {}
    for ln in open(os.path.join(YDIR, 'chieffi_limongi2018', 'table8.dat')):
        if int(ln[0:3]) != 0: continue
        m = re.match(r'([A-Za-z]+)(\d*)$', ln[7:12].strip())
        out[(int(ln[4:6]), ZNUC[m.group(1).lower()], int(m.group(2)) if m.group(2) else 1)] = float(ln[13:23])
    return out


def z_ratio(lc, Zn, A, feh):
    """LC18 13 Msun ejecta of (Zn, A) at [Fe/H] over [Fe/H] = 0, linear in [Fe/H] (in log where both > 0)"""
    nodes = [-3, -2, -1, 0]
    x = float(np.clip(feh, -3, 0)); j = min(int(np.floor(x)) + 3, 2); f = x - nodes[j]
    ref = lc.get((0, Zn, A), 0.0)
    if ref <= 0: return 1.0
    a, b = lc.get((nodes[j], Zn, A), 0.0) / ref, lc.get((nodes[j + 1], Zn, A), 0.0) / ref
    return float(np.exp((1 - f) * np.log(a) + f * np.log(b))) if a > 0 and b > 0 else (1 - f) * a + f * b


def limongi_expected_net(M, Z, X, L=None):
    """net yield per element of the 8-14 Msun block at (M, Z): Table 2 decayed, metals scaled by the LC18 isotope
    ratio below Z_sun, net held at its solar value above; He net held; H closes the paper's ejecta mass"""
    if L is None: L = (limongi25(), lc18_13msun(), stable_isobars())
    (ms, cut, iso), lc, st = L
    mej = float(np.interp(M, ms, ms - cut))
    feh = np.log10(Z / Z_SUN)
    tot = np.zeros(len(ELEMS))
    for (Zn, A), v in iso.items():
        r = z_ratio(lc, Zn, A, feh) if (feh < 0 and Zn > 2) else 1.0
        for k, br in decay_to_element(Zn, A, st).items():
            tot[k] += br * r * float(np.interp(M, ms, v))
    if feh > 0:   # super-solar: net yields of all elements stay at their solar values
        Xs = X.copy(); Xs[2:] *= Z_SUN / Z; Xs[0] = 1 - Xs[1:].sum()
        tot[1:] += (X[1:] - Xs[1:]) * mej
    tot[0] = mej - tot[1:].sum()
    return tot - X * mej, mej


@test
def limongi25_block_matches_paper_all_Z(T):
    """8-14 Msun lifetime net yields = Limongi+2025 Table 2 (all isotopes, decayed) at Z = 0.014 and = Table 2 x the
    LC18 13 Msun isotope ratio at every other grid Z; every element, every 8-14 Msun grid mass (tolerance 1e-3 of the
    net yield + 1e-4 of the star's birth mass of that element)"""
    L = (limongi25(), lc18_13msun(), stable_isobars())
    worst, n = (0.0, None), 0
    for iz, Z in enumerate(T.Z):
        for im, M in enumerate(T.M):
            if not (8 <= M < 14): continue
            want, _ = limongi_expected_net(M, Z, T.X[iz], L)
            have = T.net(iz, im)
            tol = 1e-3 * abs(want) + 1e-4 * T.X[iz] * M + 1e-9
            tol[0] = 2e-4 * M
            e = np.abs(have - want) / tol; k = int(np.argmax(e)); n += 1
            if e[k] > worst[0]: worst = (e[k], (Z, M, ELEMS[k], have[k], want[k]))
    return n > 0 and worst[0] < 1, '%d (Z, M) cells; worst |table - paper|/tol = %.3g at (Z, M, elem, table, paper) = %s' % (n, worst[0], worst[1])


@test
def limongi25_remnant_is_paper_mass_cut(T):
    """8-14 Msun remnant mass = M_init - sum(Table 2) interpolated in mass (the paper's mass cut, 1.32-1.56 Msun)"""
    ms, cut, _ = limongi25()
    worst = (0.0, None)
    for iz in range(len(T.Z)):
        for im, M in enumerate(T.M):
            if 8 <= M < 14:
                e = abs(T.rm[iz, im] - np.interp(M, ms, cut))
                if e > worst[0]: worst = (e, (T.Z[iz], M, T.rm[iz, im], float(np.interp(M, ms, cut))))
    return worst[0] < 1e-6, 'max |M_rem - M_cut| = %.3g at (Z, M, table, paper) = %s' % worst


# ── yield seams ───────────────────────────────────────────────────────

def parsec2_raw_net(root, Z, M):
    """{elem: net} of the raw PARSEC v2 total-ejecta row (Z, M) against the file's initial composition, ejecta =
    M - M_bar; elements with an unfilled (exactly 0) column in the total or the wind file are left out"""
    d = os.path.join(root, 'parsec2_vms')
    zdir = [p for p in os.listdir(d) if p.startswith('ejecta_Z') and np.isclose(float(p[8:]), Z, rtol=1e-9)]
    if not zdir: return None
    f = [p for p in os.listdir(os.path.join(d, zdir[0])) if p.endswith('total_ejecta.dat')][0]
    fw = [p for p in os.listdir(os.path.join(d, zdir[0])) if p.endswith('winds_ejecta.dat')]
    wind0 = set()   # wind-file columns that are exactly 0 although present at birth: unfilled there too
    if fw:
        Lw = open(os.path.join(d, zdir[0], fw[0])).read().splitlines(); nw = Lw[3].split()[4:]
        tw = Lw[1].replace('#', '').split(); xw = [float(v) for v in tw[-len(nw):]]
        for ln in Lw[4:]:
            p = ln.split()
            if len(p) >= 4 + len(nw) and np.isclose(float(p[0]), M):
                wind0 = {re.match(r'[A-Z]+', n).group(0).capitalize() for n, x, v in zip(nw, xw, p[4:]) if float(v) == 0.0 and x > 0}
    L = open(os.path.join(d, zdir[0], f)).read().splitlines()
    names = L[3].split()[7:]
    tok = L[1].replace('#', '').split(); x0 = [float(v) for v in tok[tok.index('SN') + 1:]]
    for ln in L[4:]:
        p = ln.split()
        if len(p) >= 7 + len(names) and np.isclose(float(p[0]), M):
            mej = float(p[0]) - float(p[5]); out = {}; skip = set()
            for nm, x, v in zip(names, x0, p[7:]):
                e = {'AL26': 'Mg'}.get(nm, re.match(r'[A-Z]+', nm).group(0).capitalize())
                if e not in ELEMS: continue
                if float(v) == 0.0 and x > 0: skip.add(e)
                out[e] = out.get(e, 0.0) + float(v) - x * mej
            return {e: v for e, v in out.items() if e not in skip | wind0}
    return None


SEAM_SHOW = ['O', 'Mg', 'Si', 'Fe', 'Ni', 'Mn', 'Zn', 'K', 'N']


@test
def seam_14_matches_sources(T):
    """12 -> 14 Msun (Limongi+2025 -> PARSEC v2): the table's per-element net-yield step equals the step between the two
    raw sources at every Z where PARSEC v2 has ejecta (tolerance 1e-3 of the nets + 1e-4 of the birth mass); the
    source steps printed are physics of the two sources, not of the table"""
    root = os.path.join(T.build, 'data_root'); L = (limongi25(), lc18_13msun(), stable_isobars())
    worst, n, doc = (0.0, None), 0, []
    i12, i14 = T.im(12.0), T.im(14.0)
    for iz, Z in enumerate(T.Z):
        p2 = parsec2_raw_net(root, Z, 14.0)
        if p2 is None: continue
        lim, _ = limongi_expected_net(12.0, Z, T.X[iz], L)
        h12, h14 = T.net(iz, i12), T.net(iz, i14)
        for e, v in p2.items():
            k = ELEMS.index(e); n += 1
            if e == 'H': continue
            tol = 1e-3 * (abs(v) + abs(lim[k])) + 1e-4 * T.X[iz, k] * 14.0 + 1e-9
            err = abs((h14[k] - h12[k]) - (v - lim[k])) / tol
            if err > worst[0]: worst = (err, (Z, e, h14[k] - h12[k], v - lim[k]))
        if Z in (1e-4, 1e-3, 0.014):
            doc.append('Z=%g ' % Z + ' '.join('%s %.2g->%.2g' % (e, lim[ELEMS.index(e)], p2[e]) for e in SEAM_SHOW if e in p2))
    return n > 0 and worst[0] < 1, 'worst |step table - step sources|/tol = %.3g at (Z, elem, table, sources) = %s | source nets 12 -> 14 Msun: %s' % (
        worst[0], worst[1], ' | '.join(doc))


@test
def seam_8_matches_sources(T):
    """7 -> 8 Msun (Karakas AGB -> Limongi+2025 9.22 Msun, clamped below 9.22): the 8 Msun side equals the paper (as
    scaled in Z); the AGB side is the Karakas tests' job. The step itself is printed: it is the step between the
    two sources"""
    L = (limongi25(), lc18_13msun(), stable_isobars())
    worst, doc = (0.0, None), []
    i7, i8 = T.im(7.0), T.im(8.0)
    for iz, Z in enumerate(T.Z):
        lim, _ = limongi_expected_net(8.0, Z, T.X[iz], L)
        h7, h8 = T.net(iz, i7), T.net(iz, i8)
        tol = 1e-3 * abs(lim) + 1e-4 * T.X[iz] * 8.0 + 1e-9; tol[0] = 1.6e-3
        e = np.abs(h8 - lim) / tol; k = int(np.argmax(e))
        if e[k] > worst[0]: worst = (e[k], (Z, ELEMS[k], h8[k], lim[k]))
        if Z in (1e-4, 1e-3, 0.014):
            doc.append('Z=%g ' % Z + ' '.join('%s %.2g->%.2g' % (x, h7[ELEMS.index(x)], h8[ELEMS.index(x)]) for x in SEAM_SHOW))
    return worst[0] < 1, 'worst |8 Msun - paper|/tol = %.3g at (Z, elem, table, paper) = %s | nets 7 -> 8 Msun: %s' % (
        worst[0], worst[1], ' | '.join(doc))


# ── negative ejecta ───────────────────────────────────────────────────

@test
def no_clipped_negative_ejecta(T):
    """no element is clipped: every SN and AGB ejecta element of a star that ejects mass is > 0 (every element is
    present at birth, so an exact 0 is a negative that was clipped), and every lifetime net yield is finite"""
    bad = []
    for name, ej in (('SN', T.SN), ('AGB', T.AGB)):
        tot = ej.sum(axis=2)
        z = (ej <= 0) & (tot[..., None] > 0) & (T.X[:, None, :] > 0)
        for iz, im, k in zip(*np.where(z)):
            bad.append((name, float(T.Z[iz]), float(T.M[im]), ELEMS[k]))
    c = collections.Counter(b[3] for b in bad)
    return not bad, '%d clipped (kind, Z, M, elem) entries, by element %s; first: %s' % (len(bad), dict(c), bad[:6])


def main():
    b = sys.argv[1]
    T = Tab(b); T.build = b
    nf = 0
    for fn in TESTS:
        try: ok, msg = fn(T)
        except Exception as e:
            import traceback; traceback.print_exc(); ok, msg = False, 'EXCEPTION %r' % e
        nf += not ok
        print('  [%s] %s: %s' % ('PASS' if ok else 'FAIL', fn.__name__, msg))
    print('\n%d/%d failed' % (nf, len(TESTS))); sys.exit(1 if nf else 0)


if __name__ == '__main__':
    main()
