# Changelog

## Reproducible build (after tables-aa87301)

- `build.sh` checks the Python environment against `requirements-build.txt` (python-waterboa 2026.07: Python 3.14.6,
  numpy 2.4.6, scipy 1.18.0, h5py 3.16.0), runs single-threaded and on numpy's baseline CPU code paths: the same
  tables on every machine with that environment (login node and compute node orion002 built identical files).
- Python float sums are correctly rounded (`math.fsum`); the energy table is installed from its canonical file, so it
  is the same in every environment (sha256 e1b03445..., as pinned in GIZMO).
- Content hashes in the pinned environment: `stellar_tables_unified.hdf5` b94a4edc596ba04af2143a5982b0777b9297adc43e12cd4c39d2709aea2a5c62,
  `stellar_wind_cumulative.hdf5` 89f40e36aa5abe3ed3846973204dba5fec0b2598f5af2d9b0a9fbd9495d49285 (unchanged).
  The stellar table differs from tables-aa87301 (built with Python 3.6, numpy 1.19) only in the last bits.
- Physics unchanged against tables-aa87301; the stellar table differs only in the last bits (per dataset, max
  |a-b|/max(|a|,|b|)): 29 of 52 datasets identical; `elem_ej_AGB_mass` 9.1e-14, `elem_ej_SN_mass` 6.2e-14, band
  luminosities and photon rates (`log_L_*`, `log_Q_*`) <= 1.9e-14, `X_init` 4.4e-16; `surface_abundances` differs
  in 235 of 3.4e7 entries, all O and F at ~1e-18 (round-off residue of fully depleted elements, absolute <= 4.3e-19,
  relative 3.7e-4). `stellar_wind_cumulative.hdf5` and all chemistry files are identical.
- `yield_sources` attribute: Limongi+2025 8-14 Msun (was "8-13"), as the builder uses it.
- Tests run with numpy 2 (`np.trapz` removed there).

## fix/wind-composition (not yet pinned in GIZMO)

Compared with the tables GIZMO pins (stella-tables a4830ab; unified 14a23087…, wind 5ad36096…). Built with
`./build.sh` (radiation source `hotstar_stack`, MIST v2.5 [a/Fe] = +0.4). Content hashes (GIZMO's
`scripts/stella_tables_hash.py`):

| file | hash |
|---|---|
| `stellar_tables_unified.hdf5` | 7593d9a85ff134441d0e417302b5fd1d40da7aba04f0d8727e2dca9ec3213613 |
| `stellar_wind_cumulative.hdf5` | 89f40e36aa5abe3ed3846973204dba5fec0b2598f5af2d9b0a9fbd9495d49285 |
| chemistry tables | unchanged |

### Yields and winds
- Winds (M >= 8): birth composition X_init plus the wind enrichment of the PARSEC v2 wind ejecta (all 27 elements,
  M >= 14) or BoOST (8-12 Msun, blended over 12-14); before: track surface, Na, Al, Si, Fe, ... at a solar-scaled
  initial value. Wind Na/Mg at Z = 0.01: 0.050 -> 0.133 (40 Msun), 0.052 -> 0.228 (100 Msun); birth 0.041.
- SN ejecta = X_init Mej + (lifetime net - wind net): no wind production counted twice. The old closure step
  multiplied net yields by up to 3.5e4 (H = 0 in 61 SN cells, e.g. Z = 1e-4, 11 Msun: H 0 -> 5.31, O 1.45 -> 0.34).
- PARSEC v2 ejecta: stable Ni no longer moved into Fe (Ni 3e-5 -> 5.3e-3 at Z = 0.014, 14 Msun); unfilled columns
  (exactly 0: P everywhere; Si29, Cl, K, Sc, V, Mn, Co, Cu in winds and failed SNe) at initial composition;
  26Al -> Mg in every isotope source; explosion type and yields from one model where the bracketing masses differ.
- Late wind: mass lost after the last alive age-grid point and before ZAMS is in the wind table; M_preSN = the raw
  track's final mass (Z = 0.014, 100 Msun: 37.64 -> 28.43 Msun; 9.2 Msun no longer locked in the black hole).
- No winds below 8 Msun (GIZMO injects none): M_current = M_init, M_preSN = M_init, the AGB event carries the whole
  Karakas lifetime budget.
- AGB: K2010 Z = 1e-4, 6 Msun VW93 model only (net N 0.071 -> 0.037); KL16 from Z = 0.007, linear in log Z between
  the Karakas metallicities; KL16 net yields against its own initial composition.
- 8-14 Msun: Z dependence as NuGrid net-yield difference; elements Limongi lacks from NuGrid.
- `surface_abundances`: X_init plus the track's change since ZAMS (ZAMS = birth composition everywhere; was off by
  1.4-2x for metals at Z = 1e-4-5e-4 and 0.84x for Mg).

### Evolution
(Lifetime and t_PMS ratios below 2.5 Msun quoted for [a/Fe] = 0; +0.4 changes them by <= 1% for Z <= 1e-3 and
2-8% for 1-2 Msun at Z = 0.01.)
- MIST v2.5 replaces v1.2 below 2.5 Msun (grid Z from the track headers, tracks to the hottest post-AGB point),
  [a/Fe] = +0.4 by default (`STELLA_MIST_ALPHA_FE=p0` for solar-scaled).
- Same ZAMS for MIST and PARSEC v2 (central H down by 1e-4): t_PMS 0.5-0.8x below 2.5 Msun, lifetimes 1-5% shorter
  there, no t_PMS kink at the 2.0-2.5 Msun blend.
- Mass scale of interpolated tracks from the birth mass (was the trimmed-ZAMS mass: +0.7% at 500 Msun).

### Radiation
- Source: `hotstar_stack` (TLUSTY, CMFGEN, PoWR OB/WR, Gotberg+18, ATLAS9/PHOENIX at each point's log g and Z); the
  pinned table was ATLAS9 at log g = 4, Z_sun, labelled "Planck".
- Band integrals include the band-edge intervals; the ATLAS9/PHOENIX LUT now really contains PHOENIX below 3500 K
  (the old one was ATLAS9 only); PoWR WR He+-ionizing rates no longer lost between sparse samples.
- Radiation stored as age-cell averages; the first point after death repeats the last value. Post-MS integral of
  L_bol / raw track: 0.16-0.92 -> 0.94-0.99.
- Lifetime-integrated Q(H) vs the pinned table: x1.0-1.4 at 20-40 Msun, x1.7-2.3 at 100-300 Msun; L_FUV x0.5-0.9
  above 40 Msun. Against the PARSEC v2 photon tables: Q(H) and Q(HeI) within 0.1 dex except 20 Msun at Z = 0.014
  (+0.24 / +0.77 dex, TLUSTY vs wind atmospheres at late O); Q(HeII) differs by up to 3 dex either way.
- New datasets `log_Q_He0`, `log_Q_He1`, `log_Q_ion_He0`; composite bands are exact sums of the atomic bands.

### Datasets unchanged
`Z`, `log_Z`, `M_init`, `log_M_init`, `log_age_yr`, `elements`, `X_init`, `remnant_type`, `type_ia_yields`,
`phase_transitions/t_MS_start`, `t_WR_start`.
