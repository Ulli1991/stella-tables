# Changelog

## fix/yields-8-14-and-tests (not pinned in GIZMO, not released)

Compared with tables-e3a707a. Content hashes (GIZMO's `scripts/stella_tables_hash.py`): `stellar_tables_unified.hdf5`
9e09a203c4d23235a5b80c3a10719160c1c021e853575a99956ace59459d9dce, `stellar_wind_cumulative.hdf5`
1fff3f620561bc2d5c7be9b08119aa0fcb95319304d65f89dd7a467eda1253eb. Chemistry files unchanged. Needs the GIZMO reader of
stella_dev fix/late-wind-release (elem_ej_wind_total, type_ia_yields per Z); older readers refuse these tables.

- 8-14 Msun yields: Limongi, Chieffi & Roberti 2025 (arXiv:2505.22030) Table 2, every isotope, decayed to the stable
  isobar (stable set from the NuGrid list; Cl36 and Cu64 branchings); source file in `stellar/data/yields/limongi2025/`.
  The hand-typed main-isotope values (Mg, Si, S, Ne, Fe 7-16 % short) and all NuGrid terms are gone: Ni, Zn, K, Mn, Co,
  Cr came whole from NuGrid 12-15 Msun (net Ni 0.021-0.036 per star, ~10x the paper), and NuGrid's mass-cut change with Z
  was added to explosive Fe (low-Z Fe 4-8x too high, N negative in 21 cells). Below Z_sun each metal isotope is scaled by
  its Chieffi & Limongi 2018 (set R, v = 0) 13 Msun ratio to [Fe/H] = 0 (same code and network); above, the solar net
  yields. Remnant: the paper's mass cut M_init - sum(Table 2), 1.32-1.56 Msun (was the Fe core, and 1.25 for M <= 10).
- Type Ia: Gronow+2021a,b model M10_03 (1.0 Msun core, 0.03 Msun He shell, double detonation), Z dependent (0.01, 0.1,
  1, 3 Z_sun, linear in log Z), from the arXiv appendix tables (`stellar/data/yields/gronow2021/`, choice in its README);
  `type_ia_yields` is now (N_Z, N_elements), 1.04 Msun per event. Replaces the hand-typed W7 vector whose 58Ni had been
  set by hand ([Ni/Fe] +0.53 per event).
- PARSEC v2 SN (CCSN, PPISN, PISN): the explosive part comes from the raw total and wind files of the same model
  (total - wind) instead of total net minus a wind net built with the track's wind mass (49 cells went negative and were
  clipped); explosive part of the model's own mass M_fin - M_bar, the extra envelope of a heavier track at birth
  composition. Elements a source destroys (AGB, PARSEC v2 wind and explosion) leave as the source's destroyed fraction of
  the star's birth content, produced mass is added; H closes the mass (the transfer is recorded per source:
  AGB KL16 1.1 %, K2010 0.17 %, PARSEC v2 SN 2.5 %, wind 0.9 % of the ejecta at most; the PARSEC v2 models are born with
  Y = 0.249-0.302, the table with solar He at every Z). Raw elements with exactly 0 in a PARSEC v2 explosion stay 0.
  The builder no longer clips: a negative element stops the build.
- Closure: KL16 ejecta fractions in the files sum to 1.005-1.006, renormalised to the expelled mass (this was the 4.5e-3
  booked to H); elements beyond Zn (KL16, K2010 g sink, Li/Be/B) go to H as X_init counts them. Mass booked to H:
  KL16 9e-16, Limongi 2e-15 Msun; Karakas 2010 1.6e-5 Msun = its own models' imbalance (sum of all 77 species' net
  yields), recorded as attributes `closure_residual_*`, `depletion_transfer_to_H_*`. K2010 species decayed to their
  stable ends (13N -> C, 55Fe -> Mn, ...).
- Late wind: the wind lost after the last alive grid age leaves between that age and death. The table holds physical
  cumulative values only: up to the last grid age before death the wind lost so far, after death all wind; new dataset
  `elem_ej_wind_total` (N_Z, N_M, N_elements). The reader interpolates linearly in log age from the last grid age to the
  total at lifetime_yr and holds v_wind there (GIZMO fix/late-wind-release). M_current at the last alive age is the
  track mass (up to 15 Msun more than before); M_preSN = M_init - all wind.
- Phase transitions: blended only where both tracks reach the phase, else the nearer track; a track stopping at the
  AGB tip starts its post-AGB at its AGB end. 0 transitions out of order (471 before).
- IMF ratios (Kroupa, M >= 8): [Ni/Fe] -0.03/-0.03/+0.04/+0.11 at Z = 1e-4/1e-3/6e-3/0.014 (was +0.56/+0.82/+0.56/+0.80),
  [Mg/Fe] +0.37/+0.36/+0.36/+0.34 (was +0.11/+0.27/+0.16/+0.34). Whole population after 13.7 Gyr with Ia: [Ni/Fe]
  -0.13..+0.02 (was +0.42..+0.44 with W7), [Mn/Fe] -0.06/+0.01, [N/O] +0.11/+0.24 at Z = 6e-3/0.014.
- Changes against e3a707a: elem_ej_SN_mass (8-14 Msun block; PARSEC v2 SN: up to 5.5 Msun H in the 200 Msun PISN at Z = 2e-4, a Z between PARSEC v2 nodes), elem_ej_AGB_mass
  (up to 0.05 Msun), remnant_mass (8-14 Msun), M_current (last alive age), M_preSN (last bits), type_ia_yields (shape),
  phase transitions (142 entries), surface_abundances and elem_ej_wind_cumulative (late wind; destroyed-fraction rule:
  at Z = 0.04, 0.06, which take the Z = 0.03 PARSEC v2 model, the wind of a 350 Msun star keeps 2.7 instead of 5.7 Msun O
  and 3 Msun more H; at PARSEC v2 metallicities the wind totals are unchanged); v_wind identical.
- Tests: `test_yield_ratios.py` (IMF ratios against halo and disc bands, Limongi+2025 against the paper at every Z,
  8 and 14 Msun seams against the raw sources, clipped ejecta), `phase_transitions_ordered`, `late_wind_follows_raw_track`,
  `ia_yields_are_gronow21_m10_03`, `parsec2_sn_is_raw_explosion`, `closure_residuals_at_roundoff`, wind total at death;
  the wind-composition tests follow the destroyed-fraction rule; `interp.py` follows GIZMO's wind lookup.

## fix/late-o-ionizing (tables unchanged)

- Late-O ionizing rates checked against Martins, Schaerer & Hillier 2005 (CMFGEN with winds, Tables 1-3, 36 Galactic
  O stars): the stack is +0.02 dex in Q(H)/L on average (max +0.14 at O9.5 V) and +0.04 dex in Q(HeI)/L (max +0.47 at
  O9.5 V). The PoWR OB grid agrees with Martins+05 to +0.03 dex in Q(H) at 28-34 kK. PARSEC v2's photon tables
  (WM-basic) are 0.25 dex (Q(H)) and 0.6 dex (Q(HeI)) below Martins+05 at the same Teff and log g along the 20 Msun
  track; Martins+05 Sect. 5.3.2 reports the same WM-basic offset (up to 0.3 dex for dwarfs below 40 kK). The +0.24 / +0.77
  dex offset at 20 Msun, Z = 0.014 is therefore PARSEC's atmosphere, not ours; no table change.
- New test `ostar_rates_vs_martins05` (`tests/stellar/test_radiation.py`, data `tests/stellar/data/`): Q(H)/L within
  0.15 dex at every Martins+05 point and 0.05 dex on average, Q(HeI)/L within 0.5 dex.

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
  (+0.24 / +0.77 dex: PARSEC's WM-basic rates are low at late O, see fix/late-o-ionizing); Q(HeII) differs by up to 3 dex either way.
- New datasets `log_Q_He0`, `log_Q_He1`, `log_Q_ion_He0`; composite bands are exact sums of the atomic bands.

### Datasets unchanged
`Z`, `log_Z`, `M_init`, `log_M_init`, `log_age_yr`, `elements`, `X_init`, `remnant_type`, `type_ia_yields`,
`phase_transitions/t_MS_start`, `t_WR_start`.
