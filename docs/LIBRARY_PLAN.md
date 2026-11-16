# stella-tables library plan

Status: 2026-11-16, main = e3a707a. A plan, not a record of done work.

## Goal

stella-tables is a stellar evolution library for single stars. It is a core deliverable in its own right, not a
GIZMO helper. It should:

- combine the best public sources per quantity and regime, and say which one it used;
- be tested against the raw data it reads, not only against itself;
- be reproducible bit for bit in the pinned environment;
- be extensible: a new source enters as one module with its own tests;
- use own MESA runs where no public grid is good enough, and as an independent cross-check.

## How to read this plan

Each module has the same five parts: sources now, quality bar, tests that enforce it, known gaps, candidates to
fold in. A quality bar is a number a test checks. If a module has no such number yet, the plan says so.

References are arXiv IDs where checked; author and year otherwise (IDs to be added when checked). Full catalogue with URLs and licences:
`reports/stellar_evolution_grids_review.md` (§6) and `reports/Hot star atmosphere grids.md` in the GIZMO tree.

Table grid today: 15 metallicities Z = 1e-4 ... 0.06, 110 masses 0.09-500 Msun, 768 ages, 27 elements H-Zn.

---

## 1. Tracks (L, Teff, R, mass, lifetimes, phases, core masses)

**Sources now.**
- MIST v2.5 below 2 Msun, [a/Fe] = +0.4 by default (Dotter+2026, 2602.22012).
- PARSEC v2.0 VMS from 2 Msun (Costa+2025, 2501.12917), blended in log M over 2.0-2.5 Msun.
- Both interpolated in log Z and in fractional lifetime. Same ZAMS definition (central H down by 1e-4).
- PARSEC v1.2S and BoOST as fallbacks only.

**Quality bar.**
- Grid nodes reproduce the raw track (exactly at a node, to interpolation error between).
- Lifetime and L_bol continuous across the source change (except the MIST He-flash feature, XFAIL).
- Post-MS and lifetime integral of L_bol within 15% of the raw track.
- Pre-SN mass equals the raw track's final mass.

**Tests.** `test_extended.py`: stitch_continuity_lifetime, stitch_continuity_ms_luminosity,
stitch_continuity_pms_luminosity, blend_edges_continuous_post_ms. `test_table_physics.py`: lifetime_monotonic_in_mass,
lifetime_anchors_solar, zams_luminosity_anchors, pms_positive_and_short. `test_radiation.py`:
lbol_time_integral_vs_raw_track. `test_yields_raw.py`: preSN_mass_is_raw_track_final_mass. `test_solar.py`:
mist_grid_z_matches_headers.

**Known gaps.**
- **Z < 1e-4.** The table stops at 1e-4. Lookups clamp there. RIGEL and Imladris go far lower. PARSEC v2 has tracks
  at Z = 1e-6 and 1e-11 (2-2000 Msun), and the builder already knows them (`PARSEC2_Z_GRID`). MIST v2.5 reaches
  [Fe/H] = -4 (Z about 1.4e-6). Below 2 Msun at Z < 1e-6 there is no MIST track.
- **Above Z = 0.03.** PARSEC v2 tops out at 0.03, MIST v2.5 at 0.047. Tracks are held at the edge, so lifetimes and
  remnants at Z = 0.04 and 0.06 equal Z = 0.03. Remnant types at Z >= 0.03 also look wrong: every star above
  35 Msun is a CCSN up to 500 Msun (starforge integration report, §5). Not yet traced.
- **Pre-main sequence.** Only the delay t_PMS is stored. The PMS part of each track (L, R, Teff vs age) is cut off
  by `trim_to_ZAMS`. Needed for non-accreting PMS stars (STARFORGE option C).
- **Brown dwarfs and M < 0.09 Msun.** Not covered. Lookups clamp to 0.09 Msun. Only mass bookkeeping matters there.
- **Rotation.** None. Rotation can raise ionizing output by up to x5 (Leitherer+14, 1403.5444).
- **Post-MS physics.** Final masses at 20-40 Msun differ by about x2 between codes (Costa+25 Tab. 6). We carry one
  code's answer and no spread.
- **MIST/PARSEC handover.** Lifetimes differed by 10% at 2.2 Msun in the 2026-11-15 audit (B5). Recheck after the
  common-ZAMS fix.

**Candidates to fold in.**
- PARSEC v2 Z = 1e-6 and 1e-11 tracks (have the format; 2501.12917).
- MIST v2.5 [Fe/H] = -4 ... -2 for < 2 Msun (have the format).
- PARSEC v1.2S Z = 0.04, 0.06 (Chen+15, 1506.01681) for lifetimes above 0.03; check against MIST v2.5 at 0.045.
- Brown dwarfs: Baraffe+2015 (1503.04107) and ATMO 2020 (Phillips+2020); import L, R, Teff vs age.
- Rotating cross-check: Geneva grids (Ekström+12, 1110.5049, VizieR J/A+A/537/A146), PARSEC v2 ROT (2207.08642).
- Lifetime cross-check at very low Z: Geneva Z = 1e-5 (2407.06739), Murphy+2021 Pop III.

---

## 2. Winds (mass loss, speed, composition)

**Sources now.**
- Mass loss from the track (PARSEC v2 above 2 Msun). GIZMO injects winds only from 8 Msun.
- Composition: birth X_init plus wind enrichment. PARSEC v2 wind ejecta from 14 Msun, BoOST 8-14 Msun (blended over
  12-14), timed by the track's surface He excess.
- Speed: Lamers+95 v_inf / v_esc scaling (2.6 / 1.3 / 0.7), not from the tracks.

**Quality bar.**
- M_current + cumulative wind = M_init while alive (2e-7).
- End-of-wind element masses equal the PARSEC v2 wind ejecta enrichment, all 27 elements.
- Cumulative wind monotonic per element, never negative.

**Tests.** `test_invariants.py`: test_cum_wind_monotonic_total, test_cum_wind_endpoint, test_cum_wind_per_elem_monotonic,
test_M_current_plus_wind_ledger. `test_wind_composition.py`: all five tests (BoOST column map, Na/Mg surface,
end-of-wind vs PARSEC v2, Na/Al, lifetime budget). `test_extended.py`: wind_grid_matches_table.

**Known gaps.**
- **Winds below 8 Msun.** By design there is no time-resolved wind below 8 Msun. The AGB event carries the whole
  budget. The data exist (MIST `star_mdot`, PARSEC v2 `RATE`). A separate dataset would let a code choose.
- **Wind speed** is a scaling, not data. AGB speeds (72-459 km/s at 1 Msun in the table) are far above the observed
  10-20 km/s. Needs a check and probably an AGB rule.
- **Mass-loss rates.** PARSEC v2 uses Vink+01 for hot stars and de Jager for RSGs. Modern rates are x2-3 lower for
  OB stars (Björklund+21, 2008.06066) and up to x20 lower for RSGs (Beasor+20, 2001.07222). We carry no spread.
- **Na/Al time profile, 2-14 Msun.** PARSEC v2 wind ejecta start at 14 Msun and give only lifetime totals.
- **BoOST download location** is not confirmed (`fetch/README.md`).

**Candidates.**
- ORFEO LC18 time-resolved isotopic wind files, 13-120 Msun, [Fe/H] = 0 ... -3, v = 0/150/300 km/s (1805.09640).
  The best public source for rotating wind composition vs time.
- BoOST v_wind and wind kinetic energy (2004.08203) as data in place of the v_esc scaling, 9-500 Msun.
- Björklund+21 v_inf fits; Vink & Sander 2021 for hot stars.
- AGB wind speed: observed CO line widths of AGB envelopes as a clamp.

---

## 3. Yields: core-collapse SN, AGB, Type Ia

**Sources now.**
- CCSN, PPISN, PISN from 14 Msun: PARSEC v2 total minus wind ejecta (explosive part from LC03/CL04, Woosley 2017,
  Heger & Woosley 2002, mapped by core mass).
- 8-14 Msun: Limongi+2025 (2505.22030, solar), Z dependence as NuGrid net-yield difference.
- AGB < 8 Msun: Karakas 2010 below Z = 0.007, Karakas & Lugaro 2016 above, linear in log Z.
- Type Ia: W7 (Iwamoto+99), 1.378 Msun.
- Net yields against each source's own initial composition. 26Al counted as Mg.

**Quality bar.**
- Per-star total: wind + SN/AGB + remnant = M_init (2e-7).
- Per-element: wind + SN budget equals the PARSEC v2 lifetime net yield.
- ZAMS surface equals birth composition for every Z.
- K2010, KL16 and Limongi+25 values equal the raw files where the builder takes them unchanged.

**Tests.** `test_invariants.py`: closure SN and AGB, FSN rules, test_elem_SN_nonneg. `test_yields_raw.py`: all seven.
`test_wind_composition.py`: lifetime_budget_wind_plus_sn. `test_table_physics.py`: ia_yields_w7, imf_integrated_yields
(factor-10 window: too loose, see below).

**Known gaps.**
- **Z < 1e-4.** No AGB or 8-14 Msun yields there. PARSEC v2 has ejecta at 1e-6 and 1e-11 (>= 14 Msun only).
- **8-14 Msun away from solar Z.** Limongi+25 is solar only. The Z dependence is borrowed from NuGrid.
- **Super-AGB, 8-9 Msun.** Not covered by any source we use.
- **Type Ia delay times.** The table holds only W7 yields. The delay-time distribution (t^-1 from 40 Myr,
  1.3e-3 per Msun) lives in GIZMO (`resolvedism_fb_shared.h`) and in a builder dict that GIZMO does not read. Two
  places for one number. No Z dependence of Ia yields.
- **26Al.** Folded into Mg at once. Fine for Mg after 1 Myr. Lost as a tracer: no 26Al for gamma-ray maps or
  radiogenic heating, and KL16 gives only elemental Al, so its 26Al cannot be split out.
- **IMF-integrated check** accepts a factor 10. It should compare with published IMF yields to about 30%.
- **CCSN yields** above 14 Msun come from 2003/2004 FRANEC tables grafted onto PARSEC cores, only 6 CCSN rows per Z.

**Candidates.**
- ORFEO LC18 set R (Prantzos+18 adopted it), plus the 2024 [Fe/H] = +0.3 set (Roberti+24, 2403.18945) for super-solar Z.
- ORFEO 2012 Z = 0 and Heger & Woosley 2010 Pop III (0803.3161, via StarFit) for Z < 1e-4.
- Doherty+14 super-AGB (1310.2614), 6.5-9 Msun, Z >= 0.004.
- Low-Z AGB: Campbell & Lattanzio 2008, Gil-Pons+2013, FRUITY to [Fe/H] about -2.2 (1507.07338).
- Type Ia: Seitenzahl+13 N100 with Z variants (1211.3015), Leung & Nomoto 2018/2020 (Chandrasekhar and
  sub-Chandrasekhar, Z-dependent), Shen+18. DTD: Maoz & Graur 2017 (1703.04540), Maoz & Mannucci 2012.
- 26Al yields as their own isotope: PARSEC v2 26Al columns, LC18, Brinkman+2019/2021.

---

## 4. Remnants and explosion outcome

**Sources now.**
- WD below 8 Msun from Karakas final masses, else a linear IFMR.
- ECSN 8-10 Msun (NS 1.25 Msun). 10-14 Msun: Limongi pre-SN compactness (xi_2.5 > 0.35 -> BH).
- From 14 Msun: PARSEC v2 fate (CCSN, failed SN, PPISN, PISN, direct BH) and remnant mass.
- He and CO core masses stored per (Z, M).
- Explosion energy: GIZMO uses its own rule (1e51; 0 for failed SN and direct BH). The table attributes
  (ECSN 1e50, PISN 1e52) are not read.

**Quality bar.** Remnant type equals the raw PARSEC v2 fate at grid masses. Remnant mass within [0, M_preSN].
IFMR anchor at solar Z.

**Tests.** `test_yields_raw.py`: remnant_type_is_raw_parsec2_fate. `test_table_physics.py`: remnant_types_by_mass,
ifmr_anchor_solar. `test_invariants.py`: test_rem_mass_bounds, test_rem_type_integer, test_FSN_rem_equals_M_preSN.

**Known gaps.**
- Explosion energy is flat. Explodability has no islands (structure-based models give direct BHs at
  M_CO 6.6-7.2 Msun at solar Z; Maltsev+25, 2503.23856).
- Two sources of truth for SN energy (table attributes vs GIZMO).
- ECSN window is a fixed 8-10 Msun at every Z.
- Remnant types above Z = 0.03 (see tracks).
- IFMR not checked against an observed IFMR (e.g. Cummings+18, 1809.01673).

**Candidates.** Maltsev+25 explodability code (Apache-2.0) on our M_CO; Sukhbold+16 energies, Ni and remnant masses
(1510.04643, Garching archive); Schneider+21 (2008.08599, Zenodo CC-BY); PUSH at low Z (1910.08958).

---

## 5. Radiation and atmospheres

**Sources now.** `hotstar_stack`: TLUSTY OSTAR2002/BSTAR2006 15-55 kK; ATLAS9 / PHOENIX below; PoWR OB for He+;
PoWR WR and Götberg+18 stripped stars for WR phases; CMFGEN (Martins & Palacios 2017-2022, Martins+2025 VMS) for
Eddington-limit holes and the hot edge; blackbody only off every grid. Nine disjoint bands 0.4 eV to > 70 eV,
photon rates Q(H), Q(HeI), Q(He+), stored as age-cell averages.

**Quality bar.**
- Band integration of a blackbody gives the Planck photon rates to 1%.
- Composite bands equal the sum of atomic bands exactly. Bands never exceed L_bol.
- O stars: Q(H)/L within 0.15 dex of Martins+05 (CMFGEN) at every tabulated point, 0.05 dex on average; Q(HeI)/L
  within 0.5 dex (branch fix/late-o-ionizing).
- Lifetime Q(H) within 0.3 dex of the PARSEC v2 photon tables at 20-100 Msun.
- Q(He+) is order-of-magnitude only and flagged as such.

**Tests.** `test_radiation.py` (all). `test_table_physics.py`: photon_rates_nested, band_composites_exact,
hotstar_* (seams, TLUSTY nodes, He+ flag, WR, CMFGEN), radiation_source_recorded, radiation_no_silent_fallback,
radiation_matches_recorded_source, ionizing_photon_anchor. `test_extended.py`: lut_fractions_physical,
band_le_bolometric, ionizing_rate_rises_with_mass(_vms).

**Known gaps.**
- **Hot-star Q at late O.** Resolved as a PARSEC-side difference: our rates agree with Martins+05, TLUSTY and PoWR;
  PARSEC v2 uses WM-basic, 0.25 dex lower in Q(H) and 0.6 dex in Q(HeI) at 28-34 kK. See branch fix/late-o-ionizing.
  Left: Q(HeI) +0.47 dex at O9.5 V (steep regime), and no wind-atmosphere grid at Z < 0.01 Zsun for O stars.
- **IR band (< 0.4 eV).** Missing. The spectra we integrate (ATLAS9, PHOENIX, TLUSTY to 300 micron) contain it.
  One more band in `integrate_bands`. Needed for RT_INFRARED.
- **Q(He+)** spread is orders of magnitude between wind models (Roy+2025, 2501.08376). No fix possible from data.
- **Pre-main sequence and accreting protostars:** no radiation before ZAMS.
- **Z < 1e-4 atmospheres:** only TLUSTY Z = 0 and Larkin+23 ZAMS models (2210.09185).
- **B/A-star FUV and LW:** blackbody vs atmosphere difference not yet quantified.

**Candidates.** zYBC WM-basic with an Mdot axis down to Z = 1e-4 (2606.31244; download not verified) as a second
O-star family; NLTE-OBGRID TLUSTY 208 (better Lyman series, LW); PoWR SMC/sub-SMC OB grids for He+ at low Z;
Lorenzo+2025 (2501.07569) and Pauli+2026 (2608.25003) low-Z O-star calibrations as test anchors; pySTARBURST99
(2505.24841) and BPASS for population-level Q checks.

---

## 6. Chemistry tables (SGChem inputs)

**Sources now.** Energy table generated from `eh2-table.txt` (seam at 20000 K fixed). WSS09 CIE metal cooling
(md5 against the Leiden release). `H2-cooling-ratios.dat` (no recorded origin). SGChem net 18 tables from RTZ
sources (MIT) and the chem tables from SGChem DATA statements, both hash-pinned.

**Quality bar.** Byte-identical regeneration where a generator exists; checksums elsewhere; gamma within
[1.25, 5/3]; temperature round trip through GIZMO's Fortran and C++ lookups to 2e-3.

**Tests.** `test_chemistry.py` (all), `test_sgchem_data.py`, `run_lookup_test.py` (with --gizmo).

**Known gaps.** Origin of `H2-cooling-ratios.dat` and of `eh2-table.txt`. No physics check of the H2 cooling ratios
against their source. No metal cooling below 1e4 K in the WSS table (GIZMO suppresses it).

**Candidates.** Regenerate the H2 internal energy from the partition function (H2 level data); trace H2 cooling ratios to Glover & Abel 2008 (0803.1768).

---

## 7. Cross-cutting

- **Raw-data tests first.** Every number taken from a source gets a test against the raw file. Self-fingerprints
  only catch change, not error.
- **One source of truth.** SN energy and Ia DTD move into the tables, or out of them, not both.
- **Per-module provenance attribute** listing source, version, sha256 of inputs (exists for radiation and yields;
  extend to all).
- **Spread datasets.** Where codes disagree by more than the test tolerance (final mass, Q(He+), RSG winds), store a
  low/high variant or document the spread. Not now; after v1.

---

## 8. MESA validation set (proposal)

**Why.** Three uses: (1) check we read and stitch tracks right, by rerunning published physics; (2) fill regimes no
public grid covers; (3) measure the spread from physics choices we cannot otherwise quantify.

**Code and inputs.** One pinned MESA release and SDK, recorded like `requirements-build.txt`. Start from the MIST v2.5
work directories (Zenodo 15232687, CC-BY), so set A reruns MIST physics exactly.

**Set A: reproduce (does our pipeline match MESA/MIST?).**
- Masses 0.5, 1, 2, 2.5, 5 Msun. Z = 1e-4, 1e-3, 0.014, 0.03 (MIST [Fe/H] nodes). 20 runs.
- Physics: MIST v2.5 inlists unchanged, [a/Fe] = 0 and +0.4.
- Compare: lifetime, t_PMS, L and Teff at ZAMS and TAMS, final mass, cumulative mass loss, surface C, N, O, Na vs
  age. Tolerance: 1% lifetime, 0.02 dex L. These runs should match MIST; a miss means our pipeline or our MESA setup
  is wrong.

**Set B: massive stars, the PARSEC cross-check.**
- Masses 9, 12, 15, 20, 25, 40, 60, 100 Msun. Z = 1e-4, 1e-3, 0.004, 0.014, 0.03. 40 runs, to core C depletion.
- Physics: Vink+01 and de Jager (PARSEC-like) as the base; core overshoot f_ov = 0.016.
- Variants on 15, 20, 40 Msun at all five Z (15 runs each): Björklund+21 OB rates with Beasor+20 RSG rates;
  larger overshoot (f_ov = 0.03).
- Network with NeNa and MgAl chains and 26Al (to be chosen from the MESA nets; must contain na23, mg24-26, al26,
  al27).
- Compare: lifetimes, M_final, M_He, M_CO, WR onset, cumulative wind per element, surface Na/Al vs age. Then run our
  atmosphere stack on the MESA tracks and compare lifetime Q(H), Q(HeI), L_FUV with the table at the same (M, Z).

**Set C: gap fillers.**
- Z = 1e-5, 1e-6 for 1-100 Msun (12 masses): fills Z < 1e-4 below 2 Msun and checks PARSEC v2 at 1e-6.
- Z = 0.04, 0.06 for 1-100 Msun (12 masses): replaces the edge clamp.
- 2-14 Msun at six Z with the Na/Al network: time-resolved wind composition where PARSEC v2 has none.
- 0.08-0.1 Msun and PMS phases from the birthline for 0.1-5 Msun (PMS tables).
- About 100 runs.

**Not in scope.** Core collapse and explosions (use Sukhbold+16, Maltsev+25); atmospheres (MESA gives none);
rotation (later).

**Cost estimate.** Rough planning numbers, to be replaced by a pilot (one 1 Msun and one 20 Msun run timed first).
- Low mass to WD, 4 threads: about 2-5 core-hours per run.
- Massive to C depletion, larger network, 4 threads: about 10-40 core-hours per run.
- Set A about 100 core-h; set B (40 + 30 variants) about 2000 core-h; set C about 2000 core-h.
- Total about 4000-5000 core-hours, i.e. about 40-50 node-hours on 112-core Orion nodes, run as sbatch arrays.
- Storage: histories plus sparse profiles, under 50 GB on ptmp. Only the reduced histories enter the repo.
- People time is the real cost: about two weeks for setup and set A, two more for B and C, plus analysis.

**Order.** Pilot -> set A (it validates everything after) -> set B base -> variants -> set C.

---

## 9. Independent audits, module by module

Each audit is done by someone (or an agent) who did not write the module. They get the raw data and the code, not our
conclusions. They recompute the key numbers from the raw files and write a short report into `docs/audits/`. Every
finding becomes a failing test first, then a fix.

Order and why:

1. **Yields (SN, AGB, Ia).** Most past bugs were here (net-yield scaling, Ni, P, K2010 double model, 26Al). Feeds the
   gas composition directly.
2. **Winds.** Rewritten in fix/wind-composition; timing rule and BoOST blending are new and have only self-consistency
   tests plus the PARSEC v2 end-point test.
3. **Remnants and fates.** Open high-Z remnant puzzle; ECSN window; two SN-energy sources.
4. **Tracks and stitching.** Mostly raw-data-tested already; check the 2.0-2.5 Msun handover, the Z interpolation
   between PARSEC grids, and the edge clamp.
5. **Radiation.** Late-O checked against Martins+05 on fix/late-o-ionizing. Still to audit: WR and VMS phases (PoWR
   transformed radius, CMFGEN fill), B/A-star FUV and LW, Gotberg branch at Z < 0.1 Zsun.
6. **Chemistry tables.** Small and checksum-locked; provenance of two inherited files.
7. **Build and reproducibility.** Pinned environment, fingerprints, reduced catalogues vs full downloads.

## 10. Order of work (suggestion; Uli decides)

1. Merge fix/late-o-ionizing if accepted (tests and docs only).
2. Audits 1-3.
3. Type Ia DTD into the tables; one SN-energy source.
4. IR band; PMS tables.
5. Z < 1e-4 from PARSEC v2 1e-6 / 1e-11, MIST v2.5, Pop III yields.
6. MESA pilot and set A.
7. Z > 0.03 (PARSEC v1.2S, MESA set C), brown dwarfs, 26Al as an isotope.
