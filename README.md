# stella-tables

Builds every input table the Stella feedback model in GIZMO reads, from versioned inputs, and tests them.

```
module load python-waterboa/2026.07          # MPCDF: the pinned environment (requirements-build.txt)
./build.sh                                   # stellar + chemistry tables into build/, then the full test suite
STELLA_ALLOW_ANY_ENV=1 ./build.sh            # any other environment (refused otherwise); last-bit differences
GIZMO_SRC=/path/to/gizmo ./build.sh          # also test the energy table through GIZMO's Fortran and C++ lookups
./build.sh --raw RAW_ROOT                    # re-derive the reduced track catalogues from the full downloads first
python3 tests/run_all.py BUILD_DIR [--gizmo GIZMO_SRC]   # tests only
```

Needs the Python environment in `requirements-build.txt` (Python 3.14.6, numpy 2.4.6, scipy 1.18.0, h5py 3.16.0;
MPCDF: `module load python-waterboa/2026.07`); `python3 check_build_env.py` lists any mismatch, and `build.sh`
stops on one unless `STELLA_ALLOW_ANY_ENV=1`. Prebuilt tables: the GitHub releases (`tables-<commit>`). The `--gizmo` lookup test also needs `mpifort` and `g++` (MPCDF:
`module load gcc openmpi gsl hwloc fftw-mpi hdf5-mpi`).

## Reproducible builds

The tables are bit-reproducible only in one environment: numpy's exp, log, pow and sums, Python's float `sum` (3.12+)
and scipy/LAPACK differ in the last bits between versions and CPU code paths. `build.sh` therefore checks the
environment against `requirements-build.txt` and refuses others (`STELLA_ALLOW_ANY_ENV=1` builds anyway; the stellar
table then differs from the published one in the last bits), runs single-threaded, and limits numpy to its baseline
CPU code paths (`NPY_DISABLE_CPU_FEATURES`), so the pinned versions give the same bits on every x86-64 machine. Python
float sums use `math.fsum`. The energy table is installed from `chemistry/reference/energy-table-x.txt.xz` (generator
output in the original environment, the version GIZMO pins); each build regenerates it and the tests compare (exact
in that environment, within 1e-13 with scipy 1.18). The table attribute `build_environment` records the versions.

## Use from GIZMO

A GIZMO build with `SGCHEM` or `GALSF_RESOLVEDISM_STELLAR_TABLES` runs `scripts/get_stella_tables.sh` before
linking. If any of the three generated tables is missing there, it clones this repo, runs `build.sh` (tests
included), and installs the outputs only if every suite passes; otherwise it does nothing. Settings:

```
make STELLA_TABLES_REPO=<git URL> [STELLA_TABLES_REF=<tag or commit>] [STELLA_TABLES_DIR=<clone dir>]
```

(or the same names as environment variables). By default it clones https://github.com/Ulli1991/stella-tables at the
commit pinned in GIZMO's `scripts/get_stella_tables.sh`, so each GIZMO version builds the tables it was tested with.

## Outputs (`build/`)

| file | read by GIZMO from | built by |
|---|---|---|
| `stellar_tables_unified.hdf5` | `galaxy_sf/stellar_tables/` | `stellar/build_stellar_tables.py` |
| `stellar_wind_cumulative.hdf5` | `galaxy_sf/stellar_tables/` | same (written next to the main table) |
| `chemistry/energy-table-x.txt` | `cooling/sgchem/` (`cooling/chemcool/` before the rename) | `chemistry/make_energy_table.py` |
| `chemistry/eh2-table.txt`, `H2-cooling-ratios.dat`, `wss_z_collis.txt` | `cooling/sgchem/` (`cooling/chemcool/` before the rename) | copied from `chemistry/tables/` (inherited, see below) |

## SGChem data (chemistry/sgchem/)

| file | built by | from |
|---|---|---|
| `sgchem_net18_tables.bin` (network 18 rates, 44 float64 arrays) | `make_net18_tables.py` | `net18_source/` (RTZ C sources + data, MIT) |
| `sgchem_chem_tables.bin` (45 arrays, former `chem_tables_data.h`) and legacy `chem_tables_data.h` | `make_chem_tables.py` | `chem_tables_source/` (SGChem Fortran DATA statements) |
| `reference_data/` + `data.sha256` | frozen, not built | one-zone / cooling-function references for nets 3 5 15 16 17 18 20 (10 MB) |

`build.sh` writes the two .bin files and the header to `build/chemistry/`; `tests/chemistry/test_sgchem_data.py` checks they equal the pinned hashes in
`chemistry/sgchem/expected_sha256.txt` (net18 = the hash GIZMO expects), the .bin round-trips to the header, the reference data match `data.sha256`;
with `--gizmo` also GIZMO's own hash entries and files. GIZMO-side changes (hash entries, required net 18, loader for the chem tables, install script):
`docs/GIZMO_SIDE_NOTES.md` and `docs/gizmo_side_installed_tables.diff`. Licence: BSD-3 for own code, see `LICENSE`; third-party parts and their terms in `NOTICE`.

## Layout

```
stellar/                      table builder and helpers
  data/reduced/*.tar.xz       track catalogues cut down to what the builder reads (see below)
  data/reduced/MANIFEST.txt   sha256 of every raw file the builder reads; names of the files it only lists
  data/yields/                Chieffi & Limongi 2018, Karakas 2010, Karakas & Lugaro 2016, NuGrid
  data/lut/spectral_band_lut.hdf5   band-energy fractions (stellar/build_spectral_lut.py)
chemistry/
  make_energy_table.py        SGCHEM energy -> (T, gamma) table from eh2-table.txt
  tables/                     inherited SGCHEM tables
  reference/energy-table-x.gizmo-original.txt.xz   energy table GIZMO shipped before the 20000 K seam fix
fetch/                        how to get the raw catalogues
tests/run_all.py              runs every suite below
```

## Reduced track catalogues

The full MIST / PARSEC v2 / PARSEC v1.2 / BoOST downloads are about 24 GB; the builder reads 1816 files
`stellar/reduce_tracks.py` keeps only those files, and in the MIST `.track.eep` and PARSEC v2 `.TAB`
tracks only the columns whose names appear in `build_stellar_tables.py`, as the original text, so the builder
parses identical numbers. Files the builder only discovers by name become empty files. Building from the reduced
catalogues gives tables identical to building from the full downloads (checked per dataset; see
`tests/stellar/reference_fingerprints.json`).

## Track stitching

MIST v2.5 (below 2 Msun; `--mist v1.2` for the old release) and PARSEC v2 (from 2 Msun) are each interpolated in
log Z between their bracketing metallicity grids, in fractional-lifetime space, and blended linearly in log M over
2.0-2.5 Msun (`MASS_STITCH_P2_LO`, `MASS_STITCH_P2_FULL`). Truncated catalogue tracks are skipped in the mass
bracketing (MIST v2.5 [Fe/H] = +0.5 ends its 0.10-0.45 Msun tracks before the main sequence; the next grid down is
used). MIST is read up to 3 Msun (`MIST_M_MAX`). `track_source` marks the dominant catalogue (blend weight >= 0.5).
Both catalogues start the table age at the same ZAMS, central H down by 1e-4, so t_PMS and lifetimes join across the
blend. Grid Z = the Zinit of each grid's track headers. MIST v2.5 [a/Fe] (`--mist-alpha-fe`,
`STELLA_MIST_ALPHA_FE`): default p4 = +0.4, p0 = solar-scaled selectable; at fixed total Z, +0.4 changes MS lifetimes
by <= 1% for Z <= 1e-3 and by 2-8% for 1-2 Msun at Z = 0.01. The reduced catalogue carries the p4 and p0 tracks the
builder reads.

## Radiation

Band luminosities and ionizing photon rates of every track point come from one explicitly chosen source
(`build_stellar_tables.py --radiation-source`, `STELLA_RADIATION_SOURCE` for `build.sh`), recorded in the table
attributes (`radiation_source*`, `atlas9_lookup_option`, `spectral_lut_sha256`).

- `planck`: blackbody at Teff.
- `hotstar_stack`: model atmospheres per band (`stellar/hotstar_radiation.py`, LUT `stellar/data/lut/hotstar_lut.hdf5`
  from `stellar/build_hotstar_lut.py`): TLUSTY BSTAR2006/OSTAR2002 for 15-55 kK (blended with ATLAS9 over 15-20 kK and
  into a blackbody over 55-65 kK), ATLAS9 below, He+-ionizing photons of O/B stars from the PoWR OB Vd3 log Q(HeII)
  tables (order of magnitude), WR phases (surface X < 0.4, log Teff > 4.45) from the PoWR WR grids through the
  transformed radius (Z >= 0.1 Zsun, lowest grid 0.07 Zsun) or from Gotberg+2018 stripped stars (Z < 0.1 Zsun),
  CMFGEN wind models (Martins & Palacios 2017, 2021, 2022; Martins+2025 very massive stars with Z-scaled winds; POLLUX;
  0.01-1 Zsun, scattered along their tracks, linear on a Delaunay mesh in log Teff and log g) take over where the stack
  above used filled TLUSTY cells (Eddington-limit holes at low log g) or the hot-edge blackbody (above 55 kK); up to
  0.25 dex in log g outside a set's models the nearest edge value, and where the bracketing metallicities have no model
  the next lower (then higher) one, e.g. LMC/SMC very massive stars at solar Z - both counted as clamped;
  blackbody only off every grid. Coverage per band (model / blended / clamped / blackbody) in `radiation_coverage_*`,
  the Q_H0 share from CMFGEN in `radiation_coverage_Q_H0_cmfgen_frac`.
- `atlas9_nonionizing`: ionizing bands and photon rates (> 13.6 eV) blackbody; non-ionizing bands from the LUT below
  at the point's log g and Z, blackbody where the LUT has no model (as `masked_planck`).
- `atlas9`: `stellar/data/lut/spectral_band_lut.hdf5` (ATLAS9 Castelli-Kurucz >= 3500 K, PHOENIX-ACES below), trilinear
  in log Teff, log g, log Z ([M/H] = log Z/Z_sun). The LUT spans 3.5-50 kK and log g 0-5 and has no models for hot
  low-gravity stars; its empty cells were filled with the nearest model when it was built (`has_data` marks the
  real ones). `--atlas9-lookup` (`STELLA_ATLAS9_LOOKUP`):
  - `legacy`: log g = 4 and Z = Z_sun for every point, linear extrapolation off the grid (every build before this option),
  - `clamp`: the point's log g (from track mass and radius) and initial Z, clamped to the grid edges,
  - `masked_planck`: the point's log g and Z; blackbody wherever an interpolation corner is a filled cell or the
    point is off the grid.

  `radiation_coverage_*` attributes give the share of track points (and of their blackbody Q_H0, L_FUV, L_LW)
  that hit only models, a filled cell, or nothing.

`build.sh` uses `hotstar_stack` (log g and Z dependence, hot and WR stars); `atlas9` / `legacy` was the default up to
the tables GIZMO pinned at a4830ab, whose attribute wrongly said "Planck".

Band integrals of every model spectrum include the partial wavelength intervals at the band edges
(`build_spectral_lut.integrate_bands`). Radiation is stored as the time average over each age-grid cell (edges at the
log-age midpoints, cut at death), so short post-main-sequence phases keep their time integral.

Every source gives each track point disjoint atomic bands from one spectrum: energy bands OPT_NIR 0.4-3.4, NUV_lo
3.4-6, NUV_hi 6-8, FUV_mid 8-11.2, LW 11.2-13.6, ion_H0 13.6-24.6, ion_He0 24.6-54.4, ion_He1 54.4-70, ion_He2 > 70 eV,
and photon rates Q_ion_H0 13.6-24.6, Q_ion_He0 24.6-54.4, Q_He1 > 54.4 eV. Only these (and L_bol) are interpolated onto
the table grid; the composites L_FUV (6-11.2), L_FUV_total (6-13.6), L_NUV (3.4-8), L_FUV_M1 (8-13.6), L_ion_tot,
Q_He0 and Q_ion are written as the exact sums of their atomic bands at every node, so e.g.
10^log_L_FUV_total - 10^log_L_FUV_M1 = 10^log_L_NUV_hi (the 6-8 eV bin).

## Wind composition

`stellar_wind_cumulative.hdf5` (`elem_ej_wind_cumulative`) is what GIZMO injects as winds (from 8 Msun up). The
wind mass per age step is the track's mass loss, down to the track's final mass (mass lost after the last alive
age-grid point goes into the last interval). For M >= 8 the composition is the birth composition `X_init` plus the
wind enrichment of the source star per unit wind mass, `ejecta_k / M_wind - X_source,init_k`, for all 27 elements:
PARSEC v2 wind ejecta (`winds_ejecta.dat`, 14-600 Msun, nearest PARSEC v2 Z, linear in M) for M >= 14; below, BoOST
(non-rotating tracks, surface H..Fe x mass loss over the whole track; the 15 elements BoOST lacks, P, S, Cl, Ar, K,
Ca, Sc, Ti, V, Cr, Mn, Co, Ni, Cu, Zn, stay at birth fraction), interpolated in log Z between its grids and blended
linearly in log M into PARSEC v2 over 12-14 Msun; above the BoOST metallicities (from 2 x its MW grid) PARSEC v2 at
14 Msun. The timing follows the track: x_k(t) = X_init,k + net_k f(t) / sum f dM, f = the track's surface He
excess scaled to peak 1 (0 at ZAMS), raised where needed so that no element goes negative; the sum over elements is 1
at every step and the end-of-wind totals are X_init W + net exactly. An element is at most fully absent from the
wind (net >= -X_init W, rest to H). Below 8 Msun the table has no wind, as GIZMO injects none there: M_current =
M_init during the star's life, M_preSN = M_init, and the AGB event carries the whole Karakas lifetime budget.
`surface_abundances` is X_init plus the track's own surface change since ZAMS; for M >= 8 it holds the wind fractions.
MIST tracks end at the hottest post-AGB point (MIST v2.5 continues through white-dwarf cooling).

## Yields

Lifetime element budget per star: wind + SN (or AGB) ejecta = X_init (M - rem) + the source's net yield; the SN
gets X_init Mej + (lifetime net - wind net). Net yields are taken against each source's own initial composition.
Per-element closure is done on absolute ejecta: a source's sum of net yields that is not 0 (Li/Be, rounding,
transcription) goes to H, negative elements are clipped (logged; the build stops if H runs out).

- PARSEC v2 (>= 14 Msun): total and wind ejecta files. An isotope column that is exactly 0 while the isotope is
  present initially is not filled by PARSEC v2 (P in every row; Si29, P, Cl, K, Sc, V, Mn, Co, Cu in all wind, FSN,
  PPISN and DBH rows) and is set to initial fraction x ejected mass (net 0). NI is stable Ni (FE already holds the
  decayed 56Ni). Where the two bracketing masses have different explosion types, the nearer model is used whole.
- Limongi+2025 (8-14 Msun, solar Z, lifetime totals), Z dependence added as NuGrid net yield per ejected mass at
  Z x 0.02/0.014 minus that at 0.02 (log-Z interpolated between NuGrid grids, linear in M over 12-25 Msun);
  elements Limongi lacks from the NuGrid net yield.
- AGB (< 8 Msun): Karakas 2010 below Z = 0.007, Karakas & Lugaro 2016 from 0.007, linear in log Z between their
  metallicities. K2010 Z = 1e-4, 6 Msun: the VW93 model (the Reimers variant is skipped). KL16 initial composition
  from the files' [X/H] (X(i) 10^-[X/H]), scaled so Ca, which AGB stars do not process, has net 0. GIZMO runs no
  winds below 8 Msun, so the AGB composition is built on M_init - rem and scaled to M_preSN - rem.
- 26Al (half-life 0.72 Myr) is counted as Mg in every source that lists it as an isotope (all but KL16, which gives Al
  as an element).

## Solar composition

`stellar/solar_abundances.py` is the single source for the tables and for GIZMO (`core/init.cc`,
`python3 stellar/solar_abundances.py --c` prints the array): Asplund+09 Table 1 photospheric ratios with standard
atomic weights, metals C..Zn normalized to Z_sun = 0.014, Y = 0.2485 (A09 helioseismic), X = 1 - Y - Z. The table's
X_init is this pattern scaled by Z/0.014, so it is exactly the composition GIZMO gives the gas.

## Chemistry tables

`energy-table-x.txt` is tabulated in (ekn, x) with x = abh2 / (1 + abhe - abh2), so one table holds for any
He abundance. It is generated from `eh2-table.txt` (the H2 internal energy per molecule). That inherited table
has two defects the generator handles:

* the row at T = 20000 K holds 1.5 T instead of about 3 T; it is dropped;
* at 20000 K it switches from a partition-function energy to analytic 3T, with a 45 K step. Splining across the
  step made gamma ring between 1.05 and 2.27 for 19998-20002 K. The generator now splines each side separately;
  the table changes only at 19972-20026 K (1336 of 1050000 rows), is byte-identical everywhere else, and gamma
  stays within [1.25, 5/3]. `--legacy-seam` reproduces the original GIZMO table byte for byte, and a copy of it
  is kept in `chemistry/reference/`; the tests compare against it.

`wss_z_collis.txt` is the collisional-ionization-equilibrium metal cooling table of Wiersma, Schaye & Smith
(2009), bit-identical to the official release (md5 checked against the Leiden distribution, 2026-07-28). It came
to GIZMO through the TNT code (Hu et al.); GIZMO uses the nine metal columns (C, N, O, Ne, Mg, Si, S, Ca, Fe) and
suppresses them below 1e4 K with exp(-(1e4/T)^2) when loading (`wss_z_collis.F`). `H2-cooling-ratios.dat` has no
recorded origin. Both are checksum-locked and sanity-checked, not regenerated.

## Tests

| suite | what it checks |
|---|---|
| `tests/stellar/test_invariants.py` | mass and wind closure under interpolation, monotonic cumulative winds, no NaN at the grid edges |
| `tests/stellar/test_table_physics.py` | every dataset GIZMO's loader opens, physics sanity, provenance, composite bands = sums of atomic bands; `--compare A B` diffs two tables |
| `tests/stellar/test_extended.py` | lifetime and L_bol continuity where the track source changes, solar constants, scaled-solar mixture, spectral LUT fractions, band luminosities below L_bol, ZAMS Q_ion rising with mass, table/wind grid agreement, per-dataset fingerprints against the reference build |
| `tests/stellar/test_wind_composition.py` | BoOST isotope columns vs `fdat-columns.txt`; end-of-wind element masses = PARSEC v2 wind ejecta enrichment for all 27 elements; Na/Al wind enrichment of massive stars; per-star wind + SN element budget = PARSEC v2 lifetime net yield (raw inputs from `BUILD_DIR/data_root`) |
| `tests/stellar/test_yields_raw.py` | pre-SN mass = raw PARSEC v2 track final mass; ZAMS surface = X_init; K2010 Z = 1e-4 6 Msun = VW93 model; KL16 unprocessed elements follow the raw production factors; Limongi+2025 lifetime budget at solar Z (raw inputs from `BUILD_DIR/data_root`) |
| `tests/stellar/test_solar.py` | A09 transcription reproduces A09's Z and Z/X, pattern sums, table X_init = pattern x Z/0.014, GIZMO `init.cc` array equals the generator (`--gizmo`) |
| `tests/chemistry/test_chemistry.py` | checksums, eh2 physics (1.5T at low T, c_v bounds, glitch row, seam), energy table regenerates byte for byte, layout, x = 0 block exact, gamma bounds, monotonicity, new vs original table |
| `tests/chemistry/run_lookup_test.py` | built and original energy tables through GIZMO's Fortran `calc_temp.F` and the C++ port: T round trip to 2e-3 for abhe 0 / 0.079 / 0.1 / 0.12, backends agree, original and new identical off the seam |

Expected failures (XFAIL) are listed in the test output with their reason and do not fail the run.

## Known limitations

* MIST v2.5's own tracks peak in post-MS duration at the He-flash transition near 2 Msun (+0.11 dex from 1.95 to
  2.0 Msun at [Fe/H] = -0.75, [a/Fe] = +0.4), which the table reproduces (XFAIL `blend_edges_continuous_post_ms` in `test_extended.py`).
* Above the catalogue metallicity ranges (PARSEC v2 Z > 0.03, MIST v2.5 Z > 0.047) tracks are held at the edge, so
  lifetimes and remnants at Z = 0.04, 0.06 equal those at the top grid; yields and winds follow each Z's X_init.
* The BoOST download location is not confirmed (`fetch/README.md`); BoOST enters only the wind composition of
  8-14 Msun stars.
* Q(He+) of massive stars differs from the PARSEC v2 photon tables by up to 3 dex in either direction (different
  atmospheres; the literature spread is orders of magnitude), printed by `tests/stellar/test_radiation.py`.
