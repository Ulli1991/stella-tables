# stella-tables

Builds every input table the Stella feedback model in GIZMO reads, from versioned inputs, and tests them.

```
./build.sh                                   # stellar + chemistry tables into build/, then the full test suite
GIZMO_SRC=/path/to/gizmo ./build.sh          # also test the energy table through GIZMO's Fortran and C++ lookups
./build.sh --raw RAW_ROOT                    # re-derive the reduced track catalogues from the full downloads first
python3 tests/run_all.py BUILD_DIR [--gizmo GIZMO_SRC]   # tests only
```

Needs python3 with numpy, scipy, h5py. The `--gizmo` lookup test also needs `mpifort` and `g++` (MPCDF:
`module load gcc openmpi gsl hwloc fftw-mpi hdf5-mpi`).

## Use from GIZMO

A GIZMO build with `CHEMCOOL` or `GALSF_RESOLVEDISM_STELLAR_TABLES` runs `scripts/get_stella_tables.sh` before
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
| `chemistry/energy-table-x.txt` | `cooling/chemcool/` | `chemistry/make_energy_table.py` |
| `chemistry/eh2-table.txt`, `H2-cooling-ratios.dat`, `wss_z_collis.txt` | `cooling/chemcool/` | copied from `chemistry/tables/` (inherited, see below) |

## Layout

```
stellar/                      table builder and helpers
  data/reduced/*.tar.xz       track catalogues cut down to what the builder reads (see below)
  data/reduced/MANIFEST.txt   sha256 of every raw file the builder reads; names of the files it only lists
  data/yields/                Chieffi & Limongi 2018, Karakas 2010, Karakas & Lugaro 2016, NuGrid
  data/lut/spectral_band_lut.hdf5   band-energy fractions (stellar/build_spectral_lut.py)
chemistry/
  make_energy_table.py        CHEMCOOL energy -> (T, gamma) table from eh2-table.txt
  tables/                     inherited CHEMCOOL tables
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

MIST (below 2 Msun) and PARSEC v2 (from 2 Msun) are each interpolated in log Z between their bracketing
metallicity grids, in fractional-lifetime space, and blended linearly in log M over 2.0-2.5 Msun
(`MASS_STITCH_P2_LO`, `MASS_STITCH_P2_FULL`). Truncated catalogue tracks are skipped in the mass bracketing.
`track_source` marks the dominant catalogue (blend weight >= 0.5). MIST grids sit at their own initial Z
(track headers: Zinit = 1.42857e-2 x 10^[Fe/H]), PARSEC v2 at the Z of its grid.

## Radiation

Band luminosities and ionizing photon rates of every track point come from one explicitly chosen source
(`build_stellar_tables.py --radiation-source`, `STELLA_RADIATION_SOURCE` for `build.sh`), recorded in the table
attributes (`radiation_source*`, `atlas9_lookup_option`, `spectral_lut_sha256`).

- `planck`: blackbody at Teff.
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
  that hit only models, a filled cell, or nothing. Default `atlas9` / `legacy`.

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
| `tests/stellar/test_table_physics.py` | every dataset GIZMO's loader opens, physics sanity, provenance; `--compare A B` diffs two tables |
| `tests/stellar/test_extended.py` | lifetime and L_bol continuity where the track source changes, solar constants, scaled-solar mixture, spectral LUT fractions, band luminosities below L_bol, ZAMS Q_ion rising with mass, table/wind grid agreement, per-dataset fingerprints against the reference build |
| `tests/stellar/test_solar.py` | A09 transcription reproduces A09's Z and Z/X, pattern sums, table X_init = pattern x Z/0.014, GIZMO `init.cc` array equals the generator (`--gizmo`) |
| `tests/chemistry/test_chemistry.py` | checksums, eh2 physics (1.5T at low T, c_v bounds, glitch row, seam), energy table regenerates byte for byte, layout, x = 0 block exact, gamma bounds, monotonicity, new vs original table |
| `tests/chemistry/run_lookup_test.py` | built and original energy tables through GIZMO's Fortran `calc_temp.F` and the C++ port: T round trip to 2e-3 for abhe 0 / 0.079 / 0.1 / 0.12, backends agree, original and new identical off the seam |

Expected failures (XFAIL) are listed in the test output with their reason and do not fail the run.

## Known limitations

* MIST v1.2's own tracks have a ~4% lifetime uptick at one mass step near 1.9 Msun ([Fe/H] = -0.5, -0.25),
  which the table reproduces (XFAIL in `test_table_physics.py`).
* Above the catalogue metallicity ranges (PARSEC v2 Z > 0.03, MIST Z > 0.0425) tracks are held at the edge.
* The BoOST download location is not confirmed (`fetch/README.md`); no BoOST track enters the current build.
