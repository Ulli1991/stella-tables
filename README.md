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

The full MIST / PARSEC v2 / PARSEC v1.2 / BoOST downloads are about 24 GB; the builder reads 1170 files
(4.8 GB). `stellar/reduce_tracks.py` keeps only those files, and in the MIST `.track.eep` and PARSEC v2 `.TAB`
tracks only the columns whose names appear in `build_stellar_tables.py`, as the original text, so the builder
parses identical numbers. Files the builder only discovers by name become empty files. Building from the reduced
catalogues gives tables identical to building from the full downloads (checked per dataset; see
`tests/stellar/reference_fingerprints.json`).

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

`H2-cooling-ratios.dat` and `wss_z_collis.txt` have no recorded origin. They are checksum-locked and
sanity-checked, not regenerated.

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

* Below 2 Msun the builder takes the nearest MIST [Fe/H] for each table Z instead of interpolating (up to
  about 0.1 dex mismatch; MIST [Fe/H] = -1.00 is never used).
* Source URLs for the PARSEC v2, PARSEC v1.2 and BoOST downloads are not recorded yet (`fetch/README.md`).
