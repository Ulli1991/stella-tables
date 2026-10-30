# GIZMO-side changes for the SGChem data in stella-tables

GIZMO is not changed by this branch. Everything below is for a GIZMO commit after this branch is pushed and has a commit id.

1. `gizmo_side_installed_tables.diff` (`patch -p1` in the GIZMO root; dry-run clean against stella_dev 6690d3b7):
   - `scripts/get_stella_tables.sh`: `sgchem_net18_tables.bin` moves from OPTIONAL to the required list, `sgchem_chem_tables.bin` is added, the
     frozen reference data are copied to `PREFIX/sgchem_reference_data` and checked with `sha256sum -c data.sha256`; `--into-tree` links both .bin files.
   - `galaxy_sf/stella_tables_check.h`: new hash entry `sgchem_chem_tables.bin` (2d623ca6...); set `STELLA_TABLES_REF` to the new stella-tables commit
     (also the pin in the script header; the placeholder in the diff is `<new stella-tables commit>`).
   - `core/begrun.cc`: startup check requires `sgchem_chem_tables.bin` for every C++ network and `sgchem_net18_tables.bin` for network 18.
   - The net18 hash is unchanged (e3c41f2c...): the new builder reproduces the installed file byte for byte.
2. `chem_tables_data.h` (1385 lines, five translation units include it as `static const double dat_*[N]`) becomes an installed file:
   - format of `sgchem_chem_tables.bin`: `"SGCT0001"`, uint64 N, then N records `{char name[48]; uint64 n; float64 data[n]}`, little-endian, 45 arrays.
   - `chemistry/sgchem/make_chem_tables.py OUT.bin --loader DIR` writes the GIZMO-side pair for it: `chem_tables_data.h` (extern declarations, same
     identifiers and sizes as today, so the .cc users compile unchanged) and `chem_tables_data.cc` (definitions + `int chem_tables_data_load(const char* dir)`
     that checks magic, count, names, sizes and order). Compile-tested standalone with g++ -Wall.
   - GIZMO to do: add those two (generated, checked-in) files in place of the 1385-line header, add `chem_tables_data.cc` to the cooling/sgchem/cxx
     sources, call `chem_tables_data_load(All.SGCHEMTablesDir)` once next to `n18_tables_load` (and in the Fortran-backend helper that includes the header),
     and delete `tables_extract.py` plus its Fortran inputs from GIZMO (they live in `chemistry/sgchem/chem_tables_source/` here).
   - not static const any more: the arrays are written once at start-up, read-only afterwards.
3. `make_net18_tables.py` (header generator) and its RTZ inputs leave GIZMO the same way; `chemistry/sgchem/net18_source/` holds them with RTZ's MIT LICENSE.
4. `fortran_reference/data.sha256` stays in GIZMO (identical copy here); point `SGCHEM_REFERENCE_DATA` at `PREFIX/sgchem_reference_data`.
