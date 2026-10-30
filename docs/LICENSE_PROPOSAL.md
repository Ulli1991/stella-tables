# LICENSE: BSD-3 for own code (decided); third-party parts keep their terms

Needs Uli's choice. What bears on it:
- own code (builders, tests, chemistry/make_energy_table.py): free choice (MIT or BSD-3 are the permissive options; GPL-3 if copyleft is wanted).
- `chemistry/sgchem/net18_source/`: RTZ by HarleyKatz, MIT (LICENSE kept in that directory; MIT must stay with it).
- `chemistry/sgchem/chem_tables_source/` and `chemistry/tables/`: SGChem Fortran/data by Simon Glover (and inherited CHEMCOOL tables), GPL as distributed in Arepo (no GPL version is named in Arepo's SGChem sources); derived files keep that licence. See the NOTICE in `chem_tables_source/`.
- `stellar/data/`: third-party track catalogues and yields (MIST, PARSEC, BoOST, Chieffi & Limongi, Karakas, NuGrid): each has its own terms/citation request; a repo licence would not cover them.
A LICENSE file would therefore cover "own code" only, with a NOTICE listing the third-party parts. 
