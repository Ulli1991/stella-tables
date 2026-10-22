#!/bin/bash
# Build every table Stella/GIZMO reads, then run the test suite.
#   ./build.sh                 build from the versioned reduced catalogues (stellar/data/reduced/*.tar.xz)
#   ./build.sh --raw RAW_ROOT  re-derive the reduced catalogues from the full downloads first (see fetch/)
#   GIZMO_SRC=/path/to/gizmo ./build.sh   also test the energy table through GIZMO's Fortran and C++ lookups
#   STELLA_RADIATION_SOURCE=atlas9|planck|atlas9_nonionizing|hotstar_stack (default atlas9), STELLA_ATLAS9_LOOKUP=legacy|clamp|masked_planck
#     (default legacy): band luminosities, see README "Radiation"
# Outputs land in build/: stellar_tables_unified.hdf5, stellar_wind_cumulative.hdf5, chemistry/energy-table-x.txt
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/build"; ROOT="$OUT/data_root"
mkdir -p "$OUT/chemistry" "$ROOT"

if [ "${1:-}" = "--raw" ]; then
    python3 "$HERE/stellar/reduce_tracks.py" "$2" "$ROOT"
else
    for t in "$HERE"/stellar/data/reduced/*.tar.xz; do tar -C "$ROOT" -xJf "$t"; done
fi
for y in chieffi_limongi2018 karakas2010 karakas_lugaro2016 nugrid; do ln -sfn "$HERE/stellar/data/yields/$y" "$ROOT/$y"; done

echo "== stellar tables"
RAD=(--radiation-source "${STELLA_RADIATION_SOURCE:-atlas9}")
[ "${STELLA_RADIATION_SOURCE:-atlas9}" = atlas9 ] && RAD+=(--atlas9-lookup "${STELLA_ATLAS9_LOOKUP:-legacy}")
echo "   radiation: ${RAD[*]}"
python3 "$HERE/stellar/build_stellar_tables.py" "$ROOT" "$OUT/stellar_tables_unified.hdf5" "${RAD[@]}" > "$OUT/build_stellar.log" 2>&1
# (the builder writes stellar_wind_cumulative.hdf5 next to the main output)

echo "== chemistry tables"
python3 "$HERE/chemistry/make_energy_table.py" --eh2 "$HERE/chemistry/tables/eh2-table.txt" --mode x --dx 0.05 \
        --out "$OUT/chemistry/energy-table-x.txt" > "$OUT/build_chemistry.log" 2>&1
for f in eh2-table.txt H2-cooling-ratios.dat wss_z_collis.txt; do cp "$HERE/chemistry/tables/$f" "$OUT/chemistry/"; done

echo "== tests"
python3 "$HERE/tests/run_all.py" "$OUT" ${GIZMO_SRC:+--gizmo "$GIZMO_SRC"}
