#!/bin/bash
# Build every table Stella/GIZMO reads, then run the test suite.
#   ./build.sh                 build from the versioned reduced catalogues (stellar/data/reduced/*.tar.xz)
#   ./build.sh --raw RAW_ROOT  re-derive the reduced catalogues from the full downloads first (see fetch/)
#   GIZMO_SRC=/path/to/gizmo ./build.sh   also test the energy table through GIZMO's Fortran and C++ lookups
#   STELLA_RADIATION_SOURCE=hotstar_stack|atlas9|planck|atlas9_nonionizing (default hotstar_stack), STELLA_ATLAS9_LOOKUP=legacy|clamp|masked_planck
#     (default legacy): band luminosities, see README "Radiation"
#   STELLA_MIST_ALPHA_FE=p4|p0|m2|p2|p6: MIST v2.5 [a/Fe] = +0.4 (default), 0, -0.2, +0.2, +0.6
# Outputs land in build/: stellar_tables_unified.hdf5, stellar_wind_cumulative.hdf5, chemistry/energy-table-x.txt
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# bit-reproducible tables: pinned environment (requirements-build.txt), one thread, numpy's baseline CPU code paths
python3 "$HERE/check_build_env.py" || { [ "${STELLA_ALLOW_ANY_ENV:-0}" = 1 ] && echo "   STELLA_ALLOW_ANY_ENV=1: building anyway; tables will differ from the published ones in the last bits"; } || {
    echo "build.sh: Python environment differs from requirements-build.txt (MPCDF: module load python-waterboa/2026.07)"; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -c 'import numpy, sys; sys.exit(int(numpy.__version__.split(".")[0]) < 2)' && export NPY_DISABLE_CPU_FEATURES="${NPY_DISABLE_CPU_FEATURES-X86_V3 X86_V4 AVX512_ICL AVX512_SPR}"
OUT="$HERE/build"; ROOT="$OUT/data_root"
mkdir -p "$OUT/chemistry" "$ROOT"

if [ "${1:-}" = "--raw" ]; then
    python3 "$HERE/stellar/reduce_tracks.py" "$2" "$ROOT"
else
    for t in "$HERE"/stellar/data/reduced/*.tar.xz; do tar -C "$ROOT" -xJf "$t"; done
fi
for y in chieffi_limongi2018 gronow2021 karakas2010 karakas_lugaro2016 limongi2025 nugrid; do ln -sfn "$HERE/stellar/data/yields/$y" "$ROOT/$y"; done

echo "== stellar tables"
RAD=(--radiation-source "${STELLA_RADIATION_SOURCE:-hotstar_stack}")
[ "${STELLA_RADIATION_SOURCE:-hotstar_stack}" = atlas9 ] && RAD+=(--atlas9-lookup "${STELLA_ATLAS9_LOOKUP:-legacy}")
RAD+=(--mist-alpha-fe "${STELLA_MIST_ALPHA_FE:-p4}")
echo "   radiation and MIST: ${RAD[*]}"
python3 "$HERE/stellar/build_stellar_tables.py" "$ROOT" "$OUT/stellar_tables_unified.hdf5" "${RAD[@]}" > "$OUT/build_stellar.log" 2>&1
# (the builder writes stellar_wind_cumulative.hdf5 next to the main output)

echo "== chemistry tables"
# energy table: the canonical file (generator output in the recorded environment) is installed; the generator runs
# every build and the tests compare its output with the canonical table (to 1e-12, exact in the recorded environment)
python3 "$HERE/chemistry/make_energy_table.py" --eh2 "$HERE/chemistry/tables/eh2-table.txt" --mode x --dx 0.05 \
        --out "$OUT/chemistry/energy-table-x.regenerated.txt" > "$OUT/build_chemistry.log" 2>&1
xz -dc "$HERE/chemistry/reference/energy-table-x.txt.xz" > "$OUT/chemistry/energy-table-x.txt"
for f in eh2-table.txt H2-cooling-ratios.dat wss_z_collis.txt; do cp "$HERE/chemistry/tables/$f" "$OUT/chemistry/"; done

echo "== SGChem data tables"
python3 "$HERE/chemistry/sgchem/make_net18_tables.py" "$OUT/chemistry/sgchem_net18_tables.bin" >> "$OUT/build_chemistry.log" 2>&1
python3 "$HERE/chemistry/sgchem/make_chem_tables.py" "$OUT/chemistry/sgchem_chem_tables.bin" --header "$OUT/chemistry/chem_tables_data.h" >> "$OUT/build_chemistry.log" 2>&1

echo "== tests"
python3 "$HERE/tests/run_all.py" "$OUT" ${GIZMO_SRC:+--gizmo "$GIZMO_SRC"}
