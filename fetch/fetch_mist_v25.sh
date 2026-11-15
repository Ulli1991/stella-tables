#!/bin/bash
# Download and unpack the MIST v2.5 non-rotating EEP tracks (Dotter+2026, https://mist.science/model_grids.html)
# into RAW_ROOT/mist_v2.5, and check them against fetch/mist_v2.5_sha256.txt where listed.
#   fetch_mist_v25.sh RAW_ROOT [AFE ...]      AFE in m2 p0 p2 p4 p6 ([a/Fe] = -0.2 ... +0.6); default p0
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
BASE="https://mist.science/data/tarballs_v2.5/eeps"
mkdir -p "$1/mist_v2.5"; cd "$1/mist_v2.5"
shift; AFES=("${@:-p0}")
for afe in "${AFES[@]}"; do
    for feh in m400 m350 m300 m275 m250 m225 m200 m175 m150 m125 m100 m075 m050 m025 p000 p025 p050; do
        f="MIST_v2.5_feh_${feh}_afe_${afe}_vvcrit0.0_EEPS.txz"
        [ -f "$f" ] || wget -q "$BASE/$f" || { echo "not available: $f"; rm -f "$f"; continue; }
        want=$(awk -v f="$f" '$2 == f {print $1}' "$HERE/mist_v2.5_sha256.txt" 2>/dev/null || true)
        [ -z "$want" ] || [ "$(sha256sum "$f" | cut -d' ' -f1)" = "$want" ] || { echo "CHECKSUM MISMATCH $f"; exit 1; }
        [ -d "${f%.txz}" ] || tar -xJf "$f"
    done
done
