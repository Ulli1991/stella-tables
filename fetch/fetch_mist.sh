#!/bin/bash
# Download and unpack the MIST v1.2 non-rotating EEP tracks into RAW_ROOT/mist (URL as used for the original build).
#   fetch_mist.sh RAW_ROOT
set -euo pipefail
BASE="http://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2"
mkdir -p "$1/mist"; cd "$1/mist"
for feh in m4.00 m3.50 m3.00 m2.50 m2.00 m1.75 m1.50 m1.25 m1.00 m0.75 m0.50 m0.25 p0.00 p0.25 p0.50; do
    f="MIST_v1.2_feh_${feh}_afe_p0.0_vvcrit0.0_EEPS.txz"
    [ -f "$f" ] || wget -q "$BASE/$f"
    [ -d "${f%.txz}" ] || tar -xJf "$f"
done
