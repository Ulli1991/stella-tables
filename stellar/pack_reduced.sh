#!/bin/bash
# Pack a reduced catalogue tree (reduce_tracks.py output) into one xz tarball per catalogue subdirectory,
# plus one per catalogue for its loose files, so every file stays well below GitHub's size limits.
#   pack_reduced.sh REDUCED_ROOT OUT_DIR
set -euo pipefail
R="$(cd "$1" && pwd)"; mkdir -p "$2"; O="$(cd "$2" && pwd)"
for c in mist parsec2_vms parsec1.2 boost; do
    loose=()
    for e in "$R/$c"/*; do
        n="$(basename "$e")"
        if [ -d "$e" ]; then tar -C "$R" --sort=name --owner=0 --group=0 --numeric-owner --mtime=@0 -cf - "$c/$n" | xz -9 -T1 > "$O/${c}__${n}.tar.xz"
        else loose+=("$c/$n"); fi
    done
    [ ${#loose[@]} -eq 0 ] || tar -C "$R" --sort=name --owner=0 --group=0 --numeric-owner --mtime=@0 -cf - "${loose[@]}" | xz -9 -T1 > "$O/${c}__files.tar.xz"
done
