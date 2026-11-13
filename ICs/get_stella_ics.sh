#!/bin/bash
# Fetch Stella initial conditions by name and check them against manifest.txt (sha256).
#   get_stella_ics.sh [--list] [--from DIR] [--prefix DIR] NAME...
#   NAME      a set (tallbox, clouds, wlm, zoom, parents), a path from the manifest (e.g. zoom/ic_1e8_c_L10_gas.dat), or "all"
#   --prefix  install directory (default: $STELLA_ICS_PREFIX, else ~/.stella-ics/ics-v1); files land at PREFIX/<set>/<file>
#   --from    take the files from this directory (same layout as PREFIX) instead of downloading
#   $STELLA_ICS_URL  release download URL (default: https://github.com/Ulli1991/stella-tables/releases/download/ics-v1)
# Already installed files are only re-fetched if their hash does not match. Sizes: see --list (total of "all": about 16 GB).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; MAN="$HERE/manifest.txt"
PREFIX="${STELLA_ICS_PREFIX:-$HOME/.stella-ics/ics-v1}"; FROM=""; LIST=0; NAMES=()
URL="${STELLA_ICS_URL:-https://github.com/Ulli1991/stella-tables/releases/download/ics-v1}"
while [ $# -gt 0 ]; do case "$1" in
    --list) LIST=1;; --from) FROM="$2"; shift;; --prefix) PREFIX="$2"; shift;; -*) echo "unknown option $1" >&2; exit 2;; *) NAMES+=("$1");; esac; shift; done
if command -v sha256sum >/dev/null; then sha() { sha256sum "$1" | cut -d' ' -f1; }; else sha() { shasum -a 256 "$1" | cut -d' ' -f1; }; fi
entries() { grep -v '^#' "$MAN" | grep -v '^$'; }
if [ $LIST = 1 ] || [ ${#NAMES[@]} = 0 ]; then
    entries | awk '{printf "%-48s %8.1f MB\n", $3, $2/1e6}'; [ ${#NAMES[@]} = 0 ] && [ $LIST = 0 ] && echo "give a name, a set or 'all'" >&2; exit 0; fi
mkdir -p "$PREFIX"; PREFIX="$(cd "$PREFIX" && pwd)"; found=0
while read -r want size path rest; do
    pick=0; for n in "${NAMES[@]}"; do case "$n" in all) pick=1;; "$path") pick=1;; */) [ "${path#$n}" != "$path" ] && pick=1;; *) [ "${path#$n/}" != "$path" ] && pick=1;; esac; done
    [ $pick = 1 ] || continue; found=1
    parts=1; [ -n "$rest" ] && parts="${rest#parts=}"; dest="$PREFIX/$path"; mkdir -p "$(dirname "$dest")"
    if [ -f "$dest" ] && [ "$(stat -c %s "$dest")" = "$size" ] && [ "$(sha "$dest")" = "$want" ]; then echo "ok (already there) $path"; continue; fi
    tmp="$dest.tmp.$$"; : > "$tmp"
    if [ -z "$rest" ]; then srcs=("$path"); else srcs=(); for ((i=0;i<parts;i++)); do srcs+=("$(printf '%s.part%03d' "$path" $i)"); done; fi
    for s in "${srcs[@]}"; do
        if [ -n "$FROM" ]; then cat "$FROM/$s" >> "$tmp"; else curl -fL --retry 3 -sS "$URL/${s//\//__}" >> "$tmp" || { rm -f "$tmp"; echo "download failed: $s" >&2; exit 1; }; fi
    done
    got="$(sha "$tmp")"
    if [ "$got" != "$want" ]; then rm -f "$tmp"; echo "HASH MISMATCH for $path: got $got, expected $want" >&2; exit 1; fi
    mv "$tmp" "$dest"; echo "installed $path ($size bytes, hash checked)"
done < <(entries)
[ $found = 1 ] || { echo "no manifest entry matches: ${NAMES[*]}" >&2; exit 1; }
echo "initial conditions are in $PREFIX"
