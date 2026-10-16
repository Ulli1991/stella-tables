# Raw catalogues

`build.sh` does not need anything from here: it builds from the reduced catalogues in
`stellar/data/reduced/*.tar.xz`. The raw downloads are only needed to re-derive those
(`./build.sh --raw RAW_ROOT`). `reduce_tracks.py` checks the sha256 of every file it reads against
`stellar/data/reduced/MANIFEST.txt`, so a raw tree that differs from the one the tables were built from is
caught.

Expected layout of `RAW_ROOT`:

| directory | content | how to get it |
|---|---|---|
| `mist/MIST_v1.2_feh_*_afe_p0.0_vvcrit0.0_EEPS/` | MIST v1.2 EEP tracks, non-rotating, 15 [Fe/H] | `fetch_mist.sh RAW_ROOT` |
| `parsec2_vms/tracks_Z*/`, `parsec2_vms/ejecta_Z*/` | PARSEC v2 very-massive-star tracks (`.TAB`) and ejecta | source URL not recorded yet (TODO) |
| `parsec1.2/Z*Y*/` | PARSEC v1.2S tracks (11 Z) | source URL not recorded yet (TODO) |
| `boost/BoOST-stellarmodels-*/` | BoOST stellar models | source URL not recorded yet (TODO); only file names are used |

The builder reads 1170 files (4.8 GB) out of about 24 GB of downloads; the other 7306 entries in the manifest
are files it finds by name only. `reduce_tracks.py` writes those as empty files.
