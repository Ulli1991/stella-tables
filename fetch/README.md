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
| `parsec2_vms/tracks_Z*/`, `parsec2_vms/ejecta_Z*/` | PARSEC v2.0 tracks (`.TAB`) and ejecta, 13 Z (Costa et al. 2025, arXiv:2501.12917) | "Get" links per Z at https://stev.oapd.inaf.it/PARSEC/tracks_v2_VMS.html |
| `parsec1.2/Z*Y*/` | PARSEC v1.2S tracks (Bressan et al. 2012) | "Get" per Z at https://stev.oapd.inaf.it/PARSEC/tracks_v12s.html |
| `hotstar/` (svo, powr, gotberg18) | hot-star atmospheres for `stellar/build_hotstar_lut.py`: TLUSTY OSTAR2002/BSTAR2006 (SVO theory server), PoWR OB Vd3 log Q tables and WR SEDs, Gotberg+2018 SEDs; checksums in `stellar/data/lut/hotstar_MANIFEST.txt` | `fetch_hotstar.py RAW_ROOT` |
| `boost/BoOST-stellarmodels-*/` | BoOST stellar models (Szécsi et al., arXiv:2004.08203) | download location not confirmed; the current build uses no BoOST track (file names only) |

The builder reads 1816 files out of about 24 GB of downloads; the other 6660 entries in the manifest
are files it finds by name only. `reduce_tracks.py` writes those as empty files.
