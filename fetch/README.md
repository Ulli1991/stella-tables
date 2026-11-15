# Raw catalogues

`build.sh` does not need anything from here: it builds from the reduced catalogues in
`stellar/data/reduced/*.tar.xz`. The raw downloads are only needed to re-derive those
(`./build.sh --raw RAW_ROOT`). `reduce_tracks.py` checks the sha256 of every file it reads against
`stellar/data/reduced/MANIFEST.txt`, so a raw tree that differs from the one the tables were built from is
caught.

Expected layout of `RAW_ROOT`:

| directory | content | how to get it |
|---|---|---|
| `mist_v2.5/feh_*_afe_<a>_vvcrit0.0/eeps/` | MIST v2.5 EEP tracks (Dotter et al. 2026, arXiv:2602.22012), non-rotating, 17 [Fe/H] from -4 to +0.5, [a/Fe] = -0.2 ... +0.6 (`<a>` = m2 p0 p2 p4 p6); the build uses `p0` (`--mist-alpha-fe`), tracks up to 3 Msun; tarball checksums in `fetch/mist_v2.5_sha256.txt` | `fetch_mist_v25.sh RAW_ROOT [p0 p4 ...]` (https://mist.science/model_grids.html) |
| `mist/MIST_v1.2_feh_*_afe_p0.0_vvcrit0.0_EEPS/` | MIST v1.2 EEP tracks, non-rotating, 15 [Fe/H] (only for `--mist v1.2`; not in the reduced catalogue) | `fetch_mist.sh RAW_ROOT` |
| `parsec2_vms/photons/Z*_photons/` | PARSEC v2.0 photon rates per track (Q_H, Q_HeI, Q_HeII, Q_OII, Q_Werner vs age), used by `tests/stellar/test_radiation.py` (Z = 0.001 and 0.014, 20/40/100 Msun in the reduced catalogue) | `all_photons.zip` (sha256 df31fd6b9fe49abb776ba1701fcc1aef19a407c3ab95768c58537f11b4175f53) from https://stev.oapd.inaf.it/PARSEC/Database/PARSECv2.0_VMS/, unzip, then unzip each `Z*_photons.zip` into its own directory. The server sends no intermediate certificate: verify with the ZeroSSL RSA DV SSL CA 2 intermediate from its AIA URL (http://crt.sectigo.com/ZeroSSLRSADVSSLCA2.crt) added to the system CAs (`curl --cacert`), not with verification off |
| `parsec2_vms/tracks_Z*/`, `parsec2_vms/ejecta_Z*/` | PARSEC v2.0 tracks (`.TAB`) and ejecta, 13 Z (Costa et al. 2025, arXiv:2501.12917) | "Get" links per Z at https://stev.oapd.inaf.it/PARSEC/tracks_v2_VMS.html |
| `parsec1.2/Z*Y*/` | PARSEC v1.2S tracks (Bressan et al. 2012) | "Get" per Z at https://stev.oapd.inaf.it/PARSEC/tracks_v12s.html |
| `hotstar/` (svo, powr, gotberg18, pollux) | hot-star atmospheres for `stellar/build_hotstar_lut.py`: TLUSTY OSTAR2002/BSTAR2006 (SVO theory server), PoWR OB Vd3 log Q tables and WR SEDs, Gotberg+2018 SEDs, CMFGEN O-star and very-massive-star SEDs from POLLUX (Martins & Palacios 2017-2022, Martins+2025, with the Martins+2025 arXiv source for its parameter tables); checksums in `stellar/data/lut/hotstar_MANIFEST.txt` | `fetch_hotstar.py RAW_ROOT` |
| `boost/BoOST-stellarmodels-*/` | BoOST stellar models (Szécsi et al., arXiv:2004.08203) | download location not confirmed; the build reads the non-rotating tracks up to 25 Msun (wind composition below 14 Msun) and `fdat-columns.txt` |

Files listed with a sha256 in the manifest are read by the builder or the tests; `placeholder` entries are files it
finds by name only, which `reduce_tracks.py` writes as empty files.
