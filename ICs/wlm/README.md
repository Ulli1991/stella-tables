**For the time being, please notify Ulrich Steinwandel (GitHub: Ulli1991) before using these initial conditions, for example by opening an issue in this repository.** They are in development. Licence of the data: CC BY 4.0 (see `LICENSE-DATA`); the scripts here are BSD 3-Clause (`LICENSE`).

# wlm: WLM-like isolated dwarf (Smith 2026 m10 specification)

The isolated dwarf `m10` of Smith 2026 ("Imladris", arXiv:2604.00100, Sec. 3.1 and Table 1; the system of Smith et al. 2021, loosely WLM), with a live halo, an old stellar disc and a gas disc. Two resolutions: **1 Msun** gas (production, finer than the paper's finest 2.5 Msun) and 10 Msun (coarse test).

| | 1 Msun (production) | 10 Msun (test) |
|---|---|---|
| file | `wlm/wlm_m10_mg1.hdf5` (5,297,409,012 bytes, 3 parts) | `wlm/wlm_m10_mg10.dat` (281,139,288 bytes) |
| format | HDF5 (`ICFormat 3`), one file, explicit `Masses` | Gadget binary format 1 (`ICFormat 1`) |
| gas (type 0) | 68,276,196 x 1.000 Msun | 6,827,620 x 10.0 Msun |
| halo (type 1, live DM) | 61,533,770 x 200.0 Msun (the paper's m10mg2.5 value) | 1,262,310 x 9,749 Msun |
| old stellar disc (type 2, gravity only) | 9,753,737 x 1.000 Msun | 975,374 x 10.0 Msun |
| total | 139,563,703 | 9,065,304 |
| IDs | uint64, unique 1..N | 32-bit, unique 1..N (read by the record length) |
| parameter file | `wlm_m10_mg1.params` | `wlm_m10.params` |

Why the 10 Msun halo particles are 9,749 Msun: the generator gives each halo particle the mass of the Hernquist halo inside its velocity grid (r < 454 kpc, 1.2307e10 Msun) divided by `N_HALO`; that file used `N_HALO` = 1,262,310 (nominally 1e4 Msun). The 1 Msun IC sets `N_HALO` = 61,533,770 = 1.2307e10 / 200, which gives 200.0 Msun.

Model (both files, same generator input apart from particle numbers): Hernquist halo matched to NFW with M200 = 1.0e10 Msun (baryons included), c = 15, V200 = 30.77 km/s, R200 = 45.43 kpc, Hernquist a = 5.80 kpc, H0 = 67.74; stellar disc 9.75e6 Msun and gas disc 6.83e7 Msun, both exponential with scale length 1.10 kpc (spin lambda = 0.0601 sets it); gas in vertical hydrostatic equilibrium at T = 1e4 K (isothermal, u = 210.5 (km/s)^2, i.e. 1e4 K for ionised gas; Stella starts the gas atomic, so the same pressure reads as about 2e4 K); stellar disc sech^2 with z0 = 0.86 kpc, which has the same half-mass height (0.472 kpc) as the paper's Gaussian with scale height 0.7 kpc read as the dispersion (the sech^2 has heavier tails, rms height 0.78 instead of 0.70 kpc); no bulge, hot halo or black hole. Units kpc, 1e10 Msun, km/s.

Metals: none in the files. On a fresh start GIZMO sets every particle's metals from `InitMetallicity` and the chemistry from `InitialMetallicity` (file_io/read_ic.cc skips IC blocks past the internal energy; core/init.cc writes the 27-element solar pattern scaled by `InitMetallicity`). Both are 0.0907: the paper's 0.1 Zsun with Zsun = 0.0127 is Z = 0.00127, and Stella's solar pattern (Asplund 2009) sums to Z = 0.014. The paper scales between a primordial and the Asplund 2021 pattern; Stella scales the 2009 solar pattern.

Softenings (paper's rule 0.75 (m / 20 Msun)^(1/3) pc for halo and initial stars; gas adaptive): 1 Msun set: halo 1.62 pc, disc stars 0.28 pc; 10 Msun set: halo 5.9 pc, disc stars 0.60 pc. Gas and new stars use Stella's adaptive softening (0.1 pc minimum).

Run length: 1.3 Gyr (`TimeMax` 1.3295). The paper relaxes its IC for 300 Myr with a Jeans pressure floor and a cap on stellar lifetimes (its App. C) and analyses 100 to 1000 Myr after that. Stella has no such relaxation mode: the run starts with the full model, and the first 300 Myr are the settling phase. Snapshots every 10 Myr.

Generator: Uli's modified MakeDiskGalaxy (https://github.com/Ulli1991/DiskModel, commit 0c1aa45 of the `bugfixes` branch, after Springel & White 1999 and Springel, Di Matteo & Hernquist 2005; build flags `-DREDSHIFT_SCALING -DV_SCALING -DMY_STRUCTURE`). Inputs `m10_mg1.txt`, `m10_mg10.txt` here:
```
./MakeDiskGalaxy m10_mg1.txt                              # m10_mg1.dat (4.2 GB), serial, 9 min, 8.7 GB memory
python3 check_ic.py m10_mg1.dat m10_mg1.vc                # masses, IDs, scale lengths, heights, Toomre Q, rotation curve (13 GB)
python3 gadget1_to_hdf5.py m10_mg1.dat wlm_m10_mg1.hdf5   # uint64 IDs, explicit masses (4 GB)
```
The generator sources are GPL-licensed and are not part of this repository.

Config: `Config_stella_wlm.sh` = the Stella default configuration (`Config_stella_default.sh`, also here) unchanged: no periodic box, no PM grid; merge/split on (the default, no `PREVENT_PARTICLE_MERGE_SPLIT`); the type-2 disc stars are not radiation or feedback sources (Stella's sources are types 4 and 5).

## Structure checks (1 Msun file; the 10 Msun file agrees to the digits shown)

| R [kpc] | 0.5 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| v_c (spherical M(<r)) [km/s] | 26.9 | 35.3 | 43.3 | 46.8 | 48.3 | 48.9 |
| gas v_rot [km/s] | 24.1 | 31.8 | 38.8 | 45.9 | 47.8 | 48.6 |
| Sigma_gas [Msun/pc^2] | 5.55 | 3.59 | 1.46 | 0.59 | 0.24 | 0.10 |
| Toomre Q gas / stars / combined | 13.8 / 182 / 12.8 | 14.3 / 160 / 13.1 | 20.5 / 167 / 18.2 | 35 / 220 / 30 | 65 / 325 / 54 | 128 / 523 / 103 |

Masses: gas 6.8276e7, disc 9.7537e6 Msun, M(<R200) = 1.000e10 Msun (DM 9.924e9). Exponential scale lengths (fit 0.5 to 5 kpc): gas 1.103, stars 1.102 kpc. Median |z|: stars 0.472 kpc, gas 0.270 kpc (R < 3 kpc). Centre of mass within 9 pc and 0.002 km/s of the origin. Very stable at 1e4 K (Q > 12); unstable once the gas cools, which is what the settling phase is for.

## Run
Build and set up as in `../README.md` with `Config_stella_wlm.sh`; link `wlm_m10_mg1.hdf5` into the run directory.

Memory (1 Msun): this build stores 896 bytes per particle plus 3080 per gas cell; with `PartAllocFactor` 2.5 that is 0.84 TB, about 0.95 TB with tree and buffers, growing with the stars formed. It needs **at least 3 nodes** of 480 GB (on 2 nodes a rank would need more than the 4.28 GB per core); 4 nodes leave room for imbalance. On Orion: `p.exclusive` (multi-node; `p.shared` is one node only):
```
srun -p p.exclusive --nodes=4 --ntasks-per-node=112 ./GIZMO_wlm wlm_m10_mg1.params     # MaxMemSize 4000 MB x 112 ranks per node; add 1 to restart
```
`TimeLimitCPU` (85000 s) and `CpuTimeBetRestartFile` (2 h) are set for 24-hour segments with restarts. 10 Msun: about 9e6 particles, one node, `mpirun -np 112 ./GIZMO_wlm wlm_m10.params`.

Output: `output/snapshot_NNN.hdf5` every 10 Myr (131 snapshots), plus the usual text logs. Tables: pin a4830ab.

## Check done
Built `Config_stella_wlm.sh` at stella_dev df1f23bfb (`make clean` first) and started it with `wlm_m10_mg1.params` on 3 login cores, using a 100 Msun version of the same model (`m10_mg100.txt` in the build area: 6.8e5 gas, 6.2e5 DM at 2e4 Msun, 1e5 disc stars; only `InitCondFile`, `MaxMemSize`, `TimeLimitCPU` and the snapshot interval changed), 2026-11-13: no ignored or missing tags; `InitMetallicity`/`InitialMetallicity` 0.0907 read; IDs and masses read from the hdf5; about 100 steps to 4.8 Myr (step 40 to 80 kyr); no NaN or Inf in any snapshot field; disc unchanged (gas half-mass radius 1.836 kpc, gas median |z| 269 pc, disc median |z| 471 to 473 pc, centre-of-mass drift < 0.2 pc); gas cooling from 2e4 K; no stars yet; 2.6 to 2.7 GB per rank (4.6e5 particles per rank). The snapshot at t = 0 shows `Temperature` 100 K: the chemistry temperature is a placeholder until the first cooling step (core/init.cc); the internal energy is the IC value. The 1 Msun IC itself does not fit a login node and was not started.
