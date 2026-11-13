**For the time being, please notify Ulrich Steinwandel (GitHub: Ulli1991) before using these initial conditions, for example by opening an issue in this repository.** They are in development. Licence of the data: CC BY 4.0 (see `LICENSE-DATA`); the scripts here are BSD 3-Clause (`LICENSE`).

# zoom: cosmological zoom-in ICs of two ~1e8 Msun halos at z = 6

Two halos from one 1.509 Mpc/h parent box (256^3 dark matter particles, see `../parents/`), re-simulated with MUSIC at higher resolution. Cosmology: Omega_m 0.31, Omega_Lambda 0.69, Omega_b 0.049, h 0.68, sigma_8 0.81, n_s 0.965, Eisenstein-Hu transfer function; start z = 50, second-order Lagrangian perturbations; comoving box 1509 ckpc/h. All files: Gadget binary format 1 (`ICFormat 1`), **64-bit particle IDs** (MUSIC `gadget_longids = yes`), unique 0..N-1; masses in 1e10 Msun/h (physical Msun = value x 1e10 / 0.68); coarse DM is type 2 (per-particle masses).

| halo | Rockstar M_vir at z = 6 | position (unit box) | picked because |
|---|---|---|---|
| `1e8_c` | 7.2e7 Msun | (0.319, 0.077, 0.760) | isolated, no comparable neighbour |
| `1e8_a` | 1.10e8 Msun | (0.770, 0.722, 0.671) | isolated, no comparable neighbour |
(match table of the parent box against an earlier one: `music/halo_match_boxA2.txt`.)

| file | size | gas | high-res DM (type 1) | coarse DM (type 2) | gas mass | high-res DM mass | total N |
|---|---|---|---|---|---|---|---|
| `ic_1e8_c_L10_gas.dat` | 635 MB | 444,032 | 444,032 | 16,793,329 | 64.0 Msun | 341 Msun | 17,681,393 |
| `ic_1e8_c_L12_gas.dat` | 2.6 GB (2 parts) | 28,418,048 | 28,418,048 | 18,603,753 | 1.0 Msun | 5.3 Msun | 75,439,849 |
| `ic_1e8_a_L12_gas.dat` | 3.5 GB (2 parts) | 40,882,176 | 40,882,176 | 19,141,474 | 1.0 Msun | 5.3 Msun | 100,905,826 |
| `ic_1e8_c_L12_dm.dat` | 1.6 GB | 0 | 28,418,048 | 18,603,753 | n/a | 6.3 Msun | 47,021,801 |
| `ic_1e8_a_L12_dm.dat` | 2.0 GB | 0 | 40,882,176 | 19,141,474 | n/a | 6.3 Msun | 60,023,650 |

Level 10 = 64 Msun gas, level 12 = 1 Msun gas. **Missing**: `1e8_a` at level 10 and the DM-only level 10 versions were not generated; make them with the MUSIC step below (`levelmax = 10`).
Metals: the gas ICs carry none; the parameter files start pristine (`InitMetallicity 0`, network 17 primordial chemistry, Pop III channel on). Gas temperature: MUSIC value (about 43 K), `InitGasTemp 0`.

## How they were made (MUSIC 2.0a, git 967651b)
Run MUSIC in a directory with `region_1e8_c.txt` (or `_a`) next to the conf; it writes white-noise files (`wnoise_00NN.bin`, deterministic from the seeds in the conf) and the IC:
```
MUSIC music/music_1e8_c_L12_gas.conf      # about an hour with 32 threads (a 200 GB, 2 h job was used); levelmax = 10 for the L10 file
```
Confs: `music/music_1e8_{c_L10,c_L12,a_L12}_gas.conf` and `music/music_1e8_{c,a}_L12_dm.conf`; region point files (`music/region_1e8_{a,c}.txt`, 4514 and 6698 points at region level 9, unit-box coordinates of the Lagrangian region of the halo). The ID rule: **never convert the IDs to 32 bit** (a 32-bit rewrite makes the snapshot IDs non-unique); check after the first snapshot that the IDs are unique and the smallest is below 2^32.

## Parameter files and Config
- `Config_stella_zoom.sh` (Stella: net 17, periodic box with PMGRID 256 and `PM_PLACEHIGHRESREGION=1+2` (gas + high-res DM on the fine PM grid), Pop III stars, Remy-Ruyer DGR, primordial helium, CR rate scaled with the local G0). `merge/split` is on (all Stella runs).
- `zoom_L10.params` for the level 10 IC, `zoom_L12.params` for level 12 (softenings scale with the spacing: `SofteningGas` 0.004 ckpc/h at level 10, 0.001 at level 12; `MaxMemSize` 4000 / 4300 MB per rank, `BufferSize` 100 / 400). The `1e8_a` level 12 IC uses `zoom_L12.params` with `InitCondFile ic_1e8_a_L12_gas.dat`. They are derived from the parameter files of the zoom production runs (level 10: the latest full run; level 12: the same with its memory settings); softenings are first guesses, untuned.
- DM only: `Config_zoom_dm.sh` (five lines) and `zoom_dm_1e8_c_L12.params`, `zoom_dm_1e8_a_L12.params`.

## Run
Build as in `../README.md` with `Config_stella_zoom.sh` (zoom branch a6fb413ef). One full node:
```
mpirun -np 112 ../GIZMO_zoom zoom_L10.params          # --mem=480G; 112 ranks x 4000 MB = 448 GB (L12: 4300 MB = 482 GB, ask for the whole node)
```
**Memory**: level 10 (17.7M particles) and level 12 (75M to 100M) both need a whole 480 GB node. Time: from the z = 50 start to z = 6 takes several 24-hour segments; the first 24 h segment of level 10 reached z = 8.5. Runs restart from the restart files (`CpuTimeBetRestartFile`; restart flag 1), hence the chain of 24 h jobs. Output: `snapshot_NNN.hdf5` about every 10% in scale factor (21 up to z = 6), `SFinfo.txt` etc. The run directory needs `ewald_spc_table_64_dbl.dat` (GIZMO writes it on the first start) and `TREECOOL`.
Tables: pin a4830ab.

## Check done
Parameter files parsed with the zoom-config build (zoom branch a6fb413ef): no ignored tags, none missing; stopped at the memory allocation on purpose. ID block sizes (8 bytes per ID) and the particle counts read from the files directly; the unique-ID check of each file was done when it was made (all PASS, ID range 0..N-1).
