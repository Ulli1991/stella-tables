**For the time being, please notify Ulrich Steinwandel (GitHub: Ulli1991) before using these initial conditions, for example by opening an issue in this repository.** They are in development. Licence of the data: CC BY 4.0 (see `LICENSE-DATA`); the scripts here are BSD 3-Clause (`LICENSE`).

# wlm: WLM-like isolated dwarf (Smith 2026 m10 specification)

Status: the **10 Msun** IC is here and its structure is checked. The **1 Msun** IC (the intended production resolution, about 1e8 particles, needs 3 or more nodes) **is to come**. The parameter file is a first draft for this IC (not yet run past start-up).

| | |
|---|---|
| file | `wlm/wlm_m10_mg10.dat` (281,139,288 bytes), Gadget binary format 1 (`ICFormat 1`), single file |
| particles | gas (type 0) 6,827,620 at 10.0 Msun; halo (type 1, live dark matter) 1,262,310 at 9,749 Msun; old stellar disc (type 2, gravity only) 975,374 at 10.0 Msun; total 9,065,304 |
| IDs | 32-bit unique 1..N (Gadget-1 default of the generator); `stella_dev` reads them by the record length (file_io/read_ic.cc) |
| metals | none in the file; set from `InitMetallicity` / `InitialMetallicity` = 0.15 Zsun in `wlm_m10.params` (Lee, Skillman & Venn 2005, 12+log O/H = 7.83) |
| model | Hernquist halo matched to NFW with M200 = 1.0e10 Msun, c = 15, V200 = 30.8 km/s (R200 = 45.4 kpc); exponential discs: stellar mass 9.75e6, gas 6.825e7 Msun, scale length 1.1 kpc, gas hydrostatic at 1e4 K (isothermal), sech^2 height 0.86 kpc, lambda = 0.0601, no bulge, no hot halo, no black hole |
| units | kpc, 1e10 Msun, km/s; H0 = 67.74, z = 0 |
| known difference | the halo particles are 9,749 Msun (the model's own matching to 10 Msun gas), not the 200 Msun of the 1 Msun target; to be revisited for the 1 Msun IC |

Generator: Uli's modified MakeDiskGalaxy (https://github.com/Ulli1991/DiskModel, commit 8b98575, after Springel & White 1999 and Springel, Di Matteo & Hernquist 2005; build flags `-DREDSHIFT_SCALING -DV_SCALING -DMY_STRUCTURE`, `make`), input file `m10_mg10.txt` here:
```
./MakeDiskGalaxy m10_mg10.txt        # writes m10_mg10.dat, serial
python3 check_ic.py m10_mg10.dat     # masses, IDs, scale length, Toomre Q, rotation curve
```
The generator sources are GPL-licensed and are not part of this repository.

Config: `Config_stella_wlm.sh` = the Stella default configuration (`Config_stella_default.sh`, also here) with no periodic box and no PM grid; disc stars are gravity-only (type 2). Parameters: `wlm_m10.params` (isolated, end time 500 Myr, snapshots every 5 Myr, softenings gas 0.1 pc minimum, halo 50 pc, disc stars 1 pc, G0 floor 0, solar-pattern metals at 0.15 Zsun).

## Run
Build as in `../README.md` with `Config_stella_wlm.sh`. Memory estimate (scaled from an older run, about 13 kB per particle): 9e6 particles need about 120 GB, so one node:
```
mpirun -np 112 ../GIZMO_wlm wlm_m10.params       # MaxMemSize 4000 MB x 112 ranks; add the restart flag 1 to continue
```
Output: `output/snapshot_NNN.hdf5` every 5 Myr, plus the usual text logs. Expected speed: not measured for this IC. Tables: pin a4830ab.

## Check done
Header, particle counts, masses and the ID block size read directly from the file (2026-11-13). Parameter file parsed by the `Config_stella_default.sh` build (stella_dev df1f23bfb): no ignored or missing tags; the start stopped at the memory allocation on purpose (9e6 particles do not fit a login node).
