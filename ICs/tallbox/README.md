**For the time being, please notify Ulrich Steinwandel (GitHub: Ulli1991) before using these initial conditions, for example by opening an issue in this repository.** They are in development. Licence of the data: CC BY 4.0 (see `LICENSE-DATA`); the scripts here are BSD 3-Clause (`LICENSE`).

# tallbox: stratified box, 1 x 1 x 20 kpc

Solar-neighbourhood-like gas slab after Hu, Sternberg & van Dishoeck 2021 (the box of the Stella tallbox showcase).

| | |
|---|---|
| file | `tallbox/tallbox.hdf5` (1,560,006,680 bytes; hash in `../manifest.txt`) |
| format | HDF5, one file, gas only (PartType0) |
| cells | 10,000,000 gas cells, every cell 1.0 Msun (total 1e7 Msun) |
| box | x, y in [0, 1] kpc (periodic); z in [8.89, 11.08] kpc inside a box of 20 kpc (needs `BOX_LONG_Z=20`, mid-plane at z = 10 kpc, outflow in z) |
| gas | Sigma_gas = 10 Msun/pc^2, sech^2 profile with scale height H_g = 250 pc, uniform in x, y, T = 1e4 K (u = 206.4 (km/s)^2), no velocities beyond the generator's random seed 42 field, no magnetic field |
| IDs | uint32, unique (10 .. 10,000,009) |
| metals | `Flag_Metals` = 27, `ElementAbundance` (27 mass fractions, scaled-solar at 0.1 Zsun: H 0.7236, He 0.275, C 2.5e-4, ...). On a fresh start GIZMO sets the metals from `InitMetallicity` (1.0, solar, in `tallbox.params`), so the IC value is not used |
| external potential | not in the IC: applied by the code (`GALSF_RESOLVEDISM_TALLBOX`: Hu+21 stellar disc Sigma_* = 40 Msun/pc^2, H_* = 0.25 kpc, plus vertical NFW) |

Generator (BSD 3-Clause), defaults reproduce the file **bit for bit** (all six datasets and the header checked, 2026-11-13):
```
python3 make_tallbox.py --out tallbox.hdf5      # = --sigma 10 --Hg 0.25 --Lxy 1.0 --zmax 10.0 --mgas 1.0 --T 1e4 --Zsun 0.1 --seed 42; needs numpy, h5py; 3 s
```
Files: `Config_stella_tallbox.sh` (Stella default plus periodic x/y, outflow z, tallbox potential, GHOSDT-style G0 and CR rate), `tallbox.params`.

## Run
Build and set up as in `../README.md` (Config `Config_stella_tallbox.sh`). Then, on one node:
```
mpirun -np 112 ../GIZMO_tallbox tallbox.params          # 1e7 cells: 112 ranks x MaxMemSize 4000 MB (480 GB node)
```
Output goes to `output/`: `snapshot_000.hdf5` ... every 4.9 Myr (about 100 over the run to 489 Myr), plus `cpu.txt`, `SFinfo.txt`, `ENERGYinfo.txt` and other text logs. The run is long (days): it was run in 24-hour segments with restarts (`CpuTimeBetRestartFile`, restart flag 1). A 1e7-cell box does not fit on a login node.
Tables: stella-tables pin a4830ab.

## Check done
Parameter file parsed by a `Config_stella_tallbox.sh` build (stella_dev df1f23bfb): no ignored tags, no missing tags. The run itself was not repeated for this release (needs a node).
