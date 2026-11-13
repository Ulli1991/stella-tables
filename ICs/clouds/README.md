**For the time being, please notify Ulrich Steinwandel (GitHub: Ulli1991) before using these initial conditions, for example by opening an issue in this repository.** They are in development. Licence of the data: CC BY 4.0 (see `LICENSE-DATA`); the scripts here are BSD 3-Clause (`LICENSE`).

# clouds: turbulent molecular cloud ICs (MakeCloud)

All files: HDF5, gas only, one file, equal-mass cells (a cloud sphere inside a periodic box filled with diffuse gas at 1/1000 of the cloud density), uniform magnetic field direction +z (set in the parameter file, `MHD_B_SET_IN_PARAMS`: the file value is ignored), IDs uint32 unique. Metals: none in the file; `InitMetallicity 1.0` (solar) from the parameter file. `Flag_Metals` = 0.

| file | cells | mass/cell | box | cloud | variant |
|---|---|---|---|---|---|
| `gmc_Z1_m1.0.hdf5` | 12,377 (10,000 cloud + 2,377 box) | 1.0 Msun | 70 pc | 1e4 Msun, R = 7.01 pc, n_H = 200 cm^-3, alpha_vir = 2.0, sigma_3D = 2.71 km/s, t_ff = 3.1 Myr, no field in the file | molecular start: x_H2 = 0.4999, x_CO = 1.44e-4, cloud 10 K, box gas u = 77.84 (km/s)^2 |
| `starforge_M2e4_m1.0_mol10K.hdf5` | 24,754 (20,000 + 4,754) | 1.0 Msun | 100 pc | STARFORGE reference cloud M2e4 (Grudic+22 fiducial; Guszejnov+22 Table 1): 2e4 Msun, R = 10 pc, alpha_vir = 2, B = 2 uG | molecular start as above |
| `starforge_M2e4_m1.0_atomic20K.hdf5` | 24,754 | 1.0 Msun | 100 pc | same | paper-literal: all H atomic, all gas 20 K (u = 0.195 (km/s)^2) |
| `starforge_M2e4_m0.1_mol10K.hdf5` | 247,550 (200,000 + 47,550) | 0.1 Msun | 100 pc | same | molecular |
| `starforge_M2e4_m0.1_atomic20K.hdf5` | 247,550 | 0.1 Msun | 100 pc | same | atomic |

Geometry, velocities, field and masses are identical between the variants of one resolution; only `InternalEnergy` and `ChemcoolTracAbund` (x_H2, x_H+, x_CO) differ. All mass checks: sum of cell masses = cloud + box mass; verification numbers are in the `*_verify.json` files (virial parameter, free-fall time, net momentum ~1e-12, field energies). `STARFORGE_M2e4_details.md` lists each parameter against the paper values, with the places where the two papers disagree (20 K vs 10 K start, mu) and what was chosen.

## Generators (BSD 3-Clause wrappers around Grudic's MakeCloud, https://github.com/mikegrudic/MakeCloud, commit cd7b4f3)
Needs python 3.10+ (MPCDF: `module load python-waterboa/2025.06`), numpy, h5py, docopt, `MAKECLOUD=<MakeCloud checkout>`.
- `gmc_Z1_m1.0.hdf5`: `python3 make_gmc_Z1_m1.0.py --mgas 1.0 --out gmc_Z1_m1.0.hdf5` (defaults M = 1e4, R = 7, alpha = 2, Z = 1, seed 42). Verification numbers: `gmc_Z1_m1.0_verify.json`. Not regenerated for this release (needs MakeCloud); the file is identical in 176 copies used in earlier tests.
- M2e4:
```
python3 $MAKECLOUD/MakeCloud.py --M=2e4 --R=10 --N=20000 --alpha_turb=2 --turb_slope=2 --turb_sol=0.5 --minmode=2 --turb_seed=42 --bfixed=2e-6 --Z=1 --ISRF=1.7 --turb_path=$PWD/turb --glass_path=$HOME/.makecloud_glass/glass_256.hdf5   # N=200000 for 0.1 Msun
python3 build_ic.py raw/<MakeCloud output>.hdf5 M2e4_m1.0_mol10K    mol10K    pairwise
python3 build_ic.py raw/<MakeCloud output>.hdf5 M2e4_m1.0_atomic20K atomic20K pairwise
```
(`build_ic.py` converts MakeCloud's pc / Msun / Gauss to kpc / 1e10 Msun / km/s / Gauss; the released files carry the prefix `starforge_`.)

## Parameter files and Config
`Config_stella_cloud.sh` (Stella default plus periodic box, isolated cloud gravity, kept IC abundances). 
- `gmc_Z1_m1.0.params`: for the 12k-cell cloud, B_z = 6.61e-6 G (mu about 3.3), end time 4.9 Myr.
- `starforge_M2e4.params`: the same with box 0.1 kpc, B_z = 2 uG, end time 9.3 Myr (Grudic+22), `InitCondFile` set to the mol10K 1 Msun file (edit the name for the other files). Written for this release from the tested cloud parameters; the physics values of Stella for the STARFORGE comparison (resolution caveats: 1 Msun cells are far coarser than the 1e-3 Msun of the paper) are still to be fixed. G0 = 1.7 Habing, cosmic-ray rate 3e-17 per s.

## Run (login node or any node)
Build as in `../README.md` with `Config_stella_cloud.sh`, then in the run directory
```
mpirun -np 4 ../GIZMO_cloud gmc_Z1_m1.0.params       # or starforge_M2e4.params
```
`MaxMemSize 1800` MB per rank. Output `output/snapshot_000.hdf5`, ... every 0.25 Myr plus text logs (`SFinfo.txt`, `IMFinfo.txt`, `ENERGYinfo.txt`, `FeedbackBudget.txt`, `mass_balance.txt`, `cpu.txt`).

## Checks done (2026-11-13, login node, 4 cores, stella_dev df1f23bfb, `make clean` + `make`, tables pin a4830ab)
| run | end time | wall time | result |
|---|---|---|---|
| `gmc_Z1_m1.0` | 0.02 Myr (3 sync points) | seconds | exit 0, no ignored-parameter warnings |
| `starforge_M2e4_m1.0_mol10K` | 0.2 Myr (18 sync points) | 28 s | exit 0, no warnings |
| `starforge_M2e4_m1.0_atomic20K` | 0.1 Myr | 28 s | exit 0, no warnings |
| `starforge_M2e4_m0.1_mol10K` | 0.03 Myr (7 sync points) | 69 s, 1.1 GB peak | exit 0, no warnings |
Expect the cost to rise steeply once stars form (the full 5 to 10 Myr is a node-hours job at 1 Msun, and more at 0.1 Msun; not measured here). The `HDF5-DIAG` lines on stderr at start are GIZMO probing optional datasets.
