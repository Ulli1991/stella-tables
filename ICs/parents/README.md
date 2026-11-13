**For the time being, please notify Ulrich Steinwandel (GitHub: Ulli1991) before using these initial conditions, for example by opening an issue in this repository.** They are in development. Licence of the data: CC BY 4.0 (see `LICENSE-DATA`); the scripts here are BSD 3-Clause (`LICENSE`).

# parents: the parent boxes of the zooms

DM-only boxes with MUSIC, Gadget binary, 64-bit IDs, `gadget_longids = yes`, unique 0..N-1 (the ID is the lattice index ix*N^2 + iy*N + iz).

| file | size | box | particles | particle mass | start | seed |
|---|---|---|---|---|---|---|
| `boxA2_z6_parent_ic.dat` | 537 MB | 1.509 Mpc/h | 256^3 = 16,777,216 | 1.762e4 Msun/h | z = 50 | seed[8] = 95064 |
| `z0_parent_ic.dat` | 4.3 GB (3 parts) | 12.072 Mpc/h | 512^3 = 134,217,728 | 1.128e6 Msun/h | z = 99 | seed[9] = 74023 |

Same cosmology as the zooms (see `../zoom/README.md`). Confs: `boxA2_z6_parent.conf`, `z0_parent.conf`.
- z = 6 box: `MUSIC boxA2_z6_parent.conf`; run to z = 6 with `Config_parent_z6.sh` (BOX_PERIODIC, PMGRID 256) and `boxA2_z6_parent.params`; the zoom halos were found in snapshot 048 (a = 1/7) with Rockstar. Both zooms in `../zoom/` come from this box.
- z = 0 box (for field-dwarf zooms): `MUSIC z0_parent.conf`, then `Config_parent_z0.sh` (PMGRID 512) with `z0_parent.params` and `z0_output_times.txt` (101 outputs); the run to z = 0 took about 14 h on one 112-core node (`MaxMemSize` 3800). **The z = 0 zoom ICs do not exist yet** (halo selection pending); `z0_zoom_template.conf` is the MUSIC template that will be used (multibox region, levelmax 13 to 15).
Run: build the DM-only config with the tested source (stella_dev), `mpirun -np 112 GIZMO_parent PARAMS`, one node (480 GB), `--mem=480G`.

Not parsed or run for this release (these are inputs for making zooms; the runs were made earlier with exactly these files).
