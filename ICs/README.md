**For the time being, please notify Ulrich Steinwandel (GitHub: Ulli1991) before using these initial conditions, for example by opening an issue in this repository.** They are in development.

# Stella initial conditions (ICs)

Standard initial conditions for the Stella model in GIZMO, each shipped as a ready-to-run set: the IC file, the GIZMO parameter file, the `Config` file it was run with, the generator and its exact command, and a short "how to run". The IC files themselves are large and are **not in git**: they are release assets, fetched and hash-checked with `get_stella_ics.sh`.

**Licence.** The data (tables and ICs, including the files this script downloads) are licensed CC BY 4.0, see `LICENSE-DATA` in the repository root. The scripts in this repository (generators, `get_stella_ics.sh`) are BSD 3-Clause, see `LICENSE`. Please credit the data as:
"Stella initial conditions, U. Steinwandel (2026), https://github.com/Ulli1991/stella-tables (ICs/), CC BY 4.0", and cite the papers named in the README of the set you use (generator and observational sources; the Stella paper once it is out).

## What is here

| set | what | IC files (size) | needs |
|---|---|---|---|
| [`tallbox/`](tallbox/README.md) | stratified 1 x 1 x 20 kpc box, 1e7 gas cells of 1 Msun, Sigma_gas = 10 Msun/pc^2 | `tallbox.hdf5` (1.6 GB) | 112 cores, 480 GB |
| [`clouds/`](clouds/README.md) | small test cloud (12k cells, 1 Msun) and the reference cloud M2e4 at 1 and 0.1 Msun (STARFORGE, Grudic+22), 2 chemistry variants each | 5 files, 1 MB to 29 MB | laptop to login node |
| [`wlm/`](wlm/README.md) | WLM-like dwarf (Smith 2026 m10 spec), 10 Msun gas, live halo and old stellar disc. The 1 Msun version is **to come** | `wlm_m10_mg10.dat` (0.3 GB) | 112 cores, 480 GB |
| [`zoom/`](zoom/README.md) | cosmological zooms of two ~1e8 Msun halos at z = 6: gas + DM at level 10 (64 Msun gas) and 12 (1 Msun gas), DM-only versions, MUSIC configs, region files | 5 files (0.6 to 3.5 GB) | 1 node, 480 GB, 64-bit IDs |
| [`parents/`](parents/README.md) | the parent boxes the zooms were cut from (z = 6 box, z = 0 box) with configs and runs | 2 files (0.5 and 4.3 GB) | to regenerate zooms |

Total of all files: 17.0 GB. The z = 0 dwarf zoom ICs do not exist yet (the z = 0 parent box is finished, the halo selection is not).

## Getting the files

```
ICs/get_stella_ics.sh --list                       # names and sizes
ICs/get_stella_ics.sh --prefix $HOME/stella-ics tallbox clouds   # sets, or single files like zoom/ic_1e8_c_L10_gas.dat
ICs/get_stella_ics.sh all
```
Files land in `PREFIX/<set>/<file>`, hash-checked against `manifest.txt` (plain sha256 of the file bytes; the tables use a dataset hash instead). Files above 2 GiB are stored in parts and joined on download. `--from DIR` takes the files from a local copy with the same layout instead of downloading; `STELLA_ICS_URL` overrides the download address (default: the `ics-v1` release of this repository).

## Using a set (common to all)

Every set is run from an empty run directory with three links: the tables, the IC, and the UV background file. The tables are the stella-tables build that GIZMO checks at start-up; **pin a4830ab** (`scripts/get_stella_tables.sh` of GIZMO installs it, and GIZMO refuses to start with other tables).

```
module purge && module load gcc openmpi gsl hwloc fftw-mpi hdf5-mpi           # MPCDF
cd GIZMO_SOURCE_DIR && cp /path/to/ICs/<set>/Config_xxx.sh . && make clean && make CONFIG=Config_xxx.sh EXEC=GIZMO_xxx -j8
mkdir run && cd run
ln -s STELLA_TABLES_INSTALL tables             # the directory holding stellar_tables_unified.hdf5 and the SGChem tables
ln -s PREFIX/<set>/<ic file> .                 # name as in InitCondFile of the parameter file
cp GIZMO_SOURCE_DIR/cooling/TREECOOL .         # UV background table, read from the run directory
cp /path/to/ICs/<set>/<set>.params .
mpirun -np N ../GIZMO_SOURCE_DIR/GIZMO_xxx <set>.params     # add 1 as last argument to restart
```
The parameter files use only relative paths (`tables/`, the IC name in the run directory), so nothing needs editing except `MaxMemSize` (MB per rank) when you change the number of ranks. `MaxMemSize` x ranks must stay below the node memory.

Tested code: Stella branch `stella_dev` at df1f23bfb for tallbox, clouds and wlm; the zoom configuration lives on a zoom branch (a6fb413ef, `merge/zoom-fixes`) and is not yet in `stella_dev`.

## IC conventions

- Units in the hdf5 ICs: kpc, 1e10 Msun, km/s, Gauss (as the parameter files); the Gadget binary ICs carry the same units via `gadget_lunit = kpc` and masses in 1e10 Msun.
- IDs: GIZMO uses 64-bit particle IDs. The zoom and parent Gadget ICs have 64-bit IDs (MUSIC `gadget_longids = yes`): **keep it**, a 32-bit ID block read as 64-bit gives garbage IDs on older code. The hdf5 ICs store unique uint32 IDs (converted on read); the WLM Gadget IC has 32-bit IDs, which `stella_dev` reads by the record length.
- Metals: on a fresh start GIZMO sets all metal fields from `InitMetallicity` (and `InitialMetallicity` for the chemistry) in the parameter file, so the metal blocks of an IC are not used; the parameter files carry the intended value of each IC.
- Chemistry abundances: the cloud ICs carry `ChemcoolTracAbund` (old name of `SGChemTracAbund`, both read); the cloud `Config` keeps them (`SGCHEM_KEEP_IC_ABUNDANCES`), the other configs reset them.

## How each set was checked

Built with the stated `Config` (stella_dev df1f23bfb, zoom branch a6fb413ef; `make clean` first) and started on the login node with the shipped parameter file, no edits other than the end time or `MaxMemSize`. Small sets (clouds): run for 0.1 to 0.2 Myr on 4 cores, exit code 0, no ignored-parameter warnings. Large sets (tallbox, wlm, zoom): parameter file parsed with no warnings and no missing tags (the start was stopped on purpose by a small `MaxMemSize`), because they need a whole node; they were run in production before in these configurations. Details in each set's README.

## Repository

`manifest.txt` hashes and sizes; `get_stella_ics.sh` fetcher. Release assets are named like the manifest paths with `/` replaced by `__` (and `.partNNN` for split files).
