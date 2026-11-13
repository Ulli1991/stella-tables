# STARFORGE reference cloud M2e4 -- initial conditions for paper B

Cloud: Grudic+22 (MNRAS 512, 216, arXiv:2201.00882) fiducial M2e4 = Guszejnov+22 (MNRAS 515, 4929, arXiv:2205.10413) Table 1 / C_M_J_R_W.
Generator: `$MAKECLOUD/MakeCloud.py` (git cd7b4f3), python-waterboa/2025.06 (python 3.13; system python 3.6 cannot run it, needs `match`).
Post-processing / verification: `build_ic.py` (raw MakeCloud pc/Msun/km/s/Gauss -> paper-B convention kpc / 1e10 Msun / km/s / Gauss).

## Files
| file | cells | variant |
|---|---|---|
| `M2e4_m1.0_atomic20K.hdf5` (+`_verify.json`) | 24754 (20000 cloud + 4754 box) | paper-literal: all H atomic, all gas 20 K |
| `M2e4_m1.0_mol10K.hdf5` | 24754 | paper-B legacy convention (as `gmc_Z1_m1.0`, `A_B/snap_000`): cloud xH2=0.4999, xCO=1.4446e-4, cloud 10 K, box gas u=77.84 (km/s)^2 |
| `M2e4_m0.1_{atomic20K,mol10K}.hdf5` | 247550 | same at 0.1 Msun |
`raw/` = untouched MakeCloud output + generated MakeCloud params files (not used).
Geometry, velocities, B and masses are IDENTICAL between variants; only InternalEnergy and ChemcoolTracAbund differ.

## MakeCloud command (both resolutions; --N=20000 or 200000)
```
module load python-waterboa/2025.06
python3 $MAKECLOUD/MakeCloud.py --M=2e4 --R=10 --N=20000 --alpha_turb=2 --turb_slope=2 --turb_sol=0.5 \
  --minmode=2 --turb_seed=42 --bfixed=2e-6 --Z=1 --ISRF=1.7 --turb_path=$PWD/turb --glass_path=$HOME/.makecloud_glass/glass_256.hdf5
python3 build_ic.py raw/<MakeCloud output>.hdf5 M2e4_m1.0_atomic20K atomic20K pairwise
python3 build_ic.py raw/<MakeCloud output>.hdf5 M2e4_m1.0_mol10K    mol10K    pairwise
```
Default units of this MakeCloud = pc, Msun, km/s, Gauss. Box = MakeCloud default 10 R = 100 pc (periodic domain in STARFORGE), diffuse box gas at rho_cloud/1000 (MakeCloud default), equal-mass cells, not derefined.

## Parameters: paper value, source, IC value
| parameter | paper value (quote) | source | IC value |
|---|---|---|---|
| mass | "uniform-density sphere with mass M0=2e4 Msun" | G22 S2.8, Table 1 | 20000.0 Msun |
| radius | R0 = 10 pc | G22 S2.8, Table 1; Gu22 Table 1 | 10.004 pc (outermost cell, N=2e4); `--R=10` |
| density profile | uniform (n_H=146 cm^-3, Sigma=64 Msun/pc^2, t_ff=3.7 Myr) | G22 S2.8 Table 1; Gu22 S2.1.3 "uniform density, surrounded by diffuse gas with a density contrast of 1000" | uniform, `--density_exponent=0`; t_ff = 3.705 Myr; n_H = 138 cm^-3 with X=0.717 (the paper's 146 implies X~0.756; mass/radius identical, definition differs) |
| turbulent virial parameter | alpha_turb = 2 (BM92, 2E_turb/|E_grav|) "marginally bound" | G22 Table 1, S2.8; Gu22 S2.1.3 | 2.001 (uniform-sphere E_g), 2.004 (pairwise) |
| velocity field | "pseudo-turbulent random velocity field constructed in Fourier space", E_k ~ k^-2, E_sol = 2 E_comp, no driving | G22 S2.8; Gu22 S2.1.3 | Gaussian random field, `--turb_slope=2`, `--turb_sol=0.5` (= natural 2:1 solenoidal:compressive mix of the projection), `--minmode=2` |
| sigma_3D | 3.2 km/s | Gu22 Table 1 | 3.213 km/s |
| turbulence seed | NOT STATED | -- | MakeCloud default 42 |
| min populated wavenumber | NOT STATED | -- | MakeCloud default minmode=2 (in pi/R) |
| magnetic field | "initially uniform in the +z direction with a strength B=2 uG"; mu0=4.2; alpha_B=0.02; E_mag ~1% of turbulence | G22 S2.8, S3.2; Gu22 Table 1 (mu=4.2, alpha_B=0.02) | uniform +z, 2.000 uG (also in the diffuse gas, MakeCloud). E_mag/|E_g|=0.00956, 2E_mag/|E_g|=0.0191, E_mag/E_turb=0.0096. mu=4.29 with the Table-1 constant 0.4 (see caveat). `--bfixed=2e-6` (MakeCloud default bturb=0.1 NOT used; that would be 10x the field energy) |
| temperature | "The gas temperature is also initialized to 20 K" (Gu22 S2.1.3: 10 K) | G22 S2.8 vs Gu22 | atomic20K: 20 K (u=0.1951 (km/s)^2, mu=1.27); mol10K: cloud 10 K (u=0.0531, as legacy ICs) |
| composition / chemistry | "All gas is initially of Solar composition"; "all hydrogen starts atomic" | G22 S2.8 | atomic20K: ChemcoolTracAbund=(xH2 0, xH+ 1e-4, xCO 0) (the 1e-4 ion floor is ours, copied from legacy ICs); mol10K: cloud (0.4999, 1e-4, 1.4446e-4) |
| metallicity | Solar, no number given | G22 S2.8; Gu22 Table 2 | Z=1 Zsun (`--Z=1`, params InitMetallicity 1.0); no per-cell metal fields in IC (same as legacy ICs) |
| ISRF | "fixed 1.7 Habing background" (FUV 6-13.6 eV); FIR 20 K blackbody 0.3 eV/cm^3 | G22 S2.4, S2.8, Table 1 | params `G0 1.7` (already in run params); MakeCloud `--ISRF=1.7` only for its params file |
| cosmic rays | no numerical rate | G22 | not in IC; run params CosmicRayIonRate 3e-17 (ours) |
| box / ambient medium | periodic L=100 pc; rest of box at rho_cloud/1000, ~5000 Msun, static (Gu22: "periodic 10 R_cloud box") | G22 S2.8; Gu22 S2.1.3 | box 100 pc = 0.1 kpc, 4754 Msun (rho ratio 1.001e-3), v=0; ambient temperature NOT STATED -> atomic20K: 20 K ("all gas"); mol10K: u=77.84 as legacy ICs. MakeCloud's own default would be u=101 (km/s)^2 (~1e4 K neutral) |
| mass resolution | 1e-3 Msun/cell (2e7 cells) | G22 S2.8 | 1 Msun (2e4 cells) and 0.1 Msun (2e5 cells): test resolution, deliberately coarser |
| run end / SFE | "run to 9.3Myr"; "integrated star formation efficiency of 8%" (~1600 Msun); Gu22: ~7-8%, 6-8% across realisations; cloud dispersed ~6 Myr | G22 S3, S3.1, S3.5 | n/a (9.3 Myr = 2.51 t_ff for this cloud) |

Not stated in either paper: turbulence seed, minmode, ambient temperature, numerical metallicity/dust-to-gas, CR rate, Sphere-box size for Gu22 (Table 1 lists L_box=16 pc, but that is the Box-run size; text says 10 R_cloud = 100 pc).

## Caveats
- Temperature: G22 says 20 K, Gu22 Table says 10 K. Chemistry equilibrates in << t_ff; the two variants bracket this.
- Atomic start: H2 then forms on a timescale ~1e9 yr/n ~ several Myr at n~150, i.e. the cloud starts mostly atomic exactly as in STARFORGE; mol10K avoids this transient but is NOT the paper setup. Config flag CHEMCOOL_KEEP_IC_ABUNDANCES is on in the cloud configs, so IC abundances are used as given (u is set for the stated composition).
- mu: Gu22 S2.1.3 text says alpha_B=0.02 "translates to mu=0.4" while Table 1/2 and G22 give mu=4.2; the Table-1 expression as extracted is "0.4 G^1/2 M/(pi R^2 B0)", which evaluates to 0.69 for B=2 uG, whereas 4.29 is obtained with G^1/2 M/(0.4 pi R^2 B). I used the latter (reproduces 4.2 to 2%) -- the printed form was read from an HTML scrape and should be checked against the PDF. B=2 uG (explicit) was used rather than tuning to mu=4.2 exactly (would be 2.04 uG).
- The IC carries a MagneticField dataset (Gauss, as MakeCloud writes). The cloud configs use MHD_B_SET_IN_PARAMS, so the file B is ignored: set `BiniX 0, BiniY 0, BiniZ 2.0e-6`.
- 1 Msun cells: the Jeans-resolution floor is far above STARFORGE (Delta m = 1e-3); results are an ablation at test resolution, not a STARFORGE reproduction.
