#!/usr/bin/env python3
"""build_ic.py: convert MakeCloud raw output (pc, Msun, km/s, Gauss) to the paper-B GIZMO IC convention
(kpc, 1e10 Msun, km/s, Gauss; Header attrs + ChemcoolTracAbund as gmc_diag/ics/gmc_*.hdf5) and write <name>_verify.json.
usage: build_ic.py raw.hdf5 out_basename variant(atomic20K|mol10K) [pairwise]"""
import sys, json, h5py, numpy as np
raw, base, variant = sys.argv[1:4]
pairwise = len(sys.argv) > 4
f = h5py.File(raw, 'r'); p = f['PartType0']
m = p['Masses'][:]; x = p['Coordinates'][:]; v = p['Velocities'][:]; B = p['MagneticField'][:]
boxpc = float(f['Header'].attrs['BoxSize']); N = len(m)
mcell = m[0]; Nc = int(round(2e4 / mcell)); assert np.allclose(m, mcell, rtol=1e-12)
cl = np.arange(N) < Nc                                   # MakeCloud: cloud cells first, then diffuse box
xc = x - boxpc / 2; r = np.sqrt((xc**2).sum(1))
# thermal / chemistry (see README): kT/mp = 8.254e-3 (km/s)^2/K ; mu from H mass fraction implied by existing IC (mu 2.33 mol)
KMP = 8.2544e-3; Xh = 0.717; mu_at = 1 / (Xh + (1 - Xh) / 4)
if variant == 'atomic20K':
    uc = ub = 1.5 * KMP * 20.0 / mu_at; xh2, xco = 0.0, 0.0
elif variant == 'mol10K':
    uc = 0.0531; ub = 77.84; xh2, xco = 0.4999, 1.4446454e-4   # numbers copied from gmc_Z1_m1.0.hdf5 / snap_000.hdf5
else: raise SystemExit('variant')
u = np.where(cl, uc, ub)
tr = np.zeros((N, 3)); tr[:, 1] = 1e-4; tr[cl, 0] = xh2; tr[cl, 2] = xco
L, MU, KMS = 1e-3, 1e-10, 1.0                           # pc->kpc, Msun->1e10 Msun
out = base + '.hdf5'
with h5py.File(out, 'w') as g:
    h = g.create_group('Header')
    for k, val in dict(BoxSize=boxpc * L, Flag_Cooling=0, Flag_Feedback=0, Flag_Metals=0, Flag_Sfr=0, Flag_StellarAge=0,
                       HubbleParam=1.0, MassTable=np.zeros(6), NumFilesPerSnapshot=1,
                       NumPart_ThisFile=np.array([N, 0, 0, 0, 0, 0], dtype=np.int32),
                       NumPart_Total=np.array([N, 0, 0, 0, 0, 0], dtype=np.uint32),
                       NumPart_Total_HighWord=np.zeros(6, dtype=np.uint32), Omega0=0.0, OmegaLambda=0.0, Redshift=0.0, Time=0.0).items():
        h.attrs[k] = val
    q = g.create_group('PartType0')
    q['ChemcoolTracAbund'] = tr; q['Coordinates'] = x * L; q['InternalEnergy'] = u; q['Masses'] = m * MU
    q['ParticleIDs'] = p['ParticleIDs'][:].astype(np.uint32); q['Velocities'] = v
    q['MagneticField'] = B                               # Gauss; ignored when MHD_B_SET_IN_PARAMS (set BiniZ=2e-6)
# verification (physical units)
G = 4.30091e-3                                           # pc (km/s)^2 / Msun
Mc = m[cl].sum(); Rm = r[cl].max(); vc = v[cl]; vcm = (vc * m[cl, None]).sum(0) / Mc
sig3 = np.sqrt(((vc - vcm)**2).sum(1).mean()); Ek = 0.5 * Mc * sig3**2
Eg_sph = 0.6 * G * Mc**2 / Rm
Eg_pair = None
if pairwise:
    xs = xc[cl]; s = 0.0
    for i in range(0, Nc, 500):
        d = np.sqrt(((xs[i:i + 500, None, :] - xs[None, :, :])**2).sum(2)); d[d == 0] = np.inf
        s += (1 / d).sum() * mcell * mcell
    Eg_pair = -0.5 * G * s
rho = Mc / (4 / 3 * np.pi * Rm**3); Msun_pc3 = 1.989e33 / 3.0857e18**3; mp = 1.6726e-24
rho_cgs = rho * Msun_pc3
Bmag = np.sqrt((B[cl]**2).sum(1)).mean(); Bg = Bmag
Emag = Bg**2 / (8 * np.pi) * (4 / 3 * np.pi * Rm**3) * 3.0857e18**3          # erg
Eg_cgs = Eg_sph * 1.989e33 * 1e10; Ek_cgs = Ek * 1.989e33 * 1e10
Msun = 1.989e33; Gc = 6.674e-8
mu_B = (Mc * Msun / (np.pi * (Rm * 3.0857e18)**2 * Bg)) * np.sqrt(Gc) / 0.4          # Grudic+22 Table 1 constant 0.4 (see README)
bgm = ~cl
j = dict(M_cloud_Msun=float(Mc), M_bg_Msun=float(m[bgm].sum()), N_cloud=int(Nc), N_bg=int(N - Nc), m_cell_Msun=float(mcell),
         R_pc=float(Rm), box_pc=boxpc, nH_mean=float(rho_cgs * Xh / mp),
         rho_bg_over_rho_cloud=float((m[bgm].sum() / (boxpc**3 - 4 / 3 * np.pi * Rm**3)) / rho),
         sigma3D_kms=float(sig3), alpha_vir_2Ek_Eg=float(2 * Ek / Eg_sph), alpha_vir_BM92=float(5 * (sig3 / np.sqrt(3))**2 * Rm / (G * Mc)),
         alpha_vir_2Ek_Eg_pairwise=(float(2 * Ek / abs(Eg_pair)) if Eg_pair else None),
         t_ff_Myr=float(np.sqrt(3 * np.pi / (32 * G * rho)) * 0.9778),
         B_uG=float(Bg * 1e6), B_dir=[float(c) for c in (B[cl].mean(0) / Bmag)], Emag_over_Egrav=float(Emag / Eg_cgs),
         alphaB_2Emag_Egrav=float(2 * Emag / Eg_cgs), Emag_over_Eturb=float(Emag / Ek_cgs), mu_norm_0p4_def=float(mu_B),
         xH2_cloud=xh2, xCO_cloud=xco, T_cloud_K=(20.0 if variant == 'atomic20K' else 10.0),
         T_bg_K=(20.0 if variant == 'atomic20K' else 'u=77.84 (km/s)^2 as gmc_Z1_*'), u_cloud_kms2=uc, u_bg_kms2=ub,
         net_momentum_cloud_Msun_kms=[float(c) for c in (vc * m[cl, None]).sum(0)],
         bg_max_speed_kms=float(np.abs(v[bgm]).max()),
         units=dict(UnitLength_in_cm=3.085678e21, UnitMass_in_g=1.989e43, UnitVelocity_in_cm_per_s=1e5, HubbleParam=1.0),
         source=dict(raw=raw, MakeCloud='MakeCloud.py (git cd7b4f3)', variant=variant, turb_seed=42))
json.dump(j, open(base + '_verify.json', 'w'), indent=1)
print(json.dumps(j, indent=1))
