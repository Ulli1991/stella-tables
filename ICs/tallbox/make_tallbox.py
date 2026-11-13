#!/usr/bin/env python3
"""Stratified tallbox IC per Hu, Sternberg & van Dishoeck 2021 (MW-like), Stella units.

Gas: rho(z) = Sigma_g/(2 H_g) sech^2(z/H_g)  [Hu+21 eq. before (3)], uniform in x/y,
T uniform, 27-element scaled-solar abundances. External potential is applied by the
code (GRAVITY_ANALYTIC hook): a_* = -2 pi G Sigma_* tanh(z/H_*) with Sigma_*=40, H_*=0.25 kpc,
plus vertical NFW (M_vir=1e12, c=12, r_s=17 kpc, rho_s=9.5e-3 Msun/pc^3, R0=8 kpc).

Usage:
  python3 make_tallbox.py --out tallbox.hdf5 [--sigma 10] [--Hg 0.25] [--Lxy 1.0]
                          [--zmax 10.0] [--mgas 1.0] [--T 1e4] [--Zsun 0.1] [--seed 42]
"""
import argparse
import numpy as np, h5py

ap = argparse.ArgumentParser()
ap.add_argument('--out', required=True)
ap.add_argument('--sigma', type=float, default=10.0, help='Msun/pc^2')
ap.add_argument('--Hg', type=float, default=0.25, help='gas scale height [kpc] (Hu+21: 0.25)')
ap.add_argument('--Lxy', type=float, default=1.0, help='horizontal box [kpc]')
ap.add_argument('--zmax', type=float, default=10.0, help='half-height [kpc] (Uli: 10; Hu+21: 5)')
ap.add_argument('--mgas', type=float, default=1.0, help='particle mass [Msun]')
ap.add_argument('--T', type=float, default=1.0e4)
ap.add_argument('--Zsun', type=float, default=0.1)
ap.add_argument('--seed', type=int, default=42)
a = ap.parse_args()

rng = np.random.RandomState(a.seed)
KB, MP = 1.3807e-16, 1.6726e-24

M_tot_msun = a.sigma * (a.Lxy * 1000.0)**2          # full Gaussian integral (all mass in box if zmax>>Hg)
N = int(round(M_tot_msun / a.mgas))
print(f"Sigma={a.sigma} Msun/pc^2, Lxy={a.Lxy} kpc, zmax=+-{a.zmax} kpc -> M={M_tot_msun:.3e} Msun, N={N}")

# x,y: jittered lattice for uniformity; z: Gaussian(0, Hg) truncated at +-zmax
nxy = int(np.ceil(np.sqrt(N)))
# sample positions: uniform x,y; Gaussian z via inverse-CDF on truncated normal
x = rng.uniform(0, a.Lxy, N)
y = rng.uniform(0, a.Lxy, N)
# sech^2 profile: CDF = (1 + tanh(z/Hg))/2  ->  z = Hg * atanh(2u - 1)
u = rng.uniform(0, 1, N)
z = a.Hg * np.arctanh(np.clip(2*u - 1, -0.999999999, 0.999999999))
bad = np.abs(z) > a.zmax
while bad.any():
    ub = rng.uniform(0, 1, bad.sum())
    z[bad] = a.Hg * np.arctanh(np.clip(2*ub - 1, -0.999999999, 0.999999999))
    bad = np.abs(z) > a.zmax
z = z + a.zmax   # shift to [0, 2 zmax] box coordinates

pos = np.stack([x, y, z], axis=1)
vel = np.zeros((N, 3), dtype=np.float32)
mass = np.full(N, a.mgas / 1.0e10, dtype=np.float32)

mu = 0.6 if a.T > 6.0e3 else 1.22
u_code = 1.5 * KB * a.T / (mu * MP) / 1.0e10
uu = np.full(N, u_code, dtype=np.float32)

ZSUN = 0.014
X_metals = a.Zsun * ZSUN; X_He = 0.275; X_H = 1 - X_He - X_metals
EA = np.zeros((N, 27), dtype=np.float32)
EA[:, 0] = X_H; EA[:, 1] = X_He
solar = {2: 2.53e-3, 3: 7.41e-4, 4: 6.13e-3, 5: 1.34e-3, 6: 7.57e-4,
         7: 7.12e-4, 8: 3.31e-4, 9: 6.87e-5, 10: 1.38e-3}
for k, xs in solar.items():
    EA[:, k] = xs * a.Zsun
EA[:, 0] = X_H + (1.0 - EA.sum(axis=1))

with h5py.File(a.out, 'w') as f:
    h = f.create_group('Header')
    h.attrs['BoxSize'] = a.Lxy                       # base box; z via BOX_LONG_Z multiplier
    h.attrs['NumPart_ThisFile'] = np.array([N, 0, 0, 0, 0, 0], dtype=np.int32)
    h.attrs['NumPart_Total'] = np.array([N, 0, 0, 0, 0, 0], dtype=np.uint32)
    h.attrs['NumPart_Total_HighWord'] = np.zeros(6, dtype=np.uint32)
    h.attrs['MassTable'] = np.zeros(6, dtype=np.float64)
    h.attrs['Time'] = 0.0; h.attrs['Redshift'] = 0.0
    h.attrs['NumFilesPerSnapshot'] = 1
    h.attrs['Flag_Sfr'] = 1; h.attrs['Flag_Cooling'] = 1; h.attrs['Flag_StellarAge'] = 1
    h.attrs['Flag_Metals'] = 27; h.attrs['Flag_Feedback'] = 1; h.attrs['Flag_DoublePrecision'] = 0
    h.attrs['HubbleParam'] = 1.0; h.attrs['ComovingIntegrationOn'] = 0
    h.attrs['UnitLength_In_CGS'] = 3.085678e+21
    h.attrs['UnitMass_In_CGS'] = 1.989e+43
    h.attrs['UnitVelocity_In_CGS'] = 1.0e5
    g = f.create_group('PartType0')
    g.create_dataset('Coordinates', data=pos)
    g.create_dataset('Velocities', data=vel)
    g.create_dataset('Masses', data=mass)
    g.create_dataset('ParticleIDs', data=np.arange(N, dtype=np.uint32) + 10)  # IDs>=10: promotion-safe
    g.create_dataset('InternalEnergy', data=uu)
    g.create_dataset('ElementAbundance', data=EA)
print(f"wrote {a.out}  (z in [0, {2*a.zmax}] kpc box coords; needs BOX_LONG_Z={int(2*a.zmax/a.Lxy)})")
