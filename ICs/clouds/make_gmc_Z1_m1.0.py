# R26 molecular GMC start (CHEMCOOL_KEEP_IC_ABUNDANCES): Grudic's MakeCloud turbulent cloud (1e4 Msun, R=7 pc,
# alpha_vir=2, no B, uniform background), converted to the Stella params units, molecular IC abundances (net 5).
# Short run of the production model at 1 Msun. Grade: IC abundances survive the fresh start (first snapshot is
# molecular, flag off gives x_H2 ~ 5e-4 after 20 kyr), cloud stays molecular, stars form (first at ~0.6 Myr),
# mass is conserved, no NaN.
# Also the IC generator: python3 r26_gmc_molecular.py --mgas 1.0 --out snap_000.hdf5 [--Z 1 --R 7 ...]
import argparse, glob, json, os, re, shutil, subprocess, sys, tempfile
import numpy as np
import h5py

PROD_PARAMS = os.environ.get('STELLA_PARAMS', 'stella_cloud.params')
MAKECLOUD = os.environ.get('MAKECLOUD', os.path.expanduser('~/MakeCloud/MakeCloud.py'))
MAKECLOUD_PY = os.environ.get('MAKECLOUD_PY', '/mpcdf/soft/SLE_15/packages/x86_64/python-waterboa/2025.06/bin/python3')
MAKECLOUD_PYPATH = os.environ.get('MAKECLOUD_PYPATH', '')   # docopt
HERE = os.path.abspath(__file__)

CASE = dict(
    doc='MakeCloud GMC, molecular IC kept on fresh start, production model at 1 Msun to ~1 Myr (first stars + HII)',
    requires_flag='CHEMCOOL_KEEP_IC_ABUNDANCES',
    base='Config_stella_cloud.sh',
    flags_on=['CHEMCOOL_KEEP_IC_ABUNDANCES', 'STELLA_SELFCHECK'],
    params_file=PROD_PARAMS,
    ic_cmd='python3 ' + HERE + ' --mgas 1.0 --out {out}',
    params=dict(TimeMax=0.001, TimeBetSnapshot=0.0001, TimeOfFirstSnapshot=0, MaxSizeTimestep=0.00005,
                InitMetallicity=1.0, InitialMetallicity=1.0, BiniX=0, BiniY=0, BiniZ=0),
    ranks=4, timeout_min=20,   # ~1 min on 4 login ranks
    selfchecks=['mass_sync', 'tree_insert'],
)

# MakeCloud defaults (passed explicitly): pc, Msun, 1e3 m/s, gauss
MC_UL, MC_UM, MC_UV = 3.085678e18, 1.989e33, 1.0e5
M_P, K_B, G_CGS, MSUN, PC = 1.6726e-24, 1.380649e-16, 6.674e-8, 1.989e33, 3.085678e18
ABHE = 0.1                              # net 5 (chemcool_consts.h)
X_H_CHEM = 1.0 / (1.0 + 4.0 * ABHE)     # H mass fraction the chemistry uses for n_H
ABUNDC_SOLAR = 2.4765351227e-03 / 12.0 / X_H_CHEM   # init.cc solar C mass fraction -> n_C/n_H (chem_interface.cc)


def stella_units(params):
    u = {}
    for l in open(params):
        p = l.split()
        if len(p) >= 2 and p[0] in ('UnitLength_in_cm', 'UnitMass_in_g', 'UnitVelocity_in_cm_per_s', 'UnitMagneticField_in_gauss', 'HubbleParam'):
            u.setdefault(p[0], float(p[1]))
    return u


def run_makecloud(M, R, N, alpha, Z, seed, wd):
    env = dict(os.environ, PYTHONPATH=MAKECLOUD_PYPATH)
    cmd = [MAKECLOUD_PY, MAKECLOUD, '--M=%g' % M, '--R=%g' % R, '--N=%d' % N, '--alpha_turb=%g' % alpha, '--bturb=0',
           '--turb_seed=%d' % seed, '--Z=%g' % Z, '--length_unit=1', '--mass_unit=1', '--v_unit=1e3', '--B_unit=1']
    subprocess.run(cmd, cwd=wd, env=env, check=True, stdout=subprocess.DEVNULL)
    pf = glob.glob(os.path.join(wd, 'params_*.txt'))[0]
    mc = {}
    for l in open(pf):
        p = l.split()
        if len(p) >= 2 and p[0].startswith('Unit'): mc[p[0]] = float(p[1])
    for k, v in (('UnitLength_in_cm', MC_UL), ('UnitMass_in_g', MC_UM), ('UnitVelocity_in_cm_per_s', MC_UV)):
        assert abs(mc[k] / v - 1) < 0.01, ('MakeCloud units changed', k, mc[k])
    return [f for f in glob.glob(os.path.join(wd, '*.hdf5'))][0]


def mu_of(xh2, xhp):
    return (1.0 + 4.0 * ABHE) / (1.0 - xh2 + 2.0 * xhp + ABHE)   # per H: HI + H2 + H+ + e- + He


def make_ic(out, mgas=1.0, M=1e4, R=7.0, alpha=2.0, Z=1.0, seed=42, fco=0.5, xh2=0.4999, T_cloud=10.0, T_bg=8000.0,
            params=PROD_PARAMS):
    N = int(round(M / mgas))
    wd = tempfile.mkdtemp(prefix='makecloud_')
    try:
        f = run_makecloud(M, R, N, alpha, Z, seed, wd)
        with h5py.File(f, 'r') as F:
            g = F['PartType0']; x = g['Coordinates'][:]; v = g['Velocities'][:]; m = g['Masses'][:]
            box = float(F['Header'].attrs['BoxSize'])
            assert 'MagneticField' not in g
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    su = stella_units(params); h = su['HubbleParam']
    UL, UM, UV = su['UnitLength_in_cm'], su['UnitMass_in_g'], su['UnitVelocity_in_cm_per_s']
    fx, fm, fv = MC_UL / UL * h, MC_UM / UM * h, MC_UV / UV        # MakeCloud -> Stella code units
    x, m, v, box = x * fx, m * fm, v * fv, box * fx
    ncl = N                                                       # MakeCloud writes the cloud cells first
    cloud = np.zeros(len(m), bool); cloud[:ncl] = True
    xco = fco * ABUNDC_SOLAR * Z
    trac = np.zeros((len(m), 3)); trac[:, 1] = 1e-4               # net 5: H2, H+, CO (per H nucleus)
    trac[cloud, 0] = xh2; trac[cloud, 2] = xco
    T = np.where(cloud, T_cloud, T_bg)
    mu = np.where(cloud, mu_of(xh2, 1e-4), mu_of(0.0, 1e-4))
    u = 1.5 * K_B * T / (mu * M_P) / UV**2                        # specific energy in Stella code units
    with h5py.File(out, 'w') as F:
        H = F.create_group('Header'); n = len(m)
        H.attrs['NumPart_ThisFile'] = np.array([n, 0, 0, 0, 0, 0], dtype=np.int32)
        H.attrs['NumPart_Total'] = np.array([n, 0, 0, 0, 0, 0], dtype=np.uint32)
        H.attrs['NumPart_Total_HighWord'] = np.zeros(6, dtype=np.uint32)
        H.attrs['MassTable'] = np.zeros(6); H.attrs['Time'] = 0.0; H.attrs['Redshift'] = 0.0
        H.attrs['BoxSize'] = box; H.attrs['NumFilesPerSnapshot'] = 1
        for k in ('Omega0', 'OmegaLambda'): H.attrs[k] = 0.0
        H.attrs['HubbleParam'] = h
        for k in ('Flag_Sfr', 'Flag_Cooling', 'Flag_StellarAge', 'Flag_Metals', 'Flag_Feedback'): H.attrs[k] = 0
        g = F.create_group('PartType0')
        g['Coordinates'] = x; g['Velocities'] = v; g['Masses'] = m
        g['ParticleIDs'] = np.arange(1, n + 1, dtype=np.uint32)
        g['InternalEnergy'] = u; g['ChemcoolTracAbund'] = trac
    return verify(out, ncl, params)


def verify(f, ncl, params=PROD_PARAMS):
    """recompute cloud properties from the written IC, in Stella code units -> physical"""
    su = stella_units(params); h = su['HubbleParam']
    UL, UM, UV = su['UnitLength_in_cm'], su['UnitMass_in_g'], su['UnitVelocity_in_cm_per_s']
    G = G_CGS * UM / UL / UV**2                                    # code G
    with h5py.File(f, 'r') as F:
        g = F['PartType0']; x = g['Coordinates'][:ncl] / h; v = g['Velocities'][:ncl]; m = g['Masses'][:ncl] / h
        mb = g['Masses'][ncl:] / h; box = F['Header'].attrs['BoxSize'] / h; xa = g['ChemcoolTracAbund'][:]
        u = g['InternalEnergy'][:]
    M = m.sum(); c = (m[:, None] * x).sum(0) / M; r = np.sqrt(((x - c)**2).sum(1)); R = r.max()
    vc = v - (m[:, None] * v).sum(0) / M; s3 = np.sqrt((m * (vc**2).sum(1)).sum() / M)
    o = np.argsort(r); Mr = np.cumsum(m[o]); Egrav = G * np.sum(Mr * m[o] / np.maximum(r[o], 1e-30))
    Ekin = 0.5 * np.sum(m * (vc**2).sum(1))
    rho = M * UM / (4 / 3 * np.pi * (R * UL)**3)
    res = dict(M_cloud_Msun=M * UM / MSUN, M_bg_Msun=mb.sum() * UM / MSUN, N_cloud=int(ncl), N_bg=int(len(mb)),
               m_cell_Msun=float(np.median(m) * UM / MSUN), R_pc=R * UL / PC, box_pc=box * UL / PC,
               nH_mean=rho * X_H_CHEM / M_P, sigma3D_kms=s3 * UV / 1e5,
               alpha_vir_2Ek_Eg=2 * Ekin / Egrav, alpha_vir_BM92=5 * (s3**2 / 3) * R / (G * M),
               t_ff_Myr=np.sqrt(3 * np.pi / (32 * G_CGS * rho)) / 3.156e13,
               xH2_cloud=float(xa[:ncl, 0].mean()), xCO_cloud=float(xa[:ncl, 2].mean()),
               T_cloud_K=float(u[0] * UV**2 * mu_of(xa[0, 0], xa[0, 1]) * M_P / 1.5 / K_B),
               units=dict(UnitLength_in_cm=UL, UnitMass_in_g=UM, UnitVelocity_in_cm_per_s=UV, HubbleParam=h))
    return {k: (round(float(v), 6) if isinstance(v, (float, np.floating)) else v) for k, v in res.items()}


def grade(w, txt):
    sn = sorted(glob.glob(os.path.join(w, 'output', 'snapshot_*.hdf5')))
    if len(sn) < 2: return {'snapshots': {'pass': False, 'n': len(sn)}}
    def load(s):
        with h5py.File(s, 'r') as F:
            g = F['PartType0']; ms = F['PartType4/Masses'][:].sum() if 'PartType4' in F else 0.0
            return g['Masses'][:].astype(float), g['Density'][:].astype(float), g['ChemcoolTracAbund'][:], g['InternalEnergy'][:], ms, F['Header'].attrs['Time']
    m0, d0, a0, u0, s0, t0 = load(sn[0]); m1, d1, a1, u1, s1, t1 = load(sn[-1])
    dense0 = d0 > np.percentile(d0, 30)
    nH1 = d1 * 1.989e43 / 3.085678e21**3 * X_H_CHEM / M_P
    fH2_0 = float(np.sum(m0[dense0] * 2 * a0[dense0, 0]) / np.sum(m0[dense0]))
    mol1 = nH1 > 100
    fH2_1 = float(np.sum(m1[mol1] * 2 * a1[mol1, 0]) / max(np.sum(m1[mol1]), 1e-30))
    mtot0, mtot1 = m0.sum() + s0, m1.sum() + s1
    finite = all(np.isfinite(q).all() for q in (m1, d1, a1, u1))
    return {'ic_molecular_kept': {'pass': bool(fH2_0 > 0.9), 'f_H2_dense_t0': round(fH2_0, 4), 't0': float(t0)},
            'stays_molecular': {'pass': bool(fH2_1 > 0.5), 'f_H2_nH>100_tend': round(fH2_1, 4), 't_end': float(t1)},
            'stars_formed': {'pass': bool(s1 > 0), 'M_star_Msun': float(s1 * 1e10)},
            'mass_conserved': {'pass': bool(abs(mtot1 / mtot0 - 1) < 1e-6), 'rel': float(mtot1 / mtot0 - 1)},
            'finite': {'pass': bool(finite)}}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True); ap.add_argument('--mgas', type=float, default=1.0)
    ap.add_argument('--M', type=float, default=1e4); ap.add_argument('--R', type=float, default=7.0)
    ap.add_argument('--alpha', type=float, default=2.0); ap.add_argument('--Z', type=float, default=1.0)
    ap.add_argument('--seed', type=int, default=42); ap.add_argument('--fco', type=float, default=0.5)
    ap.add_argument('--params', default=PROD_PARAMS)
    a = ap.parse_args()
    r = make_ic(a.out, a.mgas, a.M, a.R, a.alpha, a.Z, a.seed, a.fco, params=a.params)
    print(json.dumps(r, indent=1))
    json.dump(r, open(os.path.splitext(a.out)[0] + '_verify.json', 'w'), indent=1)
