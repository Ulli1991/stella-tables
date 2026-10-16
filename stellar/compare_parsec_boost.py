#!/usr/bin/env python3
"""
Compare PARSEC v1.2s vs BoOST tracks for massive stars (>= 20 Msun)
at MW-like metallicity (PARSEC Z=0.008, BoOST MW Z=0.00876).

Quantities compared:
  1. Current mass M(t)
  2. Mass loss rate Mdot(t)
  3. Ionizing photon rate Q_ion(t) — PI radiation
  4. FUV luminosity L_FUV(t) — PE heating (6-11.2 eV)
  5. Bolometric luminosity L_bol(t) — radiation pressure proxy (L/c)
  6. Stellar radius R(t)
"""

import numpy as np
import os
import re
import sys
from glob import glob
from scipy.integrate import quad

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ── Physical constants (CGS) ──
h_planck = 6.62607015e-27
k_boltz  = 1.380649e-16
c_light  = 2.99792458e10
sigma_sb = 5.670374419e-5
eV_to_erg = 1.602176634e-12
L_sun    = 3.828e33
R_sun    = 6.957e10

E_FUV_lo = 6.0    # eV
E_FUV_hi = 11.2   # eV
E_LW_lo  = 11.2   # eV
E_LW_hi  = 13.6   # eV


# ═══════════════════════════════════════════════════════════════
# Blackbody radiation computation
# ═══════════════════════════════════════════════════════════════

def _planck_photon_integrand(x):
    if x > 500: return 0.0
    return x**2 / (np.exp(x) - 1.0)

def _planck_energy_integrand(x):
    if x > 500: return 0.0
    return x**3 / (np.exp(x) - 1.0)

def compute_radiation(Teff, L_bol_erg):
    """Given Teff [K] and L_bol [erg/s], compute Q_ion, L_FUV, L_LW."""
    kT = k_boltz * Teff
    x_fuv_lo = E_FUV_lo * eV_to_erg / kT
    x_fuv_hi = E_FUV_hi * eV_to_erg / kT
    x_lw_lo  = E_LW_lo  * eV_to_erg / kT
    x_lw_hi  = E_LW_hi  * eV_to_erg / kT

    I_Qion = 0.0
    I_FUV  = 0.0
    I_LW   = 0.0

    if x_lw_hi < 200:
        I_Qion, _ = quad(_planck_photon_integrand, x_lw_hi, min(700, x_lw_hi+200), limit=200)
    if x_fuv_lo < 200:
        I_FUV, _ = quad(_planck_energy_integrand, x_fuv_lo, min(x_fuv_hi, 500), limit=200)
    if x_lw_lo < 200:
        I_LW, _ = quad(_planck_energy_integrand, x_lw_lo, min(x_lw_hi, 500), limit=200)

    inv_sigT4 = 1.0 / (sigma_sb * Teff**4)
    Q_ion = L_bol_erg * inv_sigT4 * (2*np.pi / c_light**2) * (kT/h_planck)**3 * I_Qion
    energy_pref = L_bol_erg * inv_sigT4 * (2*np.pi*h_planck / c_light**2) * (kT/h_planck)**4
    L_FUV = energy_pref * I_FUV
    L_LW  = energy_pref * I_LW

    return Q_ion, L_FUV, L_LW


def compute_radiation_vectorized(Teff_arr, L_bol_erg_arr):
    """Vectorized version for arrays."""
    n = len(Teff_arr)
    Q_ion = np.zeros(n)
    L_FUV = np.zeros(n)
    L_LW  = np.zeros(n)
    for i in range(n):
        if Teff_arr[i] > 3000 and L_bol_erg_arr[i] > 0:
            Q_ion[i], L_FUV[i], L_LW[i] = compute_radiation(Teff_arr[i], L_bol_erg_arr[i])
    return Q_ion, L_FUV, L_LW


# ═══════════════════════════════════════════════════════════════
# PARSEC track loading
# ═══════════════════════════════════════════════════════════════

def load_parsec_track(base_dir, Z_target, M_target):
    """Load a single PARSEC track closest to (Z_target, M_target)."""
    zdirs = sorted(glob(os.path.join(base_dir, 'Z*Y*')))
    zdirs = [d for d in zdirs if os.path.isdir(d)]

    # Find closest Z
    best_zdir = None
    best_zdist = 1e10
    for zdir in zdirs:
        m = re.match(r'Z([\d.]+)Y', os.path.basename(zdir))
        if m:
            Z = float(m.group(1))
            dist = abs(np.log10(Z) - np.log10(Z_target))
            if dist < best_zdist:
                best_zdist = dist
                best_zdir = zdir

    # Find closest mass
    dat_files = sorted(glob(os.path.join(best_zdir, '*.DAT')))
    best_file = None
    best_mdist = 1e10
    best_mass = None
    for f in dat_files:
        m = re.search(r'_M([\d.]+)\.DAT', os.path.basename(f))
        if m:
            M = float(m.group(1))
            dist = abs(M - M_target)
            if dist < best_mdist:
                best_mdist = dist
                best_file = f
                best_mass = M

    Z_actual = float(re.match(r'Z([\d.]+)Y', os.path.basename(best_zdir)).group(1))

    # Parse
    data = []
    with open(best_file) as fh:
        header = fh.readline().split()
        for line in fh:
            vals = line.split()
            if len(vals) == len(header):
                data.append([float(v) for v in vals])
    arr = np.array(data)
    d = {name: arr[:, i] for i, name in enumerate(header)}

    # Compute radiation
    Teff = 10**d['LOG_TE']
    L_bol = 10**d['LOG_L'] * L_sun
    Q_ion, L_FUV, L_LW = compute_radiation_vectorized(Teff, L_bol)

    return {
        'source': 'PARSEC',
        'Z': Z_actual, 'M_init': best_mass,
        'age_yr': d['AGE'],
        'mass': d['MASS'],
        'Teff': Teff,
        'logL': d['LOG_L'],
        'L_bol': L_bol,
        'R_Rsun': 10**d['LOG_R'],  # PARSEC LOG_R is log(R/Rsun)
        'log_Mdot': d['LOG_RAT'],  # log(Mdot in Msun/yr)
        'Q_ion': Q_ion,
        'L_FUV': L_FUV,
        'L_LW': L_LW,
    }


# ═══════════════════════════════════════════════════════════════
# BoOST track loading
# ═══════════════════════════════════════════════════════════════

def load_boost_track(boost_dir, grid_name, M_target):
    """Load a single BoOST track closest to M_target from given grid."""
    grid_dir = os.path.join(boost_dir, f'BoOST-stellarmodels-{grid_name}')
    dat_files = glob(os.path.join(grid_dir, 'f*.dat'))

    best_file = None
    best_mdist = 1e10
    best_mass = None
    for f in dat_files:
        m = re.match(r'f(\d+)-', os.path.basename(f))
        if m:
            M = int(m.group(1))
            dist = abs(M - M_target)
            if dist < best_mdist:
                best_mdist = dist
                best_file = f
                best_mass = M

    # Parse — space-separated, no header line (header is in fdat-columns.txt)
    data = np.loadtxt(best_file)

    # Compute Z from initial surface abundances (first row)
    X0 = data[0, 26]   # sH1
    He3_0 = data[0, 28]  # sHe3
    He4_0 = data[0, 29]  # sHe4
    Z_actual = 1.0 - X0 - He3_0 - He4_0

    age_yr = data[:, 0]
    mass = data[:, 1]
    Teff = data[:, 2]
    logL = data[:, 3]
    R_Rsun = data[:, 4]
    log_Mdot = data[:, 5]

    L_bol = 10**logL * L_sun
    Q_ion, L_FUV, L_LW = compute_radiation_vectorized(Teff, L_bol)

    return {
        'source': 'BoOST',
        'Z': Z_actual, 'M_init': best_mass,
        'age_yr': age_yr,
        'mass': mass,
        'Teff': Teff,
        'logL': logL,
        'L_bol': L_bol,
        'R_Rsun': R_Rsun,
        'log_Mdot': log_Mdot,
        'Q_ion': Q_ion,
        'L_FUV': L_FUV,
        'L_LW': L_LW,
    }


# ═══════════════════════════════════════════════════════════════
# Comparison plots
# ═══════════════════════════════════════════════════════════════

def compare_tracks(parsec_dir, boost_dir, mass_targets, output_file):
    """
    For each mass, load PARSEC (Z=0.008) and BoOST (MW, Z~0.00876)
    and plot side-by-side.
    """
    fig = plt.figure(figsize=(18, 5 * len(mass_targets)))
    gs = GridSpec(len(mass_targets), 6, figure=fig, hspace=0.35, wspace=0.35)

    panel_labels = ['M(t) [Msun]', r'$\dot{M}$ [Msun/yr]',
                    r'$Q_{ion}$ [phot/s]', r'$L_{FUV}$ [erg/s]',
                    r'$L_{bol}$ [erg/s]', r'R [Rsun]']

    for row, M_target in enumerate(mass_targets):
        print(f"Loading M={M_target} Msun...", flush=True)
        p = load_parsec_track(parsec_dir, 0.008, M_target)
        b = load_boost_track(boost_dir, 'MW', M_target)

        print(f"  PARSEC: Z={p['Z']}, M_init={p['M_init']}, {len(p['age_yr'])} steps, "
              f"lifetime={p['age_yr'][-1]:.3e} yr")
        print(f"  BoOST:  Z={b['Z']:.5f}, M_init={b['M_init']}, {len(b['age_yr'])} steps, "
              f"lifetime={b['age_yr'][-1]:.3e} yr")

        # Normalize age to fraction of lifetime
        p_t = p['age_yr'] / 1e6  # Myr
        b_t = b['age_yr'] / 1e6

        datasets = [
            # (PARSEC y, BoOST y, ylabel, log_scale)
            (p['mass'], b['mass'], panel_labels[0], False),
            (10**p['log_Mdot'], 10**b['log_Mdot'], panel_labels[1], True),
            (p['Q_ion'], b['Q_ion'], panel_labels[2], True),
            (p['L_FUV'], b['L_FUV'], panel_labels[3], True),
            (p['L_bol'], b['L_bol'], panel_labels[4], True),
            (p['R_Rsun'], b['R_Rsun'], panel_labels[5], True),
        ]

        for col, (py, by, ylabel, logscale) in enumerate(datasets):
            ax = fig.add_subplot(gs[row, col])

            # Filter out zeros/negatives for log plots
            if logscale:
                p_mask = py > 0
                b_mask = by > 0
                ax.semilogy(p_t[p_mask], py[p_mask], 'b-', lw=1.2, alpha=0.8,
                           label=f'PARSEC {p["M_init"]:.0f}M')
                ax.semilogy(b_t[b_mask], by[b_mask], 'r-', lw=1.2, alpha=0.8,
                           label=f'BoOST {b["M_init"]}M')
            else:
                ax.plot(p_t, py, 'b-', lw=1.2, alpha=0.8,
                       label=f'PARSEC {p["M_init"]:.0f}M')
                ax.plot(b_t, by, 'r-', lw=1.2, alpha=0.8,
                       label=f'BoOST {b["M_init"]}M')

            ax.set_ylabel(ylabel, fontsize=9)
            ax.set_xlabel('Age [Myr]', fontsize=9)
            ax.tick_params(labelsize=8)
            if row == 0 and col == 0:
                ax.legend(fontsize=7, loc='best')
            elif col == 0:
                ax.legend(fontsize=7, loc='best')

            if row == 0:
                ax.set_title(ylabel, fontsize=10)

        # Add radiation pressure panel annotation
        ax = fig.add_subplot(gs[row, 4])
        # L/c momentum rate
        p_Lc = p['L_bol'] / c_light  # dyn = g cm/s^2
        b_Lc = b['L_bol'] / c_light
        p_Lc_Msun_km_s_Myr = p_Lc * (3.156e13) / (1.989e33 * 1e5)  # Msun km/s per Myr
        b_Lc_Msun_km_s_Myr = b_Lc * (3.156e13) / (1.989e33 * 1e5)

    # Overall title
    fig.suptitle(f'PARSEC (Z={p["Z"]}) vs BoOST MW (Z={b["Z"]:.5f}) — Massive Stars',
                 fontsize=14, y=1.01)

    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {output_file}")


def print_summary_table(parsec_dir, boost_dir, mass_targets):
    """Print a text summary comparing key quantities at ZAMS and end of life."""
    print("\n" + "="*100)
    print(f"{'M_init':>6s} | {'Source':>7s} | {'Z':>8s} | {'lifetime':>10s} | "
          f"{'M_final':>8s} | {'Mdot_max':>10s} | {'Q_ion_max':>10s} | "
          f"{'L_bol_max':>10s} | {'L/c_max':>10s}")
    print(f"{'':>6s} | {'':>7s} | {'':>8s} | {'[Myr]':>10s} | "
          f"{'[Msun]':>8s} | {'[Msun/yr]':>10s} | {'[phot/s]':>10s} | "
          f"{'[erg/s]':>10s} | {'[Msun km/s/Myr]':>10s}")
    print("-"*100)

    for M_target in mass_targets:
        p = load_parsec_track(parsec_dir, 0.008, M_target)
        b = load_boost_track(boost_dir, 'MW', M_target)

        for track in [p, b]:
            t = track
            lifetime_Myr = t['age_yr'][-1] / 1e6
            M_final = t['mass'][-1]
            Mdot_max = 10**np.max(t['log_Mdot'])
            Q_ion_max = np.max(t['Q_ion'])
            L_bol_max = np.max(t['L_bol'])
            # L/c in Msun km/s per Myr
            Lc_max = L_bol_max / c_light * 3.156e13 / (1.989e33 * 1e5)

            print(f"{t['M_init']:6.0f} | {t['source']:>7s} | {t['Z']:8.5f} | "
                  f"{lifetime_Myr:10.3f} | {M_final:8.2f} | {Mdot_max:10.2e} | "
                  f"{Q_ion_max:10.2e} | {L_bol_max:10.2e} | {Lc_max:10.2f}")
        print("-"*100)


if __name__ == '__main__':
    parsec_dir = '/raven/u/uli/phil/stellar_tables/parsec1.2'
    boost_dir  = '/raven/u/uli/phil/stellar_tables/boost'

    # Compare at masses where both have tracks (BoOST MW: 25, 40, 60, 80, 120, 250, 500)
    # PARSEC closest: 24, 40, 60, 80, 120, 250, 350
    mass_targets = [25, 40, 60, 80, 120, 250]

    print("="*60)
    print("PARSEC vs BoOST comparison at MW metallicity")
    print("="*60)

    # Summary table
    print_summary_table(parsec_dir, boost_dir, mass_targets)

    # Plots
    output_file = '/raven/u/uli/phil/stellar_tables/parsec_vs_boost.png'
    compare_tracks(parsec_dir, boost_dir, mass_targets, output_file)
