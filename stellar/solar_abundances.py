#!/usr/bin/env python3
"""Solar composition used by the tables and by GIZMO: one source.

Element ratios: Asplund, Grevesse, Sauval & Scott 2009, ARA&A 47, 481 (arXiv:0909.0948), Table 1,
present-day photosphere (He, Ne, Ar: the indirect values given there). Atomic weights: IUPAC standard values.
Normalization (Stella convention): metals C..Zn sum to Z_SUN = 0.014, Y = A09 photospheric 0.2485, X = 1 - Y - Z.

  solar_abundances.py            print the 27-element pattern and A09's own X, Y, Z
  solar_abundances.py --c        print the C array for GIZMO core/init.cc
"""
import sys

ELEMENTS = ['H', 'He', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
            'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
# A09 Table 1, photosphere, log eps (log N_H = 12); Li, Be, B included for the normalization only
LOGEPS_A09 = {
    'H': 12.00, 'He': 10.93, 'Li': 1.05, 'Be': 1.38, 'B': 2.70, 'C': 8.43, 'N': 7.83, 'O': 8.69, 'F': 4.56,
    'Ne': 7.93, 'Na': 6.24, 'Mg': 7.60, 'Al': 6.45, 'Si': 7.51, 'P': 5.41, 'S': 7.12, 'Cl': 5.50, 'Ar': 6.40,
    'K': 5.03, 'Ca': 6.34, 'Sc': 3.15, 'Ti': 4.95, 'V': 3.93, 'Cr': 5.64, 'Mn': 5.43, 'Fe': 7.50, 'Co': 4.99,
    'Ni': 6.22, 'Cu': 4.19, 'Zn': 4.56,
}
ATOMIC_WEIGHT = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999,
    'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
    'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867,
    'V': 50.942, 'Cr': 51.996, 'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546,
    'Zn': 65.38,
}
Z_SUN = 0.014          # Stella / GIZMO solar metallicity (All.SolarAbundances[0])
Y_SUN = 0.2485         # A09 helioseismic surface Y_s (Sect. 3.9); log eps(He) = 10.93 is rounded


def a09_mass_fractions():
    """A09 present-day photospheric mass fractions over H..Zn (elements beyond Zn add < 1e-7)"""
    m = {e: 10 ** (LOGEPS_A09[e] - 12) * ATOMIC_WEIGHT[e] for e in LOGEPS_A09}
    tot = sum(m.values())
    return {e: v / tot for e, v in m.items()}


def solar_pattern():
    """27-element solar mass fractions in ELEMENTS order: A09 metal ratios, metals sum to Z_SUN"""
    x = a09_mass_fractions()
    metals = ELEMENTS[2:]; zsum = sum(x[e] for e in metals)
    pat = {e: x[e] / zsum * Z_SUN for e in metals}
    pat['He'] = Y_SUN; pat['H'] = 1.0 - Y_SUN - Z_SUN
    return pat


def a09_xyz():
    x = a09_mass_fractions()
    return x['H'], x['He'], 1.0 - x['H'] - x['He']


if __name__ == '__main__':
    p = solar_pattern()
    if '--c' in sys.argv:
        print('            static const double solar_elem[27] = {  /* stella-tables stellar/solar_abundances.py */')
        for i in range(0, 27, 2):
            print('                ' + '    '.join('%.10e, /* %2d: %-2s */' % (p[ELEMENTS[k]], k, ELEMENTS[k]) for k in (i, i + 1) if k < 27))
        print('            };')
    else:
        X, Y, Z = a09_xyz(); print('A09 photospheric: X=%.4f Y=%.4f Z=%.4f' % (X, Y, Z))
        for e in ELEMENTS: print('%-2s %.4e' % (e, p[e]))
        print('sum metals %.6f, sum all %.12f' % (sum(p[e] for e in ELEMENTS[2:]), sum(p.values())))
