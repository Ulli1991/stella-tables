"""
Standalone Python mirror of the C++ stellar-table interpolation in
galaxy_sf/resolvedism_stellar_tables.cc.

Bilinear in (logM, logZ) for 2D fields; trilinear (logM, logZ, log_age) for 3D.
Remnant type is nearest-corner in (logM, logZ) — matches the C++ convention
for integer fields.

Used by test_invariants.py to check closure / monotonicity / consistency
invariants under interpolation, not just at grid points.
"""

import h5py
import os
import numpy as np


REM_NAMES = {0: "WD", 1: "ECSN", 2: "CCSN", 3: "FSN",
             4: "PPISN", 5: "PISN", 6: "DBH"}
REM_WD, REM_ECSN, REM_CCSN, REM_FSN, REM_PPISN, REM_PISN, REM_DBH = range(7)
SN_TYPES = {REM_ECSN, REM_CCSN, REM_FSN, REM_PPISN, REM_PISN, REM_DBH}
NO_EJECTA_TYPES = {REM_FSN, REM_DBH}


class StellarTable:
    def __init__(self, path_unified, path_wind=None):
        with h5py.File(path_unified, "r") as f:
            self.logM       = f["log_M_init"][...]
            self.logZ       = f["log_Z"][...]
            self.log_age    = f["log_age_yr"][...]
            self.M_init     = f["M_init"][...]
            self.Z          = f["Z"][...]
            self.M_preSN    = f["M_preSN"][...]
            self.rem_mass   = f["remnant_mass"][...]
            self.rem_type   = f["remnant_type"][...]
            self.elem_ej_SN  = f["elem_ej_SN_mass"][...]
            self.elem_ej_AGB = f["elem_ej_AGB_mass"][...]
            self.X_init     = f["X_init"][...]
            self.M_current  = f["M_current"][...]
            self.lifetime   = f["lifetime_yr"][...]
            self.elements   = [e.decode() for e in f["elements"][...]]
        self.NZ, self.NM = self.M_preSN.shape
        self.NAGE        = len(self.log_age)
        self.NELEM       = len(self.elements)
        if path_wind is None:
            path_wind = path_unified.replace(
                "stellar_tables_unified", "stellar_wind_cumulative")
        try:
            with h5py.File(path_wind, "r") as f:
                # try a few common dataset names
                for k in ("elem_ej_wind_cumulative",
                          "cum_wind_per_elem",
                          "cumulative_wind_ejecta"):
                    if k in f:
                        self.cum_wind = f[k][...]
                        break
                else:
                    self.cum_wind = None
        except FileNotFoundError:
            self.cum_wind = None

    # ── factors for (logM, logZ) ───────────────────────────────────────
    def _bl_factors(self, lm, lz):
        im = int(np.clip(np.searchsorted(self.logM, lm) - 1, 0, self.NM - 2))
        iz = int(np.clip(np.searchsorted(self.logZ, lz) - 1, 0, self.NZ - 2))
        fm = (lm - self.logM[im]) / (self.logM[im + 1] - self.logM[im])
        fz = (lz - self.logZ[iz]) / (self.logZ[iz + 1] - self.logZ[iz])
        fm = max(0.0, min(1.0, fm))
        fz = max(0.0, min(1.0, fz))
        return im, iz, fm, fz

    def _bilinear(self, field2d, lm, lz):
        im, iz, fm, fz = self._bl_factors(lm, lz)
        return ((1 - fm) * (1 - fz) * field2d[iz,     im]
                + fm   * (1 - fz) * field2d[iz,     im + 1]
                + (1 - fm) * fz   * field2d[iz + 1, im]
                + fm   * fz       * field2d[iz + 1, im + 1])

    def _nearest_corner(self, field2d, lm, lz):
        im, iz, fm, fz = self._bl_factors(lm, lz)
        ii = im if fm < 0.5 else im + 1
        jj = iz if fz < 0.5 else iz + 1
        return field2d[jj, ii]

    def _bilinear_rt_aware(self, field2d, lm, lz):
        """Bilinear if all 4 corners share rem_type, else nearest-corner.
        Mirrors interp2d_rt_aware in the C++ runtime."""
        im, iz, fm, fz = self._bl_factors(lm, lz)
        rt00 = int(self.rem_type[iz,     im])
        rt01 = int(self.rem_type[iz,     im + 1])
        rt10 = int(self.rem_type[iz + 1, im])
        rt11 = int(self.rem_type[iz + 1, im + 1])
        if rt00 == rt01 == rt10 == rt11:
            return ((1 - fm) * (1 - fz) * field2d[iz,     im]
                    + fm   * (1 - fz) * field2d[iz,     im + 1]
                    + (1 - fm) * fz   * field2d[iz + 1, im]
                    + fm   * fz       * field2d[iz + 1, im + 1])
        ii = im if fm < 0.5 else im + 1
        jj = iz if fz < 0.5 else iz + 1
        return field2d[jj, ii]

    # ── factors including log_age ──────────────────────────────────────
    def _tl_factors(self, lm, lz, la):
        im, iz, fm, fz = self._bl_factors(lm, lz)
        ia = int(np.clip(np.searchsorted(self.log_age, la) - 1, 0, self.NAGE - 2))
        fa = (la - self.log_age[ia]) / (self.log_age[ia + 1] - self.log_age[ia])
        fa = max(0.0, min(1.0, fa))
        return im, iz, ia, fm, fz, fa

    def _trilinear(self, field3d, lm, lz, la):
        im, iz, ia, fm, fz, fa = self._tl_factors(lm, lz, la)
        c000 = field3d[iz,     im,     ia]
        c100 = field3d[iz,     im + 1, ia]
        c010 = field3d[iz + 1, im,     ia]
        c110 = field3d[iz + 1, im + 1, ia]
        c001 = field3d[iz,     im,     ia + 1]
        c101 = field3d[iz,     im + 1, ia + 1]
        c011 = field3d[iz + 1, im,     ia + 1]
        c111 = field3d[iz + 1, im + 1, ia + 1]
        return ((1 - fm) * (1 - fz) * (1 - fa) * c000
                + fm   * (1 - fz) * (1 - fa) * c100
                + (1 - fm) * fz   * (1 - fa) * c010
                + fm   * fz       * (1 - fa) * c110
                + (1 - fm) * (1 - fz) * fa   * c001
                + fm   * (1 - fz) * fa       * c101
                + (1 - fm) * fz   * fa       * c011
                + fm   * fz       * fa       * c111)

    # ── named accessors (mirror C++ API) ───────────────────────────────
    def remnant_mass(self, lm, lz):
        return self._bilinear_rt_aware(self.rem_mass, lm, lz)

    def remnant_type(self, lm, lz):
        return int(self._nearest_corner(self.rem_type, lm, lz))

    def m_preSN(self, lm, lz):
        return self._bilinear_rt_aware(self.M_preSN, lm, lz)

    def elem_SN(self, lm, lz, k):
        return self._bilinear_rt_aware(self.elem_ej_SN[:, :, k], lm, lz)

    def elem_SN_total(self, lm, lz):
        return self._bilinear_rt_aware(self.elem_ej_SN.sum(axis=2), lm, lz)

    def elem_AGB(self, lm, lz, k):
        return self._bilinear_rt_aware(self.elem_ej_AGB[:, :, k], lm, lz)

    def elem_AGB_total(self, lm, lz):
        return self._bilinear_rt_aware(self.elem_ej_AGB.sum(axis=2), lm, lz)

    def lifetime_yr(self, lm, lz):
        return self._bilinear(self.lifetime, lm, lz)

    def m_current(self, lm, lz, la):
        return self._trilinear(self.M_current, lm, lz, la)

    def cum_wind_per_elem(self, lm, lz, la, k):
        if self.cum_wind is None:
            return 0.0
        return self._trilinear(self.cum_wind[:, :, :, k], lm, lz, la)

    def cum_wind_total(self, lm, lz, la):
        if self.cum_wind is None:
            return 0.0
        return sum(self.cum_wind_per_elem(lm, lz, la, k)
                   for k in range(self.NELEM))


def load_default():
    """tables from $STELLA_TABLE / $STELLA_WIND, else the repo's build/ outputs"""
    build = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'build')
    return StellarTable(
        os.environ.get('STELLA_TABLE', os.path.join(build, 'stellar_tables_unified.hdf5')),
        os.environ.get('STELLA_WIND', os.path.join(build, 'stellar_wind_cumulative.hdf5')),
    )
