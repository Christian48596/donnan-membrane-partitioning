"""Physical constants and unit helpers (SI throughout, except where noted)."""

from __future__ import annotations

from scipy import constants as _c

E = _c.elementary_charge          # C
KB = _c.Boltzmann                 # J/K
EPS0 = _c.epsilon_0               # F/m
NA = _c.Avogadro                  # 1/mol
R_GAS = _c.gas_constant           # J/(mol K)
T_REF = 298.15                    # K

NM = 1e-9                         # m per nm
LITRE = 1e-3                      # m^3 per L


def kT(T: float = T_REF) -> float:
    """Thermal energy in J."""
    return KB * T


def molar_to_number_density(c_molar: float) -> float:
    """mol/L -> particles per m^3."""
    return c_molar * 1e3 * NA


def bjerrum_length(eps_r: float, T: float = T_REF) -> float:
    """Bjerrum length in m."""
    return E**2 / (4.0 * _c.pi * EPS0 * eps_r * KB * T)
