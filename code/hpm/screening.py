"""Dilute Debye reference and ion-size crossover for aqueous electrolytes.

These functions provide a molecular-scale diagnostic rather than a measurement
of screening in a heterogeneous membrane.
"""

from __future__ import annotations

import numpy as np

from .constants import E, EPS0, KB, NA, T_REF


def debye_length(I_molar, eps_r: float, T: float = T_REF):
    """Classical Debye length (m) for ionic strength I in mol/L.

    lambda_D = sqrt(eps0 eps_r kB T / (2 e^2 NA * 1000 * I))

    This point-ion expression is a dilute homogeneous reference. It does not
    determine a measured membrane correlation length.
    """
    I = np.asarray(I_molar, dtype=float) * 1e3 * NA   # ions per m^3 scale
    if np.any(I < 0) or np.any(~np.isfinite(I)) or eps_r <= 0 or T <= 0:
        raise ValueError("finite nonnegative I, positive permittivity and temperature required")
    with np.errstate(divide="ignore"):
        return np.sqrt(EPS0 * eps_r * KB * T / (2.0 * E**2 * I))


def crossover_concentration(eps_r: float, T: float = T_REF, a: float = 0.30e-9) -> float:
    """Ionic strength (mol/L) at which the dilute reference lambda_D equals ion diameter a."""
    if min(eps_r,T,a) <= 0:
        raise ValueError("positive permittivity, temperature and ion diameter required")
    I = EPS0 * eps_r * KB * T / (2.0 * E**2 * a**2)   # per m^3
    return I / (1e3 * NA)
