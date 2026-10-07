"""Fully dissociated Manning/Donnan equilibrium for ion-exchange membranes.

Concentrations are mol/L of membrane water or external solution. Fixed groups
are monovalent. The model uses positive counter-ion and co-ion charge magnitudes,
stoichiometric coefficients and the Manning parameter xi. Binary salt equations
are solved self-consistently for membrane mobile-salt concentration.

Sources: Kamcev et al. (2015), DOI 10.1021/acs.macromol.5b01654, SI Section S3;
Kamcev et al. (2016), DOI 10.1039/C5CP06747B. Model equations and assumptions
are documented in docs/MODEL_DETAILS.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

__all__ = ["SaltStoich", "STOICH", "XI_NAFION117", "XI_CR61",
           "gamma_counterion", "gamma_coion", "solve_mobile_salt",
           "coion_concentration"]

XI_NAFION117 = 3.1
XI_CR61 = 1.83


@dataclass(frozen=True)
class SaltStoich:
    """Counter-ion (g) and co-ion (c) valences and stoichiometry, for a CEM."""

    name: str
    z_g: int      # |counter-ion valence|
    nu_g: int     # counter-ion stoichiometric coefficient
    z_c: int      # |co-ion valence|
    nu_c: int     # co-ion stoichiometric coefficient
    z_p: int = 1  # fixed-charge group valence (sulfonate)

    @property
    def kind(self) -> str:
        return f"{self.z_g}:{self.z_c}"


STOICH: dict[str, SaltStoich] = {
    "NaCl":   SaltStoich("NaCl",   z_g=1, nu_g=1, z_c=1, nu_c=1),
    "LiCl":   SaltStoich("LiCl",   z_g=1, nu_g=1, z_c=1, nu_c=1),
    "MgCl2":  SaltStoich("MgCl2",  z_g=2, nu_g=1, z_c=1, nu_c=2),
    "CaCl2":  SaltStoich("CaCl2",  z_g=2, nu_g=1, z_c=1, nu_c=2),
    "Na2SO4": SaltStoich("Na2SO4", z_g=1, nu_g=2, z_c=2, nu_c=1),
    "MgSO4":  SaltStoich("MgSO4",  z_g=2, nu_g=1, z_c=2, nu_c=1),
}


def _denom(Xi: float, s: SaltStoich, xi: float) -> float:
    return Xi + s.z_g * s.z_c * xi * (s.nu_g + s.nu_c)


def _validate(Xi: float, s: SaltStoich, xi: float) -> None:
    if not np.isfinite(Xi) or Xi < 0:
        raise ValueError("Xi must be finite and nonnegative")
    if not np.isfinite(xi) or s.z_g * xi <= 1:
        raise ValueError("the condensation closure requires Z_g * xi > 1")
    if s.z_p != 1 or s.nu_g * s.z_g != s.nu_c * s.z_c:
        raise ValueError("requires monovalent fixed groups and neutral salt stoichiometry")


def gamma_counterion(Xi: float, s: SaltStoich, xi: float) -> float:
    """Counter-ion activity coefficient; Xi = C_A / C_mobile."""
    _validate(Xi, s, xi)
    pre = (Xi / (s.z_g * xi) + s.z_g * s.nu_g) / (Xi + s.z_g * s.nu_g)
    return float(pre * np.exp(-(Xi / 2.0) / _denom(Xi, s, xi)))


def gamma_coion(Xi: float, s: SaltStoich, xi: float) -> float:
    """Co-ion activity coefficient for a fully dissociated binary salt."""
    _validate(Xi, s, xi)
    return float(np.exp(-(Xi / 2.0) * (s.z_c / s.z_g) ** 2 / _denom(Xi, s, xi)))


def _closure(C_s_m: float, C_A_m: float, s: SaltStoich, xi: float,
             rhs: float) -> float:
    """LHS - RHS of the appropriate Manning/Donnan equation."""
    Xi = C_A_m / C_s_m
    gg = gamma_counterion(Xi, s, xi)
    gc = gamma_coion(Xi, s, xi)

    if (s.z_g, s.z_c) == (1, 1):                       # Eq. 2.23
        lhs = (C_A_m + C_s_m) * C_s_m * gg * gc
    elif (s.z_g, s.z_c) == (2, 1):                     # Eq. 2.24
        lhs = (C_A_m / 2.0 + C_s_m) * C_s_m**2 * gg * gc**2
    elif (s.z_g, s.z_c) == (1, 2):                     # Eq. 2.25
        lhs = (C_A_m + 2.0 * C_s_m) ** 2 * C_s_m * gg**2 * gc
    elif (s.z_g, s.z_c) == (2, 2):                     # Eq. 2.26
        lhs = (C_A_m / 2.0 + C_s_m) * C_s_m * gg * gc
    else:
        raise ValueError(f"unsupported stoichiometry {s.kind}")
    return lhs - rhs


def solve_mobile_salt(salt: str, C_A_m: float, C_s_s: float, gamma_pm_s: float,
                      xi: float = XI_NAFION117, *,
                      h_g: float = 1.0, h_c: float = 1.0) -> float:
    """Mobile salt concentration C_s^m in the membrane (mol/L sorbed water).

    C_A_m : fixed charge concentration, same basis
    C_s_s : external salt concentration (mol/L)
    gamma_pm_s : external mean activity coefficient on the CONCENTRATION scale;
                 do not pass a molal Pitzer coefficient without conversion.
    h_g, h_c : standard transfer factors (steric accessibility included once).
    """
    s = STOICH[salt]
    if not all(np.isfinite(v) for v in (C_A_m, C_s_s, gamma_pm_s, h_g, h_c)):
        raise ValueError("all model inputs must be finite")
    if C_A_m < 0 or C_s_s < 0 or min(gamma_pm_s, h_g, h_c) <= 0:
        raise ValueError("nonnegative concentrations and positive activities/transfer factors required")
    _validate(0.0, s, xi)
    if C_s_s == 0:
        return 0.0
    if (s.z_g, s.z_c) in {(1, 1), (2, 2)}:
        rhs = gamma_pm_s**2 * C_s_s**2
    elif (s.z_g, s.z_c) == (2, 1):
        rhs = gamma_pm_s**3 * C_s_s**3
    else:                                              # (1, 2)
        rhs = 4.0 * gamma_pm_s**3 * C_s_s**3
    rhs *= h_g**s.nu_g * h_c**s.nu_c

    f = lambda c: _closure(c, C_A_m, s, xi, rhs)
    lo, hi = np.finfo(float).tiny**0.2, max(1.0, C_A_m, C_s_s)
    for _ in range(30):
        if f(hi) >= 0:
            break
        hi *= 2
    else:
        raise RuntimeError(f"{salt}: Manning/Donnan root not bracketed")
    # Solve in log concentration so small sorption values retain relative accuracy.
    return float(np.exp(brentq(lambda u: f(np.exp(u)), np.log(lo), np.log(hi),
                              xtol=1e-12, rtol=1e-14)))


def coion_concentration(salt: str, C_s_m: float) -> float:
    """Membrane co-ion concentration from the mobile salt concentration.

    The stoichiometric factor matters: for MgCl2 the co-ion is 2 C_s^m, for
    Na2SO4 it is C_s^m. Getting this wrong silently shifts every comparison.
    """
    s = STOICH[salt]
    return s.nu_c * C_s_m
