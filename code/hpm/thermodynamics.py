"""Species partitioning, charge balance and association equilibrium.

Concentrations are mol/L of solution or sorbed water; activity coefficients
use the same concentration standard c_standard_M. Fixed negative charge is a
positive concentration of elementary charges. The solver counts every charged
species, including pairs, and never reduces the material charge implicitly.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class Species:
    name: str
    z: int
    gamma_solution: float = 1.0
    gamma_membrane: float = 1.0
    h: float = 1.0

    def __post_init__(self):
        if any(not np.isfinite(x) or x <= 0 for x in
               (self.gamma_solution, self.gamma_membrane, self.h)):
            raise ValueError("activity coefficients and transfer factors must be positive")

    @property
    def log_prefactor(self):
        return np.log(self.h) + np.log(self.gamma_solution / self.gamma_membrane)


def mobile_ionic_strength(composition, charges):
    """One half sum(z_i**2 c_i); exclude immobile groups and neutral species."""
    if any(not np.isfinite(c) or c < 0 for c in composition.values()):
        raise ValueError("species concentrations must be finite and nonnegative")
    return 0.5 * sum(charges[n]**2 * c for n, c in composition.items())


def solve_donnan(composition, species, fixed_charge_M):
    """Return psi_D and membrane concentrations for specified species activities.

    Coefficients supplied here are held fixed during this solve. Composition-
    dependent activity closures require an outer self-consistent solve.
    """
    if fixed_charge_M < 0 or not np.isfinite(fixed_charge_M):
        raise ValueError("fixed negative charge must be finite and nonnegative")
    names = tuple(composition)
    if not names or set(names) != set(species):
        raise ValueError("provide one Species for every external concentration")
    c = np.array([composition[n] for n in names], float)
    if np.any(c < 0) or not np.all(np.isfinite(c)):
        raise ValueError("concentrations must be finite and nonnegative")
    z = np.array([species[n].z for n in names], float)
    if not np.any((z > 0) & (c > 0)) or not np.any((z < 0) & (c > 0)):
        raise ValueError("charged reservoir requires positive and negative species")
    if abs(z @ c) > 1e-10 * max(float(np.sum(np.abs(z) * c)), 1e-12):
        raise ValueError("external composition is not electroneutral")
    pref = np.array([species[n].log_prefactor for n in names])

    def concentrations(psi):
        return c * np.exp(np.clip(pref - z * psi, -700, 700))

    def balance(psi):
        return float(z @ concentrations(psi) - fixed_charge_M)

    span = 4.0
    while span <= 512 and not (balance(-span) >= 0 >= balance(span)):
        span *= 2
    if span > 512:
        raise RuntimeError("Donnan potential not bracketed")
    psi = float(brentq(balance, -span, span, xtol=1e-13))
    return psi, dict(zip(names, concentrations(psi)))


def association_constant(c_M, c_S, c_P, gamma_M=1., gamma_S=1., gamma_P=1.,
                         c_standard_M=1.):
    """Dimensionless K_a = a_P/(a_M a_S), for species on one activity scale."""
    if min(c_M, c_S, gamma_M, gamma_S, gamma_P, c_standard_M) <= 0 or c_P < 0:
        raise ValueError("positive reactants, coefficients and standard concentration required")
    return gamma_P / (gamma_M * gamma_S) * c_P * c_standard_M / (c_M * c_S)


def membrane_association_constant(K_solution, h_M=1., h_S=1., h_P=1.):
    """Equilibrium cycle: K_m/K_s = h_P/(h_M h_S)."""
    if min(K_solution, h_M, h_S, h_P) <= 0:
        raise ValueError("association constant and transfer factors must be positive")
    return K_solution * h_P / (h_M * h_S)
