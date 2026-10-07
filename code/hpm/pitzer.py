"""Pitzer mean ionic activity coefficients for single electrolytes at 298.15 K.

The model evaluates the molal ionic strength and the Pitzer electrostatic,
binary-interaction and third-order terms. gamma_pm returns a molal coefficient;
gamma_pm_molar converts its argument and coefficient for concentration-based
equilibrium equations. Numerical parameter values and conventions are documented
in docs/MODEL_DETAILS.md.

Sources: Pitzer (1973), Pitzer and Mayorga (1973, 1974).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["PitzerSalt", "SALTS", "ln_gamma_pm", "gamma_pm", "mean_activity_product",
           "molarity_warning", "molality_from_molarity", "gamma_pm_molar",
           "molar_activity_product", "WATER_DENSITY_KG_L", "MOLAR_MASS_KG_MOL"]

A_PHI = 0.3915      # kg^1/2 mol^-1/2, water at 25 C
B_PITZER = 1.2
WATER_DENSITY_KG_L = 0.9970474  # water at 298.15 K, near 1 bar
MOLAR_MASS_KG_MOL = {"NaCl": 0.05844277, "MgCl2": 0.095211,
                     "Na2SO4": 0.14204214, "MgSO4": 0.120366,
                     "LiCl": 0.042394, "CaCl2": 0.110984}


@dataclass(frozen=True)
class PitzerSalt:
    """Single-electrolyte Pitzer parameters. beta2 is used only for 2-2 salts."""

    name: str
    z_M: int
    z_X: int
    nu_M: int
    nu_X: int
    beta0: float
    beta1: float
    C_phi: float
    beta2: float = 0.0

    @property
    def nu(self) -> int:
        return self.nu_M + self.nu_X

    @property
    def is_2_2(self) -> bool:
        return abs(self.z_M) == 2 and abs(self.z_X) == 2

    def ionic_strength(self, m: float) -> float:
        return 0.5 * (self.nu_M * m * self.z_M**2 + self.nu_X * m * self.z_X**2)


# Pitzer & Mayorga (1973) parameters at 25 C.
# Numerical values used by the calculations are listed in docs/MODEL_DETAILS.md.
SALTS: dict[str, PitzerSalt] = {
    "NaCl":   PitzerSalt("NaCl",   1, -1, 1, 1, 0.0765,  0.2664, 0.00127),
    "MgCl2":  PitzerSalt("MgCl2",  2, -1, 1, 2, 0.35235, 1.6815, 0.00519),
    "Na2SO4": PitzerSalt("Na2SO4", 1, -2, 2, 1, 0.01958, 1.113,  0.00497),
    "MgSO4":  PitzerSalt("MgSO4",  2, -2, 1, 1, 0.2210,  3.343,  0.0250,
                         beta2=-37.23),
    "LiCl":   PitzerSalt("LiCl",   1, -1, 1, 1, 0.1494,  0.3074, 0.00359),
    "CaCl2":  PitzerSalt("CaCl2",  2, -1, 1, 2, 0.3159,  1.614,  -0.00034),
}


def _B_gamma(salt: PitzerSalt, I: float) -> float:
    sq = np.sqrt(I)

    def term(beta, alpha):
        if beta == 0.0:
            return 0.0
        x = alpha * sq
        return (2.0 * beta / (alpha**2 * I)) * (
            1.0 - (1.0 + x - 0.5 * x**2) * np.exp(-x)
        )

    if salt.is_2_2:
        return 2.0 * salt.beta0 + term(salt.beta1, 1.4) + term(salt.beta2, 12.0)
    return 2.0 * salt.beta0 + term(salt.beta1, 2.0)


def ln_gamma_pm(salt: PitzerSalt | str, m: float) -> float:
    """Natural log of the mean ionic activity coefficient at molality m."""
    s = SALTS[salt] if isinstance(salt, str) else salt
    if m < 0 or not np.isfinite(m):
        raise ValueError("molality must be finite and nonnegative")
    if m == 0:
        return 0.0
    I = s.ionic_strength(m)
    sq = np.sqrt(I)

    f_gamma = -A_PHI * (sq / (1.0 + B_PITZER * sq)
                        + (2.0 / B_PITZER) * np.log(1.0 + B_PITZER * sq))

    C_gamma = 1.5 * s.C_phi
    t1 = abs(s.z_M * s.z_X) * f_gamma
    t2 = m * (2.0 * s.nu_M * s.nu_X / s.nu) * _B_gamma(s, I)
    t3 = m**2 * (2.0 * (s.nu_M * s.nu_X) ** 1.5 / s.nu) * C_gamma
    return float(t1 + t2 + t3)


def gamma_pm(salt: PitzerSalt | str, m: float) -> float:
    return float(np.exp(ln_gamma_pm(salt, m)))


def mean_activity_product(salt: PitzerSalt | str, m: float) -> float:
    """(gamma_pm)^nu * (nu_M m)^nu_M * (nu_X m)^nu_X.

    MOLAL salt activity product; use molar_activity_product for mol/L closures.
    """
    s = SALTS[salt] if isinstance(salt, str) else salt
    g = gamma_pm(s, m)
    return (g ** s.nu) * (s.nu_M * m) ** s.nu_M * (s.nu_X * m) ** s.nu_X


def molality_from_molarity(salt: str, concentration_M: float, density_kg_L: float):
    """m = c/(rho - c M), with solvent mass in kg per litre of solution."""
    solvent_mass = density_kg_L - concentration_M * MOLAR_MASS_KG_MOL[salt]
    if (not np.isfinite(concentration_M) or concentration_M < 0 or
            not np.isfinite(solvent_mass) or solvent_mass <= 0):
        raise ValueError("nonnegative molarity and positive solvent mass required")
    return concentration_M / solvent_mass


def gamma_pm_molar(salt: str, concentration_M: float, density_kg_L: float):
    """Mean coefficient normalised to gamma_c -> 1 at infinite dilution.

    gamma_c = gamma_m * rho_water/(rho_solution - c M).
    The rho_water factor is the standard-state conversion; changing the Pitzer
    concentration argument alone does not convert the activity coefficient.
    """
    solvent_mass = density_kg_L - concentration_M * MOLAR_MASS_KG_MOL[salt]
    m = molality_from_molarity(salt, concentration_M, density_kg_L)
    return gamma_pm(salt, m) * WATER_DENSITY_KG_L / solvent_mass


def molar_activity_product(salt: str, concentration_M: float, density_kg_L: float):
    """Concentration-scale neutral salt product for a mol/L Donnan closure."""
    s = SALTS[salt]
    g = gamma_pm_molar(salt, concentration_M, density_kg_L)
    return g**s.nu * (s.nu_M * concentration_M)**s.nu_M * (s.nu_X * concentration_M)**s.nu_X


def molarity_warning(m: float) -> str | None:
    if m > 1.0:
        return (f"molality/molarity distinction ignored at {m} mol/L; "
                "convert using solution density for quantitative work")
    return None
