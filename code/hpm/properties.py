"""Traceable external solution properties; no silent density extrapolation."""
from pathlib import Path
import numpy as np
import pandas as pd
from .pitzer import gamma_pm_molar

PATH = Path(__file__).resolve().parents[2] / "data" / "solution_properties.csv"


def solution_property(salt, concentration_M, column="density_kg_L", *, strict=False):
    df = pd.read_csv(PATH, comment="#")
    p = df[df.salt == salt].dropna(subset=[column]).sort_values("c_external_M")
    if p.empty or not p.c_external_M.min() <= concentration_M <= p.c_external_M.max():
        raise ValueError(f"no {column} data for {salt} at {concentration_M} M")
    if strict and p.provenance.str.contains("approximation").any():
        raise ValueError(f"{salt}: density inputs are from a constant-volume approximation rather than measured at each concentration")
    return float(np.interp(concentration_M, p.c_external_M, p[column]))


def solution_gamma(salt, concentration_M, *, strict=False):
    rho = solution_property(salt, concentration_M, strict=strict)
    return gamma_pm_molar(salt, concentration_M, rho)
