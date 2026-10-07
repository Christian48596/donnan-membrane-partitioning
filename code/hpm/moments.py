"""Closed-form moments of the oscillatory potential profile.

This supplies the quantity the whole experimental programme hinges on: the
variance of the potential across the pore, which by the cumulant expansion is the
effect size of the new electrostatic term.

    psi(s) = psi_s exp(-s/lam) cos(k0 s + theta),      s = r_p - r

Cross-sectional average in a cylinder of radius r_p:

    <f> = (2/r_p^2) INT_0^{r_p} f(r) r dr
        = (2/r_p^2) INT_0^{r_p} f(s) (r_p - s) ds

and the two moments needed are

    M1 = <psi>,        M2 = <psi^2>,       sigma_psi^2 = M2 - M1^2.

Both reduce to integrals of the form INT s^n exp(-a s) cos(b s + phi) ds, which
are elementary. The closed forms below were derived with sympy and are checked
against direct quadrature in the test suite to 1e-10.

The second cumulant is a nonnegative candidate contribution under uniform
activity prefactors and a common accessible domain. The profile below is an
illustration, not a self-consistent membrane-field prediction. Higher cumulants
matter; no separate configurational entropy correction is added.
"""

from __future__ import annotations

import numpy as np

__all__ = ["moment1", "moment2", "variance", "predicted_amplitude",
           "amplitude_profile"]


def _I0(a: float, b: float, phi: float, L: float) -> float:
    """INT_0^L exp(-a s) cos(b s + phi) ds."""
    d = a * a + b * b
    def F(s):
        return np.exp(-a * s) * (-a * np.cos(b * s + phi) + b * np.sin(b * s + phi)) / d
    return F(L) - F(0.0)


def _I1(a: float, b: float, phi: float, L: float) -> float:
    """INT_0^L s exp(-a s) cos(b s + phi) ds."""
    d = a * a + b * b
    def F(s):
        c, sn = np.cos(b * s + phi), np.sin(b * s + phi)
        # antiderivative of s e^{-as} cos(bs+phi)
        term1 = s * (-a * c + b * sn) / d
        term2 = ((b * b - a * a) * c + 2 * a * b * sn) / (d * d)
        return np.exp(-a * s) * (term1 + term2)
    return F(L) - F(0.0)


def moment1(psi_s: float, lam: float, k0: float, r_p: float,
            theta: float = 0.0) -> float:
    """<psi> over the pore cross-section."""
    a, b, L = 1.0 / lam, k0, r_p
    return (2.0 * psi_s / r_p**2) * (r_p * _I0(a, b, theta, L) - _I1(a, b, theta, L))


def moment2(psi_s: float, lam: float, k0: float, r_p: float,
            theta: float = 0.0) -> float:
    """<psi^2>, using cos^2 x = (1 + cos 2x)/2."""
    L = r_p
    a1, a2 = 2.0 / lam, 2.0 / lam
    flat = r_p * _I0(a1, 0.0, 0.0, L) - _I1(a1, 0.0, 0.0, L)
    osc = r_p * _I0(a2, 2.0 * k0, 2.0 * theta, L) - _I1(a2, 2.0 * k0, 2.0 * theta, L)
    return (psi_s**2 / r_p**2) * (flat + osc)


def variance(psi_s: float, lam: float, k0: float, r_p: float,
             theta: float = 0.0) -> float:
    """sigma_psi^2 = <psi^2> - <psi>^2. This is the effect amplitude."""
    m1 = moment1(psi_s, lam, k0, r_p, theta)
    m2 = moment2(psi_s, lam, k0, r_p, theta)
    return max(m2 - m1 * m1, 0.0)


def predicted_amplitude(I_molar, psi_s: float = 3.0, r_p_nm: float = 1.5,
                        d_nm: float = 0.65, a_ion_nm: float = 0.30,
                        theta: float = 0.0, eps_r: float = 78.3, T: float = 298.15,
                        n_exp: float = 3.0):
    """Illustrative variance using an assumed screening/profile ansatz.

    Not called by the manuscript pipeline. The amplitude and oscillatory
    profile are prescribed; these do not determine membrane electrostatics.
    """
    from .screening import crossover_concentration, effective_screening_length

    I = np.atleast_1d(np.asarray(I_molar, dtype=float))
    a_ion = a_ion_nm * 1e-9
    I_star = crossover_concentration(eps_r, T, a=a_ion)
    lam = np.atleast_1d(effective_screening_length(I, eps_r, T, a=a_ion,
                                                   exponent=n_exp))
    k0_full = 2.0 * np.pi / (d_nm * 1e-9)
    r_p = r_p_nm * 1e-9

    out = np.empty_like(I)
    for j, (Ij, lj) in enumerate(zip(I, lam)):
        # oscillation switches on above the crossover
        k0 = k0_full if Ij > I_star else 0.0
        out[j] = variance(psi_s, float(lj), k0, r_p, theta)
    return out if out.size > 1 else float(out[0])


def amplitude_profile(**kw):
    """Convenience: (I grid, sigma_psi^2) over the brine-relevant range."""
    I = np.logspace(-2, np.log10(8.0), 60)
    return I, predicted_amplitude(I, **kw)


# --------------------------------------------------------------------------- #
# Exact electrostatic factor
# --------------------------------------------------------------------------- #

def exact_ln_K(z: int, psi_s: float, lam: float, k0: float, r_p: float,
               theta: float = 0.0) -> float:
    """ln <exp(-z psi)> by direct quadrature, without truncation.

    Quadrature is performed on a dimensionless domain so the absolute tolerance
    is meaningful for nanometre pores. The charge-squared approximation is a
    cumulant truncation, not an exact scaling or a guaranteed charge asymmetry.
    """
    from scipy.integrate import quad

    if lam <= 0 or r_p <= 0:
        raise ValueError("positive decay length and pore radius required")
    def exponent(t):
        s=r_p*t
        return -z*psi_s*np.exp(-s/lam)*np.cos(k0*s+theta)
    shift=float(np.max(exponent(np.linspace(0.,1.,2001))))
    val, _ = quad(lambda t: np.exp(exponent(t)-shift)*(1-t), 0., 1.,
                  epsabs=1e-12, epsrel=1e-11, limit=400)
    return float(shift+np.log(2.*val))


def charge_asymmetry(z: int, psi_s: float, lam: float, k0: float, r_p: float,
                     theta: float = 0.0) -> float:
    """|ln K(+z) - ln K(-z)|: the counter-ion / co-ion asymmetry.

    No universal small-asymmetry claim is made for an oscillatory profile;
    its value depends on the full potential distribution, including odd cumulants.
    """
    return abs(exact_ln_K(z, psi_s, lam, k0, r_p, theta)
               - exact_ln_K(-z, psi_s, lam, k0, r_p, theta))
