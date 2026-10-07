"""Identifiability of residual mechanisms, not a fit to raw partition coefficients.

R_i = beta_i - z_i delta_psi1 + z_i^2/r_i B + z_i^2/2 psi2
      + z_i^2/Z_g^2 eta_M (direct co-ion Manning activity error).

Use I = I_solution as the shared covariate. Sharing response functions across
salts is a diagnostic assumption. Residuals conditional on charge estimated
from the same ion contents require joint uncertainty propagation. Offset
columns belong to ion identities, and a
reference salinity is fixed so constant effects are not counted twice.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy import stats


@dataclass(frozen=True)
class Measurement:
    salt: str
    ion: str
    z: int
    r_nm: float  # illustrative effective Born radius, not a measured Born cavity
    z_g: int

    def __post_init__(self):
        if self.z == 0 or self.r_nm <= 0 or self.z_g <= 0:
            raise ValueError("requires charged ions, positive radii and counter-ion valence")


R = {"Li+": .076, "Na+": .102, "Mg2+": .072, "Ca2+": .100,
     "Cl-": .181, "SO4-2": .230}


def _m(salt, ion, z, z_g):
    return Measurement(salt, ion, z, R[ion], z_g)


AVAILABLE_COIONS = (_m("NaCl", "Cl-", -1, 1), _m("MgCl2", "Cl-", -1, 2),
                    _m("Na2SO4", "SO4-2", -2, 1))
AVAILABLE = (_m("NaCl", "Na+", 1, 1), _m("NaCl", "Cl-", -1, 1),
             _m("MgCl2", "Mg2+", 2, 2), _m("MgCl2", "Cl-", -1, 2),
             _m("Na2SO4", "Na+", 1, 1), _m("Na2SO4", "SO4-2", -2, 1))
PLUS_LI_CA = AVAILABLE + (_m("LiCl", "Li+", 1, 1), _m("LiCl", "Cl-", -1, 1),
                         _m("CaCl2", "Ca2+", 2, 2), _m("CaCl2", "Cl-", -1, 2))
PLUS_MGSO4 = PLUS_LI_CA + (_m("MgSO4", "Mg2+", 2, 2), _m("MgSO4", "SO4-2", -2, 2))
PANELS = {"prospective independent counter-ion and co-ion data": AVAILABLE,
          "prospective + LiCl, CaCl2": PLUS_LI_CA,
          "prospective + MgSO4": PLUS_MGSO4}


def response_block(panel, include_psi2=True, include_mu=True):
    z = np.array([m.z for m in panel], float)
    r = np.array([m.r_nm for m in panel], float)
    zg = np.array([m.z_g for m in panel], float)
    # Born column is scaled by 1 nm; its fitted coefficient is B/(1 nm).
    cols = [-z, z*z/r]
    if include_psi2:
        cols.append(.5*z*z)
    if include_mu:
        # The derived error column is for co-ions. The counter-ion prefactor
        # has a different error structure; do not quietly apply this to it.
        cols.append(np.where(z < 0, z*z/zg**2, 0.))
    return np.column_stack(cols)


def _scaled(X):
    norms = np.linalg.norm(X, axis=0)
    return X / np.where(norms > 0, norms, 1.)


def _rank(X):
    return int(np.linalg.matrix_rank(_scaled(X)))


def _basis(X):
    X = _scaled(X)
    u, s, _ = np.linalg.svd(X, full_matrices=False)
    tol = np.finfo(float).eps * max(X.shape) * (s[0] if len(s) else 1.)
    return u[:, s > tol]


def design_matrix(panel, n_I, include_psi2, include_mu, *, basis=None):
    """Ion offsets plus salinity responses in a finite, reference-fixed basis.

    Default: independent response changes at n_I-1 external ionic strengths,
    relative to the first salinity. A supplied basis must also vanish there.
    All mechanisms use the same salinity basis for the collinearity test.
    """
    if n_I < 2 or not panel:
        raise ValueError("at least two ionic strengths and one measurement required")
    b = np.eye(n_I)[:, 1:] if basis is None else np.asarray(basis, float)
    if b.ndim != 2 or b.shape[0] != n_I or not np.allclose(b[0], 0):
        raise ValueError("salinity basis must have n_I rows and zero reference row")
    ions = tuple(dict.fromkeys(m.ion for m in panel))
    offsets = np.array([[float(m.ion == ion) for ion in ions] for m in panel])
    responses = response_block(panel, include_psi2, include_mu)
    return np.column_stack([np.tile(offsets, (n_I, 1)), np.kron(b, responses)])


def check_identifiability(panel, n_I=5):
    without = response_block(panel, True, False)
    with_error = response_block(panel, True, True)
    X0 = design_matrix(panel, n_I, False, True)
    X1 = design_matrix(panel, n_I, True, True)
    def condition(X):
        return float(np.linalg.cond(_scaled(X))) if _rank(X) == X.shape[1] else np.inf
    contrast = any(len({m.z_g for m in panel if m.ion == ion}) > 1
                   for ion in {m.ion for m in panel})
    return dict(n_measurements=len(panel), rank_without_mu=_rank(without),
                rank_with_mu=_rank(with_error), cond_without_mu=condition(without),
                cond_with_mu=condition(with_error),
                counter_ion_valence_contrast=contrast,
                full_design_rank=_rank(X1), full_design_columns=X1.shape[1],
                variance_added_rank=_rank(X1)-_rank(X0),
                psi2_identifiable_given_baseline_error=_rank(X1) > _rank(X0))


def simulate_power(panel, n_I=5, I_min=.1, I_max=6., sigma=.10,
                   amplitude=1., mu_amplitude=0., alpha=.05, n_trials=600,
                   seed=0, include_mu=True):
    """Prospective nested-model F-test, with explicit external-I and noise model.

    Undefined tests return NaN, not a measured zero power. A rank-deficient
    nuisance model is projected onto its actual span before testing variance.
    Synthetic counter-ion measurements, when included, are independent by
    assumption; these are not the charge-balanced experimental concentrations.
    """
    if sigma <= 0 or I_min <= 0 or I_max <= I_min or n_trials <= 0:
        raise ValueError("positive noise, ordered positive I range and trials required")
    rng = np.random.default_rng(seed)
    I = np.geomspace(I_min, I_max, n_I)
    X0 = design_matrix(panel, n_I, False, include_mu)
    X1 = design_matrix(panel, n_I, True, include_mu)
    Q0, Q1 = _basis(X0), _basis(X1)
    df_diff, df_res = Q1.shape[1]-Q0.shape[1], X1.shape[0]-Q1.shape[1]
    info = dict(n_obs=X1.shape[0], df_diff=df_diff, df_residual=df_res,
                n_I=n_I, sigma=sigma, I_min_external_M=I_min,
                I_max_external_M=I_max, n_trials=n_trials, seed=seed)
    if df_diff <= 0 or df_res <= 0:
        return dict(info, power=np.nan, note="variance test undefined: no extra identifiable direction")
    z = np.array([m.z for m in panel], float)
    mu_weights = response_block(panel, False, True)[:, -1]
    ramp = np.maximum(0., np.log10(I))
    truth = (amplitude*ramp[:, None]*.5*z*z +
             mu_amplitude*ramp[:, None]*mu_weights).ravel()
    y = truth[:, None] + rng.normal(0., sigma, (len(truth), n_trials))
    rss0 = np.sum((y-Q0 @ (Q0.T @ y))**2, axis=0)
    rss1 = np.sum((y-Q1 @ (Q1.T @ y))**2, axis=0)
    F = np.maximum(0., rss0-rss1)/df_diff / (rss1/df_res)
    return dict(info, power=float(np.mean(F > stats.f.ppf(1-alpha, df_diff, df_res))),
                note="synthetic residuals with independent Gaussian noise; not data validation")
