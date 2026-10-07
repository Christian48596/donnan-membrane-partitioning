"""Sulfate association with explicit species and thermodynamic cycle checks.

Free-species activities use Truesdell-Jones parameters from Sujanani et al.
(2023), DOI 10.1021/acs.jpcb.2c07900.s001, Section 9/Table S9. Parent-salt
Pitzer coefficients are NOT reused as free-species coefficients.

The default analysis uses their conductivity-derived association estimates,
linearly interpolated at the requested analytical concentration. These estimates
are not direct measurements of chemical species. K_a inferred from each state
is not asserted to be a salinity-independent thermodynamic constant.

Membrane activities extend the binary closure using SI equations S39-S40,
S58-S59 and S71-S73. Neutral pair gamma=1, h=1 are explicit assumptions.
The neutral-pair unit-partition sensitivity is distinguished from the
self-consistent sodium-pair calculation and from the published model's closure.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq, least_squares
from .properties import solution_property
from .thermodynamics import association_constant, mobile_ionic_strength

CHARGES = {"Na+": 1, "Mg2+": 2, "SO4-2": -2, "NaSO4-": -1, "MgSO4": 0}
TJ = {"Na+": (4.32, .06), "Mg2+": (5.46, .33),
      "SO4-2": (5.31, -.07), "NaSO4-": (5.40, .0)}


def truesdell_jones(composition):
    """Concentration-based free-species coefficients at 298.15 K.

    log10 gamma_i = -0.5092 z_i^2 sqrt(I)/(1 + 0.3283 a_i sqrt(I)) + b_i I.
    a_i is in Angstrom and I in mol/L. This empirical comparator is not claimed
    to be a high-salinity Pitzer replacement; its limitations remain explicit.
    """
    I = mobile_ionic_strength(composition, CHARGES)
    sq = np.sqrt(I)
    out = {}
    for n in composition:
        if CHARGES[n] == 0:
            out[n] = 1.
        else:
            a, b = TJ[n]
            out[n] = 10**(-.5092*CHARGES[n]**2*sq/(1+.3283*a*sq) + b*I)
    return out


def external_species(salt, concentration_M, alpha):
    """Element and charge balances with alpha = paired fraction of sulfate."""
    if concentration_M <= 0 or not 0 <= alpha < 1:
        raise ValueError("positive analytical concentration and 0 <= alpha < 1 required")
    c = concentration_M
    if salt == "MgSO4":
        return {"Mg2+": (1-alpha)*c, "SO4-2": (1-alpha)*c, "MgSO4": alpha*c}
    if salt == "Na2SO4":
        return {"Na+": (2-alpha)*c, "SO4-2": (1-alpha)*c, "NaSO4-": alpha*c}
    raise ValueError("association supports MgSO4 and Na2SO4 only")


def external_association(salt, concentration_M, K_a, activity_model=truesdell_jones,
                         c_standard_M=1.):
    """Solve the correct salt-dependent mass action with explicit free activities.

    Na2SO4 contains (2-alpha)c free sodium, not (1-alpha)c. An activity model
    must return coefficients for the retained species on the concentration scale.
    K_a is dimensionless on that same standard state; no undocumented default.
    """
    if K_a <= 0 or c_standard_M <= 0:
        raise ValueError("positive association constant and standard concentration required")
    M, P = ("Mg2+", "MgSO4") if salt == "MgSO4" else ("Na+", "NaSO4-")
    def f(a):
        c = external_species(salt, concentration_M, a)
        g = activity_model(c)
        return association_constant(c[M], c["SO4-2"], c[P], g[M], g["SO4-2"], g[P],
                                    c_standard_M) - K_a
    alpha = float(brentq(f, 0., 1-1e-10, xtol=1e-13))
    c = external_species(salt, concentration_M, alpha)
    return alpha, c, activity_model(c)


def manning_species(composition, fixed_charge_M, xi, counterion):
    """Generalised activity closure with one counter-ion and any co-ion species.

    All physical counter-ions remain in charge balance. The condensed population
    changes the counter-ion activity prefactor, not the material fixed charge.
    D_M = C_A + xi sum_co |z| (Z_g+|z|) c_co (concentration units here).
    The resulting log coefficients are identical to SI S39-S40.
    """
    zg = CHARGES[counterion]
    if zg*xi <= 1 or fixed_charge_M < 0:
        raise ValueError("requires condensation regime and nonnegative fixed charge")
    cg = composition[counterion]
    coions = [(n,c) for n,c in composition.items() if CHARGES[n] < 0]
    D = fixed_charge_M + xi*sum(abs(CHARGES[n])*(zg+abs(CHARGES[n]))*c
                              for n,c in coions)
    A = 0. if D == 0 else fixed_charge_M/(2*D)
    pref = (cg-fixed_charge_M/zg+fixed_charge_M/(xi*zg**2))/cg
    if pref <= 0:
        raise ValueError("counter-ion population incompatible with fixed charge")
    return {n: (pref*np.exp(-A) if n == counterion else
                np.exp(-A*(CHARGES[n]/zg)**2)) for n in composition}


def paired_partition(salt, concentration_M, fixed_charge_M, xi, *, alpha=None,
                     K_a=None, h_pair=1., gamma_pair_membrane=1.):
    """Conditional, species-consistent sulfate partitioning.

    With measured alpha, infer the external K_a at each state and propagate
    K_m/K_s=h_pair (h_M=h_S=1). With supplied K_a, solve external speciation
    instead. Sodium pairs share one Donnan potential and charge balance; neutral
    MgSO4 uptake is h_pair*c_P_solution/gamma_pair_membrane (gamma_P_s=1).
    """
    if alpha is not None and K_a is not None:
        raise ValueError("choose measured alpha or prescribed K_a, not both")
    if h_pair <= 0 or gamma_pair_membrane <= 0:
        raise ValueError("positive pair transfer factor/activity coefficient required")
    if K_a is None:
        alpha = solution_property(salt, concentration_M, "association_fraction") if alpha is None else alpha
        cs = external_species(salt, concentration_M, alpha)
        gs = truesdell_jones(cs)
        source = "conductivity-derived alpha; statewise inferred K_a"
    else:
        alpha, cs, gs = external_association(salt, concentration_M, K_a)
        source = "mass action with caller-supplied dimensionless K_a"
    M, P = ("Mg2+", "MgSO4") if salt == "MgSO4" else ("Na+", "NaSO4-")
    a = {n: gs[n]*c for n,c in cs.items()}  # c_standard=1 M
    Ka_s = association_constant(cs[M], cs["SO4-2"], cs[P], gs[M], gs["SO4-2"], gs[P])

    if salt == "MgSO4":
        pair = h_pair*a[P]/gamma_pair_membrane
        def composition(S):
            return {M: fixed_charge_M/2+S, "SO4-2": S, P: pair}
        def balance(log_S):
            cm = composition(np.exp(log_S))
            gm = manning_species(cm, fixed_charge_M, xi, M)
            return np.log(gm[M]*cm[M]*gm["SO4-2"]*cm["SO4-2"]/(a[M]*a["SO4-2"]))
        S = np.exp(brentq(balance, -100., np.log(max(10., concentration_M, fixed_charge_M)), xtol=1e-12))
        cm = composition(S)
    else:
        def composition(u):
            S, pair = np.exp(u)
            return {M: fixed_charge_M+2*S+pair, "SO4-2": S, P: pair}
        def residual(u):
            cm = composition(u)
            gm = manning_species(cm, fixed_charge_M, xi, M)
            gm[P] *= gamma_pair_membrane  # documented multiplicative pair correction
            am = {n:gm[n]*c for n,c in cm.items()}
            return [np.log(am[M]**2*am["SO4-2"]/(a[M]**2*a["SO4-2"])),
                    np.log(am[M]*am[P]/(h_pair*a[M]*a[P]))]
        if a[P] == 0:
            # A vanishing pair reduces to the fully dissociated binary salt.
            from .manning import solve_mobile_salt
            gpm = (gs[M]**2*gs["SO4-2"])**(1/3)
            S = solve_mobile_salt("Na2SO4", fixed_charge_M, concentration_M, gpm, xi)
            cm = {M:fixed_charge_M+2*S, "SO4-2":S, P:0.}
        else:
            fit = least_squares(residual, np.log([max(cs["SO4-2"]*.01,1e-12),
                                                 max(cs[P]*.1,1e-12)]),
                                bounds=(-100., 10.), xtol=1e-12, ftol=1e-12, gtol=1e-12)
            if not fit.success or max(abs(np.array(residual(fit.x)))) > 1e-8:
                raise RuntimeError("paired species equilibrium did not converge")
            cm = composition(fit.x)
    gm = manning_species(cm, fixed_charge_M, xi, M)
    if salt == "MgSO4": gm[P] = gamma_pair_membrane
    else: gm[P] *= gamma_pair_membrane
    am = {n:gm[n]*c for n,c in cm.items()}
    psi = -np.log(am[M]/a[M])/CHARGES[M]
    cycle_Km = association_constant(cm[M], cm["SO4-2"], cm[P], gm[M], gm["SO4-2"], gm[P])
    charge_error = sum(CHARGES[n]*c for n,c in cm.items())-fixed_charge_M
    return dict(alpha=float(alpha), free=cm["SO4-2"], pair=cm[P], total=cm["SO4-2"]+cm[P],
                solution=cs, membrane=cm, gamma_solution=gs, gamma_membrane=gm,
                psi_D=float(psi), K_a_solution=Ka_s, K_a_membrane=cycle_Km,
                cycle_error=cycle_Km/Ka_s/h_pair-1 if Ka_s else 0.,
                charge_error_M=charge_error, source=source)


def published_equal_association(salt, concentration_M, fixed_charge_M, xi, alpha):
    """Published SI S57/S70 closure, explicitly imposing alpha_m=alpha_s.

    This is a comparator for the equal-association case, not a reproduction of
    fitted membrane fractions or digitised published sulfate predictions. It
    does not in general enforce the individual-species equilibrium cycle.
    """
    cs = external_species(salt, concentration_M, alpha)
    gs = truesdell_jones(cs)
    M = "Mg2+" if salt == "MgSO4" else "Na+"
    stoich = {n:c/concentration_M for n,c in cs.items()}
    target = sum(stoich[n]*np.log(gs[n]*cs[n]) for n in cs if cs[n] > 0)
    def residual(log_C):
        C = np.exp(log_C)
        cm = {n:v*C for n,v in stoich.items()}
        cm[M] += fixed_charge_M/CHARGES[M]
        gm = manning_species(cm,fixed_charge_M,xi,M)
        return sum(stoich[n]*np.log(gm[n]*cm[n]) for n in cm if stoich[n] > 0)-target
    return float(np.exp(brentq(residual,-100.,5.,xtol=1e-12)))
