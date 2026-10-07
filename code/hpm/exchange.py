"""Lithium/sodium exchange equilibrium in condensed and free-ion populations.

This is a design hypothesis. The sorption data contain no lithium. A condensed
capacity theta*C_A is assumed; extra co-ion uptake adds free counter-ions and
reduces the condensed fraction of the total population. No mixed-salt Pitzer
parameters or chemically selective capacity are inferred from the data.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq

RT_KJ_MOL = 8.31446261815324e-3*298.15


def condensed_fraction(xi, z_g=1):
    if xi <= 0 or z_g <= 0 or not np.isfinite(xi):
        raise ValueError("positive finite xi and positive counter-ion valence required")
    return max(0.,1-1/(z_g*xi))


def population_weight(xi, fixed_charge_M=1., coion_M=0.):
    if fixed_charge_M <= 0 or coion_M < 0:
        raise ValueError("positive fixed charge and nonnegative monovalent co-ion uptake required")
    return condensed_fraction(xi)*fixed_charge_M/(fixed_charge_M+coion_M)


def populations(K_ex, y, w, q=1., Gamma_f=1.):
    if not 0 <= y <= 1 or not 0 <= w < 1:
        raise ValueError("0 <= y <= 1 and 0 <= w < 1 required")
    if min(K_ex,q,Gamma_f) <= 0 or not np.isfinite(q*Gamma_f):
        raise ValueError("positive exchange constant and finite positive free preference required")
    if y in (0.,1.):
        return y,y,y
    f = q*y/(1-y+q*y)
    r = K_ex*Gamma_f
    x = 1. if np.isinf(r) else r*f/(1-f+r*f)
    return f,x,w*x+(1-w)*f


def separation_factor(K_ex, xi, z_g=1, y=.5, *, fixed_charge_M=1., coion_M=0.,
                      q=1., Gamma_f=1., w=None):
    if z_g != 1:
        raise ValueError("Li/Na are monovalent; use z_g=1")
    if not 0 < y < 1:
        raise ValueError("selectivity requires 0 < y < 1")
    w = population_weight(xi,fixed_charge_M,coion_M) if w is None else w
    f,x,p = populations(K_ex,y,w,q,Gamma_f)
    return (1-y)/y*p/(w*(1-x)+(1-w)*(1-f))


def ceiling(xi, z_g=1, y=.5, *, fixed_charge_M=1., coion_M=0., q=1., w=None):
    """Reservoir strong-exchange limit; equimolar 2*xi-1 is one special case."""
    return separation_factor(np.inf,xi,z_g,y,fixed_charge_M=fixed_charge_M,
                             coion_M=coion_M,q=q,w=w)


def required_K(alpha,xi,z_g=1,y=.5,*,fixed_charge_M=1.,coion_M=0.,q=1.,Gamma_f=1.,w=None):
    if alpha <= 0 or not 0 < y < 1 or z_g != 1:
        raise ValueError("positive target selectivity, interior y and monovalent ions required")
    w = population_weight(xi,fixed_charge_M,coion_M) if w is None else w
    f,_,_ = populations(1.,y,w,q,Gamma_f)
    if w == 0:
        return 1. if np.isclose(alpha,q) else None
    p = alpha*y/(1-y+alpha*y)
    x = (p-(1-w)*f)/w
    if not 0 < x < 1:
        return None
    return x/(1-x)*(1-f)/f/Gamma_f


def required_energy(alpha,xi,z_g=1,y=.5,**kw):
    K = required_K(alpha,xi,z_g,y,**kw)
    return None if K is None else -RT_KJ_MOL*np.log(K)


def solve_finite_bath(total_Li_mol,total_Na_mol,bath_volume_L,membrane_water_L,
                      fixed_charge_M,coion_M,xi,K_ex,*,q=1.,Gamma_f=1.):
    """Element conservation for a finite bath with prescribed membrane capacity.

    Inventories include BOTH bath and membrane ions initially. Co-ion uptake and
    total membrane counter-ion concentration are supplied; their equilibrium
    calculation is a separate closure. y returned is the FINAL bath composition.
    """
    if min(total_Li_mol,total_Na_mol,bath_volume_L,membrane_water_L,fixed_charge_M) <= 0:
        raise ValueError("positive inventories, volumes and charge required")
    if not np.isfinite(K_ex):
        raise ValueError("finite bath requires a finite exchange constant; its infinite limit is singular")
    w=population_weight(xi,fixed_charge_M,coion_M)
    n_mem=(fixed_charge_M+coion_M)*membrane_water_L
    n_bath=total_Li_mol+total_Na_mol-n_mem
    if n_bath <= 0:
        raise ValueError("inventory does not fill the prescribed membrane and a positive bath")
    def balance(y):
        _,_,p=populations(K_ex,y,w,q,Gamma_f)
        return y*n_bath+p*n_mem-total_Li_mol
    # Strong exchange can leave y extremely small. An absolute 1e-14 tolerance
    # on y amplified by the binding response loses element conservation.
    y=float(brentq(balance,0.,1.,xtol=np.finfo(float).tiny,rtol=1e-13))
    f,x,p=populations(K_ex,y,w,q,Gamma_f)
    return dict(y_final=y,f=f,x=x,w=w,
                Li_bath_mol=y*n_bath,Na_bath_mol=(1-y)*n_bath,
                Li_membrane_mol=p*n_mem,Na_membrane_mol=(1-p)*n_mem,
                bath_cations_M=n_bath/bath_volume_L,
                alpha=separation_factor(K_ex,xi,y=y,w=w,q=q,Gamma_f=Gamma_f))
