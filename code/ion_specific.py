"""Generate explicitly conditional lithium exchange illustrations."""
from pathlib import Path
import numpy as np
import pandas as pd
from hpm.exchange import (RT_KJ_MOL as RT, condensed_fraction, separation_factor,
                          ceiling, required_K, required_energy, solve_finite_bath)
RES=Path(__file__).resolve().parents[1]/"results"


def main():
    RES.mkdir(exist_ok=True)
    K=np.logspace(-1,5,300)
    rows=[]
    for xi in (1.83,3.1,5.,8.):
        for y in (.01,.1,.5):
            for coion_ratio in (0.,.5,1.):
                for k in K:
                    rows.append(dict(xi=xi,y=y,coion_to_fixed_charge=coion_ratio,K_ex=k,
                                     alpha=separation_factor(k,xi,y=y,coion_M=coion_ratio),
                                     alpha_asymptote=ceiling(xi,y=y,coion_M=coion_ratio),
                                     q=1.,Gamma_f=1.))
    pd.DataFrame(rows).to_csv(RES/"ion_specific_curves.csv",index=False)
    requirements=[]
    for xi in (1.83,3.1,5.,8.):
        for y in (.01,.1,.5):
            for a in (1.5,2.,3.,5.,10.):
                dg=required_energy(a,xi,y=y)
                requirements.append(dict(xi=xi,y=y,coion_to_fixed_charge=0.,alpha=a,
                                         K_ex=required_K(a,xi,y=y),dG_kJ_per_mol=dg,
                                         favourable_energy_magnitude_kJ_per_mol=None if dg is None else -dg,
                                         reachable=dg is not None))
    pd.DataFrame(requirements).to_csv(RES/"ion_specific_requirements.csv",index=False)
    # A finite-inventory illustration: initial total inventory is supplied, not
    # inferred from single-salt observations. Bath vol 1 L, membrane water .1 L.
    finite=[dict(K_ex=k,**solve_finite_bath(.02,1.38,1.,.1,4.,0.,3.1,k))
            for k in (1.,10.,100.,1e4,1e8)]
    pd.DataFrame(finite).to_csv(RES/"ion_specific_finite_bath.csv",index=False)
    text=["CONDITIONAL Li/Na EXCHANGE ILLUSTRATIONS (no lithium validation)",
          "q=Gamma_f=1, reservoir y is fixed, condensed capacity theta*C_A.",
          "alpha_inf=1+w/[y(1-w)], w=theta*C_A/(C_A+c_coion)."]
    for xi in (1.83,3.1):
        text.append(f"xi={xi}: equimolar/negligible-coion asymptote {ceiling(xi):.3g}; "
                    f"y=.01 asymptote {ceiling(xi,y=.01):.3g}; "
                    f"equimolar, c_coion/C_A=.5 asymptote {ceiling(xi,coion_M=.5):.3g}.")
    for a in (2.,3.,5.):
        text.append(f"Nafion xi=3.1, y=.5, no coions, target {a:g}: "
                    f"DeltaG_ex^standard={required_energy(a,3.1):.3f} kJ/mol.")
    text.append("These are hypothesis-dependent reservoir limits, not universal membrane ceilings.")
    (RES/"ion_specific_summary.txt").write_text("\n".join(text)+"\n")
    print("\n".join(text))

if __name__=="__main__": main()
