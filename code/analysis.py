"""Calculate equilibrium sorption, association and identifiability diagnostics.

Sorption measurements determine the charge input. Approximate
properties and conditional speciation assumptions are reported in the outputs.
No lithium measurements are contained in the input data.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from hpm.manning import STOICH, XI_CR61, XI_NAFION117, coion_concentration, solve_mobile_salt
from hpm.pitzer import SALTS, gamma_pm, molality_from_molarity, molar_activity_product
from hpm.properties import solution_gamma, solution_property
from hpm.baseline_error import Measurement, R, check_identifiability, simulate_power
from hpm.association import paired_partition, published_equal_association, CHARGES
from hpm.screening import crossover_concentration

ROOT=Path(__file__).resolve().parents[1]
DATA,RES=ROOT/"data",ROOT/"results"
XI={"Nafion117":XI_NAFION117,"CR-61":XI_CR61}
MGSO4_GAMMA_UNCERTAINTY=.15  # illustrative sensitivity, not a confidence interval


def sorption_data():
    return pd.read_csv(DATA/"sorption_published.csv",comment="#")


def ideal_donnan_coion(salt,c_s,C_A):
    """Ideal membrane activities with concentration-scale Pitzer outside."""
    if C_A < 0 or c_s <= 0:
        raise ValueError("nonnegative fixed negative charge and positive feed required")
    s=SALTS[salt]
    rhs=molar_activity_product(salt,c_s,solution_property(salt,c_s))
    def residual(u):
        co=np.exp(u); ct=(C_A+abs(s.z_X)*co)/s.z_M
        return s.nu_M*np.log(ct)+s.nu_X*np.log(co)-np.log(rhs)
    return float(np.exp(brentq(residual,-100.,5.,xtol=1e-12)))


def manning_coion(salt,c_s,C_A,xi,gamma_scale=1.):
    return coion_concentration(salt,solve_mobile_salt(salt,C_A,c_s,
                                solution_gamma(salt,c_s)*gamma_scale,xi))


def interior_ionic_strength():
    """Mobile-ion concentration strength, with association bounds for sulfate.

    Condensed counter-ions are already represented in Manning activities and
    remain physical ions. No effective-charge or immobile-group term is added.
    Sulfate totals give a fully dissociated upper estimate, not exact speciation.
    """
    rows=[]
    for r in sorption_data().itertuples():
        s=STOICH[r.salt_system]; CA,co=r.C_fixed_M,r.c_coion_M
        ct=(CA+s.z_c*co)/s.z_g
        upper=.5*(s.z_g**2*ct+s.z_c**2*co)
        lower=CA if r.salt_system=="MgSO4" else (.5*CA+co if r.salt_system=="Na2SO4" else upper)
        rows.append(dict(membrane=r.membrane,salt=r.salt_system,c_ext=r.c_external_M,
                         C_fixed=CA,counterion_total_M=ct,z_counter=s.z_g,
                         I_mobile_dissociated_M=upper,I_mobile_min_association_M=lower,
                         sulfate_speciation_unknown=r.salt_system in ("Na2SO4","MgSO4")))
    return pd.DataFrame(rows)


def partitioning():
    ii=interior_ionic_strength()
    rows=[]
    for j,r in enumerate(sorption_data().itertuples()):
        salt=r.salt_system; xi=XI[r.membrane]
        ideal=ideal_donnan_coion(salt,r.c_external_M,r.C_fixed_M)
        man=manning_coion(salt,r.c_external_M,r.C_fixed_M,xi)
        lo=hi=np.nan
        if salt=="MgSO4":
            lo=np.log(r.c_coion_M/manning_coion(salt,r.c_external_M,r.C_fixed_M,xi,1.15))
            hi=np.log(r.c_coion_M/manning_coion(salt,r.c_external_M,r.C_fixed_M,xi,.85))
        rho=solution_property(salt,r.c_external_M)
        rows.append(dict(membrane=r.membrane,salt=salt,c_ext=r.c_external_M,C_fixed=r.C_fixed_M,
                         I_int=ii.iloc[j].I_mobile_dissociated_M,
                         I_ext=SALTS[salt].ionic_strength(r.c_external_M),
                         density_kg_L=rho,molality=molality_from_molarity(salt,r.c_external_M,rho),
                         gamma_external_molar=solution_gamma(salt,r.c_external_M),
                         density_status="constant_apparent_volume_approximation" if salt=="MgCl2" else "Sujanani2023_S2_interpolation",
                         co_obs=r.c_coion_M,co_sd=r.c_coion_sd,coion=r.coion,z_co=abs(r.z_coion),
                         z_counter=STOICH[salt].z_g,co_ideal=ideal,co_manning=man,
                         r_ideal=np.log(r.c_coion_M/ideal),r_manning=np.log(r.c_coion_M/man),
                         r_manning_lo=lo,r_manning_hi=hi))
    return pd.DataFrame(rows)


def validation():
    pub=pd.read_csv(DATA/"published_model_output.csv",comment="#")
    rows=[]
    association_predictions={"Nafion117":[.000200,.00161,.0144,.0867,.464],
             "CR-61":[.000149,.00122,.0113,.0703,.423]}
    counts={m:0 for m in association_predictions}
    for r in pub.itertuples():
        ours=solve_mobile_salt("NaCl",r.C_fixed_M,r.c_external_M,
                              solution_gamma("NaCl",r.c_external_M),XI[r.membrane])
        raw=.092 if r.membrane=="CR-61" and r.c_external_M==.1 else r.C_mobile_manning_published
        rev=association_predictions[r.membrane][counts[r.membrane]]; counts[r.membrane]+=1
        rows.append(dict(membrane=r.membrane,c_ext=r.c_external_M,experiment=r.C_mobile_experiment,
                         published_model=r.C_mobile_manning_published,published_raw=raw,
                         apparent_table_typo_corrected=raw!=r.C_mobile_manning_published,
                         published_association_model=rev,our_model=ours,
                         ratio=ours/r.C_mobile_manning_published,
                         charge_input_status="approximate_3M" if r.membrane=="CR-61" else "charge_balance"))
    return pd.DataFrame(rows)


def _panel(salts,counterions=False):
    out=[]
    for salt in salts:
        s=STOICH[salt]
        ct="Na+" if s.z_g==1 else "Mg2+"
        co="SO4-2" if s.z_c==2 else "Cl-"
        if counterions: out.append(Measurement(salt,ct,s.z_g,R[ct],s.z_g))
        out.append(Measurement(salt,co,-s.z_c,R[co],s.z_g))
    return tuple(out)


def identifiability(df=None):
    panels=[("Nafion coions: 3 salts",("NaCl","MgCl2","Na2SO4"),False),
            ("Nafion coions: + MgSO4",("NaCl","MgCl2","Na2SO4","MgSO4"),False),
            ("CR-61 coions: 2 sulfates",("Na2SO4","MgSO4"),False),
            ("Prospective both signs: 3 salts",("NaCl","MgCl2","Na2SO4"),True),
            ("Prospective both signs: + MgSO4",("NaCl","MgCl2","Na2SO4","MgSO4"),True)]
    rows=[]
    for name,salts,prospective in panels:
        p=_panel(salts,prospective); d=check_identifiability(p)
        exact=simulate_power(p,include_mu=False,n_trials=1500)
        uncertain=simulate_power(p,include_mu=True,mu_amplitude=.5,n_trials=1500)
        rows.append(dict(panel=name,prospective_independent_counterions=prospective,**d,
                         power_baseline_exact=exact["power"],power_baseline_uncertain=uncertain["power"],
                         variance_df_exact=exact["df_diff"],variance_df_uncertain=uncertain["df_diff"],
                         residual_df_uncertain=uncertain["df_residual"],
                         I_min_external_M=.1,I_max_external_M=6.,sigma_log_residual=.1,
                         n_I=5,n_trials=1500,seed=0,
                         note="hypothetical matched external-I grid; not observed-data power"))
    return pd.DataFrame(rows)


def association_fraction(salt,C_t,*,K_a=None):
    """Default conductivity-derived fraction; caller can supply an explicit K_a."""
    if K_a is None: return solution_property(salt,C_t,"association_fraction")
    from hpm.association import external_association
    return external_association(salt,C_t,K_a)[0]


def paired_model_coion(salt,c_s,C_A,xi,**kw):
    return paired_partition(salt,c_s,C_A,xi,**kw)


def pairing_check():
    rows=[]
    for r in sorption_data().itertuples():
        if r.salt_system not in ("MgSO4","Na2SO4"): continue
        p=paired_model_coion(r.salt_system,r.c_external_M,r.C_fixed_M,XI[r.membrane])
        man=manning_coion(r.salt_system,r.c_external_M,r.C_fixed_M,XI[r.membrane])
        pub=published_equal_association(r.salt_system,r.c_external_M,r.C_fixed_M,XI[r.membrane],p["alpha"])
        rows.append(dict(membrane=r.membrane,salt=r.salt_system,c_ext=r.c_external_M,alpha=p["alpha"],
                         obs=r.c_coion_M,manning=man,paired=p["total"],free_sulfate=p["free"],
                         membrane_pair=p["pair"],frac_from_pair=p["pair"]/p["total"],
                         r_manning=np.log(r.c_coion_M/man),r_paired=np.log(r.c_coion_M/p["total"]),
                         published_closure_equal_alpha=pub,r_published_equal_alpha=np.log(r.c_coion_M/pub),
                         K_a_solution_statewise=p["K_a_solution"],K_a_membrane_statewise=p["K_a_membrane"],
                         cycle_error=p["cycle_error"],charge_error_M=p["charge_error_M"],psi_D=p["psi_D"],
                         I_external_free_M=.5*sum(CHARGES[n]**2*c for n,c in p["solution"].items()),
                         model_class="neutral_pair_unit_partition_sensitivity" if r.salt_system=="MgSO4"
                                     else "conditional_charged_pair_species_equilibrium",
                         association_input="Sujanani2023_S2_conductivity_interpolation",
                         published_comparator_assumption="SI S57/S70 with alpha_m=alpha_s; not fitted cases"))
    return pd.DataFrame(rows)


def activity_basis_check(df):
    out=[]
    for r in df.itertuples():
        s=SALTS[r.salt]
        old=gamma_pm(r.salt,r.c_ext)
        argument_only=gamma_pm(r.salt,r.molality)
        out.append(dict(salt=r.salt,c_ext=r.c_ext,density_kg_L=r.density_kg_L,
                        density_status=r.density_status,
                        log_salt_product_shift_argument_only=s.nu*np.log(argument_only/old),
                        log_salt_product_shift_full_conversion=s.nu*np.log(r.gamma_external_molar/old)))
    return pd.DataFrame(out).drop_duplicates()


def main():
    RES.mkdir(parents=True,exist_ok=True)
    part,val,ii=partitioning(),validation(),interior_ionic_strength()
    ident,pair=identifiability(),pairing_check()
    tables={"partitioning":part,"validation":val,"interior_ionic_strength":ii,
            "identifiability":ident,"ion_pairing":pair,"activity_basis":activity_basis_check(part)}
    for name,df in tables.items(): df.to_csv(RES/f"{name}.csv",index=False)
    group=part.groupby(["membrane","salt"],sort=False).agg(
        ideal_mean=("r_ideal","mean"),manning_mean=("r_manning","mean"))
    pair_group=pair.groupby(["membrane","salt"],sort=False).agg(
        paired_mean=("r_paired","mean"),published_equal_alpha_mean=("r_published_equal_alpha","mean"),
        alpha_min=("alpha","min"),alpha_max=("alpha","max"))
    summary=group.join(pair_group)
    summary.to_csv(RES/"series_summary.csv")
    scalars=dict(n_measurements=len(part),rms_log_ideal=float(np.sqrt(np.mean(part.r_ideal**2))),
                 rms_log_manning=float(np.sqrt(np.mean(part.r_manning**2))),
                 largest_ideal_underprediction=float(np.exp(part.r_ideal.max())),
                 mobile_I_min_M=float(ii.I_mobile_dissociated_M.min()),
                 mobile_I_max_M=float(ii.I_mobile_dissociated_M.max()),
                 mobile_I_association_bound_min_M=float(ii.I_mobile_min_association_M.min()),
                 crossover_M=crossover_concentration(78.3),
                 max_charge_error_M=float(pair.charge_error_M.abs().max()),
                 max_association_cycle_relative_error=float(pair.cycle_error.abs().max()))
    (RES/"metrics.json").write_text(json.dumps(scalars,indent=2)+"\n")
    text=["EQUILIBRIUM MEMBRANE PARTITIONING", "="*72,
          f"30 sorption measurements. Fixed charge {part.C_fixed.min():.2f}-{part.C_fixed.max():.2f} M.",
          f"RMS log residual: ideal {scalars['rms_log_ideal']:.4f}; Manning {scalars['rms_log_manning']:.4f}.",
          f"Largest ideal underprediction: {scalars['largest_ideal_underprediction']:.1f}.",
          "\nSeries residuals (conditional association and equal-association published closure):",
          summary.to_string(float_format=lambda x:f"{x:+.4f}"),
          "\nMobile ionic strength excludes immobile polymer groups:",
          f"Fully dissociated reconstruction {scalars['mobile_I_min_M']:.3f}-{scalars['mobile_I_max_M']:.3f} M.",
          f"Conservative sulfate-association lower bound {scalars['mobile_I_association_bound_min_M']:.3f} M.",
          f"I*={scalars['crossover_M']:.5f} M is a dilute-screening diagnostic, not a decay measurement.",
          "\nThe MgCl2 density uses a constant apparent molar volume (19.77 mL/mol),",
          "anchored to Zezin 2014 SI Table 1S; other densities use Sujanani 2023 SI S2.",
          "Pitzer parameter values are defined in hpm/pitzer.py; MgSO4 +/-15%",
          "is an illustrative input sensitivity, not a statistically calibrated interval.",
          "\nIdentifiability: real co-ion measurements and hypothetical independent counter-ions:",
          ident[["panel","variance_df_uncertain","power_baseline_exact","power_baseline_uncertain"]].to_string(index=False),
          "NaN denotes an undefined variance test. The source papers report counter-ion contents.",
          "Residuals conditional on charge derived from those contents require joint uncertainty",
          "propagation rather than independent-noise assumptions. The simulation assumes matched",
          "external I grids, common response functions and independent Gaussian residual noise.",
          "\nAssociation uses free-species Truesdell-Jones activities, not parent-salt Pitzer.",
          "Default alpha is conductivity-derived; statewise K_a is not a verified universal constant.",
          "Neutral MgSO4 gamma=1 and unit transfer are assumed. Charged NaSO4- shares one",
          "Donnan potential and charge balance. Both phase equilibrium cycles are checked.",
          "Published comparator: SI S57/S70 with equal alpha in both phases; this does not",
          "reproduce fitted membrane fractions or constitute validation against sulfate model tables.",
          "Mixed residual signs show that a positive variance term ALONE is insufficient;",
          "they do not exclude variance in combination with baseline, transfer or speciation errors.",
          "\nNaCl comparison with the published fully dissociated model:",
          val.to_string(index=False,float_format=lambda x:f"{x:.5g}"),
          "CR-61 0.1 M source entry is 0.092; the adopted 0.0092 value is a possible decimal-error",
          "interpretation, not an author-confirmed correction. CR-61 charge input 3 M is approximate."]
    (RES/"summary.txt").write_text("\n".join(text)+"\n")
    print(f"Computed {len(part)} sorption observations. RMS log residual: ideal {scalars['rms_log_ideal']:.4f}; Manning {scalars['rms_log_manning']:.4f}.")
    print("Saved equilibrium, association, ionic-strength and design results to results/.")

if __name__=="__main__": main()
