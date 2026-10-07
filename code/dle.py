"""Fixed-charge transferability and equal-valence electrostatic cancellation."""
from pathlib import Path
import numpy as np
import pandas as pd
import analysis as A
from hpm.manning import XI_NAFION117, STOICH
from hpm.pitzer import molar_activity_product
from hpm.properties import solution_property

ROOT=Path(__file__).resolve().parents[1]
DATA,RES=ROOT/"data",ROOT/"results"
MEMBRANE="Nafion117"
XI=XI_NAFION117


def fit_fixed_charge(salt,c_s,c_co_obs):
    """Signed apparent charge obtained algebraically, without hiding invalid fits.

    Negative values mean no admissible negatively charged, ideal-membrane fit.
    Sulfate totals are deliberately treated as fully dissociated in this baseline
    inversion; this identifies an apparent parameter, not a measured material charge.
    """
    if c_co_obs<=0: raise ValueError("positive observed co-ion concentration required")
    s=STOICH[salt]
    rhs=molar_activity_product(salt,c_s,solution_property(salt,c_s))
    ct=(rhs/c_co_obs**s.nu_c)**(1/s.nu_g)
    return float(s.z_g*ct-s.z_c*c_co_obs)


def apparent_charge_table():
    rows=[]
    for r in A.sorption_data().itertuples():
        fitted=fit_fixed_charge(r.salt_system,r.c_external_M,r.c_coion_M)
        rows.append(dict(membrane=r.membrane,salt=r.salt_system,c_ext=r.c_external_M,
                         C_fixed=r.C_fixed_M,C_apparent=fitted,
                         admissible_negative_membrane=fitted>=0,
                         charge_underestimate_factor=r.C_fixed_M/fitted if fitted>0 else np.nan))
    return pd.DataFrame(rows)


def transferability():
    df=A.sorption_data(); df=df[df.membrane==MEMBRANE]
    cal=df[df.salt_system=="NaCl"].sort_values("c_external_M").iloc[0]
    Xfit=fit_fixed_charge("NaCl",cal.c_external_M,cal.c_coion_M)
    rows=[]
    for r in df[df.salt_system=="MgCl2"].itertuples():
        pred=A.ideal_donnan_coion("MgCl2",r.c_external_M,Xfit)
        man=A.manning_coion("MgCl2",r.c_external_M,r.C_fixed_M,XI)
        rows.append(dict(c_ext=r.c_external_M,obs=r.c_coion_M,ideal_transferred=pred,
                         manning_measured=man,r_transferred=np.log(r.c_coion_M/pred),
                         r_manning=np.log(r.c_coion_M/man),calibration_C_apparent=Xfit,
                         calibration_C_fixed=cal.C_fixed_M,calibration_feed_M=cal.c_external_M))
    out=pd.DataFrame(rows); out.attrs.update(X_fit=Xfit,X_true=cal.C_fixed_M)
    return out


def lithium_blindness():
    """Like-valence factor cancellation, not separate-salt co-ion selectivity.

    With equal activity and standard transfer prefactors, k_Li/k_Na=1 for all
    psi_D. A caller's ion-specific prefactor ratio q remains q, independent of
    the mean Donnan potential. No Li measurements or mixed Pitzer model implied.
    """
    rows=[]
    for psi in (-5.,-2.,0.,2.):
        for q in (1.,2.):
            kNa=np.exp(-psi); kLi=q*kNa
            rows.append(dict(psi_D=psi,ion_specific_prefactor_ratio=q,
                             k_Li=kLi,k_Na=kNa,alpha_Li_Na=kLi/kNa))
    return pd.DataFrame(rows)


def main():
    RES.mkdir(exist_ok=True)
    t,ca,li=transferability(),apparent_charge_table(),lithium_blindness()
    t.to_csv(RES/"dle_transferability.csv",index=False)
    ca.to_csv(RES/"apparent_fixed_charge.csv",index=False)
    li.to_csv(RES/"dle_like_valence.csv",index=False)
    text=["FIXED-CHARGE TRANSFERABILITY AND Li/Na LIMITATION", "="*72,
          f"Reference NaCl 0.01 M: measured charge {t.attrs['X_true']:.5f} M, "
          f"apparent ideal-membrane charge {t.attrs['X_fit']:.5f} M.",
          f"Underestimate factor {t.attrs['X_true']/t.attrs['X_fit']:.4f}.",
          f"MgCl2 transferred ideal mean log residual {t.r_transferred.mean():+.5f}; "
          f"Manning with measured charge {t.r_manning.mean():+.5f}.",
          f"Apparent charge range {ca.C_apparent.min():.5g} to {ca.C_apparent.max():.5g} M; "
          f"{(~ca.admissible_negative_membrane).sum()} of {len(ca)} fits require opposite-sign charge.",
          "Negative apparent charges are inadmissible for the negatively charged membrane.",
          "\nSame-charge Donnan factors cancel: alpha_Li/Na equals the ratio of transfer",
          "and activity prefactors. A charge-only model with equal prefactors gives 1.",
          "Li/Na equilibrium selectivity is the ratio of counter-ion partition coefficients."]
    (RES/"dle_summary.txt").write_text("\n".join(text)+"\n")
    print(f"NaCl reference: charge balance {t.attrs['X_true']:.5f} M; apparent charge {t.attrs['X_fit']:.5f} M.")
    print("Saved signed charge fits, MgCl2 transfer predictions and equal-valence checks to results/.")

if __name__=="__main__": main()
