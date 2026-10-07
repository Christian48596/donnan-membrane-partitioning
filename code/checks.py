"""Run invariant tests and report the two approximation diagnostics."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from hpm.moments import exact_ln_K, moment1, variance
ROOT=Path(__file__).resolve().parents[1]


def main():
    rp,lam,k0=1e-9,.5e-9,2*np.pi/.30e-9
    rows=[]
    for zpsi in (1.,2.,4.,6.,8.):
        exact=exact_ln_K(1,zpsi,lam,k0,rp)
        cumulant=-moment1(zpsi,lam,k0,rp)+.5*variance(zpsi,lam,k0,rp)
        rows.append(dict(z_psi_surface=zpsi,ln_K_exact=exact,ln_K_second_cumulant=cumulant,
                         relative_error_percent=100*(cumulant-exact)/abs(exact)))
    df=pd.DataFrame(rows)
    df.to_csv(ROOT/"results"/"cumulant_check.csv",index=False)
    activity=pd.read_csv(ROOT/"results"/"activity_basis.csv")
    text=["APPROXIMATION CHECKS (illustrative potential profile)",df.to_string(index=False),
          "Full molality-to-concentration conversion includes standard-state rescaling.",
          f"Largest absolute log salt-product shift, argument only: "
          f"{activity.log_salt_product_shift_argument_only.abs().max():.6f}.",
          f"Largest absolute shift, full conversion: "
          f"{activity.log_salt_product_shift_full_conversion.abs().max():.6f}."]
    (ROOT/"results"/"checks_summary.txt").write_text("\n".join(text)+"\n")
    print("\n".join(text))

if __name__=="__main__": main()
