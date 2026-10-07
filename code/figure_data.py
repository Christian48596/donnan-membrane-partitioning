"""Export numerical tables used by plotting; no digitisation or fitted smoothing."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from hpm.screening import debye_length, crossover_concentration

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def wilson_interval(probability, trials):
    """95% Wilson interval for the Monte Carlo detection proportion."""
    if not np.isfinite(probability):
        return np.nan, np.nan
    z = 1.959963984540054
    denom = 1 + z*z/trials
    mid = (probability + z*z/(2*trials))/denom
    half = z*np.sqrt(probability*(1-probability)/trials + z*z/(4*trials*trials))/denom
    return mid-half, mid+half


def export_figure_data(output_dir=None):
    out = (Path(output_dir) if output_dir else ROOT / "figures") / "source_data"
    out.mkdir(parents=True, exist_ok=True)
    part = pd.read_csv(RESULTS / "partitioning.csv")
    pairs = pd.read_csv(RESULTS / "ion_pairing.csv")
    interior = pd.read_csv(RESULTS / "interior_ionic_strength.csv")
    identity = ["membrane", "salt", "c_ext"]
    merged = part.merge(pairs[identity+["paired", "published_closure_equal_alpha",
                                       "r_paired", "r_published_equal_alpha"]],
                        on=identity, how="left", validate="one_to_one")
    measurements = pd.read_csv(ROOT / "data/sorption_published.csv", comment="#").rename(
        columns={"salt_system": "salt", "c_external_M": "c_ext"})
    merged = merged.merge(measurements[identity+["source"]], on=identity,
                          validate="one_to_one")
    year = merged.source.str.extract(r"Sujanani (2021|2023)", expand=False)
    if year.isna().any():
        raise ValueError("Every measured series must have a traceable author and year")
    merged["measurement_reference"] = "Sujanani et al. (" + year + ")"
    if len(merged) != 30:
        raise ValueError("Expected 30 source sorption rows")
    tables = {}
    I = np.geomspace(1e-3, 12, 500)
    tables["fig1_screening_curve.csv"] = pd.DataFrame({
        "I_M": I, "lambda_D_nm": debye_length(I, 78.3)*1e9,
        "ion_diameter_nm": .30, "I_crossover_M": crossover_concentration(78.3),
        "temperature_K": 298.15, "relative_permittivity": 78.3})
    ranges = [dict(environment="External solution", lower_M=part.I_ext.min(),
                   upper_M=part.I_ext.max(), basis="fully dissociated reference")]
    for membrane in ("Nafion117", "CR-61"):
        group = interior[interior.membrane == membrane]
        ranges.append(dict(environment=membrane,
                           lower_M=group.I_mobile_min_association_M.min(),
                           upper_M=group.I_mobile_dissociated_M.max(),
                           basis="mobile-ion reconstruction with association bounds"))
    tables["fig1_screening_ranges.csv"] = pd.DataFrame(ranges)
    tables["fig2_sorption.csv"] = merged
    tables["fig3_residuals.csv"] = merged
    summary = []
    for (membrane, salt), group in merged.groupby(["membrane", "salt"], sort=False):
        summary.append(dict(membrane=membrane, salt=salt, n=len(group),
                            manning_mean=group.r_manning.mean(),
                            manning_min=group.r_manning.min(),
                            manning_max=group.r_manning.max(),
                            conditional_mean=group.r_paired.mean(),
                            conditional_min=group.r_paired.min(),
                            conditional_max=group.r_paired.max(),
                            published_mean=group.r_published_equal_alpha.mean(),
                            published_min=group.r_published_equal_alpha.min(),
                            published_max=group.r_published_equal_alpha.max()))
    tables["fig3_series_means.csv"] = pd.DataFrame(summary)
    ident = pd.read_csv(RESULTS / "identifiability.csv")
    for basis in ("exact", "uncertain"):
        limits = [wilson_interval(p, n) for p, n in
                  zip(ident[f"power_baseline_{basis}"], ident.n_trials)]
        ident[f"power_{basis}_ci_low"] = [v[0] for v in limits]
        ident[f"power_{basis}_ci_high"] = [v[1] for v in limits]
    tables["fig4_identifiability.csv"] = ident
    tables["fig5_apparent_charge.csv"] = pd.read_csv(RESULTS / "apparent_fixed_charge.csv")
    tables["fig5_transferability.csv"] = pd.read_csv(RESULTS / "dle_transferability.csv").merge(
        part.query('membrane == "Nafion117" and salt == "MgCl2"')[["c_ext", "co_sd"]],
        on="c_ext", validate="one_to_one")
    curves = pd.read_csv(RESULTS / "ion_specific_curves.csv")
    tables["fig5_exchange_curves.csv"] = curves[np.isclose(curves.xi, 3.1)].copy()
    required = {(y, co) for y in (.01, .1, .5) for co in (0., .5, 1.)}
    present = set(zip(tables["fig5_exchange_curves.csv"].y,
                      tables["fig5_exchange_curves.csv"].coion_to_fixed_charge))
    if not required <= present:
        raise ValueError("Run ion_specific.py to generate every plotted exchange case")
    entries = []
    for name, table in tables.items():
        path = out / name
        table.to_csv(path, index=False)
        entries.append(dict(file=name, rows=len(table), columns=list(table.columns),
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    inputs = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(RESULTS.glob("*.csv"))}
    manifest = {"tables": entries, "input_sha256": inputs,
                "sorption_lines": "straight connections of calculated points; no fitted smoothing",
                "measurement_error_bars": "reported co-ion standard deviations only",
                "residual_range_bars": "min-to-max of five concentration-dependent residuals for each model; not statistical uncertainty",
                "power_error_bars": "95% Wilson intervals for Monte Carlo proportions",
                "screening_curve": "dilute Debye reference, not a measured membrane decay length",
                "lithium_curves": "conditional reservoir illustration, no lithium validation"}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(f"  exported {len(tables)} numerical figure tables to {out}")
    return out


if __name__ == "__main__":
    export_figure_data()
