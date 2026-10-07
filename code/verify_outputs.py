"""Verify numerical results, figure source data and print-size exports."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]


def main():
    measured = pd.read_csv(ROOT / "data/sorption_published.csv", comment="#")
    predictions = pd.read_csv(ROOT / "results/partitioning.csv")
    assert len(measured) == len(predictions) == 30
    reference = measured.rename(columns={"salt_system": "salt", "c_external_M": "c_ext"})
    joined = predictions.merge(reference, on=["membrane", "salt", "c_ext"], validate="one_to_one")
    np.testing.assert_allclose(joined.co_obs, joined.c_coion_M, rtol=0, atol=1e-15)
    np.testing.assert_allclose(joined.C_fixed, joined.C_fixed_M, rtol=0, atol=1e-14)
    np.testing.assert_allclose(joined.co_sd, joined.c_coion_sd, rtol=0, atol=1e-15)
    for model in ("ideal", "manning"):
        np.testing.assert_allclose(predictions[f"r_{model}"],
                                   np.log(predictions.co_obs / predictions[f"co_{model}"]), atol=1e-12)
    values = json.loads((ROOT / "results/metrics.json").read_text())
    expected = json.loads((ROOT / "tests/reference_metrics.json").read_text())
    for key in ("n_measurements", "rms_log_ideal", "rms_log_manning",
                "largest_ideal_underprediction", "mobile_I_min_M", "mobile_I_max_M",
                "mobile_I_association_bound_min_M", "crossover_M"):
        np.testing.assert_allclose(values[key], expected[key], rtol=1e-7, atol=1e-10, err_msg=key)
    assert values["max_charge_error_M"] < 1e-9
    assert values["max_association_cycle_relative_error"] < 1e-9
    np.testing.assert_allclose(values["rms_log_manning"], np.sqrt(np.mean(predictions.r_manning**2)))
    manifest = json.loads((ROOT / "figures/manifest.json").read_text())
    assert len(manifest["figures"]) == 5
    for name, digest in manifest["figure_data"]["input_sha256"].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, f"stale source: {name}"
    for entry in manifest["figure_data"]["tables"]:
        path = ROOT / "figures/source_data" / entry["file"]
        assert len(pd.read_csv(path)) == entry["rows"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
    plotdata = pd.read_csv(ROOT / "figures/source_data/fig2_sorption.csv")
    np.testing.assert_allclose(plotdata.co_obs, predictions.co_obs, rtol=0, atol=1e-15)
    assert plotdata.measurement_reference.notna().all()
    summary = pd.read_csv(ROOT / "figures/source_data/fig3_series_means.csv")
    for _, row in summary.iterrows():
        group = plotdata[(plotdata.membrane == row.membrane) & (plotdata.salt == row.salt)]
        for prefix, residual in (("manning", "r_manning"), ("conditional", "r_paired"),
                                 ("published", "r_published_equal_alpha")):
            np.testing.assert_allclose(
                [row[prefix+"_mean"], row[prefix+"_min"], row[prefix+"_max"]],
                [group[residual].mean(), group[residual].min(), group[residual].max()],
                equal_nan=True)
    for figure in manifest["figures"]:
        width = 83 if figure["stem"] in ("fig1_screening_regimes", "fig4_identifiability") else 171
        np.testing.assert_allclose(figure["width_mm"], width, atol=1e-8)
        assert figure["height_mm"] <= 233
        assert figure["minimum_text_pt"] >= manifest["style"]["minimum_text_pt"]
        assert figure["layout_bounds_check"] == "passed"
        assert {Path(f["path"]).suffix for f in figure["files"]} == {".pdf", ".svg", ".png", ".tif"}
        for entry in figure["files"]:
            path = ROOT / "figures" / entry["path"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"], f"incomplete or changed export: {path.name}"
            if path.suffix in (".png", ".tif"):
                with Image.open(path) as im:
                    dpi = im.info["dpi"]
                    assert all(abs(v-600) < .1 for v in dpi)
                    assert abs(im.width / dpi[0] * 25.4 - width) < .08
    report = {"measurement_rows_preserved": 30, "numerical_regression": "passed",
              "conservation_and_association_cycles": "passed",
              "figure_data_hashes": "passed", "figure_exports": 20,
              "print_widths_mm": {f["stem"]: f["width_mm"] for f in manifest["figures"]},
              "minimum_text_pt": manifest["style"]["minimum_text_pt"],
              "font_family": manifest["style"]["font_family"], "raster_dpi": 600,
              "layout_bounds": "passed"}
    (ROOT / "results/reproduction_verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print("Verified: 30 measurements, numerical regression, conserved charge and association cycles.")
    print(f"Verified: nine figure-data tables and 20 exports at 83/171 mm, with base text at least {manifest['style']['minimum_text_pt']} pt.")


if __name__ == "__main__":
    main()
