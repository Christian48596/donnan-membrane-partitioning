# Input data and provenance

| File | Contents |
| --- | --- |
| `sorption_published.csv` | 30 tabulated sorption observations, charge-balance inputs, source labels and reported co-ion standard deviations |
| `solution_properties.csv` | Solution densities and conductivity-derived sulfate association fractions at 298.15 K |
| `published_model_output.csv` | Published NaCl predictions for the implementation comparison |

## Sorption measurements

Measurements are transcribed from tables in:

- Sujanani et al., *Journal of Membrane Science* 638 (2021), 119687,
  [DOI 10.1016/j.memsci.2021.119687](https://doi.org/10.1016/j.memsci.2021.119687).
- Sujanani et al., *Journal of Physical Chemistry B* 127 (2023), 1842–1855,
  SI Tables S3, S4, S7 and S8,
  [DOI 10.1021/acs.jpcb.2c07900](https://doi.org/10.1021/acs.jpcb.2c07900).

Each series has external concentrations 0.01, 0.03, 0.10, 0.30 and 1.00 mol/L.
The six series are Nafion 117/NaCl, Nafion 117/MgCl2, Nafion 117/Na2SO4,
CR-61/Na2SO4, Nafion 117/MgSO4 and CR-61/MgSO4.

| Column | Meaning |
| --- | --- |
| `membrane` | Membrane identity |
| `salt_system` | Analytical salt formula |
| `c_external_M` | External analytical concentration, mol/L of solution |
| `C_fixed_M` | Fixed-charge magnitude from source ion balances, mol of elementary charges/L of sorbed water |
| `c_coion_M` | Chloride concentration or analytical total sulfate, mol/L of sorbed water |
| `c_coion_sd` | Reported standard deviation on the same membrane-water basis |
| `coion` | Analytical co-ion identity |
| `z_coion` | Positive co-ion charge magnitude for the fully dissociated reference |
| `source` | Source attribution |

The charge input is derived from the same source ion contents used for
sorption. It is not fitted to uptake and is not an independent charge
measurement. Sulfate totals include ion pairs; they are not free-sulfate
concentrations.

## Solution properties

Sulfate fractions and NaCl/sulfate densities are transcribed from Sujanani SI
Table S2, with linear concentration interpolation.
[Supporting information DOI](https://doi.org/10.1021/acs.jpcb.2c07900.s001).

MgCl2 densities use a constant apparent molar volume of 19.77 mL/mol,
anchored to Zezin et al. (2014), SI Table 1S at 298.15 K and 0.1 MPa.
[Supporting information DOI](https://doi.org/10.1021/je500371u.s001).
This is a density approximation across the analysis concentrations.

## Published model comparison

NaCl model outputs come from Sujanani SI Tables S13/S14. CR-61 NaCl is an
additional implementation-check series. Its charge input is approximately
3 mol/L. The source's CR-61 0.1 mol/L entry is 0.092 mol/L; the adopted
0.0092 mol/L interpretation is recorded explicitly in `results/validation.csv`
and has not been confirmed by the source authors.

## Scope

No lithium observations are present. Exchange curves are conditional model
predictions, and the independent-ion power panels are prospective simulations.
Calculated values are in `results/` and plotted values in `figures/source_data/`.
Full citation metadata is in `docs/references.bib`.
