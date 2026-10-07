# Models and numerical conventions

The package calculates equilibrium sorption in negatively charged membranes.
The experimental inputs cover Nafion 117 and CR-61 in NaCl, MgCl2, Na2SO4
and MgSO4. Each membrane/salt series contains five external concentrations.
[Data provenance](../data/README.md) identifies the source tables;
[References](REFERENCES.md) gives the publication metadata.

## Concentrations, charge and equilibrium

Superscripts `m` and `s` denote membrane water and external solution.
Concentrations are in mol/L of the corresponding aqueous phase. `C_A` is the
positive magnitude of the negative fixed charge, in mol of elementary charges
per litre of sorbed water. The signed charge of species `i` is `z_i` and its
concentration is `c_i`. Charge balance is

```text
sum_i z_i c_i^m = C_A
```

The fixed-charge input is calculated from the source ion contents. It is
supplied to the partition calculation, rather than fitted to co-ion uptake;
its measurement errors are consequently shared with the sorption data.
Association conserves charge, so a balance using elemental totals gives the
same fixed charge.

Let `gamma_i` denote a concentration-scale activity coefficient,
`h_i` a standard transfer factor and `psi_D = F(phi_m-phi_s)/(RT)` the
dimensionless Donnan potential. Species equilibrium gives

```text
c_i^m / c_i^s = h_i (gamma_i^s / gamma_i^m) exp(-z_i psi_D)
```

The reported calculations set the transfer factors to unity. This is a
constitutive assumption: charge neutrality alone does not determine the
transfer of a neutral species.

## Fully dissociated binary models

For a counter-ion `g` and co-ion `c`, define positive charge magnitudes
`Z_g = |z_g|` and `Z_c = |z_c|`. Their numbers per salt formula unit are
`nu_g` and `nu_c`, with `nu_g Z_g = nu_c Z_c`. The mobile formula-unit
concentration in the membrane is `C_m` and the external salt concentration is
`C_s`. Then

```text
c_g^m = C_A/Z_g + nu_g C_m
c_c^m = nu_c C_m
c_g^s = nu_g C_s
c_c^s = nu_c C_s
```

The ideal membrane model sets both membrane activity coefficients to one.
Manning/Donnan uses the dimensionless charge-to-salt ratio `Xi = C_A/C_m`
and Manning parameter `xi`. The latter is 3.1 for Nafion 117 and 1.83 for CR-61.
For monovalent fixed groups, the activity coefficients are

```text
D_M = Xi + Z_g Z_c xi (nu_g + nu_c)
gamma_g^m = [(Xi/(Z_g xi) + Z_g nu_g)/(Xi + Z_g nu_g)]
            exp[-Xi/(2 D_M)]
gamma_c^m = exp[-Xi/(2 D_M) (Z_c/Z_g)^2]
```

These expressions assume `Z_g xi > 1`. The salt equilibrium condition is

```text
(gamma_g^m c_g^m)^nu_g (gamma_c^m c_c^m)^nu_c
    = (gamma_pm^s)^(nu_g+nu_c)
      nu_g^nu_g nu_c^nu_c C_s^(nu_g+nu_c)
```

Here `gamma_pm^s` is the external mean ionic activity coefficient.
The implementation in `code/hpm/manning.py` solves this equation for `C_m`
using a logarithmic variable and Brent's bracketed method. The activity
coefficients are evaluated self-consistently at that concentration.
For MgCl2, measured chloride uptake is `2 C_m`. For the other three salts,
fully dissociated co-ion uptake is `C_m`.

The primary activity formulation is Kamcev et al. (2015), SI Section S3,
[DOI](https://doi.org/10.1021/acs.macromol.5b01654);
the membrane partition application is Kamcev et al. (2016),
[DOI](https://doi.org/10.1039/C5CP06747B).

## External activities and solution properties

The binary calculations use Pitzer activities at 298.15 K. The implementation
returns molal mean coefficients and converts them for the concentration-based
partition equations. Let `m` be molality, `rho_s` solution density,
`rho_w = 0.9970474 kg/L` water density and `M_s` salt molar mass in kg/mol:

```text
m = C_s / (rho_s - C_s M_s)
gamma_pm,c = gamma_pm,m(m) rho_w / (rho_s - C_s M_s)
```

Both the Pitzer argument and the activity coefficient are converted.
The numerical parameters are listed below; `beta2` is used for the 2:2 salt.

| Salt | beta0 | beta1 | C_phi | beta2 |
| --- | ---: | ---: | ---: | ---: |
| NaCl | 0.0765 | 0.2664 | 0.00127 | 0 |
| MgCl2 | 0.35235 | 1.6815 | 0.00519 | 0 |
| Na2SO4 | 0.01958 | 1.113 | 0.00497 | 0 |
| MgSO4 | 0.221 | 3.343 | 0.025 | -37.23 |
| LiCl | 0.1494 | 0.3074 | 0.00359 | 0 |
| CaCl2 | 0.3159 | 1.614 | -0.00034 | 0 |

The constants are `A_phi = 0.3915` and `b = 1.2`. Range parameters are
`alpha1 = 2.0` for salts with at least one monovalent ion; MgSO4 uses
`alpha1 = 1.4` and `alpha2 = 12`. The Pitzer functions operate on molality
in mol/kg and the corresponding stoichiometric molal ionic strength.
See `code/hpm/pitzer.py` and `code/hpm/properties.py`.

NaCl and sulfate densities come from Sujanani et al. (2023), SI Table S2.
MgCl2 density uses a constant apparent molar volume of 19.77 mL/mol,
anchored to Zezin et al. (2014), SI Table 1S at 298.15 K and 0.1 MPa.
Its use across concentrations is an explicit property approximation.
`solution_property(..., strict=True)` rejects approximate density inputs.
The MgSO4 mean-activity variation of +/-15% is an illustrative input sensitivity.

## Sulfate association

Two association treatments are evaluated on the four sulfate series:

1. Species equilibrium for free Mg2+ and SO4(2-) with neutral MgSO4, or
   free Na+ and SO4(2-) with charged NaSO4(-).
2. The weighted closures of Sujanani et al. (2023), SI S57/S70, evaluated
   with equal paired sulfate fractions in solution and membrane.

Both use the same external paired fractions and free-species activity model.
External fractions are conductivity-derived estimates from SI Table S2,
interpolated linearly in analytical concentration. Their interpretation
depends on the conductivity analysis's mobility assumptions.
Free-species activities use the Truesdell-Jones parameters in SI Table S9.
Parent-salt Pitzer coefficients are used in the fully dissociated calculations.

For each external state, the dimensionless association constant is inferred
from the species activities, each normalised to a 1 mol/L standard concentration:

```text
K_a = a_pair / (a_cation a_sulfate)
```

It is carried into the membrane under the unit-transfer assumption.
The inferred constant varies with the external state; it is a conditional input,
rather than an independently determined concentration-independent constant.
`paired_partition(..., K_a=...)` also accepts a consistently defined
caller-supplied association constant.

For membrane species activities, let `j` index negatively charged co-ion species.
Define a concentration-valued denominator `D_c` and factor `A_c`:

```text
D_c = C_A + xi sum_j |z_j| (Z_g + |z_j|) c_j^m
A_c = C_A / (2 D_c)
gamma_g^m = [c_g^m - C_A/Z_g + C_A/(xi Z_g^2)] / c_g^m
            exp(-A_c)
gamma_j^m = exp[-A_c (z_j/Z_g)^2]
```

These are the generalised Manning activities of Sujanani et al., SI S39/S40.
All charged species share one Donnan potential and enter the charge balance.
Neutral MgSO4 has unit activity coefficients and unit concentration partitioning.
The magnesium calculation uses a scalar root; the sodium calculation solves
free sulfate and NaSO4(-) concentrations jointly by bounded nonlinear least squares.
Logarithmic equilibrium residuals must be below `1e-8`.
Charge conservation and the association cycle are checked after solution.

Sulfate observations are analytical totals, so comparisons use the sum of
free sulfate and paired sulfate. Neutral MgSO4 predominates in the predicted
magnesium sulfate uptake under these assumptions, and unit partitioning
overpredicts both experimental series. The weighted comparator represents
the specified equal-fraction case, rather than the source authors' fitted
membrane association fractions.

## Ionic strength and charge transferability

Mobile ionic strength is `I = 0.5 sum_i z_i^2 c_i`; it excludes neutral
species and immobile polymer groups. For analytical sulfate, full dissociation
gives an upper reconstruction and charge-conserving association gives a
lower bound. Condensed counter-ions remain physical ions in this inventory.

The Debye curve uses water permittivity 78.3, temperature 298.15 K and
ion diameter 0.30 nm. It is a dilute screening reference for comparing
concentration and molecular scales.

Charge transferability is tested by inverting the ideal binary model at
each observed uptake. Signed apparent charges are retained, including
negative values that cannot describe negative polymer charge within that model.
The NaCl reference condition at 0.01 mol/L calibrates the charge transferred
to the MgCl2 series. The Manning comparison uses each condition's
charge-balance input.

## Identifiability and prospective power

The diagnostic response includes mean-potential, radius-weighted Born,
potential-variance and optional co-ion Manning-error terms. Chemical identities
share ion offsets across salts; salinity responses have zero change at the
reference ionic strength. Numerical ranks are evaluated after column normalisation.

Analytical sulfate totals are represented as free-sulfate channels only
for the favourable structural-rank check. The co-ion panel's variance
response lies in the specified nuisance span. Counter-ion measurements also
determine charge; their conditional residuals require joint uncertainty propagation.

Prospective simulations use independent measurements of both charge signs,
three salts (NaCl, MgCl2 and Na2SO4), optionally adding MgSO4.
Five external ionic strengths span 0.1 to 6 mol/L. Gaussian logarithmic
noise has standard deviation 0.10, with 1500 trials and seed 0.
The variance alternative is `max(0, log10(I/I0))` with `I0 = 1 mol/L`.
In the uncertain-baseline case, the co-ion activity-error amplitude is half
the variance amplitude. Power is a nested F-test rejection frequency at 0.05;
whiskers are 95% Wilson intervals. Undefined tests are recorded as NaN.
These results depend on the shared-response and independent-noise assumptions.

## Conditional lithium/sodium exchange

The exchange illustrations use ideal condensed mixing, unit free-ion
preference and unit free-population activity ratio. The monovalent Manning
condensed fraction is `theta = max(0, 1 - 1/xi)`. Its assumed capacity is
`theta C_A`. For monovalent co-ion uptake `c_-`, the condensed population
weight is `w = theta C_A/(C_A+c_-)`.

Let `y` be the reservoir lithium fraction, `x` the condensed lithium fraction,
`K_ex` the exchange constant and `alpha` the membrane/reservoir Li/Na ratio:

```text
x/(1-x) = K_ex y/(1-y)
alpha = ((1-y)/y) [w x + (1-w)y] / [w(1-x)+(1-w)(1-y)]
alpha_infinity = 1 + w/[y(1-w)]
```

Figures use `xi = 3.1`. Reservoir composition and co-ion uptake therefore
change the asymptote; `2 xi - 1` is the equimolar, negligible-co-ion special case.
Finite-bath calculations conserve both ion inventories.
The source dataset contains no lithium, and these curves are model predictions.

## Implementation comparison and verification

The NaCl comparison uses published fully dissociated predictions in
Sujanani SI Tables S13/S14. Additional CR-61 NaCl rows are implementation
checks outside the primary 30-point dataset. Their charge input is approximately
3 mol/L, based on the separate membrane characterisation cited in the data.
The CR-61 0.1 mol/L source entry is 0.092 mol/L; the adopted 0.0092 mol/L
interpretation is recorded alongside it and is not author-confirmed.

`tests/test_theory.py` checks independent equilibrium relations, charge and
element conservation, association cycles, concentration-scale conversions,
exchange limits, finite inventories, design rank and null-test behaviour.
`code/verify_outputs.py` checks measurements, numerical reference values,
figure-data provenance, residual ranges and all 20 figure exports.
