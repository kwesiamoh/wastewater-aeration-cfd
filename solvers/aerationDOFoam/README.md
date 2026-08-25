# aerationDOFoam

`aerationDOFoam` is a custom OpenFOAM v2412 solver for transient Euler–Euler gas–liquid aeration with dissolved-oxygen (DO) transport and local oxygen mass transfer. It extends the `twoPhaseEulerFoam` framework for fine-bubble aeration in activated-sludge tanks.

This is the first stage of the project workflow. It resolves the multiphase aeration field and provides time-averaged and derived inputs for [`aerationASM1ReducedFoam`](../aerationASM1ReducedFoam/README.md).

## Modelling scope

The solver predicts the relationship between:

- gas distribution and liquid circulation;
- gas–liquid slip velocity and interfacial area;
- the local liquid-side coefficient `kL` and volumetric coefficient `kLa`;
- pressure-dependent local oxygen saturation; and
- transient dissolved-oxygen concentration.

It models physical aeration only. Biological substrate consumption and biological oxygen demand are introduced in the second-stage reduced ASM1 solver.

## Dissolved-oxygen equation

The conservative liquid-phase DO equation is

```math
\frac{\partial(\alpha_w C_O)}{\partial t}
+ \nabla\cdot(\alpha\Phi_w C_O)
- \nabla\cdot(\alpha_w D_{\mathrm{eff}}\nabla C_O)
= \alpha_w k_La(C_O^* - C_O),
```

where $\alpha_w$ is the instantaneous water volume fraction, $C_O$ is DO concentration, $\alpha\Phi_w$ is the conservative water-phase flux, and $C_O^*$ is the local saturation concentration.

The effective diffusivity is

```math
D_{\mathrm{eff}} = D_{m,O_2} + \frac{\nu_{t,w}}{Sc_T}.
```

For monodisperse spherical bubbles, the solver evaluates

```math
a = \frac{6\alpha_g}{d_b\alpha_w},
\qquad
k_L = 2\sqrt{\frac{D_{m,O_2}|U_g - U_w|}{\pi d_b}},
\qquad
\left(k_La\right)_{\mathrm{wastewater}} = \alpha_{\mathrm{factor}} k_L a.
```

The local saturation concentration is pressure dependent:

```math
C_O^* = \beta C_{O,\mathrm{ref}}^*\frac{p}{p_{\mathrm{ref}}}.
```

Consequently, the reference surface saturation value is not a uniform upper limit throughout the tank.

## Required case data

In addition to the fields and dictionaries required by the parent Euler–Euler solver, a case must contain:

- `0/DO`;
- `constant/DOProperties`.

`DOProperties` defines:

| Entry | Meaning |
|---|---|
| `DmDO` | Molecular diffusivity of oxygen in water |
| `ScT` | Turbulent Schmidt number |
| `DOsatRef` | Reference clean-water saturation concentration |
| `pRefDO` | Reference absolute pressure |
| `alphaFactor` | Wastewater mass-transfer correction |
| `betaFactor` | Wastewater saturation correction |
| `dBubble` | Effective bubble diameter |

The present study uses `DOsatRef = 0.0108 kg/m³`, corresponding to an approximate reference surface saturation DO of 10.8 mg/L at 12 °C. The wastewater-corrected value at the reference pressure is 10.26 mg/L when `betaFactor = 0.95`.

## Output fields

The solver writes the transient DO and aeration diagnostics, including:

- `DO`;
- `slipSpeed`;
- `interfacialAreaLiquid`;
- `kL`;
- `kLa`;
- `DOsatLocal`.

The project workflow subsequently time averages and derives the fields passed to Stage 2. These include `U.waterMean`, `alpha.airMean`, `alpha.waterMean`, `alphaUwaterMean`, `nut.waterMean`, `DOMean`, `DOsatLocalMean`, and `oxygenTransferCoeffPostMean`. They are not all direct solver outputs.

## Building

From this directory in an OpenFOAM v2412 environment:

```bash
wmake
```

The executable is written to `$FOAM_USER_APPBIN`. Check it with:

```bash
which aerationDOFoam
```

## Main source files

| File | Purpose |
|---|---|
| `aerationDOFoam.C` | Advances the transient Euler–Euler solution and invokes DO transport |
| `createDOFields.H` | Reads `DO` and the oxygen-transfer parameters |
| `DOEqn.H` | Calculates local transfer quantities and solves the DO equation |
| `Make/files` | Defines the application target |
| `Make/options` | Defines OpenFOAM include paths and libraries |

## Assumptions and limitations

- A single effective bubble diameter is used; breakup, coalescence, and bubble-size distributions are not solved.
- `alphaFactor` and `betaFactor` are prescribed modelling inputs without site-specific calibration.
- Direct atmospheric free-surface reaeration is not included.
- The transient solver retains the dispersed-bubble interfacial-area expression near the resolved atmosphere interface with a numerical lower bound on water volume fraction. Atmosphere-adjacent bubble-transfer contributions are removed when preparing the Stage 2 coupling field.

## Verification

The implementation has been checked for:

- dimensional consistency of the DO transport and oxygen-transfer terms;
- the oxygen-transfer source sign;
- zero net transfer at local saturation;
- positive oxygen transfer below local saturation;
- stripping tendency above local saturation;
- agreement with the expected first-order oxygen-transfer behavior; and
- positivity and physically consistent ranges in the production DO field.

These checks cover mathematical and numerical verification of the implemented model.

## References

- Fayolle, Y., Cockx, A., Gillot, S., Roustan, M., & Héduit, A. (2007). Oxygen transfer prediction in aeration tanks using CFD. *Chemical Engineering Science, 62*(24), 7163–7171. https://doi.org/10.1016/j.ces.2007.08.082
- Sánchez, F., Rey, H., Viedma, A., Nicolás-Pérez, F., Kaiser, A. S., & Martínez, M. (2018). CFD simulation of fluid dynamic and biokinetic processes within activated sludge reactors under intermittent aeration regime. *Water Research, 139*, 47–57. https://doi.org/10.1016/j.watres.2018.03.067
- Shah, K. A., Jiao, Y., & Chen, J. (2024). CFD investigation of dissolved oxygen distribution in a full-scale aeration tank of an industrial wastewater treatment plant. *Journal of Water Process Engineering, 59*, 105078. https://doi.org/10.1016/j.jwpe.2024.105078
