# CFD Investigation of Aeration-Tank Diffuser Layouts

## Project overview

This project investigates how diffuser arrangement affects air
distribution, water circulation, oxygen transfer, dissolved oxygen,
and biodegradable COD removal in an activated-sludge aeration tank.

Three diffuser layouts are compared under identical total airflow:

- Case A: distributed diffuser rows
- Case B: central diffuser bank
- Case C: lateral diffuser banks

## Numerical approach

The aeration tank is modelled using an Euler-Euler air-water
multiphase formulation in OpenFOAM.

A custom aeration solver is used for dissolved-oxygen transport and
oxygen transfer.

Time-averaged hydrodynamic and aeration fields are subsequently used
in a one-way coupled biological transport model.

## Main results

| Case | Mean kLa (h^-1) | Gas holdup (%) | Low-velocity region (%) |
|------|----------------:|---------------:|------------------------:|
| A | 10.22 | 0.561 | 37.4 |
| B | 7.13 | 0.369 | 39.6 |
| C | 8.75 | 0.448 | 58.8 |

Biological performance:

| Case | Outlet DO (mg/L) | Outlet biodegradable COD (mg/L) | COD removal (%) |
|------|-----------------:|---------------------------------:|----------------:|
| A | 10.52 | 7.64 | 97.19 |
| B | 7.50 | 92.79 | 65.86 |
| C | 9.13 | 52.65 | 80.63 |

## Repository structure

- `cases/` OpenFOAM case setup files
- `solvers/` custom OpenFOAM solver source
- `postProcessing/` Python analysis and KPI data
- `figures/` selected CFD visualisations
- `videos/` selected animations
- `docs/` methodology and modelling notes

## Software

- OpenFOAM v2412
- ParaView
- Python
