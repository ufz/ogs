# Double diffusion experiment — simulation metadata

## Purpose

This directory is a placeholder for an OpenGeoSys benchmark of the saturated
double through-diffusion experiment in Boom Clay. It is intended for the
**DigBen** Model Hub and is based on the experiment and modelling study by
Pitz et al. (2024).

## Project context

The benchmark was developed within the European Joint Programme on Radioactive
Waste Management (EURAD), Work Package GAS (EURAD WP-GAS). The project was
funded by the Horizon 2020 Euratom programme under grant agreement No. 847593
(2019--2024).

## Scientific reference

Pitz, M., et al. (2024). *On Multi-Component Gas Migration in Single-Phase
Systems*. Rock Mechanics and Rock Engineering, 57, 4251--4264.
<https://doi.org/10.1007/s00603-024-03838-1>

## Experiment

| Field | Value |
| --- | --- |
| Scale | Laboratory; double through-diffusion cell |
| Rock | Boom Clay from borehole ON-Mol-1, Mol, Belgium |
| Specimen | Fully liquid-saturated cylinder; 30 mm long, 40 mm radius |
| Reservoirs | Two 1 L vessels, each containing water and a gas phase |
| Initial gases | Helium in one reservoir; methane in the other |
| Initial total pressure | 1.0 MPa in both reservoirs |
| Measurement duration | 72 days |
| Measured quantities | Reservoir pressure and gas composition by gas chromatography |

Helium diffuses towards the methane reservoir, while methane diffuses in the
opposite direction. The laboratory setup does not impose a total gas- or
liquid-pressure gradient across the sample, so the intended transport mechanism
is diffusion of dissolved gas rather than advection.

## OpenGeoSys representation

The model uses the OpenGeoSys `TH2M` process. Temperature and mechanical
quantities are fixed in this benchmark, so the active part of the formulation
is effectively hydraulic (`H2`).
