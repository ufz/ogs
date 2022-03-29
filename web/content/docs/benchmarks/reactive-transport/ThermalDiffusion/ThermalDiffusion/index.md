+++
author = "Jaime Garibay-Rodriguez, Renchao Lu, Vanessa Montoya"
weight = 142
project = "Parabolic/ComponentTransport/ThermalDiffusion/TemperatureField_transport.prj"
date = "2019-07-08T14:41:09+01:00"
title = "Tracer diffusion in a thermal gradient"

[menu]
  [menu.benchmarks]
    parent = "Reactive Transport"

+++

{{< data-link >}}


## Overview

This benchmark simulates the diffusion of a non-reactive tracer in clay during a thermal gradient. Here are the relevant parts of the benchmark:

1. The test follows the 'HT' process in a first model to generate a temperature field (`TemperatureField.prj`) in the subsurface. A 8 by 4 m 2D domain is selected with finer elements on the left side (borehole) for which a Dirichlet boundary condition is applied with a temperature of 353.15 K. The initial temperature of the media is 289.15 K.

2. Opalinus Clay is selected as porous media. Full saturation and instantaneous thermal equilibrium with its porewater is assumed.

3. After one year of heating, it is assumed that the temperature profile is at a quasi-steady-state. The output of this is then used to compute temperature dependent diffusion coefficients to be used in the next part of the model.

4. With a newly generated mesh containing the temperature-dependent diffusion coefficients, a second model is set up (`TemperatureField_transport.prj`) with the same geometry and the 'ComponentTransport' process. The model follows the diffusion of a non-reactive tracer on the left side.

Both models `TemperatureField.prj` and `TemperatureField_transport.prj` run independently from each other.

## Model setup and results

See [this PDF](DiffusionThermalGradient.pdf).
