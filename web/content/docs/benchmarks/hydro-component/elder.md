+++
author = "Marc Walther"
weight = 142
project = "Parabolic/ComponentTransport/elder/"
date = "2017-09-07T14:41:09+01:00"
title = "Saturated Variable-Density Flow and Mass Transport (Elder)"

[menu]
  [menu.benchmarks]
    parent = "hydro-component"

+++

{{< project-link >}}


## Overview

This benchmark is one of the classical free-convection density-driven flow and mass transport setups. It was originally published by Elder (1965) and has since then been used as a basic test case (Diersch & Kolditz, 1998; Guo & Langevin, 2002), or has been subject to various investigations concerning grid convergence (Graf & Degener, 2011) and numerical stability (Musuuza et al., 2009; Johannsen, 2003).

For the setup and parameterization, see the chapter "Density dependent flow - The Elder Problem" in Kolditz et al. (2012).


## Problem description

The Elder benchmark describes free convection of a dense fluid in mixable, single-phase environment. A high-concentration solute increases fluid density on the upper boundary and perpetrates the domain by evolving concentration fingers. Here, we compare numerical results of OGS-6 to those of OGS-5. Settings of both simulators were chosen to be as identical as possible. Simulation times were $3300 s$ and $7800 s$ for OGS-6 and OGS-5, respectively.


### Model results

A comparison of the numerical data is shown in the figure below. The numerical results of OGS-6 coincide with those of OGS-5.

{{< img src="../gif/elder.png" title="Results for numerical (OGS-5 - green, OGS-6 - white) results together with concentration distribution in the domain and mesh resolution for different time steps.">}}

[The project files are here. ](../../../../../Tests/Data/Parabolic/ComponentTransport/elder)

## Literature

Diersch, H.-J.G., Kolditz, O., 1998. Coupled groundwater flow and transport: 2. Thermohaline and 3D convection systems. Adv. Water Resour. 21, 401–425. doi:10.1016/S0309-1708(97)00003-1.

Elder, J.W., 1965. Numerical experiments with free convection in a vertical slot. J. Fluid Mech. 24, 823. doi:10.1017/S0022112066001022.

Elder, J., Simmons, C., Diersch, H.-J., Frolkovič, P., Holzbecher, E., Johannsen, K., 2017. The Elder Problem. Fluids 2, 11. doi:10.3390/fluids2010011.

Graf, T., Degener, L., 2011. Grid convergence of variable-density flow simulations in discretely-fractured porous media. Adv. Water Resour. 34, 760–769. doi:10.1016/j.advwatres.2011.04.002.

Guo, W., Langevin, C.D., 2002. User’s Guide to SEAWAT: A computer program for simulation of three-dimensional variable-density ground-water flow, USGS Techniques of Water Resources Investigations. ISBN: 0607992573.

Johannsen, K., 2003. On the Validity of the Boussinesq Approximation for the Elder Problem. Comput. Geosci. 7, 169–182. doi:10.1023/A:1025515229807.

Kolditz, O., Görke, U.-J., Shao, H., Wang, W., 2012. Thermo-Hydro-Mechanical-Chemical Processes in Porous Media: Benchmarks and Examples, Lecture notes in computational science and engineering. Springer. ISBN: 3642271766.

Musuuza, J.L., Attinger, S., Radu, F.A., 2009. An extended stability criterion for density-driven flows in homogeneous porous media. Adv. Water Resour. 32, 796–808. doi:10.1016/j.advwatres.2009.01.012.
