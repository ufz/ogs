+++
project = "https://github.com/ufz/ogs-data/blob/master/ThermoHydroMechanics/square_1e0.prj_"
author = "Tianyuan Zheng"
date = "2017-05-15T15:10:33+01:00"
title = "Composite Axisymmetric Beam Model"
weight = 158

[menu]
  [menu.benchmarks]
    parent = "thermo-hydro-mechanics"

+++

{{< project-link >}}

## Problem description

We solve a thermo-hydro-mechanical linear biphasic model (small deformation, linear elastic, Darcy flow, incompressible constituents) in sealed composite beam where on the bottom boundary a no displacement boundary in y direction is applied and on the left boundary a no displacement boundary in x direction is applied. The domain is sealed so no fluid is allowed to come in and out of the domain. The beam is heated from the left boundary and finally reaches same temperature in the whole domain. An axisymmetric domain is used in this model.

## Assumptions

In this problem, it is assumed that the biot coefficient $\alpha = 1$ and the Storage term $S_\mathrm{s}$ is neglected.

## Results and evaluation
An analytical solution is used to verify the nuemrical result.
\begin{equation}
p = -K_\mathrm{S}\phi_\mathrm{F} (\beta_\mathrm{TF} - \beta_\mathrm{TS})\delta T
\end{equation}

The result at 1000s shows that the whole domain is heated up to 353.15K. At this stage, the pressure and stress are the same everywhere and equals 0.1042 which fits very well with the analytical solution.

{{< img src="../THM_homo_temperature.png" >}}

{{< img src="../THM_homo_pressure.png" >}}

{{< img src="../THM_homo_sigma_xx_.png" >}}
