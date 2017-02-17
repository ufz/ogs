+++
project = "Parabolic/HT/ConstViscosity/square_5500x5500.prj"
author = "Thomas Fischer"
date = "2017-02-16T15:17:53+01:00"
title = "Constant viscosity"
weight = 161

[menu]
  [menu.benchmarks]
    parent = "hydro-thermal"

+++

{{< project-link >}}

## Equations

TODO: `/assets/files/SelectedBenchmarks/HT/HT-Process.pdf`

## Problem description

This is a 2d benchmark of large-scale thermal convection that tests the temperature dependent fluid density in the hydro-thermal process monolitic approach implementation. It is defined on the domain $\Omega = [0,5500]^2.$

- The initial conditions for the pressure is a gradient starting from zero at the top surface to a pressure of circa 54 mega pascal at the bottom given in the data array 'initial_pressure' in the vtu file. The initial temperature is also almost a gradient from top (293 K) to bottom (443 K) of the domain, except there is a small pertubation given by adding $\sin \left( \pi \frac{y}{5500}\right) \cdot \cos \left( \pi \frac{x}{5500}\right).$ See the following images.

TODO 3 images

- For the pressure we set at the top left point of $\Omega$ a Dirichlet-type boundary condition with a value of zero. On the top of $\Omega$ the temperature is fixed to 293 K (20 degree C), at the bottom  of $\Omega$ the temperature is set to 443 K (170 degree C) via Dirichlet-type boundary conditions.
- The further parameter specification can be found in the project file linked at the top of this page.
- The steady state temperature is shown in the following on the right figure. The left figure shows the resulting temperature minus the initial gradient. The resulting temperatures are in good accordance with the FEFLOW results and with results from the OGS version < 6.

{{< vis path="Parabolic/HT/ConstViscosity/ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000.000000.vtu" >}}

## Comparison with FEFlow solution

{{< img src="../compare.png" >}}
