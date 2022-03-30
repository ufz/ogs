+++
project = "Parabolic/HT/HeatTransportInStationaryFlow/HeatTransportInStationaryFlow.prj"
author = "Wenqing Wang"
title = "One dimensional heat transport in stationary flow"
date = 2020-12-14T08:00:32+01:00
weight = 153

[menu]
  [menu.benchmarks]
     parent = "hydro-thermal"
+++

## Problem description
We consider one dimensional heat transport in stationary flow in a porous medium.
 This benchmark was first introduced as an exercise of
 [Geoenergy Modeling I – Geothermal Processes in Fractured Porous Media](https://www.opengeosys.org/books/geoenergy-modeling-i/)
 for [OGS 5](https://github.com/ufz/ogs5).

## Numerical setting
 The size of the domain is 1 m in the horizontal direction.
  The material properties of fluid are:

| Property               | Value | Unit        |
|------------------------|-------|-------------|
| Density                | 1000  | kg/m^3      |
| Viscosity              | 1.e-3 | Pa⋅s        |
| Specific heat capacity | 4182  | J/K/kg      |
| Thermal conductivity   | 0.6   | W/(m⋅K)     |

The material properties of porous medium are:

| Property               | Value | Unit        |
|------------------------|-------|-------------|
| Density                | 2850  | kg/m^3      |
| Specific heat capacity | 0  | J/K/kg      |
| Thermal conductivity   | 0   | W/(m⋅K)     |

The intrinsic permeability is 1.e-11 m^2. Since the flow equation is steady state,
 the porosity $n$ is not applied in that equation. We set $n=1$
 in order to use the specific heat capacity and the thermal conductivity of fluid
 as the effective ones.

The initial boundary conditions are <em>T(0)=</em>0 °C and
  <em>p(0)=</em>1.e+5 Pa.

At the left boundary, a constant temperature of <em>T=</em>1 °C
 and a constant pressure of <em>p=</em>1.01e+5 Pa are prescribed.
 At the right boundary condition, there is no heat flux and the pressure is
 set as the initial one. The pressure boundary conditions lead to a
 stationary flow with a velocity of 1.e-5 m/s.

The time duration is 5.e+4 s. A fixed time step size of 250 s is used
 for the temporal discretization.

This example is also set as one of the benchmarks of
 [ThermoHydroMechanics (THM)](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/ThermoHydroMechanics/Linear/HeatTransportInStationaryFlow/HeatTransportInStationaryFlow.prj)
 and ThermoRichardsMechanics (TRM), respectively.
 In order to provide a reference result for the same benchmark of
 THM and TRM, a 2D domain of 1 m $\times$ 0.1 m is used, which
 is discretized into 3$\times$39 quadrilateral elements.

## Result

The temperature distribution  at <em>t=</em>  5.e+4 s together with the mesh is
 illustrated in the following figure:
{{< img src="heat_transport_in_stationary_flow_domain.png" >}}

The temperature  profile at <em>t=</em>  5.e+4 s along a horizontal line in the
  2D domain is given in the following figure:
{{< img src="heat_transport_in_stationary_flow_profile.png" >}}

