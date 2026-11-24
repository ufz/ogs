+++
date = "2024-12-13T12:00:00+01:00"
title = "Introduction"
author = "Maximilian Dörnbrack"
weight = 1


[menu.docs]
name = "Process information"
identifier = "processes"
weight = 5
post = "Get more insight into process implementation."
[menu.docs.params]
category = "User"
+++

<div class="flex flex-col-reverse gap-5 items-center md:flex-row">
<div class="flex-none w-1/3 min-w-[200px] not-prose mr-4">
{{< figure src="/images/coupling-icons/thmc-tet-externals.svg" alt="THMC processes diagram" class="mt-0 w-full" img-class="w-full" >}}
</div>
<div class="flex-1">
Coupled THMC processes form the core of OpenGeoSys's function as a multi-physics/chemistry code. The following pages provide examples and workflows for setting up OGS processes from scratch. For complex chemistry and material constitutive modelling, OGS can be coupled with PHREEQC and MFront.
</div>
</div>

The overview is organised according to the fundamental processes of the THMC concept. Triple-coupled THM processes are differentiated by various flow processes. We also provide an elliptic equation solver for steady-state / equilibrium processes. The short process descriptions are organised as follows: introduction and key features; model equations; process specification in the prj project file; and main parameters. Workflows for model setup and links to the benchmark gallery will help you start and customise your own OGS simulations.

- Thermal processes
  - [Heat transport]({{< relref "/docs/processes/heat-conduction/HeatConduction" >}})
  - Borehole heat exchanger
- Hydraulic processes
  - [Groundwater flow]({{< relref "/docs/processes/liquid-flow/LiquidFlow" >}})
  - Gas flow
  - [Richards flow]({{< relref "/docs/processes/richards-flow/RichardsFlow" >}})
  - [Multiphase flow]({{< relref "/docs/processes/multiphase/Multiphase_Flow_Overview" >}})
- Mechanical processes
  - [Small deformations]({{< relref "/docs/processes/small-deformation/SmallDeformation" >}})
  - [Large deformations]({{< relref "/docs/processes/large-deformation/f-bar" >}})
- Chemical processes
  - [Component transport]({{< relref "/docs/processes/component-transport/hydro-component" >}})
  - Reactive transport
- THM processes
  - [THM: Thermo-Hydro-Mechanics]({{< relref "/docs/processes/thermal-processes/THM" >}})
  - [TRM: Thermo-Richards-Mechanics]({{< relref "/docs/processes/thermal-processes/TRM" >}})
  - TH2M: Two-phase flow in THM
- [Steady-state: Elliptic equation]({{< relref "/docs/processes/steady-state-diffusion/SteadyStateDiffusion" >}})
