+++
date = "2017-02-16T14:59:01+01:00"
title = "Linear"
weight = 151
project = "HydroMechanics/Linear/Confined_Compression/square_1e2.prj"
author = "Dmitri Naumov"

[menu]
  [menu.benchmarks]
    parent = "hydro-mechanics"

+++

{{< data-link >}}

## Problem description

We solve a hydro-mechanical linear biphasic model (small deformation, linear elastic, Darcy flow, incompressible constituents) in square domain where on the top boundary a displacement boundary condition is applied. The boundary condition is a linear ramp from 0 to the final displacement of 0.05m at time t = 100s. After that the top is held fast in the final position. The fluid is allowed to escape through the top, which regulates the pressure in the domain.

## Results and evaluation

The analytical solution of the problem can be found in Mov _et.al._ (1980)"Biphasic Creep and Stress Relaxation of Articular Cartilage in Compression: Theory and Experiments."

The result showing initial pore pressure increase in the bottom of the domain due to applied displacement on the top and subsequent pressure drop after the displacement ramp finishes at time 100s:

{{< img src="../HM_confined_compression_analytical.png" >}}

Comparison with the numerical solution shows good agreement with the analytical solution:

{{< img src="../HM_confined_compression_simulation_error.png" >}}
