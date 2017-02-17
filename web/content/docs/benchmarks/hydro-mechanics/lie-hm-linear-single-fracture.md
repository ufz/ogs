+++
date = "2017-02-16T15:03:28+01:00"
title = "Linear, single fracture"
weight = 152
project = "https://github.com/ufz/ogs-data/blob/master/LIE/HydroMechanics/single_fracture.prj"
author = "Norihiro Watanabe"

[menu]

  [menu.benchmarks]
    parent = "hydro-mechanics"

+++

{{< project-link >}}

## Equations

We solve a hydromechanics problem (small deformation, linear elastic, Darcy flow) with a pre-existing fracture using the lower-dimensional interface element (LIE) approach.

See the [LIE_HM.pdf](https://docs.opengeosys.org/assets/files/Documentation/Selected-Benchmarks/LIE_HM.pdf) for detailed problem description.


## Results and evaluation

Result showing pore pressure increase due to the injection and subsequent fracture aperture increases at t=500s. A small discrepancy of the aperture near the injection is due to the interpolation method for converting cell data to point data in ParaView.

{{< img src="../x_p_t500.png" >}}
{{< img src="../x_b_t500.png" >}}
