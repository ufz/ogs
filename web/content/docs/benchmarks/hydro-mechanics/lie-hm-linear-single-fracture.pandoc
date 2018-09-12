+++
date = "2017-02-16T15:03:28+01:00"
title = "Linear, single fracture"
weight = 152
project = "LIE/HydroMechanics/single_fracture.prj"
author = "Norihiro Watanabe"

[menu]

  [menu.benchmarks]
    parent = "hydro-mechanics"

+++

{{< data-link >}}

## Equations

We solve a hydromechanics problem (small deformation, linear elastic, Darcy flow) with a pre-existing fracture using the lower-dimensional interface element (LIE) approach.

See [this PDF](../LIE_HM.pdf) for detailed problem description.


## Results and evaluation

Result showing pore pressure increase due to the injection and subsequent fracture aperture increases at t=500s. A small discrepancy of the aperture near the injection is due to the interpolation method for converting cell data to point data in ParaView.

{{< img src="../x_p_t500.png" >}}
{{< img src="../x_b_t500.png" >}}

## 3D setup

Same setup as for the given 2D case above with additional plane strain
conditions in the front and back x-y planes.
Warp of the 1000-times oversized displacement and the fracture's aperture are
shown in the following figure.
{{< img src="../single_fracture_3D.png" >}}

Comparison with 2D setup yields identical results (up to numerical differences
in order of 1e-15):
TODO: Image missing!

<!-- {{< img src="../single_fracture_3D_vs_2D.png" >}} -->
