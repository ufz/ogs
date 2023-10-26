+++
author = "Francesco Parisio"
title = "Pressure boundary conditions"
project = ["Mechanics/Linear/PressureBC/hollow_sphere.prj"]
weight = 117
image = "pipe_axisymmetric.png"
+++

## Problem description

Five different benchmarks are reported in which pressure boundary conditions are tested: an axisymmetric elastic cylinder, a plain strain elastic cylinder, an axisymmetric elastic sphere, a tri-dimensional elastic sphere and an axisymmetric elasto-plastic sphere.
See [this PDF](pressure_bc.pdf) for detailed problem description.

## Results and evaluation

Plain strain elastic cylinder comparison between numerical and analytical results:
{{< figure src="pipe_plane_strain.png" >}}

Axisymmetric elastic cylinder comparison between numerical and analytical results.
{{< figure src="pipe_axisymmetric.png" >}}

Axisymmetric elastic sphere comparison between numerical and analytical results:
{{< figure src="sphere_axisymmetric.png" >}}

Tri-dimensional elastic sphere comparison between numerical and analytical results:
{{< figure src="sphere_3d.png" >}}

Axisymmetric plastic sphere comparison between numerical and analytical results:
{{< figure src="sphere_axisymmetric_pl.png" >}}

Axisymmetric plastic sphere residuals of stress:
{{< figure src="sphere_axisymmetric_pl_residual_stress.png" >}}
