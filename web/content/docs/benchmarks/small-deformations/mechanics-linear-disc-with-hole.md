+++
date = "2017-02-15T12:12:18+01:00"
title = "Linear; Disc with hole"
project = "Mechanics/Linear/disc_with_hole.prj"
author = "Dmitri Naumov"
weight = 111

[menu]
  [menu.benchmarks]
    parent = "small-deformations"

+++

{{< project-link >}}

## Equations

We solve a linear elastic small deformation problem on a quarter of a plate with hole put under tension on the top boundary.

See the [Circular_hole.pdf](https://docs.opengeosys.org/assets/files/Documentation/Selected-Benchmarks/Circular_hole.pdf) for detailed problem description.


## Results and evaluation

Result showing the displacement field and ten-fold overscaled distortion:

{{< vis path="Mechanics/Linear/disc_with_hole_pcs_0_ts_4_t_1.000000.vtu" >}}
