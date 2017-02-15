+++
date = "2017-02-15T14:46:38+01:00"
title = "Ehlers; special case - Drucker-Prager"
weight = 115
project = "Mechanics/Ehlers/cube_1e0_dp.prj"
author = "Xing-Yuan Miao"

[menu]
  [menu.benchmarks]
    parent = "small-deformations"

+++

{{< project-link >}}

## Problem description

The Ehlers material model can be reduced to the well-known criteria, such as the Drucker-Prager and von Mises. Here, we provide a Drucker-Prager model test as an additional verification of the Ehlers material model. The interested reader is referred to the document attached in the standard Ehlers material model test for obtaining more details.

## Results and evaluation

Triaxial compression test:

{{< img src="../ss_load.png" >}}

Variations of the stress states and the plastic volumetric strain with a monotonic loading process:

{{< img src="../dp_test.png" >}}
