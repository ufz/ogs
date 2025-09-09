+++
author = "Dmitri Naumov"
weight = 113
project = ["Mechanics/Ehlers/cube_1e3.prj"]
date = "2017-02-15T14:39:39+01:00"
title = "Ehlers; Single-surface yield function"
image = "plasticity_ss.png"
models = [ "lab" ]
+++

{{< data-link >}}

## Problem description

We use a seven-parametric yield function for geomaterials to describe the plastic response. The triaxial compression test is setup.

See [this PDF](Plasticity.pdf) for detailed problem description.

## Results and evaluation

Triaxial compression test:

{{< figure src="ss_load.png" >}}

Variations of the stress states and the plastic volumetric strain with a monotonic loading process:

{{< figure src="plasticity_ss.png" >}}
