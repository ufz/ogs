+++
author = "Christian Silbermann, Thomas Nagel"
project = ["Mechanics/ModifiedCamClay/square_1e0_shear.prj",
           "Mechanics/ModifiedCamClay/square_1e0_biax.prj",
           "Mechanics/ModifiedCamClay/triaxtest.prj",
           "Mechanics/ModifiedCamClay/triaxtest_original.prj",
           "Mechanics/ModifiedCamClay/triaxtest_original_abs.prj"]
date = "2020-12-14T14:39:39+01:00"
title = "Modified Cam clay model"
image = ""
+++

## Test cases

Five tests are presented:

{{< data-link >}}

of which the last three have the same test program but use different implementations of the modified Cam clay model.
The mfront-files can be found at [here](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/MaterialLib/SolidModels/MFront).

## Problem description

We perform plane strain and axisymmetric mechanical tests using
the modified Cam clay model revealing both the features and
the limitations of this material model.
See [this PDF](ModifiedCamClay_report.pdf) for the detailed description.
