+++
author = "Christian Silbermann, Thomas Nagel"
date = "2022-05-19"
title = "Heat conduction with phase change"
project = ["Parabolic/T/1D_freezing_column_Stefan/Stefan_problem.prj"]
image = ""
+++

## Test cases

Three tests are presented:

- {{< data-link "1D homogeneous" "Parabolic/T/1D_freezing_column_Stefan/Stefan_problem_homogen.prj" >}},
- {{< data-link "1D heterogeneous" "Parabolic/T/1D_freezing_column_Stefan/Stefan_problem.prj" >}}, and
- {{< data-link "2D heterogeneous" "Parabolic/T/2D_freezing_disk/circle_disk.prj" >}}.

## Problem description

We simulate the melting of a 1D ice column including the latent heat necessary for this phase change. No solid (porous medium) is involved and the problem can either be homogeneous (volumetric heat source) or heterogeneous (boundary condition as heat input).
Further, the heating of fully saturated porous medium with the shape of a disk is simulated including the phase change of pore water to pore ice.

See [this PDF](Heatconduction_with_phase_change.pdf) for the detailed description.
