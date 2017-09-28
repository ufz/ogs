+++
author = "Marc Walther"
weight = 142
project = "Parabolic/ComponentTransport/goswami/"
date = "2017-09-07T14:41:09+01:00"
title = "Saturated Variable-Density Flow and Mass Transport (Goswami)"

[menu]
  [menu.benchmarks]
    parent = "hydro-component"

+++

{{< project-link >}}


## Overview

This benchmark features an Henry-like variable-density flow and transport setup that evaluates saltwater intrusion on laboratory scale. The original benchmark was published by Goswami & Clement (2007); OpenGeoSys5 has been verified for all steady and transient states of the original benchmark (see Kolditz et al., 2012).

For the setup and parameterization, see the chapter "Density dependent flow - The Goswami Problem" in Kolditz et al. (2012).


## Problem description

The Goswami-Clement benchmark is based on experiment observations for intruding and receding saltwater in a laboratory-scale sand tank. Here, we compare numerical results of ogs6 to those of OGS5 and the original observation data.

We had to extend the modeling domain to include the required boundary conditions on the left and right side.

### Model results

A comparison of numerical and laboratory data is shown in the figure below. The numerical results of ogs6 coincide with those of OGS5 and likewise with the laboratory observations.

{{< img src="../goswami.png" title="Results for numerical (OGS5 - green, ogs6 - white) and laboratory data (black squares) together with concentration distribution in the domain and mesh resolution for steady state 1 (see original research paper).">}}

[The project files are here. ](../../../../../Tests/Data/Parabolic/ComponentTransport/gosami)

## Literature

Goswami, R.R., Clement, T.P., 2007. Laboratory-scale investigation of saltwater intrusion dynamics. Water Resour. Res. 43, 1–11. doi:10.1029/2006WR005151.

Kolditz, O., Görke, U.-J., Shao, H., Wang, W., 2012. Thermo-Hydro-Mechanical-Chemical Processes in Porous Media: Benchmarks and Examples, Lecture notes in computational science and engineering. Springer. ISBN: 3642271766.
