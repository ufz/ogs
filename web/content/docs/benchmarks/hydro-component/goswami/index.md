+++
author = "Jasper Bathmann, Marc Walther"
weight = 142
project = "Parabolic/ComponentTransport/goswami/"
date = "2017-09-07T14:41:09+01:00"
title = "Saturated Variable-Density Flow and Mass Transport (Goswami)"

[menu]
  [menu.benchmarks]
    parent = "Hydro-Component"

+++

{{< data-link >}}

## Overview

This benchmark features an Henry-like variable-density flow and transport setup that evaluates saltwater intrusion on laboratory scale. The original benchmark was published by Goswami & Clement (2007); OpenGeoSys5 has been verified for all steady and transient states of the original benchmark (see Kolditz et al., 2012).

For the setup and parameterization, see the chapter "Density dependent flow - The Goswami Problem" in Kolditz et al. (2012).

## Problem description

The Goswami-Clement benchmark is based on experiment observations for intruding and receding saltwater in a laboratory-scale sand tank. Here, we compare numerical results of ogs6 to the original observation data.

### Model results

An example for the intruding salt front is shown below. The numerical results of OGS-6 coincide with those of OGS-5.

{{< img src="goswami.gif" title="Results for numerical experiment. The steady state SS2 from the original experimental work is well reproduced.">}}

{{< data-link >}}

A comparison of numerical and laboratory data is shown in the figure below. The numerical results of ogs6 coincide with those of OGS5 and likewise with the laboratory observations.

{{< img src="Goswami_Exp_Num_Comp.png" title="Results for numerical (colored diamonds) and laboratory data (colored straight lines) on the steady state location of the concentration front (see original research paper).">}}

{{< img src="Goswami_Transient_States.png" title="Results for numerical (colored diamonds) and laboratory data (colored straight lines) on the transient state locations of the concentration front (see original research paper).">}}

## Literature

Goswami, R.R., Clement, T.P., 2007. Laboratory-scale investigation of saltwater intrusion dynamics. Water Resour. Res. 43, 1–11. doi:10.1029/2006WR005151.

Kolditz, O., Görke, U.-J., Shao, H., Wang, W., 2012. Thermo-Hydro-Mechanical-Chemical Processes in Porous Media: Benchmarks and Examples, Lecture notes in computational science and engineering. Springer. ISBN: 3642271766.
