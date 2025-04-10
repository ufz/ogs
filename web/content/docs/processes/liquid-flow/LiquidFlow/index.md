+++
title = "Liquid Flow (groundwater flow)"
author = "Tom Fischer, Philipp Selzer"
weight = 3
+++

In the following pdf-document you find a description of the major governing equations of the Liquid-Flow-Process implemented in OpenGeoSys (OGS) as well as derivations and hints on how to set the parameterization of the process correctly:
[documentation](https://gitlab.opengeosys.org/ogs/documentation/liquidflow/-/jobs/artifacts/main/raw/main.pdf\?job\=build).

The Liquid-Flow-Process models saturated single-phase flow, which can be of variable density, in porous and fractured media. It may be used for modeling water flow but can be equally used to model flow of other liquids and also gas as long as flow follows the Darcy regime.
However, the natural application of the Liquid-Flow-Process is in modeling groundwater flow.
As such, the documentation above specifically targets hydrogeologists.
In the context of groundwater modeling, it is important to note that the pressure-based variant of the groundwater flow equation is implemented in OGS, which may not be the most intuitive for hydrogeologists.
To this end, the documentation above aims, among other things, to give conversions between the pressure-based and the hydraulic-head based formulations of the groundwater flow equation as well as it outlines how to use the Liquid-Flow-Process to natively model the hydraulic head-based formulation of the groundwater flow equation commonly used in hydrogeology.
