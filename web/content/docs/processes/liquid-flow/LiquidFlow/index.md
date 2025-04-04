+++
title = "Liquid Flow (groundwater flow)"
author = "Tom Fischer, Philipp Selzer"
weight = 3
+++

In the following pdf-document you find a description of the major governing equations of the Liquid-Flow-Process implemented in OpenGeoSys (OGS) as well as derivations and hints how to set the parameterization of the process correctly:
[documentation](https://gitlab.opengeosys.org/ogs/documentation/liquidflow/-/jobs/artifacts/main/raw/main.pdf\?job\=build).

The Liquid-Flow-Process models saturated single-phase flow, which may be of variable density, in porous and fractured media. It may be used for modelling water flow but can be equally used to model flow of other liquids and also gas as long as flow follows the Darcy regime.
However, the natural application of the Liquid-Flow-Process is in modelling groundwater flow.
As such, the documentation above targets especially towards hydrogeologists.
In the context of groundwater modelling, it is important to note that the pressure-based variant of the groundwater flow equation is implemented in OGS, which may not be the most intuitive for hydrogeologists.
To this end, the documentation above aims, among others, in giving conversions between the pressure-based and the hydraulic-head based formulations of the groundwater flow equation as well as it outlines how to use the Liquid-Flow-Process such that it models natively the hydraulic-head based formulation of the groundwater flow equation commonly used in hydrogeology.
