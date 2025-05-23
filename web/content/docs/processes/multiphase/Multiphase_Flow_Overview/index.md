+++
author = "Boyan Meng and Haibing Shao"
date = "2021-11-02T18:52:00+01:00"
title = "Overview of Multiphase Flow Processes (without mechanics)"
weight = 3
+++


This documentation gives an overview of the multiphase flow processes (except `TH2M`) in OGS. Currently all processes assume two-phase flow (The `RICHARDS_FLOW` process is not considered as a "two-phase flow" process since the gas phase is assumed static). A comparison of the process-specific features is listed in the following table.

| Feature | `TwoPhaseFlowWithPP` | `ThermalTwoPhaseFlowWithPP` |
|---|---|---|
| \# of phases | 2 | 2 |
| \# of components | 2 | 2 |
| Phase composition | Pure |  Compositional (only gas phase) |
| Primary Variables | $P_g, P_c$ | $P_g, P_c, T$ |
| Temperature effect | Isothermal  | Non-isothermal |
| Phase change | N.A. | Evaporation/condensation |
| Phase appearance/disappearance | N.A. | Yes, only liquid phase |

Nomenclature: $P_g$: gas pressure; $P_c$: capillary pressure; $T$: temperature.

Some remarks:

1. The `TwoPhaseFlowWithPP` process assumes that the two fluid phases are immiscible. Thus, it is most suitable for simulating two-phase flow under capillary effects (e.g. replacement of one phase by another due to gravity). Note that the wetting and non-wetting phases are not limited to water and gas, see the [McWhorter benchmark]({{< ref "docs/benchmarks/two-phase-flow/two-phase-flow-pp-mcwhorter" >}}) for example.
2. The `ThermalTwoPhaseFlowWithPP` process simulates the temperature-dependent two-phase flow and moisture transport. Water evaporation and recondensation can be modeled thanks to that the gas phase is compositional. This process is most favorably used for shallow geothermal applications (e.g. borehole thermal energy storage), especially in unsaturated soils (see the [heat pipe benchmark](https://www.opengeosys.org/docs/benchmarks/thermal-two-phase-flow/heatpipe/)).
