+++
author = "Boyan Meng and Haibing Shao"
date = "2021-11-2T18:52:00+01:00"
title = "Overview of Multiphase Flow Processes (without mechanics)"
weight = 43

[menu]
  [menu.userguide]
    parent = "process-dependent-configuration"

+++


This documentation gives an overview of the multiphase flow processes (except `TH2M`) in OGS6. Currently all processes assume two-phase flow (The `Richards Flow` process is not considered as a "two-phase flow" process since the gas phase is assumed static). A comparison of the process-specific features is listed in the following table.

| Feature | `TwoPhaseFlowWithPP` | `TwoPhaseFlowWithPrho` | `ThermalTwoPhaseFlowWithPP` |
|---|---|---|---|
| \# of phases | 2 | 2 | 2 |
| \# of components | 2 | 2 | 2 |
| Phase composition | Pure | Compositional (only liquid phase) | Compositional (only gas phase) |
| Primary Variables | $P_g, P_c$ | $P_g, \rho^h_{tot}$ | $P_g, P_c, T$ |
| Temperature effect | Isothermal  | Isothermal | Non-isothermal |
| Phase change | N.A. | Gas dissolution/degassing | Evaporation/condensation |
| Phase appearance/disappearance | N.A. | Yes, only gas phase | Yes, only liquid phase |

Nomenclature: $P_g$: gas pressure; $P_c$: capillary pressure; $T$: temperature; $\rho^h_{tot}$: total density of the light component.

Some remarks:

1\. The `TwoPhaseFlowWithPP` process assumes that the two fluid phases are immiscible. Thus, it is most suitable for simulating two-phase flow under capillary effects (e.g. replacement of one phase by another due to gravity). Note that the wetting and non-wetting phases are not limited to water and gas, see the [McWhorter benchmark]({{< ref "../../benchmarks/two-phase-flow-pp-form/two-phase-flow-pp-mcwhorter.md" >}}) for example.

2\. The `TwoPhaseFlowWithPrho` process assumes that the main component of the gas phase can be dissolved in the liquid phase. Water evaporation is neglected here. The appearance/disappearance of the gas phase is controlled by the solubility (given by the Henry's Law) of the gaseous component, e.g. H2. It is therefore most suitable for nuclear waste repository (see the [MoMaS benchmark]({{< ref "../../benchmarks/two-phase-flow/momas.md" >}})) or CO2 storage problems.

3\. The `ThermalTwoPhaseFlowWithPP` process simulates the temperature-dependent two-phase flow and moisture transport. Water evaporation and recondensation can be modeled thanks to that the gas phase is compositional. This process is most favorably used for shallow geothermal applications (e.g. borehole thermal energy storage), especially in unsaturated soils (see the [heat pipe benchmark]({{< ref "../../benchmarks/thermal-two-phase-flow/heat-pipe.md" >}})). 