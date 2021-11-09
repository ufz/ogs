+++
author = "Boyan Meng, Norbert Grunwald and Haibing Shao"
date = "2021-11-2T18:52:00+01:00"
title = "Overview of Multiphase Flow Processes"
weight = 43

[menu]
  [menu.userguide]
    parent = "process-dependent-configuration"

+++


As of November 2021, there are four implemented processes in OGS6 that are able to simulate the multiphase flow and transport phenomena in porous media. An overview of the process-specific features is given in the following table. Currently all processes assume two-phase flow.

| Feature | `TwoPhaseFlowWithPP` | `TwoPhaseFlowWithPrho` | `ThermalTwoPhaseFlowWithPP` | `TH2M` |
|---|---|---|---|---|
| \# of phases | 2 | 2 | 2 | 2 |
| \# of components | 2 | 2 | 2 | 2 |
| Phase composition | Pure | Compositional (only liquid phase) | Compositional (only gas phase) | Flexible |
| Primary Variables | $P_g, P_c$ | $P_g, \rho^h_{tot}$ | $P_g, P_c, T$ | $P_g, P_c, T, \sigma$ |
| Temperature effect | Isothermal  | Isothermal | Non-isothermal | Non-isothermal |
| Phase change | N.A. | Gas dissolution/degassing | Evaporation/condensation | Evaporation, etc. (WIP) |
| Phase appearance/disappearance | N.A. | Yes, only gas phase | Yes, only liquid phase | ?? (TODO: Norbert) |

Nomenclature: $P_g$: gas pressure; $P_c$: capillary pressure; $T$: temperature; $\sigma$: displacement; $\rho^h_{tot}$: total density of the light component.

Some remarks:
1. The `TwoPhaseFlowWithPP` process assumes that the two fluid phases are immiscible. Thus, it is most suitable for simulating two-phase flow under capillary effects (e.g. replacement of one phase by another due to gravity). Note that the wetting and non-wetting phases are not limited to water and gas, see the [McWhorter benchmark]({{< ref "../../benchmarks/two-phase-flow-pp-form/two-phase-flow-pp-mcwhorter.md" >}}) for example.

2. The `TwoPhaseFlowWithPrho` process assumes that the main component of the gas phase can be dissolved in the liquid phase. Water evaporation is neglected here. The appearance/disappearance of the gas phase is controlled by the solubility (given by the Henry's Law) of the gaseous component, e.g. H2. It is therefore most suitable for nuclear waste repository (see the [MoMaS benchmark]({{< ref "../../benchmarks/two-phase-flow/momas.md" >}})) or CO2 storage problems.

3. The `ThermalTwoPhaseFlowWithPP` process simulates the temperature-dependent two-phase flow and moisture transport. Water evaporation and recondensation can be modeled thanks to that the gas phase is compositional. This process is most favorably used for shallow geothermal applications (e.g. borehole thermal energy storage), especially in unsaturated soils (see the [heat pipe benchmark]({{< ref "../../benchmarks/thermal-two-phase-flow/heat-pipe.md" >}})). 

4. The `TH2M` process is very flexible. It can be used for the coupled effect of two-phase flow with mechanics or as a proxy for two-phase flow simulations without considering the mechanical or temperature effects (see the [TH2M bechmarks](<https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/Tests/Data/TH2M>). Currently the phase transition model of `TH2M` process is being extended to include more functionalities. 