+++
author = "Shuang Chen, Vinay Kumar"
project = ["Parabolic/ComponentTransport/ReactiveTransport/CO2Injection/calcite_pwp.prj"]
date = "2023-09-07T14:41:09+01:00"
title = "OGS-PHREEQC-Benchmark: CO2 injection into Opalinus Clay"
image = "results_comparision_sce1.png"
+++

{{< data-link >}}

## Overview

The present numerical work is motivated by the CO2LPIE project (shortened as CL-experiment) [1], which is an in-situ experiment that is being conducted at the Mont Terri
rock laboratory.
In the experiment, carbon dioxide (CO2) is injected into the Opalinus Clay leading to changes in its hydraulic, mechanical and chemical properties.
In general, these processes are of great interest in the evaluation of barrier intergrity.

Two scenarios are considered in this benchmark.
In the first scenario, the pure CO2 gas induced Calcite dissolution is simulated by OGS6-IPHREEQC and the results are verified with the related PHREEQC example presented by Appelo and Postma [2].
A comprehensive information regarding the computational procedure of OGS6-IPHREEQC can be found by the work from Lu et al. [3].
In the second scenario, the simulation considers a more accurate representation of chemical environments based on the CL Project.
This includes the incorporation of primary minerals typically found in Opalinus Clay and the relevant composition of species present in the porewater.

## Scenario #1: Kinetics of CO2 induced calcite dissolution

### Problem description

In this case, carbon dioxide ($CO_2$) gas with a partial pressure $10\^{-1.5}$ bar is injected into a fluid containing Calcite.
The temperature is maintained at 10 °C, the fluid’s pH value is set at 6 and the initial concentration of Calcite is 1 mol/kgw.
The dissolution of the CO2 gas in water together with the Calcite dissolution pathways can be described with the following reactions (Plummer et al. [4], Appelo and Postma [2]).

$$ \ce{CO2_{(g)} + H2O -> H2CO3^*} $$

$$ \ce{CaCO3 + H+ <=> Ca^2+ + HCO3-} $$

$$ \ce{CaCO3 + H2CO3^* <=> Ca^2+ + 2HCO3-} $$

$$ \ce{CaCO3 + H2O <=> Ca^2+ + HCO3- + OH-} $$

Two kinetic rates are adopted to describe the Calcite dissolution.
The first approach involves a simplified rate calculation derived directly from the current concentration of calcium,
employing the following formula: $rate = 10^{-6.91} - 10^{-1.52} * (c_{Ca})^2$.
And another one is the well known PWP rate which proposed by Plummer et al. [4].
The formulae of the PWP approach is available to be found in the PHREEQC database e.g. phreeqc.dat.
In the numerical experiment, the total simulation time is 30 000 s.
Details about the case study is described in the example 5.9 from Appelo and Postma [2].
The related PHREEQC script is available online to be found [5].

### Model and results

A simple 1D line-element model with one element is constructed in this work.
The coupled Hydraulic-Chemical (HC) Process is adopted for the OGS6-IPHREEQC simulation.
To match the PHREEQC example, the advection and diffusion have been set to zero in the modelling.
The initial conditions and chemical parameters are provided as outlined in the associated example.
The Fig. 1 depicts the comparison between the computed results obtained from the OGS6-IPHREEQC model and the results derived from PHREEQC, which shows a very
good agreement with each other.

{{< img src="results_comparision_sce1.png" width="200" title="pH and calcium increase with kinetic dissolution of Calcite." >}}

## Scenario #2: Modelling of the CO2 injection into Opalinus Clay

### model description

The identical 1D line-element model utilized in scenario #1 is employed for the OGS6-IPHREEQC simulation.
In this case, chemical environments based on the CL-experiment are considered.
For the porous medium in the OGS-6 model, in terms of the solid component, the following mineral composition is assumed in this work: 36% illite, 24% kaolinite, 7.5% calcite and 2.5% dolomite-dis (namely sedimentary (disordered) dolomite).
A more detailed description of the mineral composition of Opalinus Clay can be found in the work of Thury [6].
During the simulation, a constant carbon dioxide ($CO_2$) gas with a partial pressure of $10\^{1.5}$ bar is applied to the model domain as the equilibrium phase.
Consequently, the dissolved CO2 results in an alteration in the acidity of the porewater and leads to the different chemical reactions of each mineral.
The related reaction formulae are described in the PHREEQC database $llnl.dat$.
The adopted kinetic dissolution rate of each mineral are referenced from the transition state theory-based reaction mechanism following the work by Palandri and Kharaka [7]. The general equation formula reads

$${{Rate_\mathrm{mineral}}} = [ k_{\mathrm{acid_{mineral}}}^\mathrm{298.15K}{a_\mathrm{H^{+}}} + k_{\mathrm{neutral_{mineral}}}^\mathrm{298.15K} + k_{\mathrm{base_{mineral}}}^\mathrm{298.15K}{a_\mathrm{{H_{2}CO_{3}}^{*}}} ](1 - \mathrm{SR_{mineral}}),
\label{eq:transition_state}$$ with
$$\mathrm{SR_{mineral}} = \frac{\mathrm{IAP_{mineral}}}{K_\mathrm{eq,mineral}}.$$

where $a$ denotes the activity of the ion, and $\mathrm{SR}$ is the abbreviation for the saturation ratio of a phase, which describes the ion activity product $\mathrm{IAP}$ divided by equilibrium constant $K_\mathrm{eq}$.
In the simulation, the main species composition within the Opalinus Clay porewater are considered and the corresponding values are listed in table 1, following the
work of Wersin et al. [8].
In the model, the initial porosity is set to 0.15. Similarly to the scenario #1, hydraulic advection and diffusion are not considered in the simulation.

| Species | Value in [mol/kgw]      |
|---------|-------------------------|
| C(4)    | $3.89 \times 10^{-3}$   |
| Ca      | $1.89 \times 10^{-2}$   |
| Mg      | $2.197 \times 10^{-2}$  |
| Cl      | $3.2667 \times 10^{-1}$ |
| K       | $1.92 \times 10^{-3}$   |
| Na      | $2.8067 \times 10^{-1}$ |
| S(6)    | $1.679 \times 10^{-2}$  |
| Al      | $2.295 \times 10^{-10}$ |
| Si      | $1.7 \times 10^{-5}$    |
| Sr      | $4.6 \times 10^{-5}$    |

Table 1: concentration of the species in the Opalinus Clay porewater

For verification purposes, a corresponding PHREEQC model was constructed with identical parameter settings.

### Results

The left figure in Fig. 2 illustrates the evolution of calcium and magnesium molality during the CO2 injection computed from the OGS-6 and PHREEQC model.
The results of the two software programs match perfectly.
And the right figure in Fig. 2 illustrates the evolution of the dissoluted molality of each minerals during the simulation.
Only calcite and dolomite have been partially dissolved due to CO2 injection.
In contrast, clay minerals of illite and kaolinite were hardly affected.
The kinetic dissolution rate of calcite and dolomite is mostly controlled by their saturation ratio states.
When the SR value reaches 1, the dissolution process stops.

{{< img src="results_cl_comparision_sce2.png" width="200" title="Left: evolution of the calcium and magnesium molality during the CO2 injection; Right: evolution of the dissoluted minearals molality and the related saturation ratio of calcite and dolomite over the time." >}}

## Literature

<!-- vale off -->

[1] BGR, CO2LPIE project, 2023. URL: [https://www.bgr.bund.de/DE/Themen/Endlagerung/Projekte/Wirtsgesteine_geotechnische_Barrieren/laufend/Nur-Deutsch/mont_terri_experimente.html?nn=1542156](https://www.bgr.bund.de/DE/Themen/Endlagerung/Projekte/Wirtsgesteine_geotechnische_Barrieren/laufend/Nur-Deutsch/mont_terri_experimente.html?nn=1542156).

[2] C. A. J. Appelo, D. Postma, Geochemistry, groundwater and pollution, CRC press, 2004.

[3] R. Lu, T. Nagel, J. Poonoosamy, D. Naumov, T. Fischer, V. Montoya, O. Kolditz, H. Shao, A new operator-splitting finite element scheme for reactive transport modeling in saturated porous media, Computers and Geosciences 163 (2022) 105106. URL: https://doi.org/10.1016/j.cageo.2022.105106. doi:10.1016/j.cageo.2022.105106.

[4] L. N. Plummer, T. M. Wigley, D. L. Parkhurst, KINETICS OF CALCITE DISSOLUTION IN CO2-WATER SYSTEMS AT 5 degree TO 60 degree C AND 0. 0 TO 1. 0 ATM CO2., 1978.

[5] C. A. J. Appelo, D. Postma, Example 5.9: Kinetic dissolution of calcite, 2023. URL: https://www.hydrochemistry.eu/a&p/5/ex_5_9.phr.

[6] M. Thury, The characteristics of the Opalinus Clay investigated in the Mont Terri underground rock laboratory in Switzerland, 2002. doi:10.1016/S1631-0705(02)01372-5.

[7] J. L. Palandri, Y. K. Kharaka, A compilation of rate parameters of water-mineral interaction kinetics for application to geochemical modeling (2004).

[8] P. Wersin, M. Mazurek, T. Gimmi, Porewater chemistry of Opalinus Clay revisited: Findings from 25 years of data collection at the Mont Terri Rock Laboratory, Applied Geochemistry 138 (2022) 105234. URL: https://doi.org/10.1016/j.apgeochem.2022.105234. doi:10.1016/j.apgeochem.2022.105234.
