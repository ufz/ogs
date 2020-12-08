+++
author = "Jaime Garibay-Rodriguez, Renchao Lu, Vanessa Montoya"
weight = 142
project = "Parabolic/ComponentTransport/ReactiveTransport/SurfaceComplexation/RadionuclideSorption.prj"
date = "2019-07-08T14:41:09+01:00"
title = "Sorption of U(VI) in mineral phases"

[menu]
  [menu.benchmarks]
    parent = "Reactive Transport"

+++

{{< data-link >}}

## Overview and model setup

This benchmark is focused on the simulation of the migration of $U(VI)$ in a porous media. An important chemical process happening is the surface sorption, which favors the retardation of radionuclide migration. Therefore, it is necessary to quantify its overall impact on the transport.

We use a Surface Complexation approach as a mechanistic model to describe the chemical reactions happening at the solid surface. The surface reactions considered in this benchmark were added to the PSI/NAGRA chemical thermodynamic database version 12/07 (Thoenen *et al.*, 2014) for the speciation calculations. A schematic representation of the 1D model with relevant parameters is shown in Fig. 1

A permeability = $1.019368 \times 10^{-11}~m^2$, porosity = $0.2$, molecular diffusion coefficient = $1.0 \times 10^{-9}~m^2/s$, longitudinal/ horizontal/ vertical dispersivities = $0.2~m$, tortuosity = $1$, density = $1000~kg/m^3$ and viscosity = $1\times 10^{-3}~kg~m^{-1} s^{-1}$ are the input parameters of the porous medium.

The benchmark uses the `ComponentTransport` process in OGS-6 (see [HC-Process.pdf](/docs/benchmarks/hydro-component/HC-Process.pdf)) coupled with the IPhreeqc interface for the chemical speciation calculations. The porewater initial composition and injected for the first 10 000 s is shown in Table 1.

{{< img src="../domain.png" title="Spatial and temporal discretization of the 1D model. Solution concentrations with/without U(VI) are applied at the inflow boundary. Initial concentration of U(VI) in the domain is 0.">}}

-----------------------------------------

|Parameter | Unit | Initial solution | Inlet solution |
|:-------- | :------ | :---- | :--- |
| $Na^+$   |  mol/L  | $1.0 \times 10^{-3}$  | $1.0 \times 10^{-3}$ |
| $Cl^-$ | mol/L | $60.0 \times 10^{-3}$ | $60.0 \times 10^{-3}$ |
| $Ca^{2+}$    | mol/L | $25.02 \times 10^{-3}$ | $25.02 \times 10^{-3}$ |
| $DIC$  | mol/L | $1.0 \times 10^{-3}$ | $1.0 \times 10^{-3}$ |
| $U(VI)$   | mol/L | 0  | $1 \times 10^{-8}$ |
-----------------------------------------

Table 1: **Porewater solution compositions.**

-----------------------------------------

## Results

To approximate the results obtained with the ESTRAL database, we use thermodynamic data from the **RES³T - Rossendorf Expert System for Surface and Sorption Thermodynamics** (https://www.hzdr.de/db/RES3T.queryData) for the available mineral phases (see Table 2). Note that the mineral groups of the ESTRAL database contains a larger number of individual phases. Because the composition of each mineral group is not provided for the surface characterization of Table 2, we choose only one representative phase at a time for each mineral group in our simulations. Further, to keep the chemical system simple, the chosen chemical reactions and *log K* values for the surface complexation model are only those of mineral-OH sites, with the exception of gibbsite, where only Aluminol sites are available for $U(VI)$. In total, 19 chemical reactions are considered in the surface complexation model. All *log K* values are reported for 298.15 K and normalized to 2.31 sites/nm$^2$.  

-----------------------------------------

|Mineral group/phase | Representative phases | Reference binding sites (sites/nm$^2$) | Specific surface area (m$^2$/g) |  Solid mass (g/kgw) |
|:-------- | :------ | :---- | :--- | :--- |
| Quartz   |  -  | 2.31  | 0.007 | 9010 |
| Feldspars | albite, orthoclase | 2.31 | 0.21 | 1060 |
| Mica    | muscovite | 2.31     | 1.72 | 53 |
| Fe(III)-oxids/-hydroxids  | goethite, hematite| 2.31   | 0.26 | 53 |
| Al-hydroxids   | gibbsite   | 2.31  | 0.11 | 53 |
| 2-layer-clay minerals | kaolinite | 2.31 | 0.07 | 159 |
-----------------------------------------

Table 2: **Surface parameters and characterization used in the simulations.**

-----------------------------------------

Four different combinations can be simulated taking the albite and orthoclase phases of the Feldspars group and the goethite and hematite phases for the Fe(III)-oxids/-hydroxids group. Mineral combinations from 1 to 4 (see Fig. 2) is as follows: 1) albite-goethite, 2) albite-hematite, 3) orthoclase-goethite and 4) orthoclase-hematite. From the concentration profiles in Fig. 2, it is clear that the combination 2 approximates better the profile obtained with the ESTRAL database. We choose this combination for the next part of our simulations. Furthermore, this combination is written in the `RadionuclideSorption.prj` file of this benchmark.

{{< img src="../Fig1.png" title="Comparison of concentration profiles at final simulation time (115 000 s) for various representative minerals of the Feldspar and Fe(III)-oxids/-hydroxids groups. The mineral combinations profiles are obtained using the PSI/Nagra database version 12/07 and the dashed profile is obtained with the ESTRAL database.">}}

The temporal evolution of the concentration profiles of the chosen mineral combination (albite-hematite) compared to the ESTRAL database is shown in Fig. 3. In addition, a simulation of the reactive transport treating $U(VI)$ as a non-sorbing radionuclide is presented. Recall that the contaminated solution with $U(VI)$ is injected only for the first 10 000 s of simulation. On the one hand, we note the difference between the profile with the augmented PSI/Nagra database and with the ESTRAL database. This is expected, since different reactions and significant differences in *log K* values are considered for each simulation. However, we note that the trend is similar enough to capture the relevant sorption process happening at the surface.

On the other hand, the enormous difference between sorbing and non-sorbing reactive transport is evident from the resulting concentration profiles. Therefore, we highlight the importance of considering the impact of sorption in the transport of radionuclides, as this is paramount for the safety assessment in the design of nuclear waste repositories. Finally, the CPU time of the simulation taking into account surface complexation is roughly double of the simulation with only aqueous speciation. This posses the necessity of choosing a good compromise between accuracy (large number of reactions and chemical parameters) and performance.

{{< img src="../Fig2.gif" title="Time evolution of mineral combination 2 (albite/hematite) in comparison to the results obtained with the ESTRAL database. The green dotted line shows the temporal evolution of U(IV) as a non-sorbing radionuclide.">}}

{{< data-link >}}

## References

Thoenen, T., Hummel, W., Berner, U., & Curti, E. (2014). *The PSI/Nagra Chemical Thermodynamic Database 12/07*.