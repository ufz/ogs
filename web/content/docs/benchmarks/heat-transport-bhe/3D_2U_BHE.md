+++
author = "Chaofan Chen, Haibing Shao"
date = "2018-09-25T13:44:00+01:00"
title = "3D 2U BHE"
weight = 123
project = "Parabolic/T/3D_2U_BHE/3D_2U_BHE.prj"

[menu]
  [menu.benchmarks]
    parent = "heat-transport-bhe"

+++

{{< data-link >}}

## Model Setup

The numerical model was established based on dual continuum method developed by Diersch et al. (2011). In the 3D 2U type BHE benchmark, the depth of the whole domain is 18.5 m with a square cross section of 5 m per side. The BHE bottom is situated at -18 m and its total length is 18 m. Detailed parameters for the model configuration are listed in the follwoing table.

| Parameter                          | Symbol             |  Value              | Unit                        |
| ---------------------------------- |:------------------ | -------------------:| --------------------------: |
| Soil thermal conductivity          | $\lambda_{s}$      | 2.78018             | $\mathrm{W m^{-1} K^{-1}}$  |
| Soil heat capacity                 | $(\rho c)_{soil}$  | $3.2\times10^{6}$   | $\mathrm{Jm^{-3}K^{-1}}$    |
| Length of the BHE                  | $L$                | 18                  | $\mathrm{m}$                |
| Diameter of the BHE                | $D$                | 0.13                | $\mathrm{m}$                |
| Diameter of the pipeline           | $d$                | 3.78                | $\mathrm{cm}$               |
| Wall thickness of the pipeline     | $b$                | 0.29                | $\mathrm{cm}$               |
| Distance between pipelines         | $w$                | 5.3                 | $\mathrm{cm}$               |
| Pipeline wall thermal conductivity | $\lambda_{p}$      | 0.42                | $\mathrm{W m^{-1} K^{-1}}$  |
| Grout thermal conductivity         | $\lambda_{g}$      | 1                   | $\mathrm{W m^{-1} K^{-1}}$  |
| Grout heat capacity                | $(\rho c)_{grout}$ | $2.5\times10^{6}$   | $\mathrm{Jm^{-3}K^{-1}}$    |

## OGS Input Files

The detailed input parameters can be seen from the 3D_2U_BHE.prj file. The inflow temperature of the BHE, which was imposed as boundary condition of the BHE is shown in Figure 1. All the initial temperatures are set as 22 $^{\circ}$C. The flow rate within each U-pipe is set to $2.0\times10^{-4}$ $\mathrm{m^{3} s^{-1}}$ during the whole simulation time.

{{< img src="../3D_2U_BHE_figures/In_out_temperature_comparison.png" width="200">}}

Figure 1: Inflow temperature curve and outflow temperature comparison

## Results

The OGS numerical outflow temperature over time was compared against results of the FEFLOW software as shown in the Figure 1. And the vertical distributed temperature of circulating water was presented in Figure 2 after operation for 3000 s. The comparison figures demonstrate that the OGS numerical results and FEFLOW results can match very well and the biggest absolute error of outflow temperature is 0.19 $^{\circ}$C at the starting up stage, while such error decreases to 0.05 $^{\circ}$C after 3600 s' operation. The maximum relative error of vertical temperature is 0.015 \% after operation for 3000 s.

{{< img src="../3D_2U_BHE_figures/vertical_temperature_distribution.png" width="200">}}

Figure 2: Comparison of vertical temperature distribution

## References

[1] Diersch, H. J., Bauer, D., Heidemann, W., Rühaak, W., & Schätzl, P. (2011). Finite element modeling of borehole heat exchanger systems: Part 1. Fundamentals. Computers & Geosciences, 37(8), 1122-1135.
