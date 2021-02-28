+++
author = "Jasper Bathmann, Dimitri Yu. Naumov, Marc Walther"
weight = 142
project = "Parabolic/ComponentTransport/VariableNeumannBoundary/"
date = "2018-11-01T10:41:09+01:00"
title = "Variable Dependent Boundary Condition"

[menu]
  [menu.benchmarks]
    parent = "Hydro-Component"

+++

{{< data-link >}}

## Overview

The component transport process is used for the benchmark setup. Here, a analytical solution of a simple setup is derived and compared to the numerical results.
This Benchmark is described in [this PDF](../HC-VDBCTest.pdf).

For the setup and parameterization, see the chapter "Density dependent flow - The Goswami Problem" in Kolditz et al. (2012).

## Results

{{< img src="../VDBC_num_ana_comp.png" title="UPPER PART: Analytical solution on the right boundary in dependence of time $t$ of the problem indicated with red dashed line in comparison to numerical solution indicated by blue crosses; LOWER PART: development of relative error in dependence of time $t$. Grid spacing for simulations: 0.1; widest timestep 10. The relative error is below $5 \times 10^{-5}$ for all simulation times.">}}
