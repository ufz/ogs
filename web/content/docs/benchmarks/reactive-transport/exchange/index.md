+++
author = "Jaime Garibay-Rodriguez, Renchao Lu, Vanessa Montoya"
project = ["Parabolic/ComponentTransport/ReactiveTransport/CationExchange/exchange.prj"]
date = "2019-07-08T14:41:09+01:00"
title = "Transport and Cation Exchange"
image = "fig2.png"
+++

{{< data-link >}}

## Overview

This benchmark simulates the chemical composition of the effluent from a column containing a cation exchanger (Example 11 in the PHREEQC 3 documentation).
The following setup is used for the simulation:

{{< figure src="fig1.png" title="Model setup for the simulation of the column with a cation exchanger in OGS6.">}}

Full details of the model setup and parameters are given in the PHREEQC3 example (consulted MAY-2021):

<https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-73.htm#50593807_46434>

The benchmark uses the `ComponentTransport` process in OGS-6 coupled with the IPhreeqc software (Parkhurst et al. 2013). The results show good agreement between codes. More details about the implementation of the `ComponentTransport` process in OGS-6 can be found in  [`HC-Process.pdf`](/docs/benchmarks/hydro-component/saturated-mass-transport/HC-Process.pdf).

{{< figure src="fig2.png" title="Comparison between PHREEQC and OGS6 of simulated concentrations of solutes at time = 18,000 s reacting with an exchanger.">}}

{{< data-link >}}

## References

<!-- vale off -->

Parkhurst, D.L., Appelo, C.A.J., 2013. Description of Input and Examples for PHREEQC Version 3 - a Computer Program for Speciation, Batch-reaction, One-dimensional Transport, and Inverse Geochemical Calculations.
