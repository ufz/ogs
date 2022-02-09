+++
project = "ThermoRichardsMechanics/CTF1/CTF1.prj"
title = "A test based on CTF1 experiment for water vapour diffusion model for non-isothermal Richards flow"
date = "2022-02-07T16:01:49+01:00"
author = "Wenqing Wang and Chaudhry Aqeel Afzal"

weight = 70

[menu]
  [menu.benchmarks]
    parent = "thermo-richards-mechanics"

+++
This test simulates the coupled thermal hydraulic processes in
 In the CTF1 experiment carried out by Villar et al. [[1]](#1).


The description of this test can be found in
the paper by Wang et al. [[2]](#2). In the calculation, the formula of
specific heat capacity of solid phase has already taken account of
 the porosity factor.

The following figures compare the results of this test against the results
 presented in [[2]](#2):

<img src="../CTF1_results_T.jpg" alt="drawing" width="400"/>
<img src="../CTF1_results_S.jpg" alt="drawing" width="400"/>

## References
<a id="1">[1]</a>
Villar MV, Fernandez AM, Cuevas J (1997) Caracterizacio ́n
Geoquı ́mica de bentonita compactada: efectos producidos por
flujo termohidra ́ulico. Interim Report FEBEX, Informe 70-IMA-
M-0-2, CIEMAT, Madrid

<a id="2">[2]</a>
Wang, W., Rutqvist, J., Görke, U., Birkholzer, J., Kolditz, O. (2011)
 Non-isothermal flow in low permeable
 porous media: a comparison of Richards’ and two-phase flow approaches.
 Environ Earth Sci 62, 1197–1207.
 DOI: [10.1007/s12665-010-0608-1](https://doi.org/10.1007/s12665-010-0608-1)
