+++
project = "ThermoRichardsMechanics/TaskCDECOVALEX2023/Decovalex-0.prj"

title = "A test based on DECOVALEX2023 Task C for water vapour diffusion model for non-isothermal Richards flow"
date = "2022-02-07T13:19:49+01:00"
author = "Wenqing Wang and Sonja Kaiser"

weight = 70

[menu]
  [menu.benchmarks]
    parent = "thermo-richards-mechanics"

+++

{{< data-link >}}

This test is prepared exactly according to the specifications of Step 0c of
 [Task C of
 the DECOVALEX 2023 project](https://decovalex.org/D-2023/task-c.html),
 which is aimed to simulate the coupled THM processes
 in the full-scale emplacement experiment (FE experiment) at the Mont Terri
 Underground Rock Laboratory [[1]](#1).

The description of the Step 0c can be found in the task specifications.

Running the entire test takes long time (three hours on a PC with i7-8565U CPU).
 As a test for code, it is limited to run only 20 steps.

The test is used to mainly verify the implementation of water vapour diffusion model,
 heat latency model, and equation wise residuum compensation with non equilibrium
 initial state, and etc.

 The results obtained are comparable with that obtained by other teams in
 Task C. The following figure shows the distribution of temperature in the domain,
 water saturation in the vicinity of the heater, and displacement magnitude
 in the domain after nearly 3 years' heating:
{{< img src="decovalex_2023_c.png" >}}

The following two figures show the temporal variations of temperature and water
 saturation, respectively, at a node near the heater:

<img src="decovalex_2023_c_T_t.png" alt="drawing" width="450"/>
<img src="decovalex_2023_c_S_t.png" alt="drawing" width="450"/>

As shown the water saturation variation curve, the de-saturation -
 re-saturation process is well captured by the numerical simulation.

## References
<a id="1">[1]</a>
Müller, H.R., Garitte, B., Vogt, T. et al. Implementation of the full-scale
 emplacement (FE) experiment at the Mont Terri rock laboratory.
 Swiss J Geosci 110, 287–306 (2017). DOI:
[10.1007/s00015-016-0251-2](https://sjg.springeropen.com/articles/10.1007/s00015-016-0251-2)
