+++
project = "PhaseField/me2_pf_case1.prj"
author = "Keita Yoshioka"
date = "2020-07-27"
title = "Fluid driven percolation"
weight = 158

[menu]
  [menu.benchmarks]
    parent = "Phase-Field"

+++

{{< data-link >}}

## Problem description
**Note**, this project file runs only with a modified version of OGS
which you can find [here](https://github.com/KeitaYoshioka/ogs/tree/H2M_phasefield).

We simulate two different fluid percolation experiments performed on rock salt samples with a true tri-axial loading system as described in [this PDF](../Yoshioka_percolation.pdf). The experiments were performed under two different stress configurations as below.
{{< img src="../ME2_stress_state_1.pdf" >}}
{{< img src="../ME2_stress_state_2.pdf" >}}

## Results and evaluation

Simulated crack pathes (phase-field) for the two cases are shown below:

{{< img src="../VPF_ME2_case1.png" >}}
{{< img src="../VPF_ME2_case2.png" >}}
