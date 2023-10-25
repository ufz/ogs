+++
author = "Keita Yoshioka"
date = "2020-07-27"
title = "Fluid driven percolation"
image = "VPF_ME2_case1.png"
+++

{{< data-link >}}

## Problem description

**Note**, this project file runs only with a modified version of OGS
which you can find [here](https://github.com/KeitaYoshioka/ogs/tree/H2M_phasefield).

We simulate two different fluid percolation experiments performed on rock salt samples with a true tri-axial loading system as described in [this PDF](./Yoshioka_percolation.pdf). The experiments were performed under two different stress configurations as below.

{{< figure src="ME2_stress_state_1.pdf" class="w-1/2 float-left" >}}
{{< figure src="ME2_stress_state_2.pdf" class="w-1/2 float-left" >}}

## Results and evaluation

Simulated crack paths (phase-field) for the two cases are shown below:

{{< figure src="VPF_ME2_case1.png" class="w-1/2 float-left" >}}
{{< figure src="VPF_ME2_case2.png" class="w-1/2 float-left" >}}
