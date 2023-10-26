+++
author = "Keita Yoshioka"
date = "2020-07-27"
title = "Three point bending test on anisotropic rock"
image = "ME1_ext_2D_orth_result.png"
+++

{{< data-link >}}

## Problem description

**Note**, this project file runs only with a modified version of OGS
which you can find [here](https://github.com/KeitaYoshioka/ogs/tree/H2M_phasefield).

We simulate three point bending tests performed on lamination orthogonal and parallel. The layer lamination is represented through the fracture toughness $G_c$ as shown below.
{{< figure src="ME1_ext_2D_orthogonal_init.png" >}}
{{< figure src="ME1_ext_2D_parallel_init.png" >}}

## Results and evaluation

Simulated crack path (phase-field) for the lamination orthogonal and the parallel are shown below:

{{< figure src="ME1_ext_2D_orth_result.png" >}}
{{< figure src="ME1_ext_2D_para_result.png" >}}

The responses of crack mouth opening displacement (CMOD) vs. force is as follows.

{{< figure src="VPF_ME1_ex_NF_CMOD.png" >}}
