+++
project = "PhaseField/tpb_rgranite.prj"
author = "Keita Yoshioka"
date = "2020-07-27"
title = "Three point bending test"
weight = 158

[menu]
  [menu.benchmarks]
    parent = "Phase-Field"

+++

{{< data-link >}}

## Problem description
**Note**, this project file runs only with a modified version of OGS
which you can find [here](https://github.com/KeitaYoshioka/ogs/tree/H2M_phasefield).

We simulate a three point bending test performed on Rockville granite as shown below. Details of the experiment can be found in Tarok et al. 2017.
{{< img src="../TPB_exp.png" >}}

## Results and evaluation

Developed crack (phase-field) and the crack mouth opening displacement (CMOD) vs. the force are shown below.

{{< img src="../VPF_ME1_frac.png" >}}
{{< img src="../VPF_ME1_NF_CMOD_comp.png" >}}


The model is able to simulate up to the brittle elastic failure, but as cracked surfaces are currently treated as frictionless, the behavior after the failure deviates from the experiment results.
