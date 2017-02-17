+++
date = "2017-02-17T14:33:45+01:00"
title = "Theis' problem"
weight = 171
project = "Parabolic/LiquidFlow/AxiSymTheis/axisym_theis.prj"
author = "Wenqing Wang"

[menu]
  [menu.benchmarks]
    parent = "liquid-flow"

+++

{{< project-link >}}

## Problem description

Theis' problem examines the transient lowering of the water table induced by a pumping well. Theis' fundamental insight was to recognize that Darcy's law is analogous to the law of heat flow by conduction, i.e., hydraulic pressure being analogous to temperature, pressure-gradient to thermal gradient.

The assumptions required by the Theis solution are:
- the aquifer is homogeneous, isotropic, confined, infinite in radial extent,
- the aquifer has uniform thickness, horizontal piezometric surface
- the well is fully penetrating the entire aquifer thickness,
- the well storage effects can be neglected,
- the well has a constant pumping rate,
- no other wells or long term changes in regional water levels.


## Analytical solution

The analytical solution of the drawdown as a function of time and distance is expressed by
$$
\begin{eqnarray}
h_0 - h(t,x,y) = \frac{Q}{4\pi T}W(u)
\label{theis}
\end{eqnarray}
$$

$$
\begin{eqnarray}
u = \frac{(x^{2}+y^{2})S}{4Tt}
\label{theis_u}
\end{eqnarray}
$$

where $h_0$ is the constant initial hydraulic head $[L]$, $Q$ is the constant discharge rate [$L^{3}T^{-1}$], $T$ is the aquifer transmissivity [$L^{2}T^{-1}$], $t$ is time $[T]$, $x,y$ is the coordinate at any point $[L]$ and $S$ is the aquifer storage $[-]$. $W(u)$ is the well function defined by an infinite series for a confined aquifer as

$$
\begin{eqnarray}
W(u) = -\gamma -lnu + \sum^{\infty}_{k=1}{\frac{(-1)^{k+1}u^k}{k\cdot k!}}
\label{theis_wu}
\end{eqnarray}
$$

where $\gamma\approx$ 0.5772 is the Euler-Mascheroni constant. For practical purposes, the simplest approximation of $W(u)$ was proposed as $W(u)=-0.5772-lnu$  for $u <$ 0.05. Other more exact approximations of the well function were summarized by R. Srivastava and A. Guzman-Guzman

## Results and evaluation
