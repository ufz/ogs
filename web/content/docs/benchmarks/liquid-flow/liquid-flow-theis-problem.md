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

{{< data-link >}}

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

The following figure compares the analytical solution, the result by ogs5, and
 the result by ogs6 (labeled as `pressure`) within the range that satisfies
 $u <$ 0.05.
{{< img src="../theis_comparison.png" >}}
The figure shows that there is a good match between the analytical solution and
 the numerical solution obtained by using ogs5 or ogs6.

<p>Some of the data of the above curves are given in the following table.</p>
<table>
<caption>Comparison of solutions (where <em>Distance</em> means the distance
 from the position where the source term is applied).</caption>
<tbody>
<tr class="odd">
<td style="text-align: left;">Distance</td>
<td style="text-align: left;">Analytic</td>
<td style="text-align: left;">OGS5</td>
<td style="text-align: left;">OGS6</td>
<td style="text-align: left;">Error</td>
<td style="text-align: left;">Error</td>
</tr>
<tr class="even">
<td style="text-align: left;"></td>
<td style="text-align: left;">Solution (<span class="math inline"><em>h</em><sub><em>a</em></sub></span>)</td>
<td style="text-align: left;">(<span class="math inline"><em>h</em><sub>5</sub></span>)</td>
<td style="text-align: left;">(<span class="math inline"><em>h</em><sub>6</sub></span>)</td>
<td style="text-align: left;">(<span class="math inline">$$|\frac{h_5-h_a}{h_a}|$$</span>)</td>
<td style="text-align: left;">(<span class="math inline">$$|\frac{h_6-h_a}{h_a}|$$</span>)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">0</td>
<td style="text-align: left;">12.8141</td>
<td style="text-align: left;">12.474</td>
<td style="text-align: left;">12.474</td>
<td style="text-align: left;">0.0265</td>
<td style="text-align: left;">0.0272</td>
</tr>
<tr class="even">
<td style="text-align: left;">1.21799</td>
<td style="text-align: left;">8.9441</td>
<td style="text-align: left;">8.79341</td>
<td style="text-align: left;">8.79341</td>
<td style="text-align: left;">0.0168</td>
<td style="text-align: left;">0.0171</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2.43597</td>
<td style="text-align: left;">7.48878</td>
<td style="text-align: left;">7.34717</td>
<td style="text-align: left;">7.34717</td>
<td style="text-align: left;">0.0189</td>
<td style="text-align: left;">0.0192</td>
</tr>
<tr class="even">
<td style="text-align: left;">3.65396</td>
<td style="text-align: left;">6.59978</td>
<td style="text-align: left;">6.46176</td>
<td style="text-align: left;">6.46176</td>
<td style="text-align: left;">0.0209</td>
<td style="text-align: left;">0.0213</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4.87195</td>
<td style="text-align: left;">5.94548</td>
<td style="text-align: left;">5.81072</td>
<td style="text-align: left;">5.81072</td>
<td style="text-align: left;">0.0226</td>
<td style="text-align: left;">0.0231</td>
</tr>
<tr class="even">
<td style="text-align: left;">6.08994</td>
<td style="text-align: left;">5.43381</td>
<td style="text-align: left;">5.30267</td>
<td style="text-align: left;">5.30267</td>
<td style="text-align: left;">0.0241</td>
<td style="text-align: left;">0.0247</td>
</tr>
<tr class="odd">
<td style="text-align: left;">7.30792</td>
<td style="text-align: left;">5.00981</td>
<td style="text-align: left;">4.88283</td>
<td style="text-align: left;">4.88283</td>
<td style="text-align: left;">0.0253</td>
<td style="text-align: left;">0.0260</td>
</tr>
<tr class="even">
<td style="text-align: left;">8.52591</td>
<td style="text-align: left;">4.65012</td>
<td style="text-align: left;">4.52793</td>
<td style="text-align: left;">4.52793</td>
<td style="text-align: left;">0.0262</td>
<td style="text-align: left;">0.0269</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8.83041</td>
<td style="text-align: left;">4.56714</td>
<td style="text-align: left;">4.44623</td>
<td style="text-align: left;">4.44623</td>
<td style="text-align: left;">0.0264</td>
<td style="text-align: left;">0.0271</td>
</tr>
<tr class="even">
<td style="text-align: left;">9.4394</td>
<td style="text-align: left;">4.4116</td>
<td style="text-align: left;">4.29344</td>
<td style="text-align: left;">4.29344</td>
<td style="text-align: left;">0.0267</td>
<td style="text-align: left;">0.0275</td>
</tr>
<tr class="odd">
<td style="text-align: left;">10.6574</td>
<td style="text-align: left;">4.12707</td>
<td style="text-align: left;">4.01501</td>
<td style="text-align: left;">4.01501</td>
<td style="text-align: left;">0.0271</td>
<td style="text-align: left;">0.0279</td>
</tr>
<tr class="even">
<td style="text-align: left;">15.2248</td>
<td style="text-align: left;">3.28072</td>
<td style="text-align: left;">3.19698</td>
<td style="text-align: left;">3.19698</td>
<td style="text-align: left;">0.0255</td>
<td style="text-align: left;">0.0261</td>
</tr>
<tr class="odd">
<td style="text-align: left;">20.0968</td>
<td style="text-align: left;">2.61899</td>
<td style="text-align: left;">2.57517</td>
<td style="text-align: left;">2.57517</td>
<td style="text-align: left;">0.0167</td>
<td style="text-align: left;">0.0170</td>
</tr>
<tr class="even">
<td style="text-align: left;">22.8373</td>
<td style="text-align: left;">2.31338</td>
<td style="text-align: left;">2.29626</td>
<td style="text-align: left;">2.29626</td>
<td style="text-align: left;">0.0074</td>
<td style="text-align: left;">0.0074</td>
</tr>
<tr class="odd">
<td style="text-align: left;">24.0553</td>
<td style="text-align: left;">2.18892</td>
<td style="text-align: left;">2.1846</td>
<td style="text-align: left;">2.1846</td>
<td style="text-align: left;">0.0019</td>
<td style="text-align: left;">0.0019</td>
</tr>
<tr class="even">
<td style="text-align: left;">25.2732</td>
<td style="text-align: left;">2.07055</td>
<td style="text-align: left;">2.07958</td>
<td style="text-align: left;">2.07958</td>
<td style="text-align: left;">0.0043</td>
<td style="text-align: left;">0.0043</td>
</tr>
<tr class="odd">
<td style="text-align: left;">25.5777</td>
<td style="text-align: left;">2.04164</td>
<td style="text-align: left;">2.05407</td>
<td style="text-align: left;">2.05407</td>
<td style="text-align: left;">0.0060</td>
<td style="text-align: left;">0.0060</td>
</tr>
<tr class="even">
<td style="text-align: left;">29.8407</td>
<td style="text-align: left;">1.67134</td>
<td style="text-align: left;">1.73518</td>
<td style="text-align: left;">1.73518</td>
<td style="text-align: left;">0.0381</td>
<td style="text-align: left;">0.0367</td>
</tr>
</tbody>
</table>
<p>The analytical solutions are for an ideal problem with point wise pumping
 term. While for the FEM analysis, the point wise source value is distributed to
 the surface of a small hole around the source point in order to avoid the
 singularity. One can see from the table that the precisions of the FEM
 solutions are still acceptable with such transform of the point pumping
 term.</p>
