+++
date = "2021-03-17"
title = "Staggered Scheme"
weight = 151
project = "HydroMechanics/StaggeredScheme/InjectionProduction1D/InjectionProduction1D.prj"
author = "Dominik Kern and Wenqing Wang"

[menu]
  [menu.benchmarks]
    parent = "hydro-mechanics"

+++

{{< data-link >}}

---

## Injection and Production in 1D Linear Poroelastic Medium

This benchmark simulates a soil column with fluid injection at the bottom and a production well at the top.
It is taken from from Kim [[1]](#1), in detail it coincides with one of his examples (case 2, coupling strength $\tau=1.21$).
A brief description of the used staggered scheme follows at the end.

{{< img src="InjectionProduction_model.png" >}}
_Simulation model with fluid source, sink, observation point and boundary conditions_

The fluid enters and leaves only via the source and sink in the domain, there is no flow across the boundaries.
The displacements at the bottom are fixed, whereas there is a vertical traction applied on top.
Originally the problem is one-dimensional, for simulation with OpenGeoSys it is created in two dimensions with corresponding boundary conditions.
All parameters are concluded in the following tables.
<table>
<caption>Material Properties</caption>
<thead>
<tr class="header">
<th align="left">Property</th>
<th align="left">Value</th>
<th align="left">Unit</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Fluid density</td>
<td align="left">$10^3$</td>
<td align="left">kg/m$^3$</td>
</tr>
<tr class="odd">
<td align="left">Viscosity</td>
<td align="left">$10^{-3}$</td>
<td align="left">Pa$\cdot$s</td>
</tr>
<tr class="even">
<td align="left">Fluid compressibility</td>
<td align="left">$27.5\cdot 10^{-9}$</td>
<td align="left">Pa$^{-1}$</td>
</tr>
<tr class="even">
<td align="left">Porosity</td>
<td align="left">$0.3$</td>
<td align="left">-</td>
</tr>
<tr class="odd">
<td align="left">Permeability</td>
<td align="left">$493.5\cdot 10^{-16}$</td>
<td align="left">m$^2$</td>
</tr>
<tr class="even">
<td align="left">Young’s modulus (bulk)</td>
<td align="left">$300\cdot 10^6$</td>
<td align="left">Pa</td>
</tr>
<tr class="odd">
<td align="left">Poisson ratio (bulk)</td>
<td align="left">$0$</td>
<td align="left">-</td>
</tr>
<tr class="even">
<td align="left">Biot coefficient</td>
<td align="left">$1.0$</td>
<td align="left">-</td>
</tr>
<tr class="even">
<td align="left">Solid density</td>
<td align="left">$3\cdot 10^3$</td>
<td align="left">kg/m$^3$</td>
</tr>
<tr class="even">
<td align="left">Solid compressibility</td>
<td align="left">$0$</td>
<td align="left">Pa$^{-1}$</td>
</tr>
</tbody>
</table>

<table>
<caption>Dimensions and Discretization</caption>
<thead>
<tr class="header">
<th align="left">Property</th>
<th align="left">Value</th>
<th align="left">Unit</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Height ($y$)</td>
<td align="left">$150$</td>
<td align="left">m</td>
</tr>
<tr class="odd">
<td align="left">Width ($x$)</td>
<td align="left">$10$</td>
<td align="left">m</td>
</tr>
<tr class="odd">
<td align="left">Finite Elements</td>
<td align="left">$15$ Taylor-Hood quadrilateral elements</td>
<td align="left">10 m $\times$ 10 m</td>
</tr>
<tr class="odd">
<td align="left">Time step</td>
<td align="left">$86.4\cdot10^3$</td>
<td align="left">s</td>
</tr>
</tbody>
</table>

<table>
<caption>Sources/Sinks, Initial and Boundary Conditions</caption>
<thead>
<tr class="header">
<th align="left"></th>
<th align="left">Hydraulic</th>
<th align="left">Mechanical</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">injection over area $10$m$\times 10$m</td>
<td align="left">$+1.16\cdot 10^{-4}$ kg/(m$^3$s)</td>
<td align="left"> - </td>
</tr>
<tr class="even">
<td align="left">production over area $10$m$\times 10$m</td>
<td align="left">$-1.16\cdot 10^{-4}$ kg/(m$^3$s)</td>
<td align="left"> - </td>
</tr>
<tr class="odd">
<td align="left">top</td>
<td align="left">no flow</td>
<td align="left">$\sigma_{yy}=-2.125\cdot 10^6$ Pa</td>
</tr>
<tr class="even">
<td align="left">left</td>
<td align="left">no flow</td>
<td align="left">$u_x=0$</td>
</tr>
<tr class="odd">
<td align="left">right</td>
<td align="left">no flow</td>
<td align="left">$u_x=0$</td>
</tr>
<tr class="even">
<td align="left">bottom</td>
<td align="left">no flow</td>
<td align="left">$u_y=0$</td>
</tr>
<tr class="odd">
<td align="left">initial state</td>
<td align="left">$p(x,y)=2.125\cdot 10^6$ Pa</td>
<td align="left">$u_x(x,y)=u_y(x,y)=0$</td>
</tr>
</tbody>
</table>


The gravity related terms are neglected in both: the Darcy velocity and the momentum balance equation.

Note that 100 time steps were used for the following results, whereas the provided input file is set to 1 time step (1 day = 86400 s).
Kim plots his results over nondimensional time, referring to the time at which the produced fluid volume equals the pore volume of the domain (450 days).

{{< img src="InjectionProduction_results.png" >}}
_Pressure at observation point (marked by circle) versus time (t=0...100 days) and spatial pressure distribution at t=100 days_

## Staggered Scheme: Fixed-stress splitting ##
For each time step run alternating simulations of the hydraulic (H) problem and the mechanical (M) problem until a convergence criterium is met.
The fixed-stress split starts with the mass balance (H) followed by the momentum balance (M).
These coupling iterations (H,M,H,M,...) add another iteration level compared to the monolithic formulation (HM).
However, due to splitting into smaller problems this may result in a speedup.

The weak form of mass (scalar test function $v$) and momentum balance (vectorial test function $\mathbf{v}$) after discretization in time (Euler-backward) reads at the $n^\mathrm{th}$ time step for the field equations (boundary terms omitted)
$$
\eqalign{
\int_V \left( \alpha \frac{ \varepsilon_V(\mathbf{u}^{n})- \varepsilon_V(\mathbf{u}^{n-1})}{\Delta t}v
+S\frac{p^n-p^{n-1}}{\Delta t}v
+\frac{k}{\mu}\nabla p^n \cdot \nabla v \right) \mathrm{d}V
&=&\int_V \left(s v-\frac{k}{\mu}\varrho_\mathrm{f} \nabla \cdot \mathbf{g} v \right) \mathrm{d}V, \cr
\int_V \left( \boldsymbol{\sigma}^\mathrm{eff}(\mathbf{u}^n):\boldsymbol{\varepsilon}(\mathbf{v})
-\alpha p^n \varepsilon_V(\mathbf{v}) \right) \mathrm{d}V
&=& \int_V \ \varrho_\mathrm{b} \mathbf{g}\cdot \mathbf{v}  \ \mathrm{d}V,
}
$$
where $\alpha$ denotes Biot coefficient, $S$ storage coefficient, $k$ permeability, $\mu$ viscosity, $\varrho_\mathrm{f}$ fluid density and $\varrho_\mathrm{b}$ bulk density.
The fixed-stress split makes the assumption of constant volumetric stress to eliminate the displacement from the mass balance.
At the $i^\mathrm{th}$ coupling iteration this assumption leads to
$$
K_\mathrm{b} \varepsilon_V(\mathbf{u}^{n,i})-\alpha p^{n,i} = K_\mathrm{b} \varepsilon_V(\mathbf{u}^{n,i-1})-\alpha p^{n,i-1},
$$
where $K_\mathrm{b}$ denotes the drained bulk modulus.
Thus the weak form of the mass balance (H) depends only on pressure $p^{n,i}$ as unknown
$$
\int_V \left( \beta_\mathrm{FS}\frac{p^{n,i}-p^{n,i-1}}{\Delta t}v
+S\frac{p^{n,i}-p^{n-1}}{\Delta t}v
+\frac{k}{\mu}\nabla p^{n,i} \cdot \nabla v \right) \mathrm{d}V
=\int_V \left( s v-\frac{k}{\mu}\varrho_\mathrm{f} \nabla \cdot \mathbf{g} v
-\alpha\frac{ \varepsilon_V(\mathbf{u}^{n,i-1})- \varepsilon_V(\mathbf{u}^{n-1})}{\Delta t} v \right) \mathrm{d}V,
$$
where $\beta_\mathrm{FS}=\beta_\mathrm{FS}^\mathrm{ph}=\frac{\alpha^2}{K_\mathrm{b}}$ would correspond to the physically motivated volumetric stress.
However, choosing another artificial volumetric stress, i.e. using a different bulk modulus than $K_\mathrm{b}$, may accelerate the convergence of the coupling iterations.
The analysis by Mikelic and Wheeler [[2]](#2) revealed that $\beta_\mathrm{FS}^\mathrm{MW}=\frac{\alpha^2}{2K_\mathrm{b}}$ is generally close to the a-priori unknown optimal value $\beta_\mathrm{FS}^\mathrm{opt}$, for more information see [user guide - conventions]({{< ref "conventions.md#staggered-scheme" >}}).

The obtained pressure $p^{n,i}$ is then inserted into the momentum equation (M) with $\mathbf{u}^{n,i}$ as unknown
$$
\int_V \ \boldsymbol{\sigma}^\mathrm{eff}(\mathbf{u}^{n,i}):\boldsymbol{\varepsilon}(\mathbf{v}) \ \mathrm{d}V
= \int_V \left( \varrho_\mathrm{b} \mathbf{g}\cdot \mathbf{v} + \alpha p^{n,i} \varepsilon_V(\mathbf{v}) \right) \mathrm{d}V.
$$
These two steps (H,M) are repeated until convergence is reached and thus the solution of the $n^\mathrm{th}$ time step is found.


## References
<a id="1">[1]</a>
{{< bib "kimtchjua2009" >}}

<a id="2">[2]</a>
{{< bib "mikwhe2013" >}}
