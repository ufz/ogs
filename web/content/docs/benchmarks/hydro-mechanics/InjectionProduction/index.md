+++
date = "2021-03-17"
title = "Injection and Production in 1D Linear Poroelastic Medium with the Staggered Scheme"
weight = 151
project = ["HydroMechanics/StaggeredScheme/InjectionProduction1D/InjectionProduction1D.prj"]
author = "Wenqing Wang and Dominik Kern"
image = "InjectionProduction_results.png"
+++

{{< data-link >}}

---

## Injection and Production in 1D Linear Poroelastic Medium

This benchmark simulates a soil column with fluid injection at the bottom and a production well at the top.
It is taken from Kim [[1]](#1), in detail it coincides with one of his examples (case 2, coupling strength $\tau=1.21$).
A brief description of the used staggered scheme follows at the end.

{{< figure src="InjectionProduction_model.png" >}}
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
Kim plots his results over non dimensional time, referring to the time at which the produced fluid volume equals the pore volume of the domain (450 days).

{{< figure src="InjectionProduction_results.png" >}}
_Pressure at observation point (marked by circle) versus time (t=0...100 days) and spatial pressure distribution at t=100 days_

## Staggered Scheme: Fixed-stress splitting

For each time step run alternating simulations of the hydraulic (H) problem and the mechanical (M) problem until a convergence criteria is met.
The fixed-stress split starts with the mass balance (H) followed by the momentum balance (M).
These coupling iterations (H,M,H,M,...) add another iteration level compared to the monolithic formulation (HM).
However, due to splitting into smaller problems this may result in a speedup.

The balance equations of mass and momentum for the fully saturated porous medium read
$$
\eqalign{
S\varrho_\mathrm{f}{\dot p}-\frac{k}{\mu}\nabla \left(\varrho_\mathrm{f}\left(\nabla p-
 \varrho_\mathrm{f} \mathbf{g}  \right) \right) + \alpha \varrho_\mathrm{f} {\dot\varepsilon_v} &=& 0\cr
\nabla\cdot \left( \boldsymbol{\sigma}^\mathrm{eff}(\mathbf{u})
-\alpha p \mathbf{I} \right)
&=&  \varrho_\mathrm{b} \mathbf{g}.
}
$$
where $\alpha$ denotes Biot coefficient, $S$ is the storage coefficient,
 $k$ is the intrinsic permeability, $\mu$ the liquid viscosity,
 $\varrho_\mathrm{f}$ is the fluid density, and $\varrho_\mathrm{b}$ is the
  bulk density.

In the staggered scheme for solving HM coupled equations, the fixed-stress splitting
 is employed to enhance the convergence. The fixed stress splitting is based on the
 the volumetric total stress rate definition the hydro-mechanics:
$$
 \dot{\sigma}_v=K_b ({\dot \varepsilon}_v-\dot{\varepsilon}^{ne}_v)- \alpha\dot {p},
$$
with $K_b$ the drained bulk modulus of porous medium, and $\dot{\varepsilon}^{ne}_v$
the volumetric  non-elastic strain rate.

### Fixed stress rate over coupling iteration

As the first option, we consider to fix the volumetric total stress rate over coupling iteration.
 This means
$$
    \dot{\sigma}_v^{n, k} = \dot{\sigma}_v^{n, k-1},
$$
with $n$ the time step index, $k$ the coupling iteration index, and
$\dot{()}^{n, k} = \left(()^{n, k} - ()^{n}\right)/dt$.

This gives the volumetric strain rate of the current time step as
$$
\dot{\varepsilon_v}^{n, k} \approx  \dot{\epsilon}_v^{n, k-1} +
\dfrac{ \alpha}{K_b}   (\dot {p}^{n, k}-\dot {p}^{n, k-1})
    +(\dot{\epsilon}^{ne}_v|^{n, k}-
\dot{\epsilon}^{ne}_v|^{n, k-1}).
$$
Practically, we can set $\dot{\epsilon}^{ne}_v|^{n, k}-
\dot{\epsilon}^{ne}_v|^{n, k-1} = 0$.

Under that  volumetric strain rate approximation, the mass balance equation
 for coupling iteration $k$ at time step $n$ becomes
 $$
\varrho_\mathrm{f}(\beta_{\text S} +\dfrac{ \alpha^2}{K_b} ) {\dot p}^{n, k}-
\frac{k}{\mu}\nabla \left(\varrho_\mathrm{f}\left(\nabla p^{n. k}-
 \varrho_\mathrm{f} \mathbf{g}  \right) \right) +
  \varrho_\mathrm{f} \left(\alpha {\dot\varepsilon_v}^{n,k-1}-
 \dfrac{ \alpha^2}{K_b} {\dot p}^{n, k-1}\right)= 0.
$$
Denoting $\frac{ \alpha^2}{K_b} $  as  $\beta_\mathrm{FS}$, the above equation turns into
 $$
\varrho_\mathrm{f}(\beta_{\text S} +\beta_\mathrm{FS} ) {\dot p}^{n, k}-
\frac{k}{\mu}\nabla \left(\varrho_\mathrm{f}\left(\nabla p^{n. k}-
 \varrho_\mathrm{f} \mathbf{g}  \right) \right) +
  \varrho_\mathrm{f} \left(\alpha {\dot\varepsilon_v}^{n,k-1}-
 \beta_\mathrm{FS} {\dot p}^{n, k-1}\right)= 0.
$$

One can see from the above equation that $\beta_\mathrm{FS}$ can be any scalar
 number once the coupling iteration converges, i.e.
${\dot p}^{n, k}\approx{\dot p}^{n, k-1}$.
Therefore, one can choose an arbitrary value for $\beta_\mathrm{FS}$ if it can
 make the coupling iteration converge.
The analysis by Mikelic and Wheeler [2] revealed that
$\beta_\mathrm{FS}=\frac{\alpha^2}{2K_\mathrm{b}}$ is generally close to the
 a-priori unknown optimal value that enhances the convergence of the coupling
 iteration. For more information,
see the [user guide - conventions]({{< ref "conventions.md#staggered-scheme" >}}).
 For code implementation, we introduce a stabilization factor $\gamma$ to the physically
meaningful parameter $\frac{\alpha^2}{K_\mathrm{b}}$ as:
$$
\beta_\mathrm{FS}=\gamma\dfrac{\alpha^2}{K_\mathrm{b}},
$$
where $\gamma$ is treated as an input parameter.

### Fixed stress rate over time step

We  assume that the volumetric stress rate of the current time step is the same as
 that of the previous time step:

$$
  \dot{\sigma}_v^{n} = \dot{\sigma}_v^{n-1}.
$$

That means the current volumetric strain rate is approximated as

$$
\dot{\varepsilon_v}^{n} \approx  \dot{\epsilon}_v^{n-1} +
\dfrac{ \alpha}{K_b}   (\dot {p}^{n}-\dot {p}^{n-1}).
$$

Consequently, and similarly to the fixed stress over coupling iteration,
the mass balance equation at time step $n$ becomes

$$
\varrho_\mathrm{f}(S +\beta_\mathrm{FS} ) {\dot p}^{n}-
\frac{k}{\mu}\nabla \left(\varrho_\mathrm{f}\left(\nabla p^{n}-
 \varrho_\mathrm{f} \mathbf{g}  \right) \right) +
  \varrho_\mathrm{f} \left(\alpha {\dot\varepsilon_v}^{n-1}-
 \beta_\mathrm{FS} {\dot p}^{n-1}\right)= 0.
$$

In that sense, only one coupling iteration is needed, and the solution accuracy
 is dependent on the time step size. The approach of a fixed stress rate over
 the time step enables the staggered scheme to efficiently solve more HM
 problems, especially those with small strain change, e.g. hydro-mechanical
 modeling of reservoirs.

## References

<a id="1">[1]</a>
{{< bib "kimtchjua2009" >}}

<a id="2">[2]</a>
{{< bib "mikwhe2013" >}}
