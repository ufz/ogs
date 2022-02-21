+++
date = "2022-01-25"
title = "Mandel-Cryer Effect"
weight = 151
project = "HydroMechanics/StaggeredScheme/MandelCryer/MandelCryer.prj"
author = "Dominik Kern"

[menu]
  [menu.benchmarks]
    parent = "hydro-mechanics"

+++

{{< data-link >}}

---

## Mandel-Cryer Effect

This is a classical example to demonstrate the effect of hydromechanical coupling in a poroelastic medium.
For more details we refer to a textbook [[1]](#1), in which also the analytical solution is derived.
As domain we consider a sphere, by symmetry we need to simulate only an octant.

{{< img src="../MandelCryer_mesh.png" >}}
_Mesh_

The boundary conditions of hydraulics are atmospheric pressure on the sphere surface and impermeable elsewhere.
The boundary conditions of mechanics are an overburden (Neumann) applied as step load on the sphere surface at initial time $t=0$. The remaining sides are fixed in normal direction (Dirichlet).

{{< img src="../MandelCryer_model.png" >}}
_Boundary conditions_

The material is isotropic, linear elastic. Solid and fluid are assumed to be incompressible.
In its initial state the sphere is undeformed and there is ambient pressure everywhere.
A sudden load increase on the surface is instantly transferred on the pore pressure, whereas the solid needs time to deform, until it carries the load.
Since the pore fluid is squeezed out of the outer layers first, they act like a tightening belt and consequently the pressure in the center rises, it may even exceed the value of the applied load.
Finally the pore pressure approaches to ambient pressure.

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
<td align="left">Porosity</td>
<td align="left">$0.2$</td>
<td align="left">-</td>
</tr>
<tr class="odd">
<td align="left">Permeability</td>
<td align="left">$10\cdot 10^{-12}$</td>
<td align="left">m$^2$</td>
</tr>
<tr class="even">
<td align="left">Young’s modulus (bulk)</td>
<td align="left">$10\cdot 10^6$</td>
<td align="left">Pa</td>
</tr>
<tr class="odd">
<td align="left">Poisson ratio (bulk)</td>
<td align="left">$0.1$</td>
<td align="left">-</td>
</tr>
<tr class="even">
<td align="left">Biot coefficient</td>
<td align="left">$1.0$</td>
<td align="left">-</td>
</tr>
<tr class="even">
<td align="left">Solid density</td>
<td align="left">$2.5\cdot 10^3$</td>
<td align="left">kg/m$^3$</td>
</tr>
<tr class="even">
<td align="left">Overburden</td>
<td align="left">$1000$</td>
<td align="left">Pa</td>
</tr>
<tr class="even">
<td align="left">Atmospheric pressure</td>
<td align="left">$0$</td>
<td align="left">Pa</td>
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
<td align="left">Radius</td>
<td align="left">$0.4$</td>
<td align="left">m</td>
</tr>
<tr class="odd">
<td align="left">Finite Elements</td>
<td align="left">$8741$ </td>
<td align="left">Taylor-Hood tetrahedral elements</td>
</tr>
<tr class="odd">
<td align="left">Time step</td>
<td align="left">$10^{-2}$</td>
<td align="left">s</td>
<tr class="odd">
<td align="left">Coupling scheme parameter</td>
<td align="left">$0.7774$</td>
<td align="left">-</td>
</tr>
</tbody>
</table>

As predicted, the pressure in the center exceeds the applied load and then levels out.

{{< img src="../MandelCryer_results.png" >}}
_Pressure at center of sphere_

For more details about the staggered scheme we refer to the [user guide - conventions]({{< ref "conventions.md#staggered-scheme" >}}).

## References
<a id="1">[1]</a>
{{< bib "verruijt2009" >}}








