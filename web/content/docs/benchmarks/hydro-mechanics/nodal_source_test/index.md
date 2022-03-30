+++
title = "A test with nodal source term"
date = 2021-04-14T15:31:38+02:00
weight = 151
project = "HydroMechanics/NodalSourceTerm/nodal_source_test.prj"
author = "Wenqing Wang"

[menu]
  [menu.benchmarks]
    parent = "hydro-mechanics"

+++

{{< data-link >}}

---

The purpose of this benchmark is to test the source term type of NodalSourceTerm
 being applied to the hydraulic equation in the computations with the
 Taylor-Hood elements,
 e.g., the pressure interpolation using linear shape function,
 and the displacement interpolation using quadratic shape function.

The problem is solved on a rectangular domain 2 m $\times$ 1 m.
Nodal fluid flux is applied on the left boundary nodes with a value of
 0.001 kg/s via the property of NodalSourceTerm.
 The pore pressure on the right boundary is equal to the initial pore pressure
  of $10^5$ Pa.
 The displacement in the normal direction of all boundaries except
the left one is fixed.

Fluid and solid are incompressible.

All involved parameters are listed in the following tables:
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
<tr class="odd">
<td align="left">Permeability</td>
<td align="left">$10^{-14}$</td>
<td align="left">m$^2$</td>
</tr>
<tr class="even">
<td align="left">Young’s modulus (bulk)</td>
<td align="left">$2.5\cdot 10^9$</td>
<td align="left">Pa</td>
</tr>
<tr class="odd">
<td align="left">Poisson ratio (bulk)</td>
<td align="left">$0.27$</td>
<td align="left">-</td>
</tr>
<tr class="even">
<td align="left">Biot coefficient</td>
<td align="left">$0.9$</td>
<td align="left">-</td>
</tr>
</tbody>
</table>

The time period of 86400 is discretised into 100 steps.

The distributions of pressure and displacement at the end time are shown in the
 following figure:
{{< img src="nodal_source_test.png" >}}
