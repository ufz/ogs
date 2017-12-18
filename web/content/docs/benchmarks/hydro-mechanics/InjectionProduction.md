+++
date = "2017-12-29"
title = "StaggeredScheme"
weight = 151
project = "HydroMechanics/StaggeredScheme/InjectionProduction1D/InjectionProduction1D.prj"
author = "Wenqing Wang"

[menu]
  [menu.benchmarks]
    parent = "hydro-mechanics"

+++

{{< data-link >}}

---
#     Consolidation example based on fluid injection and production    application

{{< img src="../InjectionProduction.png" >}}

This benchmark is based on an example about injection and
production well that is present by Kim et al. \cite kimTchJua2011, and it is used to verified the monolithic scheme
that is implemented in ogs 6 for modelling the coupled hydraulic
mechanical (HM) processes in the porous media. The deformation is
assumed as plane strain, and the problem is defined in a 10 x 150 m <span class="math inline"><em></em><sup>2</sup></span>
domain (as shown in the figure). Initially, the pore pressure and the
initial stress are set to zero, respectively. All boundaries are
assigned with no fluid flux for the hydraulic process. On the lateral
and the bottom boundaries, no vertical displacement is prescribed. On
the top surface, vertical traction boundary condition of 2.125 MP is
applied. The gravity is neglected in both of the Darcy velocity and the
momentum balance equation.

The material properties are shown in the following table.

<p>The material properties are shown in the following table.</p>
<table>
<caption>Material properties of fluid injection and production example</caption>
<thead>
<tr class="header">
<th align="left">Property</th>
<th align="left">Value</th>
<th align="left">Unit</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Water density</td>
<td align="left">1,000</td>
<td align="left">kg/m<span class="math inline"><em></em><sup>3</sup></span></td>
</tr>
<tr class="even">
<td align="left">Porosity</td>
<td align="left">0.3</td>
<td align="left">–</td>
</tr>
<tr class="odd">
<td align="left">Viscosity</td>
<td align="left"><span class="math inline">10<sup>−3</sup></span></td>
<td align="left">Pa<span class="math inline">⋅</span>s</td>
</tr>
<tr class="even">
<td align="left">Specific storage</td>
<td align="left"><span class="math inline">10<sup>−4</sup></span></td>
<td align="left">m <span class="math inline"><em></em><sup>−1</sup></span></td>
</tr>
<tr class="odd">
<td align="left">Intrinsic permeability</td>
<td align="left"><span class="math inline">4.9346165<em>e</em> × 10<sup>−11</sup></span></td>
<td align="left">m<span class="math inline"><em></em><sup>2</sup></span></td>
</tr>
<tr class="even">
<td align="left">Young’s modulus</td>
<td align="left"><span class="math inline">5 × 10<sup>8</sup></span></td>
<td align="left">Pa</td>
</tr>
<tr class="odd">
<td align="left">Posson ratio</td>
<td align="left">0.3</td>
<td align="left">–</td>
</tr>
</tbody>
</table>

The time duration is 8.64 8.64 <span class="math inline">⋅10<sup>6</sup></span> s, and the time step size is
8.64 <span class="math inline">⋅10<sup>4</sup></span> s. The example is also solved by the monolithic
scheme for a reference solutions. The displacement solution at the last
time step is shown in the enclosed figure too.


## References

{{< bib id="kimTchJua2011" >}}
