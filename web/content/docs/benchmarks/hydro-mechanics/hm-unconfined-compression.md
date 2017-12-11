+++
project = "https://github.com/ufz/ogs-data/blob/master/HydroMechanics/Linear/Unconfined_Compression_early/square_1e2_UC_early.prj"
author = "Tianyuan Zheng"
date = "2017-02-16T15:10:33+01:00"
title = "Unconfined compression"
weight = 153

[menu]
  [menu.benchmarks]
    parent = "hydro-mechanics"

+++

{{< data-link >}}

## Equations

We consider coupled hydromechanical equation here:
$$
\begin{equation}
S_\mathrm{s}\frac{\partial p}{\partial t} + \alpha \nabla\cdot \frac{\partial \mathbf{u}}{\partial t} + \nabla \cdot\left(\frac{\mathbf{k}}{\mu}(-\nabla p + \varrho \mathbf{g}) \right) = 0 \quad \text{in }\Omega
\end{equation}$$

\begin{equation}
\nabla\cdot (\sigma^\mathrm{E} - \alpha p\mathbf{I}) + \varrho\mathbf{g} = 0 \quad \text{in }\Omega
\end{equation}
w.r.t boundary conditions
$$
\eqalign{
p = \bar{p} &\quad \text{on }\Gamma^p_D,\cr
\tilde{\mathbf{w}}\cdot\mathbf{n} = \bar{\mathbf{w}} &\quad \text{on }\Gamma^p_N, \cr
\mathbf{u} = \mathbf{\bar{u}} &\quad \text{on }\Gamma^\mathbf{u}_D, \cr
\mathbf{\sigma} \cdot \mathbf{n} = \mathbf{\bar{t}}  &\quad \text{on }\Gamma^\mathbf{u}_N,
}
$$

where $p$ could be pressure, $\mathbf{u}$ could be displacement, $\tilde{\mathbf{w}}$ coulde be the water flow ($\tilde{\mathbf{w}} = -\frac{\mathbf{k}}{\mu}(\nabla p - \varrho_\mathrm{LR} \mathbf{g})$)the subscripts $D$ and $N$ denote the Dirichlet- and Neumann-type boundary conditions, $n$ is the normal vector pointing outside of $\Omega$, and $\Gamma = \Gamma_D \cup \Gamma_N$ and $\Gamma_D \cap \Gamma_N = \emptyset$.

## Problem description

We solve a hydro-mechanical linear biphasic model (small deformation, linear elastic, Darcy flow, incompressible constituents) in square domain where on the top boundary a constant displacement boundary is applied. On the right boundary a constant pressure boundary equals zero and zeros traction boundary are applied. All other boundaries are constrained in their normal direction and all boundaries except for outer radius are sealed. The fluid is allowed to escape through the right boundary. The drainage process can be concluded into two stages. During drainage, the total stress is the sum of effective stresses in the solid and the pore pressure. Once the material is fully drained, the pore pressure is zero, so that stress- and displacement fields are determined exclusively by the properties of the solid skeleton. An axisymmetric domain is used in this model. The mesh is refined based on the distance to the outer radius.

{{< img src="../mesh_UC_final.png" >}}

## Assumptions

In this problem, it is assumed that the biot coefficient $\alpha = 1$ and the Storage term $S_\mathrm{s}$ is neglected.

## Results and evaluation

The analytical solution of the problem can be found in Armstrong _et.al._ (1984)"An Analysis of the Unconfined
Compression of Articular Cartilage."

The result after 1s shows that due to the direct displacement on the boundary, the displacement on the x direction is quite large. Where the location is close enough to the outer radius, some water has already flowed out, the pore pressure is decreased and the gradient of displacement curve is modified.

{{< img src="../num_ana_1s_refined.png" >}}

The result after 1000s shows that when the drainage process is finished. The displacement has bounced back and the curve is almost straight due to the linear elastic behavior of the solid.

{{< img src="../verification_1000s_new-1.png" >}}

In order to capture the transient process at certain location, a point at the outer boundary is chosen to show the displacement during the whole compression process.

{{< img src="../transient_validation_unconfined_compression.png" >}}
