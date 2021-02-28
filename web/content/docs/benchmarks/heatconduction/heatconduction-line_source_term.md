+++
author = "Thomas Fischer"
date = "2019-10-01T11:54:02+01:00"
title = "Heatconduction (Line Source Term)"
weight = 121
project = "Parabolic/T/1D_line_source_term_tests/line_source_term.prj"

[menu]
  [menu.benchmarks]
    parent = "heatconduction"

+++

{{< data-link >}}

## Equations

We consider the Poisson equation:
$$
\begin{equation}
\nabla\cdot(\nabla T) + Q_T = 0 \quad \text{in }\Omega
\end{equation}$$
w.r.t Dirichlet-type boundary conditions
$$
\eqalign{
T(x) = 0 &\quad \text{on }\Gamma_D,\cr
}
$$
where $T$ could be temperature, the subscripts $D$ denotes the Dirichlet-type
boundary conditions. Here, the temperature distribution under the impact of a
line shaped source term should be studied.

## Problem Specifications and Analytical Solution

In OGS there are several benchmarks for line source terms in 2d and 3d domains
available. Here, some of the 3d benchmarks are described.

### Cylindrical domain

The Poisson equation on cylindrical domain of height $1$ and radius
$r=1$ is solved. In the following figure the geometry, partly semi-transparent,
is sketched. Furthermore, the mesh resolution is shown in the cylindrical domain
within the first quadrant of the coordinate system. In the second quadrant the
simulated temperature distribution is depicted.

{{< img src="../LineSourceTermFigures/temperature_distribution_line_source_term_in_cylinder.png" >}}

The source term is defined along the line in the center of the cylinder:
$$
\begin{equation}
Q(x) = 1 \quad \text{at } x=0, y=0.
\end{equation}
$$
In the above figure the source term is the red vertical line in the origin of
the coordinate system.

The analytical solution for a line source in the cylinder is
$$
\begin{equation}
T(x) = - \frac{1}{2 \pi} \ln \sqrt{x^2 + y^2}.
\end{equation}
$$

#### Analytical solution in paraview

Since the analytical solution has a singularity at $(x, y) = (0, 0)$ the
analytical solution in paraview is generated as follow:

```none
if (coordsX^2<0.0001 & coordsY^2<0.0001, temperature, -1/(4*asin(1))*ln(sqrt(coordsX^2+coordsY^2))
```

#### Results and evaluation

The following plot shows the temperature along the white line in the figure above.

{{< img src="../LineSourceTermFigures/temperature_profile_line_source_term_in_cylinder.png" >}}

- Comparison with analytical solution:

The differences of analytical and computed solutions for two different domain
discretizations are small outside of the center. In the finer mesh the error
outside of the middle region is smaller than in the coarser mesh.
{{< img
src="../LineSourceTermFigures/comparison_plot_over_line_diff_analytical_solution_temperature_and_simulated_temperature_line_source_term_in_cylinder.png" >}}

Due to the numerical evaluation of the relative error of the computed solution
the error grows in the vicinity of the boundary and in the center.
{{< img
src="../LineSourceTermFigures/comparison_plot_over_line_rel_diff_analytical_solution_temperature_and_simulated_temperature_line_source_term_in_cylinder.png" >}}

#### Input files

The project files for the described models are
[49k.prj](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/T/3D_line_source_term_tests/3D_line_source_term_in_cylinder/49k_prisms/line_source_term_in_cylinder.prj)
and
[286k.prj](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/T/3D_line_source_term_tests/3D_line_source_term_in_cylinder/286k_prisms/line_source_term_in_cylinder.prj).
The project files describe the processes to be solved and the related process variables
together with their initial and boundary conditions as well as the source terms.

The input meshes are stored in the VTK file format and can be directly visualized in Paraview for example.

### Cylindrical domain - axisymmetric example

The Poisson equation on cylindrical domain of height $1$ and radius
$r=1$ is solved. The cylindrical domain is defined as axisymmetric.

#### Results and evaluation

{{< img
src="../LineSourceTermFigures/simulated_temperature_distribution_line_source_term_in_axisymmetric_cylinder.png" >}}
The above figure shows the computed temperature distribution.

The following plot shows the temperature along the white line in the figure above.
{{< img
src="../LineSourceTermFigures/temperature_profile_line_source_term_in_axisymmetric_cylinder.png" >}}

The error and relative error shows the same behaviour like in the simulation
models above. Outside of the center, that has a singularity in the analytical
solution, the errors decreases very fast.
{{< img
src="../LineSourceTermFigures/plot_over_line_diff_and_rel_diff_analytical_solution_temperature_and_simulated_temperature_line_source_term_in_axisymmetric_cylinder.png" >}}

#### Input files

The project file for the described model is
[line_source_term_in_cylinder.prj](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/T/3D_line_source_term_tests/3D_line_source_term_in_cylinder_axisymmetric/line_source_term_in_cylinder.prj).
