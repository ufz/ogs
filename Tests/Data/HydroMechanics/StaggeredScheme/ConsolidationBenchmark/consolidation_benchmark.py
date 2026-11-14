# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "One dimensional consolidation problem"
# date = "2026-10-21"
# author = "Wenqing Wang"
# image = "figures/consolidation_1d.png"
# web_subsection = "hydro-mechanics"
# weight = 151
# +++

# %% [markdown]
# ## 1. Introduction
#
# We consider a porous column bounded by rigid and impermeable walls, except on
# its top where a mechanical pressure load $\sigma_0$ is prescribed and which is free to
# drain. Suppose that the length of the column is $H$. The boundary and initial
# conditions for the problem are
#
# $$
# \eqalign{
# &\sigma =\sigma_0, p = 0\, \text{on}\, x=0,\cr
# &u=0, \dfrac{\partial p}{\partial x}=0\, \text{on}\, x=H,\cr
# &\sigma = p = 0\, \text{at}\, t=0.
# }
# $$
#
# where $x$ is the depth measured from the top of the column.
#
# The analytical solution for this problem can be found in [1] as
# $$
# \eqalign{
# &\sigma_D=-1+\sum_{n=0}^{\infty} \dfrac{2}{M}\sin(Mx_D)\exp(-M^2t_D), \cr
# &u_D=1-x_D-\sum_{n=0}^{\infty} \dfrac{2}{M^2}\cos(Mx_D)\exp(-M^2t_D), \cr
# &p_D =\sigma_D+1.
# }
# $$
#
# with $x_D = x/H$, $t_D = (\lambda + 2G)kt/\eta H^2$, $M = \pi(2n + 1)/2$ are
#  non-dimensional quantities and $\sigma_D = \sigma/p_0$ ,
#  $u_D = u(\lambda + 2 G)/p_0 H$, $p_D = p/p_0$ are dimensionless effective
#  stress, displacement and pore pressure, respectively, where $\lambda$ is
#  the Lame constant, $G$ the shear modulus, and $p_0=|\sigma_0|$.
#
# Since the storage is zero, the pressure change in this example is purely
#  driven by the strain change.
#
# In this example, the densities can be arbitrary non-zero values. The values
#  of the required non-zero material parameters are given in the following
#  table:
#
# | Parameter | Symbol | Value | Unit |
# |---|---|---|---|
# | Young's modulus | $E$ | $3\times10^{4}$ | Pa |
# | Poisson ratio | $\nu$ | 0.2 | – |
# | Viscosity | $\eta$ | 0.001 | Pa s |
# | Intrinsic permeability | $k$ | $10^{-10}$ | m$^2$ |
#
# The column is $H=1$ m long. The mechanical load on the top boundary is
# 1000 Pa, i.e. $\sigma_0=-1000$ Pa.
#
# The problem was also solved using OGS 5 and the mixed FEM approach in [2].
#
# **Reference**
#
# [1] M. Murad and A. Loula. [Improved accuracy in finite element analysis of Biot's
# consolidation problem](https://doi.org/10.1016/0045-7825(92)90193-N). Comput. Methods Appl. Mech. Engrg., 95:359–382, 1992.
#
# [2] J. Korsawe, G. Starke, W. Wang and O. Kolditz. [Finite element analysis of poro-elastic consolidation in porous media: Standard and mixed approaches](https://doi.org/10.1016/j.cma.2005.04.011). Comput. Methods Appl. Mech. Engrg., 195:1096–1115, 2006.

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pyvista as pv

ot.plot.setup.dpi = 320
ot.plot.setup.num_levels = 10
ot.plot.setup.fontsize = 8


# %%
def contour_plot(mesh_data, variables):
    fig, axs = plt.subplots(1, len(variables), figsize=(10, 6), sharey=True)
    for ax, variable in zip(axs, variables, strict=True):
        ot.plot.contourf(
            mesh_data,
            variable,
            fig=fig,
            ax=ax,
            cmap="turbo",
        )

    fig.tight_layout()


# %% [markdown]
# ## 2. Simulation results
#
# The simulation runs to 10 s. The same problem is solved once with the
# monolithic scheme (time step size 1 s, project file
# `consolidation_benchmark_mono.prj`) and once with the staggered scheme
# (fixed stress splitting, time step size 0.5 s, project file
# `consolidation_benchmark.prj`). Both runs are registered as CTests, which
# validate their output at $t = 6$ s and at the end time against the
# committed reference results for the same five output fields as the former
# CTests before this benchmark became a notebook.
#
# Run from the test suite, where the `OGS_TESTRUNNER` environment variable is
# set by CTest on CI and locally alike, the notebook reuses the pre-computed
# output of those two tests and analyses it. Run interactively instead, for
# example in Jupyter or on Binder, it runs the two simulations itself and
# does not rely on any pre-computed output:

# %%
if "OGS_TESTRUNNER" in os.environ:
    ms_mono = ot.MeshSeries("ConsolidationBenchmark_mono.pvd")
    ms = ot.MeshSeries("ConsolidationBenchmark.pvd")
else:
    out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in ["consolidation_benchmark_mono", "consolidation_benchmark"]:
        model = ot.Project(
            input_file=f"{name}.prj", output_file=out_dir / f"{name}.prj"
        )
        model.write_input()
        model.run_model(logfile=out_dir / f"{name}.log", args=f"-o {out_dir} -m .")
    ms_mono = ot.MeshSeries(out_dir / "ConsolidationBenchmark_mono.pvd")
    ms = ot.MeshSeries(out_dir / "ConsolidationBenchmark.pvd")

# %% [markdown]
# The staggered solution is compared with the monolithic one and both with
# the analytical solution further below. The comparisons cover only $p$,
# $u_y$ and $\sigma_{yy}$, because the problem is one-dimensional.
#
# Each comparison has its own absolute tolerances, in the units of the
# compared variables. The scales are the load $p_0 = 1000$ Pa and the
# largest settlement at the end time, about $6\times10^{-3}$ m.

# %%
output_variable_names = ["pressure", "displacement", "sigma"]

# Staggered (time step 0.5 s) vs monolithic (time step 1 s): 1 % of the
# scales. The two differ by the time discretization error, measured as
# 0.69 % for p and sigma_yy and 0.62 % for u_y.
monolithic_atols = {"pressure": 10.0, "displacement": 6.0e-5, "sigma": 10.0}
# Numerical vs analytical solution: 2 % of the scales. Measured
# discretization errors: 0.7 % (staggered, 0.5 s) and 1.4 % (monolithic,
# 1 s), dominated by the first-order time discretization.
analytical_atols = {"pressure": 20.0, "displacement": 1.2e-4, "sigma": 20.0}

# %% [markdown]
# **The time step size controls the accuracy.** Once the coupling iteration
# has converged, the staggered solution does not depend on the fixed stress
# stabilization parameter $\beta_{\text{FS}}$, the
# solutions at the end time differ by less than $10^{-6}$ Pa in pressure and
# $10^{-11}$ m in displacement, which is the level of the coupling tolerance.
# $\beta_{\text{FS}}$ only changes the number of coupling iterations, which
# is smallest for $\beta_{\text{FS}} = 0.5$ (7 iterations in the worst time
# step, compared with 40 for $\beta_{\text{FS}} = 1$). These differences and
# iteration counts were obtained by running the same setup with the other
# values of $\beta_{\text{FS}}$ outside this notebook. The error with respect
# to the analytical solution is instead set by the time discretization. With
# the first-order backward Euler scheme it halves with the time step size:
#
# | Time step size / s | Pressure error / Pa | Displacement error / m |
# |---|---|---|
# | 1.0 | 13.8 | $7.6\times10^{-5}$ |
# | 0.5 | 6.9 | $3.8\times10^{-5}$ |
# | 0.25 | 3.5 | $1.8\times10^{-5}$ |
# | 0.125 | 1.8 | $8.6\times10^{-6}$ |
#
# The errors are the largest differences from the analytical solution at
# $t = 10$ s. The rows for 1 s and 0.5 s are what the analysis below reports
# for the two recorded runs; the two finer time steps were obtained by
# refining the staggered run further in the same setup. This is why the
# staggered solution with the time step size 0.5 s is closer to the
# analytical solution than the monolithic one with 1 s, and why the
# tolerances above are based on the time discretization error.


# %%
def get_data_1d_component(data):
    """Return the y component of vector or tensor data as a 1D array.

    :param data: Scalar point data, or vector or tensor point data with one
        row per point.
    :returns: Scalar data unchanged, otherwise the second column, which is
        $u_y$ for the displacement and $v_y$ for the velocity, and the yy
        entry for the tensor fields: $\\varepsilon_{yy}$ for the strain and
        $\\sigma_{yy}$ for the stress (Kelvin order xx, yy, zz, xy).
    """
    data = np.asarray(data)
    return data if data.ndim == 1 else data[:, 1]


def check_results(mesh, expected_data, atols):
    """Assert that the results on a mesh agree with the expected data.

    :param mesh: The mesh with the numerical results as point data.
    :param expected_data: The expected 1D data per variable name.
    :param atols: The absolute tolerance per variable name. Its keys are the
        variables that are compared.
    """
    for variable_name, atol in atols.items():
        np.testing.assert_allclose(
            get_data_1d_component(mesh[variable_name]),
            expected_data[variable_name],
            rtol=0.0,
            atol=atol,
            err_msg=variable_name,
        )

    print("Test passed.")


# %%
class AnalyticalSolution:
    """Analytical solution of the 1D consolidation problem [1].

    The depth is measured from the top of the column, i.e. x = y_max - y.

    :param E: Young's modulus in Pa.
    :param nu: Poisson ratio.
    :param k: Intrinsic permeability in m^2.
    :param eta: Viscosity in Pa s.
    :param H: Length of the column in m.
    :param p0: Magnitude of the load on the top boundary in Pa.
    """

    def __init__(self, E, nu, k, eta, H, p0):
        lame_lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        G = E / (2.0 * (1.0 + nu))
        self.oedometric_modulus = lame_lambda + 2.0 * G
        self.consolidation_coefficient = self.oedometric_modulus * k / eta
        self.H = H
        self.p0 = p0
        self.M = 0.5 * np.pi * (2 * np.arange(2000) + 1)

    def _series_terms(self, depth, t):
        x_D = np.asarray(depth)[..., np.newaxis] / self.H
        t_D = (
            self.consolidation_coefficient * np.asarray(t)[..., np.newaxis] / self.H**2
        )
        return x_D, np.exp(-(self.M**2) * t_D)

    def pressure(self, depth, t):
        x_D, decay = self._series_terms(depth, t)
        return self.p0 * np.sum(2.0 / self.M * np.sin(self.M * x_D) * decay, axis=-1)

    def sigma_yy(self, depth, t):
        return self.pressure(depth, t) - self.p0

    def displacement_y(self, depth, t):
        x_D, decay = self._series_terms(depth, t)
        u_D = (
            1.0
            - x_D[..., 0]
            - np.sum(2.0 / self.M**2 * np.cos(self.M * x_D) * decay, axis=-1)
        )
        # u_D is positive downwards.
        return -u_D * self.p0 * self.H / self.oedometric_modulus


# %% [markdown]
# The staggered solution agrees with the monolithic one up to the time
# discretization error.

# %%
check_results(
    ms[-1],
    {name: get_data_1d_component(ms_mono[-1][name]) for name in output_variable_names},
    monolithic_atols,
)


# %%
def with_mid_node_pressure(mesh):
    # The mesh consists of six-node triangles, whose nodes 3, 4 and 5 lie on
    # the edges (0, 1), (1, 2) and (2, 0), respectively.
    cells = mesh.cells_dict[pv.CellType.QUADRATIC_TRIANGLE]
    pressure = np.array(mesh["pressure"])
    for mid, (a, b) in zip([3, 4, 5], [(0, 1), (1, 2), (2, 0)], strict=True):
        pressure[cells[:, mid]] = 0.5 * (pressure[cells[:, a]] + pressure[cells[:, b]])
    mesh = mesh.copy()
    mesh["pressure"] = pressure
    return mesh


# %% [markdown]
# The pressure is linear on the quadratic mesh, and its values on the
# mid-edge nodes are not written to the output (they are zero). For the plots
# below and for the comparison with the analytical solution, the pressure is
# therefore interpolated linearly to the mid-edge nodes.


# %%
ms = ms.transform(with_mid_node_pressure)
ms_mono = ms_mono.transform(with_mid_node_pressure)


# %% [markdown]
# The parameters of the analytical solution are read from the project file
# of the monolithic scheme, so that they cannot differ from those of the
# simulations.


# %%
def read_value(prj, xpath):
    return float(prj.tree.findtext(xpath))


prj_mono = ot.Project(input_file="consolidation_benchmark_mono.prj")

x_min, x_max, y_min, y_max = ms[-1].bounds[0:4]
H = y_max - y_min
x_centre = 0.5 * (x_min + x_max)
ana = AnalyticalSolution(
    E=read_value(prj_mono, "./parameters/parameter[name='E']/value"),
    nu=read_value(prj_mono, "./parameters/parameter[name='nu']/value"),
    k=read_value(
        prj_mono,
        "./media/medium/properties/property[name='permeability']/value",
    ),
    eta=read_value(
        prj_mono,
        "./media/medium/phases/phase/properties/property[name='viscosity']/value",
    ),
    H=H,
    p0=-read_value(prj_mono, "./parameters/parameter[name='load']/value"),
)

# %% [markdown]
# Both solutions agree with the analytical one up to the discretization
# error.

# %%
for series in [ms, ms_mono]:
    mesh = series[-1]
    depth = y_max - mesh.points[:, 1]
    t = series.timevalues[-1]
    expected_data = {
        "pressure": ana.pressure(depth, t),
        "displacement": ana.displacement_y(depth, t),
        "sigma": ana.sigma_yy(depth, t),
    }
    check_results(mesh, expected_data, analytical_atols)
    differences = {
        name: np.abs(get_data_1d_component(mesh[name]) - expected).max()
        for name, expected in expected_data.items()
    }
    print(
        f"largest differences from the analytical solution at t = {t} s: "
        + ", ".join(f"{name}: {value:.4g}" for name, value in differences.items())
    )

# %%
contour_variables = [
    ot.variables.pressure.replace(output_unit="Pa"),
    ot.variables.displacement["y"].replace(output_unit="mm"),
    ot.variables.stress["yy"].replace(output_unit="Pa"),
]
# The analytical solution in the same order, in SI units.
analytical_functions = [ana.pressure, ana.displacement_y, ana.sigma_yy]


# %%
def to_output_unit(variable, values):
    """Convert values from the data unit to the output unit of a variable.

    :param variable: The ogstools variable that defines both units.
    :param values: The values in the data unit of `variable`.
    :returns: The values in the output unit of `variable`.
    """
    return (
        ot.variables.u_reg.Quantity(values, variable.data_unit)
        .to(variable.output_unit)
        .magnitude
    )


# %% [markdown]
# ## 3. Results
#
# ### 3.1 Contour plots of the staggered solution
#
# The figures below show the distributions of pore pressure, displacement and
# stress at the end of the simulation.

# %%
contour_plot(ms[-1], contour_variables)

# %% [markdown]
# ### 3.2 Contour plots of the monolithic solution
#
# The same quantities obtained with the monolithic scheme are shown for
# comparison.

# %%
contour_plot(ms_mono[-1], contour_variables)

# %% [markdown]
# ### 3.3 Vertical profiles
#
# The following plots show the vertical distributions of pore pressure,
# vertical displacement $u_y$ and vertical effective stress $\sigma_{yy}$ at
# the end of the simulation, sampled along the vertical centre line of the column
# from its top to its bottom.

# %%
point_top = (x_centre, y_max, 0.0)
point_bottom = (x_centre, y_min, 0.0)

profile = ms[-1].sample_over_line(point_top, point_bottom, resolution=100)
profile_mono = ms_mono[-1].sample_over_line(point_top, point_bottom, resolution=100)

end_time = ms.timevalues[-1]
depth = np.linspace(0.0, H, 201)


# %%
fig, axs = plt.subplots(1, 3, figsize=(14, 4))
for ax, variable, analytical_function in zip(
    axs, contour_variables, analytical_functions, strict=True
):
    ax.plot(
        depth,
        to_output_unit(variable, analytical_function(depth, end_time)),
        color="r",
        lw=0.7,
        label="analytical",
    )
    ax.plot(
        profile["Distance"],
        variable.transform(profile),
        color="C2",
        marker="o",
        markevery=5,
        label="staggered",
    )
    ax.plot(
        profile_mono["Distance"],
        variable.transform(profile_mono),
        color="C0",
        linestyle="dashed",
        label="monolithic",
    )
    ax.set_xlabel("depth / m")
    ax.set_ylabel(variable.get_label())
    ax.grid()
    ax.legend()
fig.tight_layout()
plt.show()

# %% [markdown]
# ### 3.4 Temporal evolution
#
# The following plots show how the pore pressure, the vertical displacement
# and the vertical effective stress evolve over time at three depths.

# %%
depths = np.array([0.1, 0.2, 0.4])
points = np.array([[x_centre, y_max - d, 0.0] for d in depths])

times = np.linspace(0.0, end_time, 201)


# The pressure on the mid-edge nodes is filled in before probing, see above.
probes_stg = ms.probe(points)
probes_mono = ms_mono.probe(points)

fig, axs = plt.subplots(1, 3, figsize=(14, 4))
for ax, variable, analytical_function in zip(
    axs, contour_variables, analytical_functions, strict=True
):
    # Shape (number of time steps, number of points), in the output unit.
    data_stg = probes_stg.values(variable)
    data_mono = probes_mono.values(variable)
    data_ana = to_output_unit(
        variable,
        analytical_function(depths[np.newaxis, :], times[:, np.newaxis]),
    )
    for i, d in enumerate(depths):
        color = f"C{i}"
        ax.plot(
            times,
            data_ana[:, i],
            color=color,
            lw=0.7,
            label=f"depth = {d:.1f} m (analytical)",
        )
        ax.plot(
            probes_stg.timevalues,
            data_stg[:, i],
            color=color,
            marker="o",
            linestyle="none",
            label=f"depth = {d:.1f} m (staggered)",
        )
        ax.plot(
            probes_mono.timevalues,
            data_mono[:, i],
            color=color,
            marker="x",
            linestyle="none",
            label=f"depth = {d:.1f} m (monolithic)",
        )
    ax.set_xlabel("time / s")
    ax.set_ylabel(variable.get_label())
    ax.grid()
axs[0].legend(fontsize=7)
fig.tight_layout()
plt.show()

# %% [markdown]
# The following plots show the differences of the staggered and monolithic
# solutions from the analytical solution at the same three depths. The
# differences are of the order of the time discretization error, which is
# smaller for the staggered solution with the time step size 0.5 s than for
# the monolithic solution with 1 s. The initial state at $t = 0$ is left
# out, because the analytical solution jumps to the undrained response
# immediately after the load is applied.


# %%
def difference_from_analytical(probes, variable, analytical_function, depths):
    """Return the times after t = 0 and the differences there.

    :param probes: The probed numerical results.
    :param variable: The ogstools variable to compare.
    :param analytical_function: The analytical solution of `variable`, in SI
        units.
    :param depths: The depths of the probed points in m.
    :returns: The output times without t = 0, and the differences with shape
        (number of times, number of points), in the output unit.
    """
    t = probes.timevalues[1:]
    data_ana = to_output_unit(
        variable,
        analytical_function(depths[np.newaxis, :], t[:, np.newaxis]),
    )
    return t, probes.values(variable)[1:] - data_ana


fig, axs = plt.subplots(1, 3, figsize=(14, 4))
for ax, variable, analytical_function in zip(
    axs, contour_variables, analytical_functions, strict=True
):
    t_stg, diff_stg = difference_from_analytical(
        probes_stg, variable, analytical_function, depths
    )
    t_mono, diff_mono = difference_from_analytical(
        probes_mono, variable, analytical_function, depths
    )
    for i, d in enumerate(depths):
        color = f"C{i}"
        ax.plot(
            t_stg,
            diff_stg[:, i],
            color=color,
            marker="o",
            label=f"depth = {d:.1f} m (staggered)",
        )
        ax.plot(
            t_mono,
            diff_mono[:, i],
            color=color,
            marker="x",
            linestyle="dashed",
            label=f"depth = {d:.1f} m (monolithic)",
        )
    ax.set_xlabel("time / s")
    ax.set_ylabel(f"difference in {variable.get_label()}")
    ax.grid()
axs[0].legend(fontsize=7)
fig.tight_layout()
plt.show()


# %% [markdown]
# ## 4. Conclusions
#
# - Both schemes reproduce the analytical solution. However, the remaining
#   differences are small and are dominated by the first-order time
#   discretization. This is why the staggered solution with a time step of
#   0.5 s is about twice as accurate as the monolithic solution with 1 s.
#   The differences are largest shortly after the load is applied and near
#   the drained top boundary; for example, at a depth of 0.1 m and
#   $t = 2$ s they are about 31 Pa (staggered) and 56 Pa (monolithic),
#   and they decrease as consolidation proceeds. The error is set by the time
#   discretization, so refining the mesh near the top boundary cannot reduce
#   it below that.
# - Once the coupling iteration has converged, the staggered solution is
#   essentially independent of the fixed-stress stabilization parameter
#   $\beta_{\text{FS}}$; it only affects the number of coupling iterations.
#   With $\beta_{\text{FS}} = 0.5$, the fewest coupling iterations are required.
#   To improve the accuracy of the staggered solution, reduce the time step
#   size, not $\beta_{\text{FS}}$.
