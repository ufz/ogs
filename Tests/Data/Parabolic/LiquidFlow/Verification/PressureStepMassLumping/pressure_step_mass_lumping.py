# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "Pressure step: consistent and lumped mass matrices"
# date = "2026-11-10"
# author = "Thomas Fischer"
# web_subsection = "liquid-flow"
# +++
#

# %% [markdown]
# # Pressure step: consistent and lumped mass matrices
#
# ## Problem description
#
# The benchmark simulates pressure diffusion in a saturated, semi-infinite porous medium
# after a sudden pressure increase of 1 MPa at the boundary $x = 0$ at time $t = 0$.
# In this case, the mass balance equation of the liquid reduces to a diffusion equation
# for the pressure $p(x, t)$:
#
# $$S \frac{\partial p}{\partial t} = \frac{k}{\mu} \, \frac{\partial^2 p}{\partial x^2},$$
#
# where $t$ is the time, $S$ is the storage coefficient, $k$ is the permeability,
# and $\mu$ is the dynamic viscosity.
#
# With the boundary condition $p(0, t) = p_0$ and the initial condition $p(x, 0) = 0$
# for $x > 0$, the analytical solution is
#
# $$p(x, t) = p_0 \, \mathrm{erfc}\left(\frac{x}{2\sqrt{D t}}\right),$$
#
# where $p_0$ is the pressure applied at the boundary and $D = k/(\mu S)$ is
# the diffusivity.
#
# For the numerical solution the domain is truncated to $[0, 1]$ m. The left boundary at
# $x = 0$ is subject to a Dirichlet boundary condition with $p = p_0 = 1$ MPa.
# All other boundaries are impermeable.
#
# The material properties are:
#
# | Property                  | Value  | Unit   |
# |---------------------------|-------:|--------|
# | Permeability              | 1e-12  | m$^2$  |
# | Porosity                  | 0.2    | --     |
# | Storage                   | 1e-9   | Pa$^{-1}$ |
# | Liquid viscosity          | 1e-3   | Pa s   |
# | Reference temperature     | 293.15 | K      |
#
# These properties give a diffusivity of
#  $D = k / (\mu S) = 10^{-12} / (10^{-3} \cdot 10^{-9}) = 1$ m$^2$/s.
#
# The initial condition is $p = 0$ everywhere.
#
# The simulation ends at $t = 0.01$ s, before the front reaches the right end at
# $x = 1$ m.
#
# The first time step $\Delta t = 10^{-6}$ s is small compared to
# $h^2 / D = 10^{-4}$ s (element size $h = 0.01$ m). The pressure front then
# stays inside the first element. In this case the consistent mass matrix
# gives non-physical, sign-alternating pressures ahead of the front. Mass
# lumping (`<mass_lumping>true</mass_lumping>`) is meant to avoid them.
#
# Each benchmark comes in four variants:
#
# | Elements  | Mass matrix | 1D (`PressureStep1D`)                           | 2D (`PressureStep2D`)                           |
# |-----------|-------------|-------------------------------------------------|-------------------------------------------------|
# | linear    | consistent  | `pressure_step_1D.prj`                          | `pressure_step_2D.prj`                          |
# | linear    | lumped      | `pressure_step_1D_mass_lumping.xml`             | `pressure_step_2D_mass_lumping.xml`             |
# | quadratic | consistent  | `pressure_step_1D_quadratic.xml`                | `pressure_step_2D_quadratic.xml`                |
# | quadratic | lumped      | `pressure_step_1D_quadratic_mass_lumping.xml`   | `pressure_step_2D_quadratic_mass_lumping.xml`   |
#
# The 2D domain is the rectangle $[0, 1] \times [0, 0.02]$ m$^2$, split into
# $100 \times 2$ square elements (Quad4 or Quad8). For quadratic elements
# OGS uses HRZ lumping. It scales the diagonal of the consistent mass matrix
# so that the total mass stays the same. For Quad8 the HRZ matrix differs from
# the simple row-sum lumped matrix used for linear elements. Row-sum lumping
# would give negative masses at the corner nodes. OGS has no option for this,
# so the last part of the notebook computes the row-sum case in Python.

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pyvista as pv
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import erfc

# %%
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

p_0 = 1e6  # Pa
D = 1.0  # m^2/s

variants = {
    "linear, consistent": "",
    "linear, lumped": "_mass_lumping",
    "quadratic, consistent": "_quadratic",
    "quadratic, lumped": "_quadratic_mass_lumping",
}
styles = {
    "linear, consistent": {"color": "C0", "marker": "o", "linestyle": "-"},
    "linear, lumped": {"color": "C1", "marker": "s", "linestyle": "-"},
    "quadratic, consistent": {"color": "C2", "marker": "^", "linestyle": "--"},
    "quadratic, lumped": {"color": "C3", "marker": "v", "linestyle": "--"},
}
# Output times: after the first time step and at the end time.
times = [1e-6, 0.01]
t_first_step, t_end = times
dims = ["1D", "2D"]


def make_project(dim, suffix):
    """Variant of the base project file, same edits as the ctest `.xml` patches."""
    mesh_dir = Path("..") / f"PressureStep{dim}"
    prefix = f"pressure_step_{dim}{suffix}"
    prj = ot.Project(
        input_file=mesh_dir / f"pressure_step_{dim}.prj",
        # Not directly in out_dir: model.run() links the project folder into
        # out_dir/<prefix>/model/project, which would then point to its own
        # ancestor and form a directory loop.
        output_file=out_dir / "projects" / f"{prefix}.prj",
    )
    if "quadratic" in suffix:
        if dim == "1D":
            prj.replace_text("line_quadratic.vtu", "./meshes/mesh[1]")
        else:
            prj.replace_text("rectangle_quadratic.vtu", "./meshes/mesh[1]")
            prj.replace_text("rectangle_quadratic_left.vtu", "./meshes/mesh[2]")
            prj.replace_text(
                "rectangle_quadratic_left",
                "./process_variables/process_variable/boundary_conditions"
                "/boundary_condition/mesh",
            )
        prj.replace_text(3, "./processes/process/integration_order")
        prj.replace_text(2, "./process_variables/process_variable/order")
    if "mass_lumping" in suffix:
        prj.add_element("./processes/process", "mass_lumping", "true")
    prj.replace_text(prefix, "./time_loop/output/prefix")
    return prj, mesh_dir


# %% [markdown]
# ## Run the simulations


# %%
def run_simulation(dim, suffix):
    prj, mesh_dir = make_project(dim, suffix)
    prj.write_input()
    execution = ot.Execution(omp_num_threads=1, ogs_asm_threads=1)
    model = ot.Model(prj, meshes=mesh_dir, execution=execution)
    sim = model.run(out_dir / f"pressure_step_{dim}{suffix}", overwrite=True)
    # run() does not raise if OGS fails, which would let the code below
    # silently use stale or missing output.
    assert (
        sim.status == ot.Simulation.Status.done
    ), f"OGS failed for {dim}{suffix}, see {sim.log_file}."
    return sim.meshseries


series = {
    (dim, suffix): run_simulation(dim, suffix)
    for dim in dims
    for suffix in variants.values()
}


# %% [markdown]
# ## Read the results
#
# The 2D results are taken along the bottom edge $y = 0$. For quadratic
# elements with mass lumping the solution varies with $y$, because corner and
# mid-side nodes get different masses. The variation is shown separately
# below.


# %%
def output_index(dim, suffix, t):
    ms = series[(dim, suffix)]
    index = ms.closest_timestep(t)
    assert np.isclose(ms.timevalues[index], t), f"no output at t = {t:g} s"
    return index


def read_mesh(dim, suffix, t):
    return series[(dim, suffix)][output_index(dim, suffix, t)]


def nodal_pressure(dim, suffix, t):
    return series[(dim, suffix)].values("pressure")[output_index(dim, suffix, t)]


def read_pressure(dim, suffix, t):
    points = read_mesh(dim, suffix, t).points
    p = nodal_pressure(dim, suffix, t)
    if dim == "2D":
        on_bottom = np.isclose(points[:, 1], 0.0)
        points, p = points[on_bottom], p[on_bottom]
    order = np.argsort(points[:, 0])
    return points[order, 0], p[order]


def interpolate_along_edges(x_nodes, p_nodes, order, n_sub=20):
    """Evaluate the finite element field between the nodes of a line of
    element edges.

    Along an edge the field of Line2/Quad4 (order 1) or Line3/Quad8 (order 2)
    elements is the polynomial of that order through the nodes of the edge.

    :param x_nodes: Sorted coordinates of all nodes on the line, i.e. for
        order 2 corner, mid-edge, corner, ...
    :param p_nodes: The nodal values.
    :param order: Polynomial order of the elements.
    :param n_sub: Number of points per element.
    :returns: Coordinates and values of the interpolated field.
    """
    assert (len(x_nodes) - 1) % order == 0, "incomplete element on the line"
    s = np.linspace(0.0, 1.0, n_sub, endpoint=False)
    x, p = [], []
    for first in range(0, len(x_nodes) - order, order):
        x_e = x_nodes[first : first + order + 1]
        p_e = p_nodes[first : first + order + 1]
        length = x_e[-1] - x_e[0]
        # Fit in the local coordinate s in [0, 1] of the edge.
        coefficients = np.polyfit((x_e - x_e[0]) / length, p_e, order)
        x.append(x_e[0] + s * length)
        p.append(np.polyval(coefficients, s))
    x.append(x_nodes[-1:])
    p.append(p_nodes[-1:])
    return np.concatenate(x), np.concatenate(p)


def profile(mesh, y, x_max, order):
    """Pressure along the element edges at height y from x = 0 to x_max.

    Returns the pressure interpolated with the shape functions of the
    elements, i.e. piecewise quadratic for quadratic elements, and the nodal
    values on the line.
    """
    on_line = np.isclose(mesh.points[:, 1], y) & (mesh.points[:, 0] <= x_max + 1e-12)
    sorting = np.argsort(mesh.points[on_line, 0])
    x_nodes = mesh.points[on_line, 0][sorting]
    p_nodes = np.asarray(mesh["pressure"])[on_line][sorting]
    x, p = interpolate_along_edges(x_nodes, p_nodes, order)
    return x, p, x_nodes, p_nodes


def plot_profile(ax, mesh, y, x_max, order, offset=None, **style):
    """Plot the interpolated profile as line and the nodal values as markers.

    `offset(x)` is subtracted from both, e.g. an analytical solution.
    """
    marker = style.pop("marker", "o")
    x, p, x_nodes, p_nodes = profile(mesh, y, x_max, order)
    if offset is not None:
        p, p_nodes = p - offset(x), p_nodes - offset(x_nodes)
    (line,) = ax.plot(x, p, linewidth=1, **style)
    ax.plot(
        x_nodes,
        p_nodes,
        linestyle="none",
        marker=marker,
        markersize=4,
        color=line.get_color(),
    )


results = {
    (dim, name, t): read_pressure(dim, suffix, t)
    for dim in dims
    for name, suffix in variants.items()
    for t in times
}

# %% [markdown]
# ## Pressure profiles
#
# Top row: after the first time step, $t = 10^{-6}$ s. The analytical front
# is only about $2\sqrt{D t} = 2$ mm wide, which is a fifth of an element.
# Bottom row: at the end time, $t = 0.01$ s, all variants are close to the
# analytical solution, so their difference to it is shown. The lines show the
# pressure interpolated with the shape functions of the elements, the markers
# the nodal values.


# %%
def analytical(x, t):
    return p_0 * erfc(x / (2 * np.sqrt(D * t)))


x_max = {t_first_step: 0.05, t_end: 0.5}
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharey="row")
for row, t in enumerate(times):
    for col, dim in enumerate(dims):
        ax = axes[row, col]
        if row == 0:
            x_exact = np.linspace(0.0, x_max[t], 500)
            ax.plot(
                x_exact,
                analytical(x_exact, t),
                "k-",
                linewidth=1,
                label="analytical",
            )
        for name, suffix in variants.items():
            plot_profile(
                ax,
                read_mesh(dim, suffix, t),
                0.0,
                x_max[t],
                order=2 if "quadratic" in name else 1,
                offset=(lambda x, t=t: analytical(x, t)) if row == 1 else None,
                label=name,
                **styles[name],
            )
        ax.axhline(0.0, color="grey", linewidth=0.5)
        ax.set_title(f"{dim}, t = {t:g} s")
        ax.set_xlabel("x / m")
        ax.grid(True, alpha=0.3)
axes[0, 0].set_ylabel("pressure / Pa")
axes[1, 0].set_ylabel("p - p_analytical / Pa")
axes[0, 0].legend()
fig.tight_layout()

# %% [markdown]
# ## Undershoot after the first time step
#
# The smallest nodal pressure on the whole mesh shows how strong the
# non-physical oscillations are. Mass lumping removes them for linear
# elements. For quadratic elements it only makes them smaller: the
# stiffness matrix of quadratic elements has positive off-diagonal entries,
# so the scheme does not become monotone. Between the nodes the quadratic
# pressure field dips even further below zero, see the profiles above.

# %%
min_p = {
    (dim, name): nodal_pressure(dim, suffix, t_first_step).min()
    for dim in dims
    for name, suffix in variants.items()
}
print(f"{'variant':<24}" + "".join(f"{'min p, ' + dim + ' / Pa':>20}" for dim in dims))
for name in variants:
    print(f"{name:<24}" + "".join(f"{min_p[(dim, name)]:>20.4g}" for dim in dims))

for dim in dims:
    assert min_p[(dim, "linear, consistent")] < 0
    assert min_p[(dim, "linear, lumped")] >= 0
    assert min_p[(dim, "quadratic, consistent")] < min_p[(dim, "quadratic, lumped")] < 0

# %% [markdown]
# ## Variation along y in 2D (quadratic elements)
#
# Corner and mid-side nodes of Quad8 elements get different lumped masses.
# So after the first time step the pressure is not constant along $y$. The
# plot shows the nodal pressures on the element edge $x = 0.01$ m, where
# corner and mid-side nodes alternate along $y$. They are read with
# `MeshSeries.probe` at the node positions. The lines only join the nodal
# values.

# %%
x_probe = 0.01
y_probe = np.linspace(0.0, 0.02, 5)  # corner, mid-side, corner, ...
points = np.column_stack(
    [np.full_like(y_probe, x_probe), y_probe, np.zeros_like(y_probe)]
)
p_along_y = {}
fig, ax = plt.subplots(figsize=(8, 4))
for name in ["quadratic, consistent", "quadratic, lumped"]:
    ms = series[("2D", variants[name])]
    probes = ms.probe(points, "pressure")
    p_along_y[name] = probes.values("pressure")[ms.closest_timestep(t_first_step)]
    ax.plot(y_probe, p_along_y[name], linewidth=1, label=name, **styles[name])
ax.axhline(0.0, color="grey", linewidth=0.5)
ax.set_xlabel("y / m")
ax.set_ylabel("pressure / Pa")
ax.set_title(f"2D, x = {x_probe:g} m, t = 1e-06 s")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()

assert np.ptp(p_along_y["quadratic, consistent"]) < 1e-9 * p_0
assert np.ptp(p_along_y["quadratic, lumped"]) > 1e-2 * p_0

# %% [markdown]
# ## Check against the analytical solution at the end time
#
# At $t = 0.01$ s the front spans many elements. All variants should then
# agree with the analytical solution.

# %%
for dim in dims:
    for name in variants:
        x, p = results[(dim, name, t_end)]
        exact = analytical(x, t_end)
        error = np.max(np.abs(p - exact)) / p_0
        print(f"{dim} {name:<24} max relative error: {error:.2e}")
        assert error < 2e-3, f"{dim} {name}: error {error:.2e} too large"

# %% [markdown]
# The 1D and 2D benchmarks give the same solution, except for HRZ lumping on
# Quad8, where the solution varies along $y$.

# %%
for name in ["linear, consistent", "linear, lumped", "quadratic, consistent"]:
    for t in times:
        x_1D, p_1D = results[("1D", name, t)]
        x_2D, p_2D = results[("2D", name, t)]
        assert np.allclose(x_1D, x_2D)
        assert np.max(np.abs(p_1D - p_2D)) < 1e-6, f"{name}, t = {t:g} s: 1D != 2D"


# %% [markdown]
# ## Python implementation of the benchmarks
#
# The next two sections need quantities that OGS does not write: matrices
# with row-sum lumping on Quad8 elements, and the discrete flux through the
# Dirichlet boundary. A small finite element solver in Python provides them.
# It has the same discretisation as the OGS benchmarks: the same elements and
# Gauss points (2 per direction for linear, 3 for quadratic elements),
# backward Euler with the same time steps, and the same Dirichlet condition at
# $x = 0$. The equation is scaled by $1/S$, i.e. mass matrix
# $M = \int N^T N \, \mathrm{d}V$ and stiffness matrix
# $K = D \int \nabla N^T \nabla N \, \mathrm{d}V$.


# %%
def line2_shape(xi):
    return np.array([(1 - xi) / 2, (1 + xi) / 2]), np.array([[-0.5, 0.5]])


def line3_shape(xi):
    # VTK node order: end nodes first, then the middle node.
    N = np.array([0.5 * xi * (xi - 1), 0.5 * xi * (xi + 1), 1 - xi**2])
    return N, np.array([[xi - 0.5, xi + 0.5, -2 * xi]])


# Reference coordinates of the Quad4 and Quad8 nodes in VTK order: corners
# first, then the mid-side nodes of the edges 0-1, 1-2, 2-3, 3-0.
quad_nodes = np.array(
    [[-1, -1], [1, -1], [1, 1], [-1, 1], [0, -1], [1, 0], [0, 1], [-1, 0]],
    dtype=float,
)


def quad4_shape(xi, eta):
    a, b = quad_nodes[:4, 0], quad_nodes[:4, 1]
    N = 0.25 * (1 + a * xi) * (1 + b * eta)
    dN = np.array([0.25 * a * (1 + b * eta), 0.25 * b * (1 + a * xi)])
    return N, dN


def quad8_shape(xi, eta):
    """Serendipity shape functions and their derivatives w.r.t. (xi, eta)."""
    N = np.empty(8)
    dN = np.empty((2, 8))
    for i, (a, b) in enumerate(quad_nodes):
        if i < 4:
            N[i] = 0.25 * (1 + a * xi) * (1 + b * eta) * (a * xi + b * eta - 1)
            dN[0, i] = 0.25 * a * (1 + b * eta) * (2 * a * xi + b * eta)
            dN[1, i] = 0.25 * b * (1 + a * xi) * (a * xi + 2 * b * eta)
        elif a == 0:
            N[i] = 0.5 * (1 - xi**2) * (1 + b * eta)
            dN[0, i] = -xi * (1 + b * eta)
            dN[1, i] = 0.5 * b * (1 - xi**2)
        else:
            N[i] = 0.5 * (1 + a * xi) * (1 - eta**2)
            dN[0, i] = 0.5 * a * (1 - eta**2)
            dN[1, i] = -eta * (1 + a * xi)
    return N, dN


# cell type: (dimension, shape functions, Gauss points per direction)
element_types = {
    pv.CellType.LINE: (1, line2_shape, 2),
    pv.CellType.QUADRATIC_EDGE: (1, line3_shape, 3),
    pv.CellType.QUAD: (2, quad4_shape, 2),
    pv.CellType.QUADRATIC_QUAD: (2, quad8_shape, 3),
}


def element_matrices(cell_type, coordinates):
    dim, shape, n_gauss = element_types[cell_type]
    points, weights = np.polynomial.legendre.leggauss(n_gauss)
    n = coordinates.shape[0]
    M = np.zeros((n, n))
    K = np.zeros((n, n))
    for index in np.ndindex(*(dim * [n_gauss])):
        N, dN = shape(*points[list(index)])
        J = dN @ coordinates[:, :dim]
        dN_dx = np.linalg.solve(J, dN)
        dV = np.linalg.det(J) * np.prod(weights[list(index)])
        M += np.outer(N, N) * dV
        K += D * dN_dx.T @ dN_dx * dV
    return M, K


def lump(M, method):
    if method == "consistent":
        return M
    if method == "HRZ":
        diagonal = np.diag(M)
        return np.diag(M.sum() / diagonal.sum() * diagonal)
    if method == "row sum":
        return np.diag(M.sum(axis=1))
    msg = f"Unknown lumping method '{method}'."
    raise ValueError(msg)


def assemble(mesh, method):
    """Global consistent, lumped and stiffness matrices."""
    n = mesh.n_points
    rows, cols, values = [], [], {"M": [], "M_lumped": [], "K": []}
    for cell_type, cells in mesh.cells_dict.items():
        for cell in cells:
            M, K = element_matrices(cell_type, mesh.points[cell])
            rows.append(np.repeat(cell, len(cell)))
            cols.append(np.tile(cell, len(cell)))
            values["M"].append(M.ravel())
            values["M_lumped"].append(lump(M, method).ravel())
            values["K"].append(K.ravel())
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    return {
        name: sp.csr_matrix((np.concatenate(v), (rows, cols)), shape=(n, n))
        for name, v in values.items()
    }


def solve(mesh, method, time_steps):
    """Backward Euler.

    Returns the nodal pressures after every time step and the flux into the
    domain through the Dirichlet boundary during every time step.
    """
    matrices = assemble(mesh, method)
    M, K = matrices["M_lumped"], matrices["K"]
    fixed = np.isclose(mesh.points[:, 0], 0.0)
    free = ~fixed
    # As in OGS, the initial condition p = 0 also holds on the Dirichlet
    # nodes; the boundary value is imposed from the first time step on.
    p = np.zeros(mesh.n_points)
    p_fixed = np.full(fixed.sum(), p_0)
    history, inflow = [], []
    for dt, repeat in time_steps:
        A = (M / dt + K).tocsr()
        solver = spla.factorized(A[free][:, free].tocsc())
        A_free_fixed = A[free][:, fixed]
        for _ in range(repeat):
            p_old = p.copy()
            p[free] = solver((M / dt @ p_old)[free] - A_free_fixed @ p_fixed)
            p[fixed] = p_0
            # The residual of the Dirichlet rows is the boundary flux.
            residual = M @ (p - p_old) / dt + K @ p
            inflow.append(residual[fixed].sum())
            history.append(p.copy())
    return {"p": history, "inflow": np.array(inflow), **matrices}


# Time steps of the benchmarks: (delta_t, repeat).
time_steps = [(1e-6, 10), (1e-5, 99), (1e-4, 90)]
dt_steps = np.array([dt for dt, repeat in time_steps for _ in range(repeat)])
t_steps = np.cumsum(dt_steps)

meshes = {
    ("1D", "linear"): Path("..") / "PressureStep1D" / "line.vtu",
    ("1D", "quadratic"): Path("..") / "PressureStep1D" / "line_quadratic.vtu",
    ("2D", "linear"): Path("..") / "PressureStep2D" / "rectangle.vtu",
    ("2D", "quadratic"): Path("..") / "PressureStep2D" / "rectangle_quadratic.vtu",
}
# OGS lumps linear elements by row sums and higher order elements by HRZ.
python_variants = {
    "linear, consistent": ("linear", "consistent"),
    "linear, lumped": ("linear", "row sum"),
    "quadratic, consistent": ("quadratic", "consistent"),
    "quadratic, lumped": ("quadratic", "HRZ"),
}
python_results = {}
for dim in dims:
    for name, (order, method) in python_variants.items():
        mesh = ot.mesh.read(meshes[(dim, order)])
        python_results[(dim, name)] = (mesh, solve(mesh, method, time_steps))

# %% [markdown]
# The Python solver agrees with OGS for all eight benchmark variants up to
# round-off and the solver tolerance:

# %%
for (dim, name), (_mesh, result) in python_results.items():
    for step, t in [(0, t_first_step), (-1, t_end)]:
        ogs_p = nodal_pressure(dim, variants[name], t)
        difference = np.max(np.abs(result["p"][step] - ogs_p))
        print(
            f"{dim} {name:<24} t = {t:<8g} s max |p_Python - p_OGS| = {difference:.1e} Pa"
        )
        assert difference < 1e-3, f"{dim} {name}, t = {t:g} s: Python and OGS differ"

# %% [markdown]
# ## Row-sum lumping on Quad8 elements
#
# For linear elements OGS lumps the mass matrix by summing each row (the
# matrix is symmetric, so row and column sums are equal). Higher order
# elements use HRZ lumping instead. To show why, this section applies
# row-sum lumping to the 2D Quad8 mesh, which OGS cannot do.
#
# The row sums of the consistent Quad8 mass matrix are $-A/12$ for the corner
# nodes and $A/3$ for the mid-side nodes, with $A$ the element area. So
# row-sum lumping gives negative masses at the corner nodes:

# %%
quad8_mesh = ot.mesh.read(meshes[("2D", "quadratic")])
M_element, _ = element_matrices(
    pv.CellType.QUADRATIC_QUAD,
    quad8_mesh.points[quad8_mesh.cells_dict[pv.CellType.QUADRATIC_QUAD][0]],
)
element_area = M_element.sum()
print(
    "row sums / element area:",
    np.round(M_element.sum(axis=1) / element_area, 4),
)
assert np.allclose(M_element.sum(axis=1)[:4], -element_area / 12)
assert np.allclose(M_element.sum(axis=1)[4:], element_area / 3)

# For Line3 the row sums are positive and HRZ gives the same matrix as
# row-sum lumping.
line3_mesh = ot.mesh.read(meshes[("1D", "quadratic")])
M_line3, _ = element_matrices(
    pv.CellType.QUADRATIC_EDGE,
    line3_mesh.points[line3_mesh.cells_dict[pv.CellType.QUADRATIC_EDGE][0]],
)
assert np.allclose(lump(M_line3, "HRZ"), lump(M_line3, "row sum"))

quad8_results = {
    "consistent": python_results[("2D", "quadratic, consistent")][1],
    "HRZ": python_results[("2D", "quadratic, lumped")][1],
    "row sum": solve(quad8_mesh, "row sum", time_steps),
}

# %% [markdown]
# With negative masses the system matrix $M/\Delta t + K$ is no longer
# positive definite and backward Euler is no longer stable. In a scalar
# model with mass $m < 0$ and stiffness $k > 0$, one time step multiplies the
# solution by $m / (m + k \Delta t)$. Its magnitude is larger than 1 for
# $k \Delta t < 2 |m|$ and becomes very large near $k \Delta t = |m|$.
#
# The first time step looks fine, with no negative nodal pressures along
# $y = 0$ (left plot). Between the nodes the quadratic field still dips below
# zero, about as much as with HRZ lumping. The growth per step is small while
# $\Delta t = 10^{-6}$ s.
# After the switch to $\Delta t = 10^{-5}$ s the solution grows by about
# three orders of magnitude per step and reaches about $10^{283}$ Pa (right
# plot). With $\Delta t = 10^{-4}$ s, $k \Delta t > 2 |m|$ for the dominant
# mode and the solution decays again, to about $10^{266}$ Pa at $t = 0.01$ s.
# HRZ lumping keeps all masses positive and avoids this.

# %%
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
x_exact = np.linspace(0.0, 0.05, 500)
ax.plot(
    x_exact,
    analytical(x_exact, t_first_step),
    "k-",
    linewidth=1,
    label="analytical",
)
for method, color, marker in [
    ("consistent", "C2", "^"),
    ("HRZ", "C3", "v"),
    ("row sum", "C4", "D"),
]:
    mesh = quad8_mesh.copy()
    mesh.point_data["pressure"] = quad8_results[method]["p"][0]
    plot_profile(ax, mesh, 0.0, 0.05, order=2, color=color, marker=marker, label=method)
ax.axhline(0.0, color="grey", linewidth=0.5)
ax.set_xlabel("x / m")
ax.set_ylabel("pressure / Pa")
ax.set_title("2D Quad8, y = 0, t = 1e-06 s")
ax.grid(True, alpha=0.3)
ax.legend()

ax = axes[1]
# The row-sum values reach about 1e283 Pa; matplotlib's log scale cannot
# label such an axis, so the exponent is plotted on a linear axis.
for method, color in [("consistent", "C2"), ("HRZ", "C3"), ("row sum", "C4")]:
    ax.semilogx(
        t_steps,
        [np.log10(np.abs(p).max()) for p in quad8_results[method]["p"]],
        color=color,
        label=method,
    )
ax.set_xlabel("t / s")
ax.set_ylabel("log10(max |p| / Pa)")
ax.set_title("2D Quad8, largest nodal pressure")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()

# %%
max_p_row_sum = np.array([np.abs(p).max() for p in quad8_results["row sum"]["p"]])
peak = max_p_row_sum.argmax()
print(
    f"row sum lumping: max |p| is {max_p_row_sum[peak]:.2e} Pa at"
    f" t = {t_steps[peak]:g} s and {max_p_row_sum[-1]:.2e} Pa at t = {t_end:g} s"
)
assert max_p_row_sum[peak] > 1e280, "row-sum lumping on Quad8 was expected to blow up"
assert max_p_row_sum[-1] < 1e-10 * max_p_row_sum[peak], "expected decay for large dt"
# First time step: no negative nodal pressures along y = 0, but the quadratic
# field between the nodes dips below zero, about as much as with HRZ lumping.
first_step_min = {}
for method in ["HRZ", "row sum"]:
    mesh = quad8_mesh.copy()
    mesh.point_data["pressure"] = quad8_results[method]["p"][0]
    _, p, _, p_nodes = profile(mesh, 0.0, 1.0, order=2)
    first_step_min[method] = p.min(), p_nodes.min()
assert first_step_min["row sum"][1] >= 0
assert first_step_min["row sum"][0] < 0
assert 0.5 < first_step_min["row sum"][0] / first_step_min["HRZ"][0] < 2

# %% [markdown]
# ## Mass balance
#
# The fluid volume stored in the domain must equal the volume that flowed in
# through the left boundary. With the equation scaled by $1/S$, both are
# measured in Pa m$^3$ (Pa m in 1D); multiplying by $S$ gives m$^3$. Two
# quantities are compared:
#
# - Inflow $V_\text{in}(t) = \sum_n \Delta t_n \, q_n$, where $q_n$ is the
#   discrete flux through the Dirichlet boundary in step $n$, i.e. the
#   residual $M (p^n - p^{n-1}) / \Delta t_n + K p^n$ summed over the
#   Dirichlet nodes. These rows are not solved, because the pressure is
#   prescribed there. Their residual is the reaction flux that the boundary
#   condition has to supply to hold $p = p_0$. On all other nodes the
#   residual is zero, as the remaining boundaries are impermeable.
# - Stored volume of the finite element pressure field,
#   $V = \int p \, \mathrm{d}V = \mathbf{1}^T M_\text{consistent} \, p$.
#
# Summing all rows of the discrete equations gives
# $\mathbf{1}^T M (p^n - p^{n-1}) = \Delta t_n \, q_n$, because the rows of $K$
# sum to zero. With $p = 0$ at $t = 0$, also on the Dirichlet nodes, the sum
# over the time steps gives $V_\text{in} = \mathbf{1}^T M p$. So every scheme
# conserves the volume measured by its own mass matrix (lumped $M$ where
# lumping is used).
#
# Whether the pressure field itself conserves mass, $V = V_\text{in}$, depends
# on the mass matrix. Row-sum lumping keeps the row sums, so
# $\mathbf{1}^T M_\text{lumped} = \mathbf{1}^T M_\text{consistent}$ and the
# two volumes are the same. HRZ lumping changes them wherever HRZ and row-sum
# lumping differ, here for Quad8 but not for Line3.

# %%
V_in = {}
field_error = {}
for (dim, name), (mesh, result) in python_results.items():
    V_in[(dim, name)] = np.cumsum(dt_steps * result["inflow"])
    ones = np.ones(mesh.n_points)
    V = np.array([ones @ (result["M"] @ p) for p in result["p"]])
    field_error[(dim, name)] = np.abs(V - V_in[(dim, name)]) / V_in[(dim, name)]
    # Every scheme conserves the volume measured by its own mass matrix.
    V_own = np.array([ones @ (result["M_lumped"] @ p) for p in result["p"]])
    own_error = np.abs(V_own - V_in[(dim, name)]) / V_in[(dim, name)]
    assert own_error.max() < 1e-10, f"{dim} {name}: not conservative"


def print_errors(errors, header):
    print(f"{'':<27}{header:>36}")
    first, last = f"t = {t_first_step:g} s", f"t = {t_end:g} s"
    print(f"{'':<27}{'max':>12}{first:>12}{last:>12}")
    for (dim, name), e in errors.items():
        print(f"{dim} {name:<24}{e.max():>12.1e}{e[0]:>12.1e}{e[-1]:>12.1e}")


def plot_errors(errors, ylabel, title):
    fig, ax = plt.subplots(figsize=(8, 4))
    for (dim, name), e in errors.items():
        style = dict(styles[name], marker=None)
        if dim == "2D":
            style["linestyle"] = ":"
        # Floor at machine precision so that exact zeros show on the log axis.
        ax.loglog(
            t_steps,
            np.maximum(e, 1e-17),
            label=f"{dim} {name}",
            linewidth=1,
            **style,
        )
    ax.set_xlabel("t / s")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small")
    fig.tight_layout()


print_errors(field_error, "|V - V_in| / V_in")
plot_errors(field_error, "|V - V_in| / V_in", "Mass balance")

# %% [markdown]
# The mass balance of the pressure field, $|V - V_\text{in}| / V_\text{in}$,
# is at round-off level for all variants except HRZ lumping on Quad8.
# There, about 29 % of the inflow is missing from the pressure field after the
# first time step, when all of it sits in the first column of elements. The
# error shrinks to about $4 \cdot 10^{-5}$ at $t = 0.01$ s, as the
# pressure varies less within each element and HRZ moves less mass between
# corner and mid-side nodes.

# %%
for (dim, name), e in field_error.items():
    if name != "quadratic, lumped" or dim == "1D":
        assert e.max() < 1e-10, f"{dim} {name}: field not conservative"
hrz_quad8_error = field_error[("2D", "quadratic, lumped")]
assert 0.2 < hrz_quad8_error[0] < 0.4
assert 1e-5 < hrz_quad8_error[-1] < 1e-4

# %% [markdown]
# ## Boundary flux
#
# The mass balance only says that the stored volume matches the inflow, not
# that the inflow is right. Here the inflow $V_\text{in}$ is compared to the
# inflow of the semi-infinite domain. The analytical flux at $x = 0$ is
# $q(t) = -D \, \partial p / \partial x |_{x = 0} = p_0 \sqrt{D / (\pi t)}$.
# Integrated over time and multiplied by the cross section $A$ (1 in 1D, the
# height 0.02 m of the rectangle in 2D), it gives
# $V_\text{exact} = A \, p_0 \, 2 \sqrt{D t / \pi}$. This equals the stored
# volume $A \int_0^\infty p \, \mathrm{d}x$ of the analytical solution.

# %%
inflow_error = {}
for (dim, name), (mesh, _result) in python_results.items():
    cross_section = 1.0 if dim == "1D" else np.ptp(mesh.points[:, 1])
    V_exact = cross_section * p_0 * 2 * np.sqrt(D * t_steps / np.pi)
    inflow_error[(dim, name)] = np.abs(V_in[(dim, name)] - V_exact) / V_exact

print_errors(inflow_error, "|V_in - V_exact| / V_exact")
plot_errors(inflow_error, "|V_in - V_exact| / V_exact", "Boundary flux")

# %% [markdown]
# The difference is the discretisation error of the boundary flux. After the
# first step the analytical front is only about 2 mm wide, a fifth of an
# element. The discrete solution cannot store so little fluid: raising the
# Dirichlet nodes to $p_0$ alone already fills part of the first element. So
# the inflow is off by 33 % (quadratic, consistent) to 350 % (linear,
# lumped). Mass lumping does not help here. For linear elements it gives the
# largest error: all nodal pressures are non-negative, so at least the
# volume $p_0 h / 2$ of the Dirichlet node is stored. With the consistent
# mass matrix the negative pressures ahead of the front cancel part of it.
# At $t = 0.01$ s all variants are within about 0.1 %. The 1D and 2D curves
# coincide except for HRZ lumping on Quad8.

# %%
for dim in dims:
    first_step = {name: inflow_error[(dim, name)][0] for name in variants}
    # From 33 % (quadratic, consistent) to 350 % (linear, lumped).
    assert min(first_step, key=first_step.get) == "quadratic, consistent"
    assert max(first_step, key=first_step.get) == "linear, lumped"
    assert 0.3 < first_step["quadratic, consistent"] < 0.4
    assert 3 < first_step["linear, lumped"] < 4
    for name in variants:
        assert inflow_error[(dim, name)][-1] < 2e-3, f"{dim} {name}: inflow error"

# The 1D and 2D inflows coincide except for HRZ lumping on Quad8.
for name in variants:
    difference = np.max(np.abs(inflow_error[("1D", name)] - inflow_error[("2D", name)]))
    if name == "quadratic, lumped":
        assert difference > 0.1
    else:
        assert difference < 1e-10, f"{name}: 1D and 2D inflow differ"
