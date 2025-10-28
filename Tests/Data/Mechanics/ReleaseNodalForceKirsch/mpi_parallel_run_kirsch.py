# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python (.venv)
#     language: python
#     name: venv
# ---

# %% [raw]
# +++
# title = "Solving Kirsch's problem using the release nodal force approach: run with MPI parallel computing"
# date = "2025-08-27"
# author = "Wenqing Wang"
# image = "figures/kirsch_figure.png"
# web_subsection = "small-deformations"
# weight = 3
# +++

# %% editable=true slideshow={"slide_type": ""}
import os
from pathlib import Path
from subprocess import run

import gmsh
import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
from mesh import MeshGenerator

# %% [markdown]
# ### Run *"Solving Kirsch's problem using the release nodal force approach"* with MPI parallel computing

# %%
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
if not out_dir.exists():
    out_dir.mkdir(parents=True)
mesh_dir = Path(out_dir, "Mesh")

# %% [markdown]
# #### 1. Mesh generation

# %%
if not gmsh.isInitialized():
    gmsh.initialize()

gmsh.model.add("Mesh")

mesh_generator = MeshGenerator(gmsh_model=gmsh.model)
mesh_generator.generate_meshes(out_dir=mesh_dir, order=2)

gmsh.finalize()

# %% [markdown]
# #### 2. Partitioning meshes

# %%
bulk_mesh = Path(mesh_dir, "domain.vtu")
num_parts = 3
ot.cli().partmesh("-s", "-i", bulk_mesh, "-o", mesh_dir)
ot.cli().partmesh(
    "-m",
    "-n",
    str(num_parts),
    "-i",
    bulk_mesh,
    "-o",
    mesh_dir,
    "--",
    Path(mesh_dir, "left.vtu"),
    Path(mesh_dir, "right.vtu"),
    Path(mesh_dir, "top.vtu"),
    Path(mesh_dir, "bottom.vtu"),
    Path(mesh_dir, "arc.vtu"),
)

# %% [markdown]
# #### 3. Parallel computing

# %%
output_prefix = "kirsch_mpi"
temporary_project = Path(out_dir, f"{output_prefix}.prj")
prj = ot.Project(
    input_file="kirsch.prj",
    output_file=temporary_project,
)
prj.replace_text(
    output_prefix,
    xpath="./time_loop/output/prefix",
)
prj.write_input()

run(
    f"mpirun -np {num_parts} ogs {temporary_project} -m {mesh_dir} -o {out_dir}",
    shell=True,
    check=True,
)

# %% [markdown]
# #### 4. Contour plotting of the simulation results

# %%
pvd = Path(out_dir, f"{output_prefix}.pvd")
ms = ot.MeshSeries(pvd).scale(time=("s", "d"))
print(pvd)
mesh_last = ms[-1]
fig = mesh_last.plot_contourf(
    ot.variables.displacement["x"],
    figsize=(6, 4),
    fontsize=8,
    cmap="jet",
)

# %%
fig = mesh_last.plot_contourf(
    ot.variables.displacement["y"], figsize=(6, 4), fontsize=8, cmap="jet"
)

# %%
fig = mesh_last.plot_contourf(
    ot.variables.stress["xx"], figsize=(6, 4), fontsize=8, cmap="jet"
)

# %%
fig = mesh_last.plot_contourf(
    ot.variables.stress["yy"], figsize=(6, 4), fontsize=8, cmap="jet"
)

# %% [markdown]
# #### 5. Result comparison

# %%
point_a = (6.5, -857.0, 0)
point_b = (70.0, -857.0, 0)
profile = ms[-1].sample_over_line(point_a, point_b, resolution=30)
xs = profile["Distance"] + 6.5
plt.plot(
    xs, profile["sigma"][:, 0] * 1e-6, color="C2", label="Release nodal force (MPI)"
)

a = 6.5
sigma_t = -20  # MPa
sigma_x_a = np.asarray(
    [0.5 * (3.0 * a * a / (r * r) - 3 * a**4 / (r**4)) * sigma_t for r in xs],
)
plt.plot(xs, sigma_x_a, linestyle="dotted", color="C0", label="Analytical solution")
plt.xlabel("x / m")
plt.ylabel(r"stress $\sigma_{xx}$ / MPa")
plt.legend()
plt.grid()
plt.show()

# %% [markdown]
# In the above figure, the radial stress $\sigma_r$ profiles along the $\theta = 0^\circ$ axis, obtained using the release nodal force approach,
# and the analytical solution, are compared.

# %%
plt.plot(
    xs, profile["sigma"][:, 1] * 1e-6, color="C2", label="Release nodal force (MPI)"
)

a = 6.5
sigma_t = -20  # MPa
sigma_y_a = np.asarray(
    [0.5 * (2 + a * a / (r * r) + 3 * a**4 / (r**4)) * sigma_t for r in xs],
)
plt.plot(xs, sigma_y_a, linestyle="dotted", color="C0", label="Analytical solution")
plt.xlabel("x / m")
plt.ylabel(r"stress $\sigma_{yy}$ / MPa")
plt.legend()
plt.grid()
plt.show()

# %% [markdown]
# In the above two figures, the tangential stress $\sigma_{\theta}$ profiles along the $\theta = 0^\circ$ axis,
# obtained using the release nodal force approach, and the analytical solution, are compared.

# %%
expected_sigma = np.asarray(
    [
        [-5.94322287e05, -6.19135039e07, -1.27523478e07, 2.88857972e04],
        [-7.47350139e06, -3.61722449e07, -7.09372389e06, -4.31847097e03],
        [-7.02632921e06, -2.82291119e07, -4.57663234e06, -3.04634647e03],
        [-5.69811254e06, -2.49548130e07, -3.19587766e06, -2.33539317e03],
        [-4.52129471e06, -2.33397674e07, -2.35831864e06, -1.87907695e03],
        [-3.60618350e06, -2.24352157e07, -1.81241976e06, -1.49642711e03],
        [-2.90810489e06, -2.18773279e07, -1.43562983e06, -1.21999451e03],
        [-2.37358359e06, -2.15068967e07, -1.16414410e06, -1.00559463e03],
        [-1.95790335e06, -2.12450510e07, -9.60886316e05, -8.50316451e02],
        [-1.63102314e06, -2.10507797e07, -8.04540863e05, -7.21601922e02],
        [-1.36956366e06, -2.08999790e07, -6.80862806e05, -6.21075254e02],
        [-1.15718662e06, -2.07782330e07, -5.80625884e05, -5.47100152e02],
        [-9.83183398e05, -2.06769008e07, -4.98025268e05, -4.73813630e02],
        [-8.38590438e05, -2.05896257e07, -4.28464851e05, -2.91454098e02],
        [-7.18847960e05, -2.05193954e07, -3.71473022e05, -5.83689809e02],
        [-6.15597155e05, -2.04493120e07, -3.19472761e05, -2.60734404e02],
        [-5.26222722e05, -2.03832868e07, -2.72852858e05, -1.22193725e02],
        [-4.48460991e05, -2.03225868e07, -2.31314329e05, -1.40384160e02],
        [-3.81052951e05, -2.02633176e07, -1.93311164e05, 8.27105976e00],
        [-3.21267935e05, -2.02058999e07, -1.58150365e05, 1.25229397e02],
        [-2.68380912e05, -2.01489856e07, -1.25209966e05, -2.03819778e02],
        [-2.21121750e05, -2.00919675e07, -9.39267720e04, -5.98779759e02],
        [-1.79147229e05, -2.00335377e07, -6.38054888e04, 6.67763887e02],
        [-1.41291420e05, -1.99727511e07, -3.42127668e04, -2.09343375e02],
        [-1.07444543e05, -1.99093000e07, -5.02336097e03, -6.12146069e02],
        [-7.76109665e04, -1.98407209e07, 2.45004314e04, -3.90010177e02],
        [-5.18790616e04, -1.97699451e07, 5.34527661e04, -3.13328250e02],
        [-3.10211876e04, -1.96869485e07, 8.46090962e04, -3.15443054e01],
        [-1.35909887e04, -1.96029252e07, 1.15045144e05, 3.61051256e02],
        [-4.03024845e03, -1.94999513e07, 1.48805539e05, 7.61397500e02],
        [2.04509604e03, -1.93933612e07, 1.82605155e05, 2.68937637e03],
    ]
)

computed_sigma = np.asarray(profile["sigma"])

np.testing.assert_allclose(actual=computed_sigma, desired=expected_sigma, atol=1e-10)

# %%
profile_1 = ms[1].sample_over_line(point_a, point_b, resolution=30)

# %%
expected_sigma_at_time_step1 = np.array(
    [
        [-3.43936509e04, -2.24255500e07, -7.37983093e05, 1.67163179e03],
        [-4.32494293e05, -2.09358938e07, -4.10516429e05, -2.49911515e02],
        [-4.06616274e05, -2.04762218e07, -2.64851408e05, -1.76293199e02],
        [-3.29751883e05, -2.02867369e07, -1.84946624e05, -1.35150068e02],
        [-2.61649000e05, -2.01932736e07, -1.36476773e05, -1.08742879e02],
        [-2.08691175e05, -2.01409268e07, -1.04885403e05, -8.65987912e01],
        [-1.68293107e05, -2.01086417e07, -8.30804300e04, -7.06015339e01],
        [-1.37360161e05, -2.00872047e07, -6.73694502e04, -5.81941335e01],
        [-1.13304592e05, -2.00720516e07, -5.56068470e04, -4.92081280e01],
        [-9.43879134e04, -2.00608090e07, -4.65590777e04, -4.17593705e01],
        [-7.92571563e04, -2.00520821e07, -3.94017827e04, -3.59418550e01],
        [-6.69668185e04, -2.00450366e07, -3.36010349e04, -3.16608884e01],
        [-5.68971874e04, -2.00391725e07, -2.88209067e04, -2.74197703e01],
        [-4.85295392e04, -2.00341219e07, -2.47954196e04, -1.68665566e01],
        [-4.15999977e04, -2.00300576e07, -2.14972814e04, -3.37783455e01],
        [-3.56248353e04, -2.00260019e07, -1.84880070e04, -1.50887965e01],
        [-3.04527038e04, -2.00221809e07, -1.57900959e04, -7.07139610e00],
        [-2.59526036e04, -2.00186682e07, -1.33862459e04, -8.12408336e00],
        [-2.20516754e04, -2.00152383e07, -1.11869887e04, 4.78649290e-01],
        [-1.85918944e04, -2.00119155e07, -9.15222018e03, 7.24707160e00],
        [-1.55313028e04, -2.00086219e07, -7.24594709e03, -1.17951260e01],
        [-1.27963975e04, -2.00053222e07, -5.43557708e03, -3.46516064e01],
        [-1.03673165e04, -2.00019408e07, -3.69244727e03, 3.86437434e01],
        [-8.17658683e03, -1.99984231e07, -1.97990549e03, -1.21147786e01],
        [-6.21785551e03, -1.99947512e07, -2.90703760e02, -3.54251197e01],
        [-4.49137538e03, -1.99907825e07, 1.41784904e03, -2.25700334e01],
        [-3.00226051e03, -1.99866866e07, 3.09333137e03, -1.81324219e01],
        [-1.79520762e03, -1.99818836e07, 4.89635973e03, -1.82548063e00],
        [-7.86515552e02, -1.99770211e07, 6.65770509e03, 2.08941699e01],
        [-2.33231970e02, -1.99710620e07, 8.61143168e03, 4.40623553e01],
        [1.18350465e02, -1.99648936e07, 1.05674280e04, 1.55635206e02],
    ]
)


computed_sigma_t1 = np.asarray(profile_1["sigma"])

np.testing.assert_allclose(
    actual=computed_sigma_t1, desired=expected_sigma_at_time_step1, atol=1e-10
)

# %%
