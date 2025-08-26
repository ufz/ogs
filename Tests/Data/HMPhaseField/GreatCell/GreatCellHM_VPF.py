# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.1
#   kernelspec:
#     display_name: Python (.venv)
#     language: python
#     name: venv
# ---

# %% [raw] magic_args="[raw]"
# +++
# title = "A 2D GREAT cell benchmark suite simulated using hydro-mechanical variational phase-field process"
# date = "2025-04-03"
# author = "Mostafa Mollaali"
# image = "figures/hm_lie_bbar_stress_trace.png"
# web_subsection = "phase-field"
# weight = 3
# +++

# %% vscode={"languageId": "python"}
import os
import sys
from pathlib import Path
from subprocess import run

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pyvista as pv
from matplotlib import colormaps

mechanics_path = Path("../../LIE/Mechanics/GreatCelljupyterNotebook").resolve()
sys.path.insert(0, str(mechanics_path))
# Local modules
from mesh_generator import (  # noqa: E402
    mesh_GreatCell_intact,
    mesh_GreatCell_VPF,
)
from ogs_model_runner import SingleOGSModel  # noqa: E402
from Plotter import Plotter  # noqa: E402


# %% vscode={"languageId": "python"}
def truncated_cmap(name, minval=0.3, maxval=0.9, n=6):
    base = colormaps.get_cmap(name)
    return lambda i: base(minval + (maxval - minval) * i / (n - 1))


mpl.rcdefaults()
mpl.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "axes.labelsize": 26,
        "axes.titlesize": 24,
        "legend.fontsize": 18,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
    }
)

# %% vscode={"languageId": "python"}
ot.plot.setup.show_region_bounds = False

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
if not out_dir.exists():
    out_dir.mkdir(parents=True)

# %% [markdown]
# # Great cell
#
# The GREAT cell is a  poly-axial rock-testing device that reproduces subsurface conditions down to 3.5 km depth on 200 mm-diameter samples. It imposes a rotating stress field, injects fluid through a central borehole, and records both fiber-optic strain and pore-pressure data—providing a rich dataset for validating coupled hydro-mechanical models. For full details, see the GREAT cell benchmark docs:
#  [www.opengeosys.org/docs/benchmarks/small-deformations/greatcell/](https://www.opengeosys.org/docs/benchmarks/small-deformations/greatcell/)
#
#
#
# We distinguish two main types of benchmarks: **intact** and **fractured rock tests**. In the first type, an intact rock sample is subjected to rotating external stress conditions ($\texttt{M}_1$) and fluid is injected from a central borehole ($\texttt{HM}_1$):
#
# - **$\texttt{M}_1$**: Intact rock sample under rotating boundary conditions.
# - **$\texttt{HM}_1$**: Intact rock sample with fluid injection.
#
# For the fractured rock samples, two subtypes are also considered: static and propagating fractures. Static fractures are studied as full ($\texttt{M}_{2a}$, $\texttt{HM}_{2a}$) and half-fractured ($\texttt{M}_{2b}$, $\texttt{HM}_{2b}$) versions to have a symmetric and non-symmetric case. Again, mechanical and hydro-mechanical versions are studied, the latter with fluid injection into the rock fracture:
#
# - **$\texttt{M}_{2a}$**: Fully fractured sample.
# - **$\texttt{M}_{2b}$**: Half-fractured sample.
# - **$\texttt{HM}_{2a}$**: Fully fractured sample with inflow-outflow.
# - **$\texttt{HM}_{2b}$**: Half-fractured sample with inflow-outflow.


# %% [markdown]
# To run this benchmark, you need to have OGS built with PETSc and PIP support. The following steps outline the process:
# 1. **Configure & build OGS with PETSc & PIP**
#    ```bash
#    cmake -S ogs-source -B build-folder --preset release-petsc \
#      -DOGS_USE_PIP=ON
#    ```

# 2. **Run the benchmark**

#    ```bash
#    cd build-folder/release-petsc
#    ninja # OR make -j
#    ctest -R nb-HMPhaseField/GreatCell
#    ```

# 3. **Verify output files**

#    ```bash
#    ls build-folder/release-petsc/Tests/Data/HMPhaseField/GreatCell/GreatCell
#    ```
# %% [markdown]
# ---
# ## Material Properties

# %% [markdown]
#
# The material properties are provided in the following dictionary. The computational model incorporates two distinct elastic materials within its domain: a central circle ($r=0.097$ m) of rock surrounded by a rubber sheath in a 2D configuration.
#

# %% vscode={"languageId": "python"}
materials = {
    "Gneiss": {
        "young_sample": 83.9e9,  # Young's modulus (Pa)
        "nu_sample": 0.21,  # Poisson's ratio
        "biot": 0.6,  # Biot coefficient
        "porosity": 0.001,  # Porosity
        "permeability": 1e-19,  # Permeability (m²)
        "density_solid": 2750,  # Solid density (kg/m³)
        "k_n": 200e9,  # Normal stiffness (Pa/m)
        "k_t": 100e9,  # Tangential stiffness (Pa/m)
        "c_f": 4.4e-10,  # Fluid compressibility (Pa⁻¹)
        "k_s": 4.82e10,  # Solid bulk modulus (Pa)
        "S_f": 4.4e-10,  # Specific storage (Pa⁻¹)
        "t_np": 10e6,  # Peak normal traction (Pa)
        "Gc": 50,  # Fracture toughness (J/m²)
        "w_init": 1e-6,  # initial fracture width (m)
        "fluid": {
            "density": 1000.0,  # Fluid density (kg/m³)
            "viscosity": 1e-3,  # Fluid viscosity (Pa·s)
            "injectionFlowRate_Inlet": 4.167e-7,  # Injection flow rate (m³/s)
            "p_outlet": 3.45e6,  # Outlet pressure (Pa)
        },
        "rubber_sheath": {
            "young_modulus": 0.1e9,  # Young's modulus (Pa)
            "poisson_ratio": 0.4,  # Poisson's ratio
            "porosity": 0.001,  # Porosity
            "permeability": 1e-17,  # Permeability (m²)
            "density": 1500,  # Density (kg/m³)
            "biot": 0.0,  # Biot coefficient
        },
    },
    "Greywacke": {
        "young_sample": 26.87e9,  # Young's modulus (Pa)
        "nu_sample": 0.27,  # Poisson's ratio
        "biot": 0.8,  # Biot coefficient
        "porosity": 0.005,  # Porosity
        "permeability": 2.58e-19,  # Permeability (m²)
        "density_solid": 2650,  # Solid density (kg/m³)
        "k_n": 100e9,  # Normal stiffness (Pa/m)
        "k_t": 50e9,  # Tangential stiffness (Pa/m)
        "c_f": 4.4e-10,  # Fluid compressibility (Pa⁻¹)
        "k_s": 1.95e10,  # Solid bulk modulus (Pa)
        "S_f": 4.4e-10,  # Specific storage (Pa⁻¹)
        "t_np": 10e6,  # Peak normal traction (Pa)
        "Gc": 30,  # Fracture toughness (J/m²)
        "w_init": 1e-6,  # initial fracture width (m)
        "fluid": {
            "density": 1000.0,  # Fluid density (kg/m³)
            "viscosity": 1.0e-3,  # Fluid viscosity (Pa·s)
            "injectionFlowRate_Inlet": 4.167e-7,  # Injection flow rate (m³/s)
            "p_outlet": 3.45e6,  # Outlet pressure (Pa)
        },
        "rubber_sheath": {
            "young_modulus": 0.1e9,  # Young's modulus (Pa)
            "poisson_ratio": 0.4,  # Poisson's ratio
            "porosity": 0.001,  # Porosity
            "permeability": 1e-17,  # Permeability (m²)
            "density": 1500,  # Density (kg/m³)
            "biot": 0.0,  # Biot coefficient
        },
    },
}

material_names = list(materials.keys())

# %% [markdown]
# ---
# # Loads

# %% [markdown]
#
#
# ### Table: Load Conditions
#
# | Marker | PEE1 Angle to $\sigma_2$ | PEE1 & 1a | PEE2 & 2a | PEE3 & 3a | PEE4 & 4a | PEE5 & 5a | PEE6 & 6a | PEE7 & 7a | PEE8 & 8a |
# |--------|--------------------------|----------|----------|----------|----------|----------|----------|----------|----------|
# | A      | 0°                       | 10.0     | 6.64     | 4.46     | 1.17     | 1.0      | 3.82     | 7.80     | 9.95     |
# | E      | 22.5°                    | 9.95     | 10.0     | 6.64     | 4.46     | 1.17     | 1.0      | 3.82     | 7.80     |
# | B      | 45.0°                    | 7.80     | 9.95     | 10.0     | 6.64     | 4.46     | 1.17     | 1.0      | 3.82     |
# | F      | 67.5°                    | 3.82     | 7.80     | 9.95     | 10.0     | 6.64     | 4.46     | 1.17     | 1.0      |
# | C      | 90°                      | 1.0      | 3.82     | 7.80     | 9.95     | 10.0     | 6.64     | 4.46     | 1.17     |
#
# *All loads are in MPa. DSS loads are averages of adjacent PEEs.*

# %% [markdown]
# ### Schematic loading conditions

# %% vscode={"languageId": "python"}
loads = {
    "A": [
        10.0e6,
        6.64e6,
        4.46e6,
        1.17e6,
        1.0e6,
        3.82e6,
        7.80e6,
        9.95e6,
        10.0e6,
        6.64e6,
        4.46e6,
        1.17e6,
        1.0e6,
        3.82e6,
        7.80e6,
        9.95e6,
    ],
    "B": [
        7.80e6,
        9.95e6,
        10.0e6,
        6.64e6,
        4.46e6,
        1.17e6,
        1.0e6,
        3.82e6,
        7.80e6,
        9.95e6,
        10.0e6,
        6.64e6,
        4.46e6,
        1.17e6,
        1.0e6,
        3.82e6,
    ],
    "C": [
        1.0e6,
        3.82e6,
        7.80e6,
        9.95e6,
        10.0e6,
        6.64e6,
        4.46e6,
        1.17e6,
        1.0e6,
        3.82e6,
        7.80e6,
        9.95e6,
        10.0e6,
        6.64e6,
        4.46e6,
        1.17e6,
    ],
}

# %% vscode={"languageId": "python"}
angles_degrees = [
    0,
    -22.5,
    -45,
    -67.5,
    -90,
    -112.5,
    -135,
    -157.5,
    -180,
    -202.5,
    -225,
    -247.5,
    -270,
    -292.5,
    -315,
    -337.5,
]
angles_radians = np.deg2rad(angles_degrees)
circle_radius = 6
circle_x = circle_radius * np.cos(angles_radians)
circle_y = circle_radius * np.sin(angles_radians)


fig, axs = plt.subplots(1, 3, figsize=(21, 7), facecolor="none")

for idx, (label, values) in enumerate(loads.items()):
    ax = axs[idx]
    ax.set_aspect("equal")
    ax.axis("off")

    formatted_values = [rf"${v / 1e6:.1f}$" for v in values]
    scaled_values = [v / 2 for v in values]

    circle = plt.Circle(
        (0, 0),
        circle_radius,
        color="black",
        fill=False,
        linestyle="--",
        linewidth=2,
    )
    ax.add_artist(circle)

    top_points_x, top_points_y = [], []

    for i, (x, y, value) in enumerate(zip(circle_x, circle_y, scaled_values)):
        unit_vector = np.array([x, y]) / circle_radius
        line_end = np.array([x, y]) + unit_vector * value / 1e6

        top_points_x.append(line_end[0])
        top_points_y.append(line_end[1])

        ax.annotate(
            "",
            xytext=line_end,
            xy=(x, y),
            arrowprops={
                "arrowstyle": "-|>",
                "color": "blue",
                "lw": 2,
                "mutation_scale": 15,
                "fill": True,
            },
        )

        angle_offset = 1.2 * unit_vector
        angle_label = f"{-angles_degrees[i]}°"
        ax.text(
            x - angle_offset[0],
            y - angle_offset[1],
            angle_label,
            fontsize=12,
            ha="center",
            va="center",
            color="m",
        )

        value_offset = 0.5 * unit_vector
        ax.text(
            line_end[0] + 3.0 * value_offset[0],
            line_end[1] + 1.8 * value_offset[1],
            formatted_values[i],
            fontsize=18,
            ha="center",
            weight="bold",
        )

    top_points_x.append(top_points_x[0])
    top_points_y.append(top_points_y[0])
    ax.plot(
        top_points_x,
        top_points_y,
        color="green",
        linestyle="-.",
        lw=2,
        marker="o",
        markersize=6,
        markerfacecolor="lightgreen",
        markeredgewidth=1.0,
        markeredgecolor="black",
    )

    ax.set_xlim([-12, 12])
    ax.set_ylim([-12, 12])

    ax.text(
        0,
        0,
        f"Load {label}",
        fontsize=32,
        ha="center",
        va="center",
        family="serif",
    )

plt.tight_layout()
output_path = Path(out_dir, "loads_A_B_C_schematic_with_angles.png")
plt.savefig(output_path, dpi=350, bbox_inches="tight", transparent=True)
plt.show()

# %% [markdown]
# ---
#
# ## Boundary Conditions
#
# The boundary conditions applied in the simulation include both Dirichlet and Neumann conditions.
#
# - **Dirichlet conditions** (displacement control):
# \begin{equation*}
# \begin{cases}
# u_x(0, y) = 0, \quad u_y(0, y) = 0 & \quad \text{for } y = -0.09894 \text{ m}, \\
# u_y(x, 0) = 0,  & \quad \text{for } x = -0.09894 \text{ m}.
# \end{cases}
# \end{equation*}
#
# - **Neumann conditions**:
#   Normal stress $\sigma_{rr}$ is applied on PEEs and DSSs. The DSS stress is calculated as:
#   $$\sigma_\text{DSS}^i = \frac{\sigma_\text{PEE}^i + \sigma_\text{PEE}^{i+1}}{2}$$

# %% [markdown]
# ---
# # Mesh generation of intact samples

# %% [markdown]
# ### Input

# %% vscode={"languageId": "python"}
h = 0.005
meshname = "GreatCell"
mesh_path = Path(out_dir, "mesh_intact").resolve()
print(mesh_path)

# %% [markdown]
# ### Boundary meshes

# %% [markdown]
# #### Gmsh

# %% vscode={"languageId": "python"}
msh_file = mesh_GreatCell_intact(
    lc=0.005,
    lc2=h,
    r0=0.097,
    r1=0.094,
    r2=0.090,
    r3=0.065,
    out_dir=mesh_path,
    meshname=meshname,
    mode="BC",
)

# %% [markdown]
# #### Convert .msh to an OGS-compatible mesh

# %% vscode={"languageId": "python"}
msh_path = Path(mesh_path, f"{meshname}.msh")
meshes = ot.meshes_from_gmsh(filename=msh_path, dim=[1], reindex=True, log=False)

for name, mesh in meshes.items():
    print(f"{name}: {mesh.n_cells} cells")
    pv.save_meshio(Path(mesh_path, f"{name}.vtu"), mesh)

# %% [markdown]
# #### Visualization of boundary meshes

# %% vscode={"languageId": "python"}
plotter = pv.Plotter()
for name, mesh in meshes.items():
    if mesh.active_scalars is not None:
        plotter.add_mesh(
            mesh,
            scalars=mesh.active_scalars_name,
            cmap="tab20",
            show_edges=False,
            opacity=0.7,
        )
    else:
        plotter.add_mesh(mesh, show_edges=False, opacity=0.7, label=name)

    clean_name = name.replace("physical_group_", "")

    center = mesh.center
    direction = np.array(center) - np.array([0, 0, 0])
    direction[:2] = direction[:2] / (np.linalg.norm(direction[:2]) + 1e-8)
    offset = center + 0.025 * direction

    plotter.add_point_labels(
        [offset], [clean_name], font_size=12, point_size=0, text_color="black"
    )

plotter.view_xy()
plotter.enable_parallel_projection()
plotter.show()

# %% [markdown]
# ### Computational domain mesh
# #### Gmsh

# %% vscode={"languageId": "python"}
msh_file = mesh_GreatCell_intact(
    lc=0.005,
    lc2=h,
    r0=0.097,
    r1=0.094,
    r2=0.090,
    r3=0.065,
    out_dir=mesh_path,
    meshname=meshname,
    mode="domain",
)

# %% [markdown]
# #### Convert .msh to an OGS-compatible mesh

# %% vscode={"languageId": "python"}
msh_path = Path(mesh_path, f"{meshname}.msh")
meshes_volume = ot.meshes_from_gmsh(
    filename=msh_path, dim=[1, 2], reindex=True, log=False
)

for name, mesh in meshes_volume.items():
    print(f"{name}: {mesh.n_cells} cells")
    pv.save_meshio(Path(mesh_path, f"{name}.vtu"), mesh)

# %% [markdown]
# #### Visualization of computational domain mesh

# %% vscode={"languageId": "python"}
plotter = pv.Plotter()
for name, mesh in meshes_volume.items():
    if mesh.active_scalars is not None:
        plotter.add_mesh(
            mesh,
            scalars=mesh.active_scalars_name,
            cmap="Set1",
            show_edges=False,
            opacity=0.7,
        )
    else:
        plotter.add_mesh(mesh, show_edges=False, opacity=0.7, label=name)

    clean_name = name.replace("physical_group_", "")

    center = mesh.center
    direction = np.array(center) - np.array([0, 0, 0])
    direction[:2] = direction[:2] / (np.linalg.norm(direction[:2]) + 1e-8)
    offset = center + 0.025 * direction

    plotter.add_point_labels(
        [offset], [clean_name], font_size=12, point_size=0, text_color="black"
    )

plotter.view_xy()
plotter.enable_parallel_projection()
plotter.show()

# %% vscode={"languageId": "python"}
mesh_dir = Path(mesh_path).resolve()

run(
    [
        "NodeReordering",
        "-i",
        str(mesh_dir.joinpath("domain.vtu")),
        "-o",
        str(mesh_dir.joinpath("domain.vtu")),
    ],
    cwd=mesh_dir,
    check=True,
)

physical_groups = [
    "physical_group_DSS1.vtu",
    "physical_group_DSS1a.vtu",
    "physical_group_DSS2.vtu",
    "physical_group_DSS2a.vtu",
    "physical_group_DSS3.vtu",
    "physical_group_DSS3a.vtu",
    "physical_group_DSS4.vtu",
    "physical_group_DSS4a.vtu",
    "physical_group_DSS5.vtu",
    "physical_group_DSS5a.vtu",
    "physical_group_DSS6.vtu",
    "physical_group_DSS6a.vtu",
    "physical_group_DSS7.vtu",
    "physical_group_DSS7a.vtu",
    "physical_group_DSS8.vtu",
    "physical_group_DSS8a.vtu",
    "physical_group_PEE1.vtu",
    "physical_group_PEE1a.vtu",
    "physical_group_PEE2.vtu",
    "physical_group_PEE2a.vtu",
    "physical_group_PEE3.vtu",
    "physical_group_PEE3a.vtu",
    "physical_group_PEE4.vtu",
    "physical_group_PEE4a.vtu",
    "physical_group_PEE5.vtu",
    "physical_group_PEE5a.vtu",
    "physical_group_PEE6.vtu",
    "physical_group_PEE6a.vtu",
    "physical_group_PEE7.vtu",
    "physical_group_PEE7a.vtu",
    "physical_group_PEE8.vtu",
    "physical_group_PEE8a.vtu",
    "physical_group_p_bottom.vtu",
    "physical_group_p_left.vtu",
    "physical_group_p_right.vtu",
    "physical_group_p_top.vtu",
    "physical_group_Inlet.vtu",
]

group_paths = [str(mesh_dir.joinpath(name)) for name in physical_groups]

run(
    [
        "identifySubdomains",
        "-f",
        "-m",
        str(mesh_dir.joinpath("domain.vtu")),
        "--",
        *group_paths,
    ],
    cwd=mesh_dir,
    check=True,
)


# %% [markdown]
# ---
# # Intact Samples ($\texttt{M}_1$)
#
# The primary objective of the first benchmark exercise is to evaluate the **mechanical deformation** of different rock samples in the GREAT cell under various boundary conditions in a **2D plane strain setup**.
#
# This benchmark considers a **plane strain triaxial stress condition**, where the in-plane principal stresses are defined as:
#
# $$
# \sigma_2 = 10\ \text{MPa},\quad \sigma_3 = 1\ \text{MPa}
# $$
#
# The out-of-plane stress component, $\sigma_1$, is governed by the plane strain constraint:
#
# $$
# \sigma_1 = \nu (\sigma_2 + \sigma_3)
# $$
#
# This loading condition is applied to both **Greywacke** and **Gneiss** samples

# %% [markdown]
# ## Run the simulation

# %% [markdown]
# ### Inputs

# %% vscode={"languageId": "python"}
# Times for load curves
times = "0.0  1000. 3500"
simulation_end_time = 3500.0
n_fracture_p_ncs = 0
model_type = "M1"
output_prefix = "M1_VPF"
# Load
PEE_load_values = {
    "A": [10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6, 7.80e6, 9.95e6],
    "B": [7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6],
    "C": [1.0e6, 3.82e6, 7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6],
}

prj_file = Path("M1_VPF.prj")
prj = ot.Project(input_file=prj_file, output_file=Path(out_dir, f"{output_prefix}.prj"))

# %% [markdown]
# ### Run OGS

# %% vscode={"languageId": "python"}
# Create SingleOGSModel
sing_ogs_model = SingleOGSModel(
    model=prj,
    out_dir=out_dir,
    mesh_path=mesh_path,
    output_prefix=output_prefix,
    method="VPF",
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    materials=materials,
)

# Run simulations
vtu_files_dict = sing_ogs_model.run_simulations_with_fracture(
    times=times,
    base_project_file=prj_file,
    mesh_path=mesh_path,
    load_cases=PEE_load_values,
    material_names=material_names,
    materials=materials,
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    output_prefix=output_prefix,
    out_dir=out_dir,
    method="VPF",
    fracture_model_type="VolumetricDeviatoric",
    mesh_size=h,
)

# %% [markdown]
# ## Post-processing

# %% [markdown]
# #### Volumetric strain vs angle at probe circle

# %% vscode={"languageId": "python"}
json_path = Path("./external_data_dict.py").resolve()
print(f"[DEBUG] Trying path: {json_path}")

if json_path.exists():
    from external_data_dict import external_data
else:
    print("[WARNING] External data dict not found! Skipping...")
    external_data = None

plotter = Plotter(
    output_dir=out_dir,
    save_extracted_data=True,
)

data_dir = Path("external_data")
external_data = Plotter.load_external_data(data_dir, benchmark_tag="M1")

plotter.plot_volumetric_strain_vs_angle(
    vtu_files_dict,
    model_type="M1",
    ylim_range=[-7.5, 2.5],
    layout="subplots",
    external_data=external_data["strain"],
)

# %% [markdown]
# #### Profiles

# %% vscode={"languageId": "python"}
plotter.plot_field_variables(vtu_files_dict)

# %% [markdown]
# ---
# # Intact rock sample with fluid injection ($\texttt{HM}_1$)
#
#  This benchmark, which does not consider any fracture, is designed to verify the basic computational setting for hydro-mechanical simulations. In addition to the mechanical loads, a zero constant pore pressure is prescribed at the outer boundary. At the center of sample, fluid is injected at a rate of $Q_0^{\text{v}} = 2.085 \times 10^{-9}$ m$^3$/s
#
#  The hydro-mechanical simulations follows a two-stage process: a 3000~s equilibrium phase under mechanical loading to stabilize initial conditions, followed by a 500~s fluid injection phase to model fluid flow. This loading condition is applied to both **Greywacke** and **Gneiss** samples
#

# %% [markdown]
# ## Run the simulation
# ### Inputs

# %% vscode={"languageId": "python"}
# Times for load curves
times = "0.0  1000. 3500"
simulation_end_time = 3500.0
n_fracture_p_ncs = 0
model_type = "HM1"
output_prefix = "HM1_VPF"
# Load
PEE_load_values = {
    "A": [10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6, 7.80e6, 9.95e6],
    "B": [7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6],
    "C": [1.0e6, 3.82e6, 7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6],
}

for p in materials.values():
    p["fluid"][
        "injectionFlowRate_Inlet"
    ] = 2.085e-6  # Injection flow rate m³/s (125 ml/min)

prj_file = Path("HM1_VPF.prj")
prj = ot.Project(input_file=prj_file, output_file=Path(out_dir, f"{output_prefix}.prj"))

# %% [markdown]
# ### Run OGS

# %% vscode={"languageId": "python"}
# Create SingleOGSModel
sing_ogs_model = SingleOGSModel(
    model=prj,
    out_dir=out_dir,
    mesh_path=mesh_path,
    output_prefix=output_prefix,
    method="VPF",
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    materials=materials,
)

# Run simulations
vtu_files_dict_HM = sing_ogs_model.run_simulations_with_fracture(
    times=times,
    base_project_file=prj_file,
    mesh_path=mesh_path,
    load_cases=PEE_load_values,
    material_names=material_names,
    materials=materials,
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    output_prefix=output_prefix,
    out_dir=out_dir,
    method="VPF",
    fracture_model_type="VolumetricDeviatoric",
    mesh_size=h,
)

# %% [markdown]
# ## Post-processing

# %% [markdown]
# ### Volumetric strain vs angle at probe circle

# %% vscode={"languageId": "python"}
json_path = Path("./external_data_dict.py").resolve()
print(f"[DEBUG] Trying path: {json_path}")

if json_path.exists():
    from external_data_dict import external_data
else:
    print("[WARNING] External data dict not found! Skipping...")
    external_data = None

data_dir = Path("external_data")
external_data = Plotter.load_external_data(data_dir, benchmark_tag="HM1")

plotter = Plotter(
    output_dir=out_dir,
    save_extracted_data=True,
)

plotter.plot_volumetric_strain_vs_angle(
    vtu_files_dict_HM,
    model_type="HM1",
    ylim_range=[-7.5, 2.5],
    layout="subplots",
    external_data=external_data["strain"],
)

# %% [markdown]
# ### Profiles
# %% vscode={"languageId": "python"}
plotter.plot_field_variables(vtu_files_dict_HM)

# %% [markdown]
# ---
# # Mesh generation of fractured samples

# %% [markdown]
# ### Input

# %% vscode={"languageId": "python"}
h = 0.001
meshname = "GreatCell"
mesh_path_VPF = Path(out_dir, "mesh_GreatCell_VPF").resolve()

# %% [markdown]
# ### Boundary meshes

# %% [markdown]
# #### Gmsh

# %% vscode={"languageId": "python"}
msh_file_VPF = mesh_GreatCell_VPF(
    lc2=h,
    lc=20 * h,
    r0=0.097,
    r1=0.094,
    r2=0.090,
    r3=0.065,
    out_dir=mesh_path_VPF,
    meshname=meshname,
    mode="BC",
)

# %% [markdown]
# #### Convert .msh to an OGS-compatible mesh

# %% vscode={"languageId": "python"}
msh_path_VPF = mesh_path_VPF / f"{meshname}.msh"
meshes_volume_VPF = ot.meshes_from_gmsh(
    filename=msh_path_VPF, dim=[0, 1], reindex=True, log=False
)

for name, mesh in meshes_volume_VPF.items():
    print(f"{name}: {mesh.n_cells} cells")
    pv.save_meshio(mesh_path_VPF / f"{name}.vtu", mesh)

# %% [markdown]
# #### Visualization of boundary meshes

# %% vscode={"languageId": "python"}
plotter = pv.Plotter()
for name, mesh in meshes_volume_VPF.items():
    if mesh.active_scalars is not None:
        plotter.add_mesh(
            mesh,
            scalars=mesh.active_scalars_name,
            cmap="tab20",
            show_edges=False,
            opacity=0.7,
        )
    else:
        plotter.add_mesh(mesh, show_edges=False, opacity=0.7, label=name)

    clean_name = name.replace("physical_group_", "")

    center = mesh.center
    direction = np.array(center) - np.array([0, 0, 0])
    direction[:2] = direction[:2] / (np.linalg.norm(direction[:2]) + 1e-8)
    offset = center + 0.025 * direction

    plotter.add_point_labels(
        [offset], [clean_name], font_size=12, point_size=0, text_color="black"
    )

plotter.view_xy()
plotter.enable_parallel_projection()
plotter.show()


# %% [markdown]
# ### Computational domain mesh
# #### Gmsh

# %% vscode={"languageId": "python"}
msh_file_VPF = mesh_GreatCell_VPF(
    lc2=h,
    lc=20 * h,
    r0=0.097,
    r1=0.094,
    r2=0.090,
    r3=0.065,
    out_dir=mesh_path_VPF,
    meshname=meshname,
    mode="domain",
)

# %% [markdown]
# #### Convert .msh to an OGS-compatible mesh

# %% vscode={"languageId": "python"}
msh_path_VPF = mesh_path_VPF / f"{meshname}.msh"
meshes_volume_VPF = ot.meshes_from_gmsh(
    filename=msh_path_VPF, dim=[2], reindex=True, log=False
)

for name, mesh in meshes_volume_VPF.items():
    print(f"{name}: {mesh.n_cells} cells")
    pv.save_meshio(mesh_path_VPF / f"{name}.vtu", mesh)

# %% [markdown]
# #### Visualization of computational domain mesh

# %% vscode={"languageId": "python"}
plotter = pv.Plotter()
for name, mesh in meshes_volume_VPF.items():
    if mesh.active_scalars is not None:
        plotter.add_mesh(
            mesh,
            scalars=mesh.active_scalars_name,
            cmap="Set1",
            show_edges=False,
            opacity=0.7,
        )
    else:
        plotter.add_mesh(mesh, show_edges=False, opacity=0.7, label=name)

    clean_name = name.replace("physical_group_", "")

    center = mesh.center
    direction = np.array(center) - np.array([0, 0, 0])
    direction[:2] = direction[:2] / (np.linalg.norm(direction[:2]) + 1e-8)
    offset = center + 0.025 * direction

    plotter.add_point_labels(
        [offset], [clean_name], font_size=12, point_size=0, text_color="black"
    )

plotter.view_xy()
plotter.enable_parallel_projection()
plotter.show()


# %% vscode={"languageId": "python"}
mesh_dir = Path(mesh_path_VPF).resolve()

run(
    [
        "NodeReordering",
        "-i",
        str(mesh_dir.joinpath("domain.vtu")),
        "-o",
        str(mesh_dir.joinpath("domain.vtu")),
    ],
    cwd=mesh_dir,
    check=True,
)

physical_groups = [
    "physical_group_DSS1.vtu",
    "physical_group_DSS1a.vtu",
    "physical_group_DSS2.vtu",
    "physical_group_DSS2a.vtu",
    "physical_group_DSS3.vtu",
    "physical_group_DSS3a.vtu",
    "physical_group_DSS4.vtu",
    "physical_group_DSS4a.vtu",
    "physical_group_DSS5.vtu",
    "physical_group_DSS5a.vtu",
    "physical_group_DSS6.vtu",
    "physical_group_DSS6a.vtu",
    "physical_group_DSS7.vtu",
    "physical_group_DSS7a.vtu",
    "physical_group_DSS8.vtu",
    "physical_group_DSS8a.vtu",
    "physical_group_PEE1.vtu",
    "physical_group_PEE1a.vtu",
    "physical_group_PEE2.vtu",
    "physical_group_PEE2a.vtu",
    "physical_group_PEE3.vtu",
    "physical_group_PEE3a.vtu",
    "physical_group_PEE4.vtu",
    "physical_group_PEE4a.vtu",
    "physical_group_PEE5.vtu",
    "physical_group_PEE5a.vtu",
    "physical_group_PEE6.vtu",
    "physical_group_PEE6a.vtu",
    "physical_group_PEE7.vtu",
    "physical_group_PEE7a.vtu",
    "physical_group_PEE8.vtu",
    "physical_group_PEE8a.vtu",
    "physical_group_p_bottom.vtu",
    "physical_group_p_left.vtu",
    "physical_group_p_right.vtu",
    "physical_group_p_top.vtu",
    "physical_group_Inlet.vtu",
    "physical_group_Outlet_R_embeddedFracture.vtu",
    "physical_group_Outlet_R_fullFracture.vtu",
    "physical_group_Outlet_L_fullFracture.vtu",
]

group_paths = [str(mesh_dir.joinpath(name)) for name in physical_groups]

run(
    [
        "identifySubdomains",
        "-f",
        "-m",
        str(mesh_dir.joinpath("domain.vtu")),
        "--",
        *group_paths,
    ],
    cwd=mesh_dir,
    check=True,
)

# %% [markdown]
# ---
# # Half fractured sample ($\texttt{M}_{2b}$)
#
# In this section, all input data and boundary conditions provided in $\texttt{M}_{2b}$, except  a single wing fracture, defined as $\Gamma =  \left[0, 0.04\right] \times \{0\}$, is considered.

# %% [markdown]
# ## Run the simulation
# ### Inputs

# %% vscode={"languageId": "python"}
# Times for load curves
times = "0.0  1000. 3500"
simulation_end_time = 3500.0
n_fracture_p_ncs = 3
model_type = "M2b"
output_prefix = "M2b_VPF"

# Load
PEE_load_values = {
    "A": [10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6, 7.80e6, 9.95e6],
    "B": [7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6],
    "C": [1.0e6, 3.82e6, 7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6],
}

# Project file
prj_file = Path("M2_VPF.prj")
prj = ot.Project(input_file=prj_file, output_file=Path(out_dir, f"{output_prefix}.prj"))

# %% [markdown]
# ### Run OGS

# %% vscode={"languageId": "python"}

sing_ogs_model = SingleOGSModel(
    model=prj,
    out_dir=out_dir,
    mesh_path=mesh_path_VPF,
    output_prefix=output_prefix,
    method="VPF",
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    materials=materials,
)

vtu_files_dict_embedded = sing_ogs_model.run_simulations_with_fracture(
    times=times,
    base_project_file=prj_file,
    mesh_path=mesh_path_VPF,
    load_cases=PEE_load_values,
    material_names=material_names,
    materials=materials,
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    output_prefix=output_prefix,
    out_dir=out_dir,
    method="VPF",
    crack_type="half",
    fracture_model_type="VolumetricDeviatoric",
    mesh_size=h,
)

# %% [markdown]
# ## Post-processing

# %% vscode={"languageId": "python"}
plotter = Plotter(
    output_dir=out_dir,
    save_extracted_data=True,
)

data_dir = Path("external_data")
external_data = Plotter.load_external_data(data_dir, benchmark_tag="M2b")

plotter.plot_volumetric_strain_vs_angle(
    vtu_files_dict_embedded,
    model_type="M2b",
    ylim_range=[-7.5, 2.5],
    layout="single",
    external_data=external_data["strain"],
)

plotter.plot_field_variables(vtu_files_dict_embedded)


# %% [markdown]
# ---
# # Fully fractured sample ($\texttt{M}_{2a}$)
#
# In this benchmark,  2D plane strain numerical simulations are performed  to establish a baseline model for assessing the impact of fracture orientation on strain distribution under poly-axial loading. A planar fracture within the specimen, $\Gamma =  \left[-0.094, 0.094\right] \times \{0\}$, is considered, under poly-axial loading applied on PEE's and DSS's


# %% [markdown]
# ## Run the simulation
# ### Inputs

# %% vscode={"languageId": "python"}
# Times for load curves
times = "0.0  1000. 3500"
simulation_end_time = 3500.0
n_fracture_p_ncs = 3
model_type = "M2a"
output_prefix = "M2a_VPF"
# Load
PEE_load_values = {
    "A": [10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6, 7.80e6, 9.95e6],
    "B": [7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6],
    "C": [1.0e6, 3.82e6, 7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6],
}

# Project file
prj_file = Path("M2_VPF.prj")
prj = ot.Project(input_file=prj_file, output_file=Path(out_dir, f"{output_prefix}.prj"))

# %% [markdown]
# ### Run OGS

# %% vscode={"languageId": "python"}
# Create SingleOGSModel
sing_ogs_model = SingleOGSModel(
    model=prj,
    out_dir=out_dir,
    mesh_path=mesh_path_VPF,
    output_prefix=output_prefix,
    method="VPF",
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    materials=materials,
)

# Run simulations
vtu_files_dict_full = sing_ogs_model.run_simulations_with_fracture(
    times=times,
    base_project_file=prj_file,
    mesh_path=mesh_path_VPF,
    load_cases=PEE_load_values,
    material_names=material_names,
    materials=materials,
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    output_prefix=output_prefix,
    out_dir=out_dir,
    method="VPF",
    crack_type="full",
    fracture_model_type="VolumetricDeviatoric",
    mesh_size=h,
)


# %% [markdown]
# ## Post-processing

# %% vscode={"languageId": "python"}
plotter = Plotter(
    output_dir=out_dir,
    save_extracted_data=True,
)

data_dir = Path("external_data")
external_data = Plotter.load_external_data(data_dir, benchmark_tag="M2a")

plotter.plot_volumetric_strain_vs_angle(
    vtu_files_dict_full,
    model_type="M2a",
    ylim_range=[-7.5, 2.5],
    layout="subplots",
    external_data=external_data["strain"],
)

# %% vscode={"languageId": "python"}
plotter.plot_field_variables(vtu_files_dict_full)

# %% [markdown]
# ---
# # Half fractured sample with fluid injection ($\texttt{HM}_{2b}$)
#
#  We use all inputs and BCs from $\texttt{M}_{2b}$ for a single wing fracture $\Gamma = [0,\,0.04]\times\{0\}.$
#  A Dirichlet BC of $p(0.04,0)=3.45\text{ MPa}$
#  is applied at the fracture's right tip. After a 3000 s equilibration, fluid is injected at the left tip with
#  $ Q_0^{\mathrm{v}}(0,0)=4.167\times10^{-7}\,\mathrm{m}^3/\mathrm{s}\;(25\text{ ml/min}).$

# %% [markdown]
# ## Run the simulation
# ### Inputs

# %% vscode={"languageId": "python"}
# Times for load curves
times = "0.0  1000. 3500"
simulation_end_time = 3500.0
n_fracture_p_ncs = 3
model_type = "HM2b"
output_prefix = "HM2b_VPF"

# Load
PEE_load_values = {
    "A": [10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6, 7.80e6, 9.95e6],
    "E": [9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6, 7.80e6],
    "B": [7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6],
    "F": [3.82e6, 7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6],
    "C": [1.0e6, 3.82e6, 7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6],
}

for p in materials.values():
    p["fluid"][
        "injectionFlowRate_Inlet"
    ] = 4.167e-7  # Injection flow rate m³/s (25 ml/min)

# Project file
prj_file = Path("HM2b_VPF.prj")
prj = ot.Project(input_file=prj_file, output_file=Path(out_dir, f"{output_prefix}.prj"))

# %% [markdown]
# ### Run OGS

# %% vscode={"languageId": "python"}
sing_ogs_model = SingleOGSModel(
    model=prj,
    out_dir=out_dir,
    mesh_path=mesh_path_VPF,
    output_prefix=output_prefix,
    method="VPF",
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    materials=materials,
)

vtu_files_dict_embedded_HM = sing_ogs_model.run_simulations_with_fracture(
    times=times,
    base_project_file=prj_file,
    mesh_path=mesh_path_VPF,
    load_cases=PEE_load_values,
    material_names=material_names,
    materials=materials,
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    output_prefix=output_prefix,
    out_dir=out_dir,
    method="VPF",
    crack_type="half",
    fracture_model_type="Isotropic",
    mesh_size=h,
)

# %% [markdown]
# ## Post-processing

# %% vscode={"languageId": "python"}
plotter = Plotter(
    output_dir=out_dir,
    save_extracted_data=True,
)

data_dir = Path("external_data")
external_data = Plotter.load_external_data(data_dir, benchmark_tag="HM2b")

plotter.plot_volumetric_strain_vs_angle(
    vtu_files_dict_embedded_HM,
    model_type="HM2b",
    ylim_range=[-7.5, 2.5],
    layout="subplots",
    external_data=external_data["strain"],
)

# %% vscode={"languageId": "python"}
plotter.plot_field_variables(vtu_files_dict_embedded_HM)


# %% vscode={"languageId": "python"}
plotter.material_names = ["Gneiss", "Greywacke"]
plotter.vtu_file_names = {"VPF": vtu_files_dict_embedded_HM}

plotter.plot_avg_width_vs_stress(
    benchmark_tag="HM2b",
    metric="width",
    methods_to_include=["VPF"],
    pee_load_values=PEE_load_values,
    external_data=external_data["average"],
    ylim_range=(0, 10e-6),
)

plotter.plot_avg_width_vs_stress(
    benchmark_tag="HM2b",
    metric="permeability",
    methods_to_include=["VPF"],
    pee_load_values=PEE_load_values,
    external_data=external_data["average"],
    ylim_range=(0, 5e-12),
)

# %% vscode={"languageId": "python"}
vpf_profiles = {}
for load_case, file_list in plotter.vtu_file_names["VPF"].items():
    vpf_profiles[load_case] = plotter.extract_vpf_width_from_list(
        file_list,
        plotter.material_names,
    )

plotter.plot_fracture_aperture_profiles(
    widthProfile=vpf_profiles,
    benchmark_tag="HM2b",
    downsample=1,
    ylim=(0, 10e-6),
    method_label="VPF",
    external_data=external_data["widthProfile"],
)

# %% [markdown]
# ---
# # Fully fractured sample with fluid injection ($\texttt{HM}_{2a}$)
#
#  We performed five hydro-mechanical (HMca) simulations with different angles between the second principal stress direction (σ₂ = 10 MPa) and the PEE 5 fracture plane
#  $\Gamma = [-0.094,\,0.094]\times\{0\}.$
#  The initial pore pressure was uniform at $p_0 = 0.1$ MPa, with Dirichlet conditions
#  $p(-0.094,0)=p(0.094,0)=3.45\text{ MPa}.$
#  After equilibrating for 3000 s, we injected fluid at the center at a rate
#  $
#  Q_0^{\text{v}}(0,0)=4.167\times10^{-7}\,\mathrm{m}^3/\mathrm{s}\;(25\,\mathrm{ml/min})$
#  for 500 s, starting from an initial fracture aperture $w_{\text{ini}}=1\times10^{-6}$ m.

# %% [markdown]
# ## Run the simulation
# ### Inputs

# %% vscode={"languageId": "python"}
# Times for load curves
times = "0.0  1000. 3500"
simulation_end_time = 3500.0
n_fracture_p_ncs = 3
model_type = "HM2a"
output_prefix = "HM2a_VPF"

# Load
PEE_load_values = {
    "A": [10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6, 7.80e6, 9.95e6],
    # "E": [9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6, 7.80e6],
    # "B": [7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6],
    # "F": [3.82e6, 7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6],
    # "C": [1.0e6, 3.82e6, 7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6],
}

# To save time, we only run the first load case A.
# If you want to run all load cases, uncomment the lines above.

for p in materials.values():
    p["fluid"][
        "injectionFlowRate_Inlet"
    ] = 4.167e-7  # Injection flow rate m³/s (25 ml/min)

# Project file
prj_file = Path("HM2a_VPF.prj")
prj = ot.Project(input_file=prj_file, output_file=Path(out_dir, f"{output_prefix}.prj"))

# %% [markdown]
# ### Run OGS

# %% vscode={"languageId": "python"}
sing_ogs_model = SingleOGSModel(
    model=prj,
    out_dir=out_dir,
    mesh_path=mesh_path_VPF,
    output_prefix=output_prefix,
    method="VPF",
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    materials=materials,
)

vtu_files_dict_full_HM = sing_ogs_model.run_simulations_with_fracture(
    times=times,
    base_project_file=prj_file,
    mesh_path=mesh_path_VPF,
    load_cases=PEE_load_values,
    material_names=material_names,
    materials=materials,
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    output_prefix=output_prefix,
    out_dir=out_dir,
    method="VPF",
    crack_type="full",
    fracture_model_type="Isotropic",
    mesh_size=h,
)


# %% [markdown]
# ## Post-processing

# %% vscode={"languageId": "python"}
plotter = Plotter(
    output_dir=out_dir,
    save_extracted_data=True,
)

data_dir = Path("external_data")
external_data = Plotter.load_external_data(data_dir, benchmark_tag="HM2a")

plotter.plot_volumetric_strain_vs_angle(
    vtu_files_dict_full_HM,
    model_type="HM2a",
    ylim_range=[-7.5, 2.5],
    layout="subplots",
    external_data=external_data["strain"],
)

# %% vscode={"languageId": "python"}
plotter.plot_field_variables(vtu_files_dict_full_HM)


# %% vscode={"languageId": "python"}
plotter.material_names = ["Gneiss", "Greywacke"]
plotter.vtu_file_names = {"VPF": vtu_files_dict_full_HM}

plotter.plot_avg_width_vs_stress(
    benchmark_tag="HM2a",
    metric="width",
    methods_to_include=["VPF"],
    pee_load_values=PEE_load_values,
    external_data=external_data["average"],
    ylim_range=(0, 5.0e-5),
)

plotter.plot_avg_width_vs_stress(
    benchmark_tag="HM2a",
    metric="permeability",
    methods_to_include=["VPF"],
    pee_load_values=PEE_load_values,
    external_data=external_data["average"],
    ylim_range=(0, 2.5e-10),
)

# %% vscode={"languageId": "python"}
vpf_profiles = {}
for load_case, file_list in plotter.vtu_file_names["VPF"].items():
    vpf_profiles[load_case] = plotter.extract_vpf_width_from_list(
        file_list,
        plotter.material_names,
    )

plotter.plot_fracture_aperture_profiles(
    widthProfile=vpf_profiles,
    benchmark_tag="HM2a",
    downsample=1,
    ylim=(0, 5.0e-5),
    method_label="VPF",
    external_data=external_data["widthProfile"],
)


# %% [markdown]
# ---

# %% vscode={"languageId": "python"}
pairs_to_check = {
    "M1_VPF_A_Greywacke": "M1",
    "M2a_VPF_A_Greywacke": "M2a",
    "M2b_VPF_A_Greywacke": "M2b",
}

for case, label in pairs_to_check.items():
    print(f"\n===== {label} case =====")
    new_result = np.load(Path(out_dir, f"extracted_{case}_volStrain.npz"))
    expected_result = np.load(Path("expected", f"extracted_{case}_volStrain.npz"))

    eps_v_new = new_result["eps_v"]
    eps_v_expected = expected_result["eps_v"]
    phi_new = new_result["phi"]
    phi_expected = expected_result["phi"]

    np.testing.assert_allclose(eps_v_new, eps_v_expected, atol=5e-4)
    np.testing.assert_allclose(phi_new, phi_expected, atol=1e-8)
    print(f"\n{label} case passed.")


# %% [markdown]
# ---
# # Reference
# 1. McDermott, C.I., Fraser-Harris, A., Sauter, M., Couples, G.D., Edlmann, K., Kolditz, O., Lightbody, A., Somerville, J. and Wang, W., 2018. New experimental equipment recreating geo-reservoir conditions in large, fractured, porous samples to investigate coupled thermal, hydraulic and polyaxial stress processes. *Scientific reports*, 8(1), p.14549.
#
# 2. Mollaali, M., Kolditz, O., Hu, M., Park, C.H., Park, J.W., McDermott, C.I., Chittenden, N., Bond, A., Yoon, J.S., Zhou, J. and Pan, P.Z., Liu H., Hou W.,  Lei H., Zhang L., Nagel T., Barsch M., Wang W., Nguyen S., Kwon S. and Yoshioka K., 2023. Comparative verification of hydro-mechanical fracture behavior: Task G of international research project DECOVALEX–2023. *International Journal of Rock Mechanics and Mining Sciences*, 170, p.105530.
#
# 3. Mollaali, M., Wang, W., You, T., Nagel, T., Fraser-Harris, A., McDermott, C., and Kolditz, O., 2025. Numerical benchmarking of GREAT cell experiments: Poly-axial stress effects on fluid flow in fractured rock using smeared and discrete methods.
#
