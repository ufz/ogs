# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [raw] magic_args="[raw]"
# +++
# title = "A 2D GREAT cell benchmark suite simulated using small deformation and LIE mechanical processes"
# date = "2025-04-03"
# author = "Mostafa Mollaali, Wenqing Wang"
# image = "figures/hm_lie_bbar_stress_trace.png"
# web_subsection = "small deformation"
# weight = 3
# +++

# %%
import os
from pathlib import Path
from subprocess import run

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pyvista as pv
from external_data_dict import external_data
from matplotlib import colormaps

# Local modules
from mesh_generator import (
    mesh_GreatCell_embeddedFracture,
    mesh_GreatCell_fullFracture,
    mesh_GreatCell_intact,
)
from ogs_model_runner import SingleOGSModel
from Plotter import Plotter


# %%
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

# %%
ot.plot.setup.show_region_bounds = False

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
if not out_dir.exists():
    out_dir.mkdir(parents=True)

# %% [markdown]
#
#
# ## Great cell
#
# <div style="float: left; width: 50%; margin: 0 1em 1em 0;">
#     <img src="figs/great-cell.png" alt="GREAT cell facility" style="width: 100%;" />
#     <p style="text-align: center; font-size: 90%;"> <strong>Figure:</strong>  The GREAT cell facility.</p>
# </div>
#
# The GREAT cell is designed to reproduce subsurface conditions in the laboratory down to a depth of **3.5 km**, simulating both geomechanical stress and geothermal reservoir conditions on **200 mm diameter** rock samples that contain fracture networks. This enables validation of numerical predictions under realistic conditions.
#
# The GREAT cell represents a major advancement in experimental geomechanics technology. It uniquely:
#
# - Generates a truly **poly-axial rotating stress field**
# - Facilitates **fluid flow** through the rock samples
# - Employs **state-of-the-art fibre-optic strain sensing** for high-resolution circumferential strain measurements
# - Records **fluid pressure and strain** at thousands of data points per hour
#
# This facility enhances understanding of fracture propagation and its effect on fluid flow—key factors in **geo-energy applications** involving fluid storage or extraction in the subsurface.
#
# ---
#
# ### Figure: GREAT Cell Benchmarking Concept
#
#
# <div style="display: flex; justify-content: space-between; align-items: center;">
#   <div style="width: 49%; text-align: center;">
#     <img src="figs/step2b1.png" alt="Typical rock sample" style="width: 60%; height: auto;" />
#     <p>(a) Typical rock sample</p>
#   </div>
#   <div style="width: 49%; text-align: center;">
#     <img src="figs/schematic_3d_crossSection.jpg" alt="Schematic cross-section" style="width: 50%; height: auto;" />
#     <p>(b) Simplification for GREAT cell benchmarking</p>
#   </div>
# </div>
#
# <p><strong>Figure:</strong> (a) Typical rock sample and (b) its simplified 3D schematic used for GREAT cell benchmarks.
# <em>The cross-section (b) aligns with the 2D simulations discussed earlier.</em></p>
#
#

# %% [markdown]
# ## Benchmarking Strategy
#
# A comprehensive benchmarking exercise is performed to establish the models' capabilities and required settings, such as the required discretization for comparable accuracy. The exercise benchmarks are simplified versions of the GREAT cell experiments, focusing on general features like fracture patterns, rotating boundary conditions, and coupled processes.
#
# These benchmark exercises are conducted in 2D (plane strain), representing a horizontal section through the middle of the rock samples tested in the GREAT cell (see Figure 1). This 2D plane-strain approach allows the main characteristics of hydro-mechanical (HM) fracture mechanics to be studied efficiently. The benchmarking concept involves a stepwise increase in complexity, including both mechanical (M) and hydro-mechanical (HM) approaches.
#
#
#
# We distinguish two main types of benchmarks: **intact** and **fractured rock tests** (see Figure 2). In the first type, an intact rock sample is subjected to rotating external stress conditions ($\texttt{M}_1$) and fluid is injected from a central borehole ($\texttt{HM}_1$):
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
#
# Finally, fracture propagation conditions are studied by injecting fluid into a pre-existing fracture ($\texttt{HM}_{3c}$) and by initiating a fracture from a borehole ($\texttt{HM}_{3d}$):
#
# - **$\texttt{HM}_{3c}$**: Propagation of pre-existing fracture by fluid injection.
# - **$\texttt{HM}_{3d}$**: Fracture nucleation from borehole by fluid injection.
#
#
#
# <div style="display: flex; justify-content: space-between; align-items: flex-start;">
#   <div style="width: 70%; text-align: center;">
#     <img src="figs/Schematic_all3.png" alt="GREAT benchmarks definition" style="width: 100%; height: auto;" />
#     <p><strong>Figure 1:</strong><br><em>Three sample types of GREAT cell benchmarks: intact, fractured (static), and propagating. Rotating boundary conditions are applied.</em></p>
#   </div>
#   <div style="width: 29%; text-align: center;">
#     <img src="figs/schematic_pee_dss.png" alt="Schematic loading 2D UFZ" style="width: 100%; height: auto;" />
#     <p><strong>Figure 2:</strong><br><em>Geometry and location of PEEs and DSSs. Fracture aligned with PEE1/1a direction.</em></p>
#   </div>
# </div>
#
#
#

# %% [markdown]
# ## Material Properties

# %% [markdown]
#
# The material properties are provided in Table 1. The computational model incorporates two distinct elastic materials within its domain: a central circle ($r=0.097$ m) of rock surrounded by a rubber sheath in a 2D configuration.
#
# ### Table 1: Material Properties
#
# | Parameter | Unit | Rubber Sheath | Greywacke | Gneiss (Freiberg) |
# |----------|------|--------------------|-----------|-------------------|
# | Young's modulus, $E$ | GPa | 0.1 | 26.85 | 83.9 |
# | Poisson's ratio, $\nu$ | - | 0.4 | 0.27 | 0.21 |
# | Tensile strength | MPa | - | 17.01–16.67 | 16.8 |
# | Permeability, $K_m$ | m² | $1 \times 10^{-19}$ | $2.58 \times 10^{-19}$ | $1 \times 10^{-19}$ |
# | Solid density, $\rho_s$ | kg/m³ | 1200 | 2650 | 2750 |
# | Porosity, $\phi$ | - | - | 0.005 | 0.001 |
# | Biot coefficient, $\alpha_m$ | - | 0.0 | 0.8 | 0.6 |
# | Fracture energy, $G_c$ (VPF-FEM) | N/m | 0.1 | 30 | 50 |
# | Normal stiffness, $k_{nn}$ (LIE-FEM) | GPa/m | - | 100.0 | 200.0 |
# | Tangential stiffness, $k_{tt}$ (LIE-FEM) | GPa/m | - | 50.0 | 100.0 |
#
# ---
# <!--
# ### Table 2: Fluid Properties
#
# | Property | Value | Unit |
# |----------|-------|------|
# | Fluid density, $\rho_f$ | $1000$ | kg/m³ |
# | Viscosity, $\mu$ | $1 \times 10^{-3}$ | Pa·s |
# | Fluid compressibility, $\kappa_f^{-1}$ | $4.4 \times 10^{-10}$ | Pa⁻¹ | -->
#
# <!-- --- -->

# %%
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
# # Loads

# %% [markdown]
#
#
# ### Table 3: Load Conditions
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

# %%
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

# %%
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

    formatted_values = [rf"${v/1e6:.1f}$" for v in values]
    scaled_values = [v / 2 for v in values]

    circle = plt.Circle(
        (0, 0), circle_radius, color="black", fill=False, linestyle="--", linewidth=2
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
        0, 0, f"Load {label}", fontsize=32, ha="center", va="center", family="serif"
    )

plt.tight_layout()
output_path = Path(out_dir, "loads_A_B_C_schematic_with_angles.png")
plt.savefig(output_path, dpi=350, bbox_inches="tight", transparent=True)
plt.show()

# %% [markdown]
# ---
#
# ### Boundary Conditions
#
# The boundary conditions applied in the simulation include both Dirichlet and Neumann conditions.
#
# - **Dirichlet conditions** (displacement control):
#   $$
# \begin{equation*}
# \begin{cases}
# u_x(0, y) = 0, \quad u_y(0, y) = 0 & \quad \text{for } y = -0.09894 \text{ m}, \\
# u_y(x, 0) = 0,  & \quad \text{for } x = -0.09894 \text{ m}.
# \end{cases}
# \end{equation*}
#  $$
#
# - **Neumann conditions**:
#   Normal stress $\sigma_{rr}$ is applied on PEEs and DSSs. The DSS stress is calculated as:
#   $$\sigma_\text{DSS}^i = \frac{\sigma_\text{PEE}^i + \sigma_\text{PEE}^{i+1}}{2}$$

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
# (*see also loading Table referenced in earlier section*).

# %% [markdown]
# # Mesh Generation

# %% [markdown]
# ### Input

# %%
h = 0.0025
meshname = "GreatCell"
mesh_path = Path(out_dir, "mesh_intact").resolve()
print(mesh_path)

# %% [markdown]
# ### Boundary meshes

# %% [markdown]
# ### Gmsh (boundary meshes)

# %%
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
# ### Convert .msh to an OGS-compatible mesh

# %%
msh_path = Path(mesh_path, f"{meshname}.msh")
meshes = ot.meshes_from_gmsh(filename=msh_path, dim=[1], reindex=True, log=False)

for name, mesh in meshes.items():
    print(f"{name}: {mesh.n_cells} cells")
    pv.save_meshio(Path(mesh_path, f"{name}.vtu"), mesh)

# %% [markdown]
# ### Visualization of boundary meshes

# %%
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
# ### Gmsh (computational domain mesh)

# %%
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
# ### Convert .msh to an OGS-compatible mesh

# %%
msh_path = Path(mesh_path, f"{meshname}.msh")
meshes_volume = ot.meshes_from_gmsh(
    filename=msh_path, dim=[1, 2], reindex=True, log=False
)

for name, mesh in meshes_volume.items():
    print(f"{name}: {mesh.n_cells} cells")
    pv.save_meshio(Path(mesh_path, f"{name}.vtu"), mesh)

# %%

# %% [markdown]
# ### Visualization of computational domain mesh

# %%
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

# %%
# %cd {mesh_path}
run(
    "identifySubdomains -f -m domain.vtu -- physical_group_*.vtu",
    shell=True,
    check=True,
)

# %cd -

# %% [markdown]
# # Run the simulation

# %% [markdown]
# ## Inputs

# %%
# Times for load curves
times = "0.0  1000. 3500"
simulation_end_time = 3500.0
n_fracture_p_ncs = 0
model_type = "M1"
output_prefix = "M1_SD"
# Load
PEE_load_values = {
    "A": [10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6, 7.80e6, 9.95e6],
    "B": [7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6],
    "C": [1.0e6, 3.82e6, 7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6],
}

prj_file = Path("M1_SD.prj")
# Now create a Project object
prj = ot.Project(input_file=prj_file, output_file=Path(out_dir, f"{output_prefix}.prj"))

# %% [markdown]
# ## Run OGS

# %%
# Now create SingleOGSModel
sing_ogs_model = SingleOGSModel(
    model=prj,
    out_dir=out_dir,
    mesh_path=mesh_path,
    output_prefix=output_prefix,
    method="LIE",
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
    use_b_bar_value="true",
    tension_cutoff="true",
)

# %% [markdown]
# # Post-processing

# %% [markdown]
# ### Volumetric strain vs angle at probe circle

# %%
json_path = Path("./external_data_dict.py").resolve()
print(f"[DEBUG] Trying path: {json_path}")

if json_path.exists():
    from external_data_dict import external_data
else:
    print("[WARNING] External data dict not found! Skipping...")
    external_data = None

plotter = Plotter(
    output_dir=out_dir,
    markers=["o", "^"],
    material_cmaps={
        "Gneiss": truncated_cmap("Blues"),
        "Greywacke": truncated_cmap("Oranges"),
    },
)

plotter.plot_volumetric_strain_vs_angle(
    vtu_files_dict,
    model_type="M1",
    # ylim_range=[-0.035, -0.005],
    layout="subplots",
    external_data=external_data,
)

# %% [markdown]
# ## Profiles

# %% [markdown]
#
# **Figure** below shows the results of this benchmark:
#
# - The **trace of effective stress**, defined as $\sum_{i=x,y,z} \sigma'_{ii}$,
# - and the corresponding **volumetric strain** plotted versus loading angle,
# - for both Greywacke and Gneiss,
#
# These results are shown for the intact rock samples under $\texttt{M}_1$ loading.

# %%
plotter.plot_field_variables(vtu_files_dict)

# %% [markdown]
# **Note**: The stress distribution and **volumetric strain** are shown only in the central zone defined by:
#
# $$
# \sqrt{x^2 + y^2} < 0.065\ \text{m}
# $$
#
# This is done because **strain is relatively high in the PEE regions** (due to their softer material properties), and we want to isolate the deformation behavior within the **core of the sample**.
#
# <!-- ---
#
# ### Figure 4: Effective Stress Trace and Volumetric Strain (VPF-FEM)
#
# | (a) Greywacke Sample | (b) Gneiss Sample |
# |----------------------|------------------|
# | ![Greywacke Stress-Strain](figs/M1_VPF_profiles_greywacke.png) | ![Gneiss Stress-Strain](figs/M1_VPF_profiles_gneiss.png) |
#
# **Figure 4**: Trace of effective stress and volumetric strain vs angle for intact Greywacke (a) and Gneiss (b) using VPF-FEM method. Stress is visualized over the full sample domain, while strain is restricted to the central region to avoid high values near PEE materials. -->

# %%
# Plot only inner mesh (within r=0.065 m):
plotter.plot_field_variables(vtu_files_dict, inner=True, r=0.065)

# %% [markdown]
# ---
# # Half fractured sample ($\texttt{M}_{2b}$)
#
# In this section, all input data and boundary conditions provided in $\texttt{M}_{2b}$, except  a single wing fracture, defined as $\Gamma =  \left[0, 0.04\right] \times \{0\}$, is considered.

# %% [markdown]
# # Mesh Generation

# %% [markdown]
# ### Input

# %%
h = 0.0025
meshname = "GreatCell"
mesh_path_embedded = Path(out_dir, "mesh_GreatCell_embeddedFracture").resolve()

# %% [markdown]
# ### Boundary meshes

# %% [markdown]
# ### Gmsh (boundary meshes)

# %%
msh_file_embedded = mesh_GreatCell_embeddedFracture(
    lc=0.005,
    lc2=h,
    r0=0.097,
    r1=0.094,
    r2=0.090,
    r3=0.065,
    out_dir=mesh_path_embedded,
    meshname=meshname,
    mode="BC",
)

# %% [markdown]
# ### Convert .msh to an OGS-compatible mesh

# %%
msh_path_embedded = Path(mesh_path_embedded, f"{meshname}.msh")
meshes_embedded = ot.meshes_from_gmsh(
    filename=msh_path_embedded, dim=[1], reindex=True, log=False
)

for name, mesh in meshes_embedded.items():
    print(f"{name}: {mesh.n_cells} cells")
    pv.save_meshio(Path(mesh_path_embedded, f"{name}.vtu"), mesh)

# %% [markdown]
# ### Visualization of boundary meshes

# %%
plotter = pv.Plotter()
for name, mesh in meshes_embedded.items():
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
# ### Gmsh (computational domain mesh)

# %%
msh_file_volume_embedded = mesh_GreatCell_embeddedFracture(
    lc=0.005,
    lc2=h,
    r0=0.097,
    r1=0.094,
    r2=0.090,
    r3=0.065,
    out_dir=mesh_path_embedded,
    meshname=meshname,
    mode="domain",
)

# %% [markdown]
# ### Convert .msh to an OGS-compatible mesh

# %%
msh_path_embedded = Path(mesh_path_embedded, f"{meshname}.msh")
meshes_volume_embedded = ot.meshes_from_gmsh(
    filename=msh_path_embedded, dim=[1, 2], reindex=True, log=False
)

for name, mesh in meshes_volume_embedded.items():
    print(f"{name}: {mesh.n_cells} cells")
    pv.save_meshio(Path(mesh_path_embedded, f"{name}.vtu"), mesh)

# %% [markdown]
# ### Visualization of computational domain mesh

# %%
plotter = pv.Plotter()
for name, mesh in meshes_volume_embedded.items():
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

# %%
# %cd {mesh_path_embedded}
# !pwd
run(
    "identifySubdomains -f -m domain.vtu -- physical_group_*.vtu",
    shell=True,
    check=True,
)

# %cd -

# %%
# Times for load curves
times = "0.0  1000. 3500"
simulation_end_time = 3500.0
n_fracture_p_ncs = 3
model_type = "M2b"
output_prefix = "M2b_LIE"

# Load
PEE_load_values = {
    "A": [10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6, 7.80e6, 9.95e6],
    "B": [7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6],
    "C": [1.0e6, 3.82e6, 7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6],
}

# Project file
prj_file = Path("M2_LIE.prj")


# Now create a Project object
prj = ot.Project(input_file=prj_file, output_file=Path(out_dir, f"{output_prefix}.prj"))

sing_ogs_model = SingleOGSModel(
    model=prj,
    out_dir=out_dir,
    mesh_path=mesh_path_embedded,
    output_prefix=output_prefix,
    method="LIE",
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    materials=materials,
)

vtu_files_dict_embedded = sing_ogs_model.run_simulations_with_fracture(
    times=times,
    base_project_file=prj_file,
    mesh_path=mesh_path_embedded,
    load_cases=PEE_load_values,
    material_names=material_names,
    materials=materials,
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    output_prefix=output_prefix,
    out_dir=out_dir,
    use_b_bar_value="true",
    tension_cutoff="true",
)

# %%
plotter = Plotter(
    output_dir=out_dir,
    markers=["o", "s", "D", "^", "v", "P", "*", "X"],
    material_cmaps={
        "Gneiss": truncated_cmap("Blues"),
        "Greywacke": truncated_cmap("Oranges"),
    },
)

plotter.plot_volumetric_strain_vs_angle(
    vtu_files_dict_embedded, model_type="M2b", ylim_range=[-0.1, 0.005], layout="single"
)

# %%
plotter.plot_field_variables(vtu_files_dict_embedded)

# %% [markdown]
# ---
# # Fully fractured sample ($\texttt{M}_{2a}$)
#
# In this benchmark,  2D plane strain numerical simulations are performed  to establish a baseline model for assessing the impact of fracture orientation on strain distribution under poly-axial loading. A planar fracture within the specimen, $\Gamma =  \left[-0.094, 0.094\right] \times \{0\}$, is considered, under poly-axial loading applied on PEE's and DSS's

# %% [markdown]
# # Mesh Generation

# %% [markdown]
# ### Input

# %%
h = 0.0025
meshname = "GreatCell"
mesh_path_full = Path(".", out_dir, "mesh_GreatCell_fullFracture").resolve()

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# ### Boundary meshes

# %% [markdown]
# ### Gmsh (boundary meshes)

# %%
msh_file_full = mesh_GreatCell_fullFracture(
    lc=0.005,
    lc2=h,
    r0=0.097,
    r1=0.094,
    r2=0.090,
    r3=0.065,
    out_dir=mesh_path_full,
    meshname=meshname,
    mode="BC",
)

# %%
msh_path_full = Path(mesh_path_full, f"{meshname}.msh")
meshes_volume_full = ot.meshes_from_gmsh(
    filename=msh_path_full, dim=[1], reindex=True, log=False
)

for name, mesh in meshes_volume_full.items():
    print(f"{name}: {mesh.n_cells} cells")
    pv.save_meshio(Path(mesh_path_full, f"{name}.vtu"), mesh)

# %%
plotter = pv.Plotter()
for name, mesh in meshes_volume_full.items():
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

# %%
msh_file_full = mesh_GreatCell_fullFracture(
    lc=0.005,
    lc2=h,
    r0=0.097,
    r1=0.094,
    r2=0.090,
    r3=0.065,
    out_dir=mesh_path_full,
    meshname=meshname,
    mode="domain",
)

# %%
msh_path_full = Path(mesh_path_full, f"{meshname}.msh")
meshes_volume_full = ot.meshes_from_gmsh(
    filename=msh_path_full, dim=[1, 2], reindex=True, log=False
)

for name, mesh in meshes_volume_full.items():
    print(f"{name}: {mesh.n_cells} cells")
    pv.save_meshio(Path(mesh_path_full, f"{name}.vtu"), mesh)

# %%
plotter = pv.Plotter()
for name, mesh in meshes_volume_full.items():
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

# %%
# %cd {mesh_path_full}
run(
    "identifySubdomains -f -m domain.vtu -- physical_group_*.vtu",
    shell=True,
    check=True,
)

# %cd -

# %%
# Times for load curves
times = "0.0  1000. 3500"
simulation_end_time = 3500.0
n_fracture_p_ncs = 3
model_type = "M2a"
output_prefix = "M2a_LIE"
# Load
PEE_load_values = {
    "A": [10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6, 7.80e6, 9.95e6],
    "B": [7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6],
    "C": [1.0e6, 3.82e6, 7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6],
}

# Project file
prj_file = Path("M2_LIE.prj")

# Now create a Project object
prj = ot.Project(input_file=prj_file, output_file=Path(out_dir, f"{output_prefix}.prj"))

# Now create SingleOGSModel
sing_ogs_model = SingleOGSModel(
    model=prj,
    out_dir=out_dir,
    mesh_path=mesh_path_full,
    output_prefix=output_prefix,
    method="LIE",
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    materials=materials,
)

# Run simulations
vtu_files_dict_full = sing_ogs_model.run_simulations_with_fracture(
    times=times,
    base_project_file=prj_file,
    mesh_path=mesh_path_full,
    load_cases=PEE_load_values,
    material_names=material_names,
    materials=materials,
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    output_prefix=output_prefix,
    out_dir=out_dir,
    use_b_bar_value="true",
    tension_cutoff="true",
)


# %%
plotter = Plotter(
    output_dir=out_dir,
    markers=["o", "s", "D", "^", "v", "P", "*", "X"],
    material_cmaps={
        "Gneiss": truncated_cmap("Blues"),
        "Greywacke": truncated_cmap("Oranges"),
    },
)

plotter.plot_volumetric_strain_vs_angle(
    vtu_files_dict_full, model_type="M2a", ylim_range=[-0.1, 0.005], layout="subplots"
)

# %%
plotter.plot_field_variables(vtu_files_dict_full)

# %%
pairs_to_check = {
    "M1_SD_A_Greywacke": "M1",
    "M2a_LIE_A_Greywacke": "M2a",
    "M2b_LIE_A_Greywacke": "M2b",
}

for case, label in pairs_to_check.items():
    new_result = np.load(Path(out_dir, f"extracted_{case}.npz"))
    expected_result = np.load(Path("expected", f"expected_{case}.npz"))

    eps_v_new = new_result["eps_v"]
    eps_v_expected = expected_result["eps_v"]

    np.testing.assert_allclose(actual=eps_v_new, desired=eps_v_expected, atol=1e-10)
    print(f"{label} case passed.")


# %% [markdown]
# ### Reference
# 1. McDermott, C.I., Fraser-Harris, A., Sauter, M., Couples, G.D., Edlmann, K., Kolditz, O., Lightbody, A., Somerville, J. and Wang, W., 2018. New experimental equipment recreating geo-reservoir conditions in large, fractured, porous samples to investigate coupled thermal, hydraulic and polyaxial stress processes. *Scientific reports*, 8(1), p.14549.
#
# 2. Mollaali, M., Kolditz, O., Hu, M., Park, C.H., Park, J.W., McDermott, C.I., Chittenden, N., Bond, A., Yoon, J.S., Zhou, J. and Pan, P.Z., Liu H., Hou W.,  Lei H., Zhang L., Nagel T., Barsch M., Wang W., Nguyen S., Kwon S. and Yoshioka K., 2023. Comparative verification of hydro-mechanical fracture behavior: Task G of international research project DECOVALEX–2023. *International Journal of Rock Mechanics and Mining Sciences*, 170, p.105530.
#
# 3. Mollaali, M., Wang, W., You, T., Nagel, T., Fraser-Harris, A., McDermott, C., and Kolditz, O., 2025. Numerical benchmarking of GREAT cell experiments: Poly-axial stress effects on fluid flow in fractured rock using smeared and discrete methods.
#
