# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [raw] magic_args="[raw]"
# +++
# title = "A 2D GREAT cell benchmark suite simulated using hydro-mechanical and LIE hydro-mechanical processes"
# date = "2025-07-16"
# author = "Mostafa Mollaali, Wenqing Wang"
# image = "figures/hm_lie_bbar_stress_trace.png"
# web_subsection = "hydro-mechanics"
# weight = 3
# projects = [ "decovalex" ]
# models = [ "lab" ]
# +++

# %% vscode={"languageId": "python"}
import os
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot

mechanics_path = Path("..", "..", "Mechanics", "GreatCelljupyterNotebook").resolve()
sys.path.insert(0, str(mechanics_path))

from mesh_generator import (  # noqa: E402
    mesh_GreatCell_embeddedFracture_meshes,
    mesh_GreatCell_fullFracture_meshes,
    mesh_GreatCell_intact_meshes,
)
from ogs_model_runner import SingleOGSModel  # noqa: E402
from Plotter import Plotter  # noqa: E402


# %% vscode={"languageId": "python"}
def truncated_cmap(name, minval=0.3, maxval=0.9, n=6):
    base = mpl.colormaps.get_cmap(name)
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
out_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Great cell
#
# The GREAT cell is a  poly-axial rock-testing device that reproduces subsurface conditions down to 3.5 km depth on 200 mm-diameter samples. It imposes a rotating stress field, injects fluid through a central borehole, and records both fiber-optic strain and pore-pressure data—providing a rich dataset for validating coupled hydro-mechanical models. For full details, see the GREAT cell benchmark docs:
# [www.opengeosys.org/docs/benchmarks/small-deformations/greatcellm/](https://www.opengeosys.org/docs/benchmarks/small-deformations/greatcellm/)
#
#
# Here, in our 2D, small-deformation benchmark suite, we employ the LIE (lower-dimensional interface element) hydro-mechanical process to model fractures as interfaces. The cases are:
#
# - **$\texttt{HM}_1$**: Intact rock sample with fluid injection.
# - **$\texttt{HM}_{2a}$**: Fully fractured sample with inflow-outflow.
# - **$\texttt{HM}_{2b}$**: Half-fractured sample with inflow-outflow.


# %% [markdown]
# ---
# ## Material Properties
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

    for i, (x, y, value) in enumerate(
        zip(circle_x, circle_y, scaled_values, strict=True)
    ):
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
# # Boundary Conditions
#
# The boundary conditions applied in the simulation include both Dirichlet and Neumann conditions.
#
# - **Dirichlet conditions** (displacement control):
# \begin{equation*}
# \begin{cases}
# u_x(0, y) = 0, \quad u_y(0, y) = 0 & \quad \text{for } y = -0.09894 \text{ m}, \\
# u_y(x, 0) = 0,  & \quad \text{for} x = -0.09894 \text{ m}.
# \end{cases}
# \end{equation*}
#
# - **Neumann conditions**:
#
#   Normal stress $\sigma_{rr}$ is applied on PEEs and DSSs. The DSS stress is calculated as:
#   $$\begin{equation*}\sigma_\text{DSS}^i = \frac{\sigma_\text{PEE}^i + \sigma_\text{PEE}^{i+1}}{2}\end{equation*}
# $$
# ---

# %% [markdown]
# # Intact rock sample with fluid injection ($\texttt{HM}_1$)
#
# This benchmark, which does not consider any fracture, is designed to verify the basic computational setting for hydro-mechanical simulations. In addition to the mechanical loads, a zero constant pore pressure is prescribed at the outer boundary. {At the center of sample, fluid is injected at a rate of $Q_0^{\text{v}} = 2.085 \times 10^{-9}$~m$^3$/s,}$
#
# The hydro-mechanical simulations follows a two-stage process: a 3000~s _equilibrium phase_ under mechanical loading to stabilize initial conditions, followed by a 500~s _fluid injection phase_ to model fluid flow. This loading condition is applied to both **Greywacke** and **Gneiss** samples
#

# %% [markdown]
# ## Mesh Generation

# %% [markdown]
# ### Input

# %% vscode={"languageId": "python"}
h = 0.005
meshname = "GreatCell"

mesh_id = "mesh_intact"
meshes_intact = mesh_GreatCell_intact_meshes(
    lc=2 * h,
    lc2=h,
    r0=0.097,
    r1=0.094,
    r2=0.090,
    r3=0.065,
    out_dir=(out_dir / "gmsh" / mesh_id).resolve(),
    meshname=meshname,
)
meshes_intact.save(out_dir / f"Meshes/{mesh_id}", overwrite=True)

# %% [markdown]
# ## Run the simulation

# %% [markdown]
# ### Inputs

# %% vscode={"languageId": "python"}
# Times for load curves
times = "0.0  1000. 3500"
simulation_end_time = 3500.0
n_fracture_p_ncs = 0
model_type = "HM1"
output_prefix = "HM1_HM"
# Load
PEE_load_values = {
    "A": [10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6, 7.80e6, 9.95e6],
    "B": [7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6, 1.0e6, 3.82e6],
    "C": [1.0e6, 3.82e6, 7.80e6, 9.95e6, 10.0e6, 6.64e6, 4.46e6, 1.17e6],
}

for p in materials.values():
    p["fluid"][
        "injectionFlowRate_Inlet"
    ] = 2.085e-6  # Injection mass rate kg/s (m_dot = Q_v * rho_f)

prj_file = Path("HM1_HM.prj")
prj = ot.Project(input_file=prj_file, output_file=Path(out_dir, f"{output_prefix}.prj"))

# %% [markdown]
# ### Run OGS

# %% vscode={"languageId": "python"}
sing_ogs_model = SingleOGSModel(
    project=prj,
    meshes=meshes_intact.active_target.resolve(),
    out_dir=out_dir,
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
    mesh_path=meshes_intact.active_target.resolve(),
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
# ## Post-processing

# %% [markdown]
# ### Volumetric strain vs angle at probe circle

# %% vscode={"languageId": "python"}
plotter = Plotter(
    output_dir=out_dir,
    save_extracted_data=True,
)

data_dir = Path("external_data")
external_data = Plotter.load_external_data(data_dir, benchmark_tag="HM1")

plotter.plot_volumetric_strain_vs_angle(
    vtu_files_dict,
    model_type="HM1",
    ylim_range=[-7.5, 2.5],
    layout="subplots",
    external_data=external_data["strain"],
    show_error=True,
)

# %% [markdown]
# ### Profiles

# %% [markdown]
#
# **Figure** below shows the results of this benchmark:
#
# - The **trace of effective stress**, defined as $\sum_{i=x,y,z} \sigma'_{ii}$,
# - and the corresponding **volumetric strain** plotted versus loading angle for both Greywacke and Gneiss,
#
# These results are shown for the intact rock samples under $\texttt{HM}_1$ loading.

# %% vscode={"languageId": "python"}
plotter.plot_field_variables(vtu_files_dict)

# %% [markdown]
# ---
# # Half fractured sample with fluid injection ($\texttt{HM}_{2b}$)
#
# We use all inputs and BCs from $\texttt{M}_{2b}$ for a single wing fracture $\Gamma = [0,\,0.04]\times\{0\}.$
# A Dirichlet BC of $p(0.04,0)=3.45\text{ MPa}$
# is applied at the fracture's right tip. After a 3000 s equilibration, fluid is injected at the left tip with
# $ Q_0^{\mathrm{v}}(0,0)=4.167\times10^{-7}\,\mathrm{m}^3/\mathrm{s}\;(25\text{ ml/min}).$

# %% [markdown]
# ## Mesh Generation

# %% [markdown]
# ### Input

# %% vscode={"languageId": "python"}
h = 0.005
meshname = "GreatCell"

mesh_id = "embedded_fracture"
meshes_embedded = mesh_GreatCell_embeddedFracture_meshes(
    lc=2 * h,
    lc2=h,
    r0=0.097,
    r1=0.094,
    r2=0.090,
    r3=0.065,
    out_dir=(out_dir / "gmsh" / mesh_id).resolve(),
    meshname=meshname,
)
meshes_embedded.save(out_dir / f"Meshes/{mesh_id}", overwrite=True)

# %% [markdown]
# ## Run the simulation
# ### Input

# %% vscode={"languageId": "python"}
# Times for load curves
times = "0.0  1000. 3500"
simulation_end_time = 3500.0
n_fracture_p_ncs = 3
model_type = "HM2b"
output_prefix = "HM2b_LIE"

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
    ] = 4.167e-7  # Injection flow rate m³/s (125 ml/min)

# Project file
prj_file = Path("HM2b_LIE.prj")
prj = ot.Project(input_file=prj_file, output_file=Path(out_dir, f"{output_prefix}.prj"))

# %% [markdown]
# ### Run OGS

# %% vscode={"languageId": "python"}
sing_ogs_model = SingleOGSModel(
    project=prj,
    meshes=meshes_embedded.active_target.resolve(),
    out_dir=out_dir,
    output_prefix=output_prefix,
    method="LIE",
    n_fracture_p_ncs=n_fracture_p_ncs,
    model_type=model_type,
    materials=materials,
)

# Run simulations
vtu_files_dict_embedded = sing_ogs_model.run_simulations_with_fracture(
    times=times,
    base_project_file=prj_file,
    mesh_path=meshes_embedded.active_target.resolve(),
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
# ## Post-processing

# %% vscode={"languageId": "python"}
data_dir = Path("external_data")
external_data = Plotter.load_external_data(data_dir, benchmark_tag="HM2b")
custom_cb = {
    "u": {"vmin": 0, "vmax": 0.25},
    "stress": {"vmin": -20, "vmax": 10},
    "strain": {"vmin": -0.05, "vmax": 0},
    "pressure": {"vmin": 0.1, "vmax": 8},
}
plotter = Plotter(
    output_dir=out_dir,
    colorbar_opts=custom_cb,
    save_extracted_data=True,
)

plotter.plot_volumetric_strain_vs_angle(
    vtu_files_dict_embedded,
    model_type="HM2b",
    ylim_range=[-7.5, 2.5],
    layout="subplots",  # "single" or "subplots"
    external_data=external_data["strain"],
)

# %% vscode={"languageId": "python"}
plotter.plot_field_variables(vtu_files_dict_embedded)

# %% vscode={"languageId": "python"}
plotter.material_names = ["Gneiss", "Greywacke"]
plotter.vtu_file_names = {"LIE": vtu_files_dict_embedded}

plotter.plot_avg_width_vs_stress(
    benchmark_tag="HM2b",
    metric="width",
    methods_to_include=["LIE"],
    pee_load_values=PEE_load_values,
    external_data=external_data["average"],
    w0=next(iter(materials.values()))["w_init"],
)

plotter.plot_avg_width_vs_stress(
    benchmark_tag="HM2b",
    metric="permeability",
    methods_to_include=["LIE"],
    pee_load_values=PEE_load_values,
    external_data=external_data["average"],
    w0=next(iter(materials.values()))["w_init"],
)
# %% vscode={"languageId": "python"}
lie_profiles = {}
for load_case, file_list in plotter.vtu_file_names["LIE"].items():
    lie_profiles[load_case] = plotter.extract_lie_aperture_from_list(
        file_list,
        plotter.material_names,
    )

plotter.plot_fracture_aperture_profiles(
    widthProfile=lie_profiles,
    benchmark_tag="HM2b",
    downsample=1,
    ylim=(0, 10e-6),
    method_label="LIE",
    external_data=external_data["widthProfile"],
)
# %% [markdown]
# ---
# # Fully fractured sample with fluid injection ($\texttt{HM}_{2a}$)
#
# We performed five hydro-mechanical (HMca) simulations with different angles between the second principal stress direction (σ₂ = 10 MPa) and the PEE 5 fracture plane
# $\Gamma = [-0.094,\,0.094]\times\{0\}.$
# The initial pore pressure was uniform at $p_0 = 0.1$ MPa, with Dirichlet conditions
# $p(-0.094,0)=p(0.094,0)=3.45\text{ MPa}.$
# After equilibrating for 3000 s, we injected fluid at the center at a rate
# $
# Q_0^{\text{v}}(0,0)=4.167\times10^{-7}\,\mathrm{m}^3/\mathrm{s}\;(25\,\mathrm{ml/min})$
# for 500 s, starting from an initial fracture aperture $w_{\text{ini}}=1\times10^{-6}$ m.

# %% [markdown]
# ## Mesh Generation

# %% [markdown]
# ### Input

# %% vscode={"languageId": "python"}
h = 0.005
meshname = "GreatCell"

mesh_id = "fullFracture"
meshes_full = mesh_GreatCell_fullFracture_meshes(
    lc=2 * h,
    lc2=h,
    r0=0.097,
    r1=0.094,
    r2=0.090,
    r3=0.065,
    out_dir=(out_dir / "gmsh" / mesh_id).resolve(),
    meshname=meshname,
)
meshes_full.save(out_dir / f"Meshes/{mesh_id}", overwrite=True)

# %% [markdown]
# ## Run the simulation
# ### Input

# %% vscode={"languageId": "python"}
# Times for load curves
times = "0.0  1000. 3500"
simulation_end_time = 3500.0
n_fracture_p_ncs = 3
model_type = "HM2a"
output_prefix = "HM2a_LIE"
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
    ] = 4.167e-7  # Injection flow rate m³/s (125 ml/min)

# Project file
prj_file = Path("HM2a_LIE.prj")
prj = ot.Project(input_file=prj_file, output_file=Path(out_dir, f"{output_prefix}.prj"))

# %% [markdown]
# ### Run OGS

# %% vscode={"languageId": "python"}

sing_ogs_model = SingleOGSModel(
    project=prj,
    meshes=meshes_full.active_target.resolve(),
    out_dir=out_dir,
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
    mesh_path=meshes_full.active_target.resolve(),
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
# ## Post-processing

# %% vscode={"languageId": "python"}
plotter = Plotter(
    output_dir=out_dir,
    save_extracted_data=True,
)
data_dir = Path("external_data")
external_data = Plotter.load_external_data(data_dir, benchmark_tag="HM2a")

plotter.plot_volumetric_strain_vs_angle(
    vtu_files_dict_full,
    model_type="HM2a",
    ylim_range=[-7.5, 2.5],
    layout="subplots",
    external_data=external_data["strain"],
)

# %% vscode={"languageId": "python"}
plotter.plot_field_variables(vtu_files_dict_full)

# %% vscode={"languageId": "python"}
plotter.material_names = ["Gneiss", "Greywacke"]
plotter.vtu_file_names = {"LIE": vtu_files_dict_full}

plotter.plot_avg_width_vs_stress(
    benchmark_tag="HM2a",
    metric="width",
    methods_to_include=["LIE"],
    pee_load_values=PEE_load_values,
    external_data=external_data["average"],
    w0=next(iter(materials.values()))["w_init"],
)

plotter.plot_avg_width_vs_stress(
    benchmark_tag="HM2a",
    metric="permeability",
    methods_to_include=["LIE"],
    pee_load_values=PEE_load_values,
    external_data=external_data["average"],
    w0=next(iter(materials.values()))["w_init"],
)

# %% vscode={"languageId": "python"}
lie_profiles = {}
for load_case, file_list in plotter.vtu_file_names["LIE"].items():
    lie_profiles[load_case] = plotter.extract_lie_aperture_from_list(
        file_list,
        plotter.material_names,
    )

plotter.plot_fracture_aperture_profiles(
    widthProfile=lie_profiles,
    benchmark_tag="HM2a",
    downsample=1,
    # ylim=(0, 1e-5),
    method_label="LIE",
    external_data=external_data["widthProfile"],
)


# %% [markdown]
# ---
#

# %% vscode={"languageId": "python"}
pairs_to_check = {
    "HM1_HM_A_Greywacke": "HM1",
    "HM2a_LIE_A_Greywacke": "HM2a",
    "HM2b_LIE_A_Greywacke": "HM2b",
}

for case, label in pairs_to_check.items():
    print(f"\n===== {label} case =====")
    new_result = np.load(Path(out_dir, f"extracted_{case}_volStrain.npz"))
    expected_result = np.load(Path("expected", f"extracted_{case}_volStrain.npz"))

    eps_v_new = new_result["eps_v"]
    eps_v_expected = expected_result["eps_v"]
    phi_new = new_result["phi"]
    phi_expected = expected_result["phi"]

    np.testing.assert_allclose(eps_v_new, eps_v_expected, atol=7e-4)
    np.testing.assert_allclose(phi_new, phi_expected, atol=1e-6)
    print(f"\n{label} case passed.")

# %% [markdown]
# ---
#
# ### Reference
# 1. McDermott, C.I., Fraser-Harris, A., Sauter, M., Couples, G.D., Edlmann, K., Kolditz, O., Lightbody, A., Somerville, J. and Wang, W., 2018. New experimental equipment recreating geo-reservoir conditions in large, fractured, porous samples to investigate coupled thermal, hydraulic and polyaxial stress processes. *Scientific reports*, 8(1), p.14549.
#
# 2. Mollaali, M., Kolditz, O., Hu, M., Park, C.H., Park, J.W., McDermott, C.I., Chittenden, N., Bond, A., Yoon, J.S., Zhou, J. and Pan, P.Z., Liu H., Hou W.,  Lei H., Zhang L., Nagel T., Barsch M., Wang W., Nguyen S., Kwon S. and Yoshioka K., 2023. Comparative verification of hydro-mechanical fracture behavior: Task G of international research project DECOVALEX–2023. *International Journal of Rock Mechanics and Mining Sciences*, 170, p.105530.
#
# 3. Mollaali, M., Wang, W., You, T., Nagel, T., Fraser-Harris, A., McDermott, C., and Kolditz, O., 2025. Numerical benchmarking of GREAT cell experiments: Poly-axial stress effects on fluid flow in fractured rock using smeared and discrete methods.
#
