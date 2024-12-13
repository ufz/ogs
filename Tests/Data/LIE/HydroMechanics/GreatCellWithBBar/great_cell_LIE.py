# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python (.venv)
#     language: python
#     name: venv
# ---

# %% [raw]
# +++
# title = "A 2D GREAT cell experiment benchmark simulated using the LIE and the B-bar methods"
# date = "2024-12-13"
# author = "Wenqing Wang & Florian Zill"
# image = "figures/hm_lie_bbar_stress_trace.png"
# web_subsection = "hydro-mechanics"
# weight = 3
# +++

# %%
import math
import os
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import ogstools as ot
import pyvista as pv
import vtuIO

# %% [markdown]
# ## A 2D GREAT cell experiment benchmark simulated using the LIE and the B-bar methods
#
# This benchmark is based on the SAFENET task of the DECOVALEX 2027 project, which aims to better understand fracture nucleation and evolution processes in crystalline rocks with applications in nuclear waste management and geothermal reservoir engineering, and it simulates a 2D ideal model of  the GREAT Cell experiment [1][2].
#
# As shown in the figure below, paired loading conditions are applied to the paired boundary segments (denoted as PEE and DSS) of a cross section of the experiment sample, enabling a momentum balance.
#
# <img src="./figures/greatcell_loading_schematic_xy.png" alt="Schematic view of Great Cell BC." style="width:400px;">
# Schematic of 2D GREAT Cell benchmark (Courstey of Mostafa Mollaali).
#
#
# The load conditions or the pressure applied on the segment `PEE`s are given in the table below:
# | PEE 1\&1a | PEE 2\&2a | PEE 3\&3a | PEE 4\&4a | PEE 5\&5a | PEE 6\&6a | PEE 7\&7a | PEE 8\&8a |
# | :-------- | :-------- | :-------- | :-------- | :-------- | :-------- | :-------- | :-------- |
# | 7.73 MPa | 5.7 MPa | 4.39 MPa | 2.4 MPa | 2.3 MPa | 4.0 MPa | 6.4 MPa | 7.7 MPa|
#
# The load condition of each `DSS` is the average of the values of its two PEE neighbours. Rigid movement is avoided by fixing the displacement at the bottom point.
#
# An existing fracture is assumed along the $x$ axis.
#
# The material properties of the sample are given in the following table:
#
# | Parameter              |      Value      |    Unit |
# | :--------------------- | :-------------: | ------: |
# | Young's modulus        |      3.85       |     GPa |
# | Poisson ratio          |      0.4        |  -      |
# | Biot's constant        |     0.005       |  -      |
# | Intrinsic permeability |     2.58$\cdot 10^{-19}$        | m$^2$ |
# | Porosity               |     $10^{-3}$   |  -  |
# | Specific storage       |     0.0   |  Pa$^{-1}$  |
#
# For the strip portion (in the blue-colored area in the figure above), only the Young's modulus is reduced to 0.1 GPa.
#
# As for the fracture, the cubic law is applied to the intrinsic permeability. The other fracture material models are listed in the following table:
# | Parameter              |      Value      |    Unit |
# | :--------------------- | :-------------: | ------: |
# | Stiffness ($k_{nn}, k_{tt}$) |      (10, 4)   |     GPa |
# | Biot's constant        |     1.0       |  -      |
# | Intrinsic permeability |     Cubic law       | m$^2$ |
# | Specific storage       |    $10^{-11}$   |  Pa$^{-1}$ |
#
# For the HM coupled simulation, the viscosity is $10^{-3}$ Pa$\cdot$s.
#
# The initial stresses in the sample are $\sigma_{xx}$ = $\sigma_{yy}$ = $\sigma_{zz}$ =-3 MPa,
# $\sigma_{xy}$ =0 MPa, while in the fracture are $\sigma_{nn}=0$ MPa and $\sigma_{tt}=-2$ MPa.
#
# This benchmark uses the LIE and the B-bar methods.
#
# ## Reference
#
# 1. McDermott, C.I., Fraser-Harris, A., Sauter, M., Couples, G.D., Edlmann, K., Kolditz, O., Lightbody, A., Somerville, J. and Wang, W., 2018. New experimental equipment recreating geo-reservoir conditions in large, fractured, porous samples to investigate coupled thermal, hydraulic and polyaxial stress processes. *Scientific reports*, 8(1), p.14549.
#
# 2. Mollaali, M., Kolditz, O., Hu, M., Park, C.H., Park, J.W., McDermott, C.I., Chittenden, N., Bond, A., Yoon, J.S., Zhou, J. and Pan, P.Z., Liu H., Hou W.,  Lei H., Zhang L., Nagel T., Barsch M., Wang W., Nguyen S., Kwon S. and Yoshioka K., 2023. Comparative verification of hydro-mechanical fracture behavior: Task G of international research project DECOVALEX–2023. *International Journal of Rock Mechanics and Mining Sciences*, 170, p.105530.

# %% [markdown]
# ## Post-processing functions


# %%
def get_last_vtu_file_name(pvd_file_name):
    tree = ET.parse(pvd_file_name)
    root = tree.getroot()
    # Get the last DataSet tag
    last_dataset = root.findall(".//DataSet")[-1]

    # Get the 'file' attribute of the last DataSet tag
    file_attribute = last_dataset.attrib["file"]
    return Path(out_dir, file_attribute)


def get_extracted_mesh_name(vtu_file_name):
    return str(vtu_file_name).split(".")[0] + "_extracted.vtu"


# %%
def get_vol_strain_of_mesh_at_circle(vtu_mesh, radius=0.065):
    eps = vtu_mesh.point_data["epsilon"]
    eps_v = []
    phi = []

    for node_id, x in enumerate(vtu_mesh.points):
        if abs(x[0] ** 2 + x[1] ** 2 - radius**2) < 1e-5:
            eps_vol = eps[node_id][0] + eps[node_id][1] + eps[node_id][2]
            eps_v.append(eps_vol)
            theta = math.atan2(x[1], x[0])
            theta_p = theta
            if theta_p < 0:
                theta_p += 2 * math.pi
            phi.append(theta_p)

    sort_idx = np.argsort(phi)
    phi_sorted = [phi[i] for i in sort_idx]
    eps_v_sorted = [eps_v[i] for i in sort_idx]

    return np.array(phi_sorted), np.array(eps_v_sorted)


# %%
def get_sub_mesh(mesh, radius):
    try:
        center = np.array([0.0, 0.0, 0.0])

        def is_cell_within_circle(cell, center, radius):
            return np.linalg.norm(cell.center - center) < radius

        selected_cells = np.array(
            [
                is_cell_within_circle(mesh.get_cell(i), center, radius)
                for i in range(mesh.number_of_cells)
            ]
        )

        return mesh.extract_cells(selected_cells)

    except Exception as e:
        print(f"An error occurred: {e}")


# %%
def plot_subfigure_mechanics(
    fig,
    ax_i,
    mesh_data,
    triangular_mesh,
    variable_name,
    variable_scalar,
    label,
    plot_trace=True,
    region_limit=True,
):
    data = mesh_data.get_point_field(variable_name) * variable_scalar
    plot_data = (
        data[:, 0] + data[:, 1] + data[:, 2]
        if plot_trace
        else np.sqrt(data[:, 0] ** 2 + data[:, 1] ** 2)
    )

    contour_data = ax_i.tricontourf(triangular_mesh, plot_data, cmap="jet")
    fig.colorbar(contour_data, ax=ax_i, label=label)
    if region_limit:
        ax_i.set_xlim([-0.09894, 0.09894])
        ax_i.set_ylim([-0.09894, 0.09894])


def contour_plot(vtu_file_name, title):
    fig, ax = plt.subplots(2, 2, figsize=(8, 6))
    ax[0, 0].set_title(title, loc="left", y=1.12)
    plt.subplots_adjust(wspace=0.5)

    m_plot = vtuIO.VTUIO(vtu_file_name, dim=2)
    triang = tri.Triangulation(m_plot.points[:, 0], m_plot.points[:, 1])
    p_plot_data = m_plot.get_point_field("pressure") * 1.0e-6
    contour_p = ax[0, 0].tricontourf(triang, p_plot_data, cmap="jet")
    fig.colorbar(contour_p, ax=ax[0, 0], label="Pore pressure / MPa")

    plot_subfigure_mechanics(
        fig,
        ax[0, 1],
        m_plot,
        triang,
        "displacement",
        1e3,
        "Displacememt norm / mm",
        plot_trace=False,
        region_limit=False,
    )
    plot_subfigure_mechanics(
        fig,
        ax[1, 0],
        m_plot,
        triang,
        "sigma",
        1e-6,
        "Stress trace /MPa",
        plot_trace=True,
        region_limit=False,
    )
    plot_subfigure_mechanics(
        fig,
        ax[1, 1],
        m_plot,
        triang,
        "epsilon",
        1.0,
        "Strain trace",
        plot_trace=True,
        region_limit=False,
    )

    fig.tight_layout()
    # plt.savefig(vtu_file_name + ".png")
    plt.show()


# %%
def contour_plot_within_circle_mechanics(vtu_file_name, title):
    fig, ax = plt.subplots(ncols=3, figsize=(10, 2.8))
    ax[0].set_title(title, loc="left", y=1.12)
    plt.subplots_adjust(wspace=0.5)

    m_plot_extracted = vtuIO.VTUIO(get_extracted_mesh_name(vtu_file_name), dim=2)
    triang_extracted = tri.Triangulation(
        m_plot_extracted.points[:, 0], m_plot_extracted.points[:, 1]
    )

    plot_subfigure_mechanics(
        fig,
        ax[0],
        m_plot_extracted,
        triang_extracted,
        "sigma",
        1e-6,
        "Stress trace /MPa",
        plot_trace=True,
    )

    plot_subfigure_mechanics(
        fig,
        ax[1],
        m_plot_extracted,
        triang_extracted,
        "displacement",
        1e3,
        "Displacememt norm / mm",
        plot_trace=False,
    )

    plot_subfigure_mechanics(
        fig,
        ax[2],
        m_plot_extracted,
        triang_extracted,
        "epsilon",
        1.0,
        "Strain trace",
        plot_trace=True,
    )

    fig.tight_layout()
    # plt.savefig(vtu_file_name + "_m.png")
    plt.show()


# %% [markdown]
# ## Simulation class


# %%
class SingleOGSModel:
    """An OGS run model"""

    def __init__(self, project_file, output_prefix, out_dir):
        self.model = ot.Project(
            input_file=project_file, output_file=Path(out_dir, "modified.prj")
        )

        self.model.replace_text(output_prefix, xpath="./time_loop/output/prefix")

        self.out_dir = out_dir
        self.output_prefix = output_prefix
        self.pvd_file_name = Path(self.out_dir, self.output_prefix + ".pvd")
        self.extracted_mesh_radius = 0.065

    def run(self, use_b_bar=True):
        b_bar_status = "true" if use_b_bar else "false"
        self.model.replace_text(b_bar_status, xpath="./processes/process/use_b_bar")
        self.model.write_input()

        self.model.run_model(
            logfile=Path(self.out_dir, "out.txt"),
            args=f"-o {self.out_dir} -m .",
        )

    def write_extracted_mesh(self):
        last_vtu_name = self.get_last_vtu_name()
        self.resulted_mesh = pv.read(last_vtu_name)
        extracted_mesh = get_sub_mesh(self.resulted_mesh, self.extracted_mesh_radius)
        print(get_extracted_mesh_name(last_vtu_name))
        extracted_mesh.save(Path(get_extracted_mesh_name(last_vtu_name)))
        return get_extracted_mesh_name(last_vtu_name)

    def get_pvd_name(self):
        return self.pvd_file_name

    def get_last_vtu_name(self):
        return get_last_vtu_file_name(self.pvd_file_name)

    def get_vol_strain_of_mesh_at_circle(self):
        mesh = pv.read(self.get_last_vtu_name())
        return get_vol_strain_of_mesh_at_circle(mesh)


# %%
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
if not out_dir.exists():
    out_dir.mkdir(parents=True)


# %% [markdown]
# ## Run the model without the B-bar method

# %%
project_file = Path("point_injection_embedded_fracture_F.prj")
output_prefix = "hm_lie_bbar"

hm_LIE_model = SingleOGSModel(project_file, output_prefix, out_dir)
hm_LIE_model.run(use_b_bar=False)


# %% [markdown]
# ## Run the model with the B-bar method and extract the volume strain at the specified circle

# %%
(angles, eps_v_non_bbar) = hm_LIE_model.get_vol_strain_of_mesh_at_circle()
hm_LIE_model.run(use_b_bar=True)
(angles, eps_v_bbar) = hm_LIE_model.get_vol_strain_of_mesh_at_circle()
hm_LIE_model.write_extracted_mesh()

# %% [markdown]
# ## Check the the volume strain at the specified circle

# %%
expected_eps_v_bbar = np.array(
    [
        1.16335759e-04,
        1.27973166e-04,
        1.38837125e-04,
        1.49406935e-04,
        1.58911804e-04,
        1.67097558e-04,
        1.79090654e-04,
        1.80998160e-04,
        1.83571467e-04,
        1.85088157e-04,
        1.85298417e-04,
        1.83743323e-04,
        1.79559685e-04,
        1.77823124e-04,
        1.69285340e-04,
        1.58809634e-04,
        1.46132589e-04,
        1.32426076e-04,
        1.18260099e-04,
        1.04200550e-04,
        9.00669579e-05,
        7.48560170e-05,
        6.09281828e-05,
        4.58225558e-05,
        3.02608126e-05,
        1.46828752e-05,
        -6.53200135e-08,
        -1.44661280e-05,
        -2.89649675e-05,
        -4.28481973e-05,
        -5.58959887e-05,
        -6.75842675e-05,
        -8.13201602e-05,
        -8.53576417e-05,
        -8.98598843e-05,
        -9.27046572e-05,
        -9.52159031e-05,
        -9.62522997e-05,
        -9.34095126e-05,
        -9.26342277e-05,
        -8.42189512e-05,
        -7.36674808e-05,
        -6.16629574e-05,
        -4.87917552e-05,
        -3.57839890e-05,
        -2.31545708e-05,
        -5.51570610e-06,
        1.12625451e-05,
        3.04753264e-05,
        5.06382581e-05,
        7.08909206e-05,
        9.02398991e-05,
        1.08850194e-04,
        1.26692165e-04,
        1.42565666e-04,
        1.56895763e-04,
        1.69267995e-04,
        1.79639671e-04,
        1.87669304e-04,
        1.93339451e-04,
        1.96577904e-04,
        1.97482336e-04,
        1.96604188e-04,
        1.93682271e-04,
        1.88620312e-04,
        1.81381029e-04,
        1.72382726e-04,
        1.61626061e-04,
        1.48629371e-04,
        1.34537579e-04,
        1.19775201e-04,
        1.04703485e-04,
        8.94327916e-05,
        7.35543718e-05,
        5.77951985e-05,
        4.18858899e-05,
        2.58258442e-05,
        9.72674909e-06,
        -6.16231412e-06,
        -2.16070260e-05,
        -3.65785375e-05,
        -5.08321997e-05,
        -6.41058419e-05,
        -7.57508147e-05,
        -8.89250430e-05,
        -9.25773818e-05,
        -9.63208222e-05,
        -9.82829846e-05,
        -9.97005218e-05,
        -9.90579844e-05,
        -9.48959796e-05,
        -9.26119896e-05,
        -8.21899006e-05,
        -6.92310269e-05,
        -5.43801462e-05,
        -3.85970752e-05,
        -2.27052663e-05,
        -7.32999995e-06,
        1.42009461e-05,
        3.30283996e-05,
        5.20950635e-05,
        7.05415028e-05,
        8.76076018e-05,
        1.02890497e-04,
    ]
)

np.testing.assert_allclose(actual=eps_v_bbar, desired=expected_eps_v_bbar, atol=1e-10)


# %% [markdown]
# ## Compare the volume strains obtained with and without the B-bar method
#
# The figure below shows that the B-bar method improves the strain result significantly.

# %%
plt.rcParams["figure.figsize"] = [8, 4]
angles_in_degree = angles * 180 / math.pi
plt.plot(angles_in_degree, eps_v_non_bbar, "C1", linestyle="dotted", label="Non B-bar")
plt.plot(angles_in_degree, eps_v_bbar, "C0", label="With B-bar")

plt.title(
    "Volume strain variation at r = 0.065 m of the point injection simulation at t 2500 s"
)
plt.xlabel("Angle to the normal to the PEE 1 (clock wise) [°]")
plt.ylabel("Volume strain")
plt.legend()

# plt.savefig("ufz_HM1_point_injection_volume_strain_at_r0.06.png")
plt.show()

# %% [markdown]
# ## Contour plotting in the entire domain

# %%
contour_plot(str(hm_LIE_model.get_last_vtu_name()), "Entire domain")

# %% [markdown]
# ## Contour plotting within a circle $r<=$0.065 m
#
# The plotting highlights the variable changes in the vicinity of fracture.

# %%
contour_plot_within_circle_mechanics(
    str(hm_LIE_model.get_last_vtu_name()), r"Within r<=0.065 m"
)

# %%
