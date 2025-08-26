# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "Cook's membrane example"
# date = "2024-07-03"
# author = "Wenqing Wang & Thomas Nagel"
# image = "figures/cooks_membrane.png"
# web_subsection = "large-deformations"
# weight = 3
# +++

# %% [markdown]
# $$
# \newcommand{\J}{\text{J}}
# \newcommand{\F}{\text{F}}
# $$
#
# # Cook's membrane example for nearly incompressible solid under large deformation
#
# The example has been analyzed in many references, for example with the F-bar method in [1] and [2].
# This example is also used as benchmark for [SmallDeformation](https://www.opengeosys.org/docs/benchmarks/small-deformations/cooksmembranebbar/). Hereby we analyze it again under finite strain assumption using the F-bar method for the total Lagrange formulation (see the attached [PDF](figures/F-bar-ogs.pdf) for its theory).
# For the finite strain assumption, the constitutive law is replaced with a hyperelasticity one, the Neo-Hookean. The Neo-Hookean model defines an energy function as the sum of volumetric component $U_{\text{dil}}(\F)$ and deviatoric component $U_{\text{dev}}(\F)$ as
# $$
# \begin{align}
# W(\F) = U_{\text{dil}}(\F) + U_{\text{dev}}(\F),
# \end{align}
# $$
# where
# $$
# \begin{align}
# & U_{\text{dil}}(\F) = \dfrac{1}{2} K (\det(\F)-1)^2\\
# & U_{\text{dev}}(\F) = \dfrac{1}{2} G \left(\text{tr} (\det(\F)^{-\frac{2}{3}}\F\F^{\text{T}})-3\right)
# \end{align}
# $$
# with $K$ the bulk modulus, and $G$ the shear modulus. The values of $K$ and $G$ are taken from references [1] and [2], which are $40.0942\cdot 10^{4}$ MPa and $80.1938$ MPa, respectively. For OGS input, the corresponding Young's modulus and Poisson ratio are 240.565 MPa, and 0.499, respectively.
#
#
# ## Reference
#
# 1. E.A. de Souza Neto, D. Perić, M. Dutko, D.R.J. Owen, [Design of simple low order finite elements for large strain analysis of nearly incompressible solids](https://doi.org/10.1016/0020-7683(95)00259-6), International Journal of Solids and Structures, Volume 33, Issues 20–22, 1996, Pages 3277-3296.
#
# 2. T. Elguedj, Y. Bazilevs, V.M. Calo, T.J.R. Hughes (2008),
#  $\bar\B$ and $\bar{\text F}$ projection methods for nearly incompressible linear and non-linear elasticity and plasticity using higher-order NURBS elements, Computer Methods in Applied Mechanics and Engineering, 197(33–40), 2732-2762.
#

# %%
# Standard library imports
import os
from pathlib import Path

# Third-party imports
import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)


# %%
def get_top_uy(pvd_file_name):
    top_point = (48.0e-3, 60.0e-3, 0)
    mesh = ot.MeshSeries(out_dir / pvd_file_name)[-1]
    p_id = mesh.find_closest_point(top_point)
    return mesh.point_data["displacement"][p_id, 1]


# %%
def run_single_test(mesh_name, output_prefix, use_fbar=False, use_load_increment=False):
    prj = ot.Project(
        input_file="CooksMembrane.prj",
        output_file=out_dir / "modified.prj",
    )

    prj.replace_text(mesh_name, xpath="./mesh")
    if not use_fbar:
        prj.replace_text("none", xpath="./processes/process/f_bar")
    prj.replace_text(output_prefix, xpath="./time_loop/output/prefix")
    vtu_file_name = output_prefix + "_ts_1_t_1.000000.vtu"
    prj.replace_text(vtu_file_name, xpath="./test_definition/vtkdiff[1]/file")
    prj.replace_text(vtu_file_name, xpath="./test_definition/vtkdiff[2]/file")
    prj.replace_text(vtu_file_name, xpath="./test_definition/vtkdiff[3]/file")

    if use_load_increment:
        prj.replace_text(
            0.5,
            xpath="./time_loop/processes/process[1]/time_stepping/timesteps/pair/delta_t",
        )
        prj.replace_text(
            "FRamp",
            xpath="./process_variables/process_variable/boundary_conditions/boundary_condition[3]/parameter",
        )

    prj.write_input()
    prj.run_model(logfile=f"{out_dir}/out.txt", args=f"-o {out_dir} -m .")

    return get_top_uy(output_prefix + ".pvd")


# %%
nedges = [4, 10, 15, 20, 25, 30]
mesh_names = [f"mesh{n}.vtu" for n in [""] + [f"_n{ne}" for ne in nedges[1:]]]
load_increment_labels = [False, False, False, True, True, False]
output_prefices = [
    f"cooks_membrane_ld_{suffix}"
    for suffix in ["edge_div_4"] + [f"_refined_mesh_{ne}" for ne in nedges[1:]]
]

uys_at_top_fbar = []
for mesh_name, load_increment_label, output_prefix in zip(
    mesh_names, load_increment_labels, output_prefices
):
    uy_at_top = run_single_test(
        mesh_name,
        output_prefix,
        use_fbar=True,
        use_load_increment=load_increment_label,
    )
    uys_at_top_fbar.append(uy_at_top)

print(uys_at_top_fbar)


# %%
expected_uys_at_top_fbar = np.array(
    [
        0.006141383357250432,
        0.006746283955378773,
        0.0068241828268382505,
        0.0068589364838052315,
        0.006873125862208623,
        0.006891409184641067,
    ]
)
np.testing.assert_allclose(
    actual=uys_at_top_fbar, desired=expected_uys_at_top_fbar, atol=2e-4
)

# %%
output_prefices_non_fbar = [prefix + "_non_fbar" for prefix in output_prefices]

uys_at_top_non_fbar = []
for mesh_name, load_increment_label, output_prefix in zip(
    mesh_names, load_increment_labels, output_prefices_non_fbar
):
    uy_at_top = run_single_test(
        mesh_name,
        output_prefix,
        use_fbar=False,
        use_load_increment=load_increment_label,
    )
    uys_at_top_non_fbar.append(uy_at_top)

print(uys_at_top_non_fbar)


# %%
expected_uys_at_top_non_fbar = np.array(
    [
        0.0022867221436878916,
        0.002840692165858716,
        0.0033745807476269606,
        0.003881199717293749,
        0.0042429289476765735,
        0.004702757454266369,
    ]
)
np.testing.assert_allclose(
    actual=uys_at_top_non_fbar, desired=expected_uys_at_top_non_fbar, atol=1e-10
)


# %%
def plot_data(u_y_fbar, uy_non_fbar, file_name=""):
    # Plotting
    plt.rcParams["figure.figsize"] = [5, 5]

    if len(u_y_fbar) != 0:
        plt.plot(nedges, np.array(u_y_fbar) * 1e3, "o--", label="F bar")
    if len(uy_non_fbar) != 0:
        plt.plot(nedges, np.array(uy_non_fbar) * 1e3, "x--", label="non F bar")

    plt.xlabel("Number of elements per side")
    plt.ylabel("Top right corner displacement /mm")
    plt.legend()

    plt.tight_layout()
    if file_name != "":
        plt.savefig(file_name)
    plt.show()


# %% [markdown]
# ## Result
#
# ### 1. Vertical diplacement at the top point
#
# The following figure shows that the convergence of the solutions obtained by using the F bar method follows the one presented in the paper by T. Elguedj et al [1]. However, the results obtained without the F bar method are quit far from the converged solution with the finest mesh.

# %%
plot_data(uys_at_top_fbar, uys_at_top_non_fbar, "f_bar_linear.png")

# %% [markdown]
# ### 2. Contour plot of the results


# %%
def contour_plot(pvd_file_name, title, fig_name=None):
    last_mesh = ot.MeshSeries(Path(out_dir, pvd_file_name))[-1]
    fig, axs = plt.subplots(1, 2, figsize=(6, 4))

    last_mesh.plot_contourf(ot.variables.displacement[1], ax=axs[0], fontsize=8)
    last_mesh.plot_contourf(ot.variables.stress["yy"], ax=axs[1], fontsize=8)
    plt.title(title)
    if fig_name:
        plt.savefig(fig_name, bbox_inches="tight", dpi=300)


# %% [markdown]
# #### 2.1. Non F-bar: Vertical displacement (left column) and vertical stress (right column):

# %%
for nedge, output_prefix in zip(nedges, output_prefices_non_fbar):
    contour_plot(
        output_prefix + ".pvd",
        f"Number of elements per side: {nedge}",
        f"non_fbar_elements_per_side_{nedge}.png",
    )

# %% [markdown]
# #### 2.1. F-bar: Vertical displacement (left column) and vertical stress (right column):

# %%
for nedge, output_prefix in zip(nedges, output_prefices):
    contour_plot(
        output_prefix + ".pvd",
        f"Number of elements per side: {nedge}",
        f"fbar_elements_per_side_{nedge}.png",
    )

# %% [markdown]
# The contour plots show that even with the coarsest mesh, the F bar method still gives reasonable stress result.

# %%

# %%
