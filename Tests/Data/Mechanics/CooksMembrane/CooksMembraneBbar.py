# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "Cook's membrane example"
# date = "2024-06-11"
# author = "Wenqing Wang"
# image = "figures/cooks_membrane.png"
# web_subsection = "small-deformations"
# weight = 3
# +++

# %% [markdown]
# $$
# \newcommand{\B}{\text{B}}
# \newcommand{\F}{\text{F}}
# \newcommand{\I}{\mathbf I}
# \newcommand{\intD}[1]{\int_{\Omega_e}#1\mathrm{d}\Omega}
# $$
#
# # Cook's membrane example for nearly icompressible solid
#
# ## B bar method
# Considering a strain decomposition: $\mathbf\epsilon = \underbrace{\mathbf\epsilon- \frac{1}{3}(\epsilon:\mathbf I)}_{\text{deviatoric}}\I +  \underbrace{\frac{1}{3}(\epsilon:\mathbf I)}_{\text{dilatational}} \I$.
# The idea of the B bar method is to use another quadrature rule to interpolate the dilatational part, which leads to a modified B matrix [1]:
# $$
# 	     \bar\B = \underbrace{\B - \B^{\text{dil}}}_{\text{original B elements}}+ \underbrace{{\bar\B}^{\text{dil}}}_{\text{by another quadrature rule} }
# $$
# There are several methods to form ${\bar\B}^{\text{dil}}$ such as selective integration, generalization of the mean-dilatation formulation.
# In the current OGS, we use the latter, which reads
# $$
# 		  {\bar\B}^{\text{dil}} = \frac{\intD{\B^{\text{dil}}(\xi)}}{\intD{}}
# $$
#
# ## Example
# To verify the implementation of the B bar method, the so called Cook's membrane is used as a benchmark.
# Illustrated in the following figure, this example simulates a tapered and swept panel of unit thickness.
# The left edge is clamped and the right edge is applied with a distributed shearing load $F$ = 100 N/mm.
# The plane strain condition is considered.
# This numerical model is exactly the same as that is presented in the paper by T. Elguedj et al [1,2].
#
# <img src="figures/cooks_membrane.png" alt="Cook's membrane" width="320" height="320" />
#

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)


# %%
def run_single_test(n: int, use_bbar: bool) -> ot.MeshSeries:
    model = ot.Project(
        input_file="CooksMembrane.prj", output_file=out_dir / "modified.prj"
    )
    model.replace_text(f"mesh_n{n:02d}.vtu", xpath="./mesh")
    model.replace_text(str(use_bbar).lower(), xpath="./processes/process/use_b_bar")
    prefix = f"cooks_membrane_n_{n}_bbar_{str(use_bbar).lower()}"
    model.replace_text(prefix, xpath="./time_loop/output/prefix")
    model.replace_text("BiCGSTAB", xpath=".//solver_type")
    model.replace_text("ILUT", xpath=".//precon_type")
    model.write_input()
    model.run_model(logfile=out_dir / "out.txt", args=f"-o {out_dir} -m .")
    return ot.MeshSeries(out_dir / (prefix + ".pvd"))


def get_top_uy(mesh: ot.Mesh) -> np.ndarray:
    top_point = (48.0e-3, 60.0e-3, 0)
    p_id = mesh.find_closest_point(top_point)
    return mesh.point_data["displacement"][p_id, 1]


def compare(mesh_a: ot.Mesh, mesh_b: ot.Mesh) -> plt.Figure:
    fig, axs = plt.subplots(2, 2, figsize=[8, 6], sharex=True, sharey=True)
    u = ot.variables.displacement["y"].replace(output_unit="mm")
    sig_tr = ot.variables.stress.trace
    ot.plot.contourf(mesh_a, u, fig=fig, ax=axs[0, 0], show_edges=True)
    ot.plot.contourf(mesh_b, u, fig=fig, ax=axs[0, 1], show_edges=True)
    ot.plot.contourf(mesh_a, sig_tr, fig=fig, ax=axs[1, 0], show_edges=True)
    ot.plot.contourf(mesh_b, sig_tr, fig=fig, ax=axs[1, 1], show_edges=True)
    axs[0, 0].set_title("bbar = false")
    axs[0, 1].set_title("bbar = true")
    ot.plot.utils.update_font_sizes(fig.axes, 10)
    return fig


n_range = [4, 10, 15, 20, 25, 30]
# %%
results_bbar_false = {n: run_single_test(n, use_bbar=False) for n in n_range}
uy_top_bbar_false = [get_top_uy(ms[-1]) for ms in results_bbar_false.values()]
uy_top_bbar_false_ref = np.array(
    [0.002164586784123102, 0.0022603329644579383, 0.002375295856067169,
     0.002519725590136146, 0.0026515294133790837, 0.002868289617025223]
)  # fmt: skip
np.testing.assert_allclose(uy_top_bbar_false, uy_top_bbar_false_ref, atol=1e-10)

# %%
results_bbar_true = {n: run_single_test(n, use_bbar=True) for n in n_range}
uy_top_bbar_true = [get_top_uy(ms[-1]) for ms in results_bbar_true.values()]
uy_top_bbar_true_ref = np.array(
    [0.006957471385697900, 0.007772616910217863, 0.007897597955618913,
     0.007951479575082158, 0.007976349858390623, 0.007999718483861992]
)  # fmt: skip
np.testing.assert_allclose(uy_top_bbar_true, uy_top_bbar_true_ref, atol=2e-4)

# %% [markdown]
# ## Result
#
# ### Vertical diplacement at the top point
#
# The following figure shows that the convergence of the solutions obtained by using the B bar method follows the one presented in the paper by T. Elguedj et al [1]. However, the results obtained without the B bar method are quit far from the converged solution with the finest mesh.

# %%
fig, ax = plt.subplots(figsize=(4, 3))
ax.plot(n_range, np.asarray(uy_top_bbar_true) * 1e3, "--o", label="B bar")
ax.plot(n_range, np.asarray(uy_top_bbar_false) * 1e3, "--x", label="non B bar")
ax.set_xlabel("Number of elements per side")
ax.set_ylabel("Top right corner displacement /mm")
ax.legend()
fig.tight_layout()

# %% [markdown]
# ### Comparison of results

# %%
for n in n_range:
    compare(results_bbar_false[n][-1], results_bbar_true[n][-1])

# %% [markdown]
# The contour plots show that even with the coarsest mesh, the B bar method still gives reasonable stress results.
#
# ## Reference
#
# [1] T.J.R. Hughes (1980). Generalization of selective integration procedures to anisotropic and nonlinear media. International Journal for Numerical Methods in Engineering, 15(9), 1413-1418.
#
# [2] T. Elguedj, Y. Bazilevs, V.M. Calo, T.J.R. Hughes (2008),
#  $\bar\B$ and $\bar\F$ projection methods for nearly incompressible linear and non-linear elasticity and plasticity using higher-order NURBS elements, Computer Methods in Applied Mechanics and Engineering, 197(33--40), 2732-2762.
#
