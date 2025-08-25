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
# title = "Two-layer diffusion problem"
# date = "2022-03-09"
# author = "Renchao Lu, Dmitri Naumov, Lars Bilke, Christoph Lehmann, Haibing Shao"
# web_subsection = "hydro-component"
# +++
#

# %% [markdown]
# # Two-layer diffusion problem

# %% [markdown]
# ## Problem description

# %% [markdown]
# In waste repositories, radionuclide release can be expected after rupture of waste canisters to occur in the engineered barrier system, which contains multiple layers of materials and host rocks. In this benchamrk, a tracer (HTO) diffusion process through a two-layer barrier is simulated. The barrier is comprised of a bentonite buffer layer and an opalinus clay (OPA) layer.
#
# Over the one-dimensional model domain, the diffusion process of HTO can be described with the following governing equation:
#
# $$
# \frac{\partial \left( \phi c\right)}{\partial t} = \frac{\partial}{\partial x} \left( \phi \mathrm{D_p} \frac{\partial c}{\partial x} \right),
# $$
#
# where $c$ [mol/m$^3$] represents the HTO concentration. $\mathrm{D_p}$ [m$^2$/s] is the pore diffusion coefficient for HTO, and $\phi$ [-] is the porosity of the media.

# %% [markdown]
# ## Model setups

# %% [markdown]
# The computational domain is assumed to be 20 meters long. It consists of a 0.625 meter thick layer of bentonite buffer, and the rest is filled with OPA. The simulation time is one million years. Initially, the entire domain is assumed to be solute free, i.e. $c_{\mathrm{ini}}(x, 0) = 0$. The inlet concentration is held at 1 mol/L throughout the simulation, i.e. $c(0, t) = 1$ mol/L. In the numerical model, the spatial domain is discretized by linear line elements with a length of 0.005 meter each. The time step size of 1000 years is used in the simulation. The discretized governing equation is iteratively solved using the Newton-Raphson method.
#
# The table below summarizes the parameters used in the simulation.
#
# |                         Parameter                          |  Value   |  Unit   |
# | :--------------------------------------------------------: | :------: | :-----: |
# |         Porosity of bentonite $\phi_{\mathrm{b}}$          |   0.36   |    -    |
# |           Porosity of OPA $\phi_{\mathrm{OPA}}$            |   0.12   |    -    |
# | Pore diffusion coefficient in bentonite $\mathrm{D_{p,b}}$ | 5.55e-10 | m$^2$/s |
# |   Pore diffusion coefficient in OPA $\mathrm{D_{p,OPA}}$   | 8.33e-11 | m$^2$/s |
# |                 Time step size $\Delta t$                  |   1e3    |  year   |
# |                    Grid size $\Delta x$                    |   0.01   |    m    |
#
# Notes: The parameter values are sourced from Nagra (2002).

# %% [markdown]
# ## Results

# %% [markdown]
# **Analytical solution**
#
# For a two-layer diffusion problem, it is difficult to obtain the exact analytical solution. Instead, Carr and Turner (2016) presented a semi-analytical solution for this type of problem and released the corresponding <a href="https://github.com/elliotcarr/MultDiff">MATLAB script</a> for public usage.
#
# Here we provide a <a href="https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/ComponentTransport/MultiLayerDiffusion/AnalyticalSolutionInput.m">MATLAB script</a> through which we can compute the concentration profiles along the two-layer domain. The following figure illustrates the concentration profiles at $t$ = 10$^3$, 10$^4$, 10$^5$, and 10$^6$ years (see calculated values in <a href="https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/ComponentTransport/MultiLayerDiffusion/SemiAnalyticalSolutionResults.csv">SemiAnalyticalSolutionResults.csv</a>).
#
#

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pandas as pd
from IPython.display import Image, display

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

# %%
ref_data_file = "./SemiAnalyticalSolutionResults.csv"
soln = pd.read_csv(
    ref_data_file, sep=",", header=None, skiprows=0,
    names=["x", "1e3", "1e4", "1e5", "1e6"], index_col=False,
)  # fmt: skip


def plot_results(ms: ot.MeshSeries, var: ot.variables.Variable) -> plt.Figure:
    "Plot numerical results against analytical solution if available"
    n_cols = 2 if var.data_name == "HTO" else 1
    figsize = (11, 4) if n_cols == 2 else (8, 6)
    fig_res, axs_ = plt.subplots(1, n_cols, figsize=figsize)
    axs: list[plt.Axes] = [axs_] if n_cols == 1 else axs_
    labels = [f"{tv:.0e} a" for tv in ms.timevalues]
    ot.plot.line(ms, var, ax=axs[0], labels=labels)
    axs[0].plot([], [], "-k", label="OGS-numerical")

    if n_cols == 2:
        axs[0].plot([], [], "--k", label="semi-analytical")
        data = soln[["1e3", "1e4", "1e5", "1e6"]].to_numpy().T
        x = ms[0].points[:, 0]
        ref_vals = np.asarray([np.interp(x, soln["x"], vals) for vals in data])
        abs_err = ms[var.data_name] - ref_vals

        ms.point_data["ref_vals"] = ref_vals
        ms.point_data[var.abs_error.data_name] = abs_err
        ot.plot.line(ms, "ref_vals", ax=axs[0], ls="--")
        ot.plot.line(ms, var.abs_error, ax=axs[1])
        axs[1].set_yscale("symlog", linthresh=0.001)

        max_errors = [2e-1, 2e-2, 4e-3, 5e-4]  # per timevalue
        assert np.all(np.abs(abs_err.T) <= max_errors)

    axs[0].legend(loc="best", fontsize=10)
    for ax in axs:
        ax.axvspan(0, 0.625, color="royalblue", alpha=0.2)  # bentonite layer
        ax.axvspan(0.625, 20, color="orange", alpha=0.1)  # OPA host rock
        ax.margins(x=0)
    ot.plot.utils.update_font_sizes(axs, 10)

    fig_res.tight_layout()
    return fig_res


# %% [markdown]
# **Numerical solution**

# %% [markdown]
# Correspondingly, the OGS input files of this 1D mass transport benchmark can be found <a href="https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/ComponentTransport/MultiLayerDiffusion/1D_MultiLayerDiffusion.prj">here</a>.
#
# Then, the numerical solution by OpenGeoSys is plotted against the semi-analytical solution for comparison.

# %%
# Run OGS simulation
name = "1D_MultiLayerDiffusion"
model = ot.Project(
    input_file=f"{name}.prj", output_file=out_dir / f"{name}_modified.prj"
)
model.write_input()
model.run_model(logfile=out_dir / "out.txt", args=f"-o {out_dir} -m .")

# %% Read simulation results
ms_MLD = ot.MeshSeries(out_dir / f"{name}.pvd").scale(time=("s", "a"))[1:]
HTO_var = ot.variables.Scalar("HTO", "mol/L", "mol/L", "HTO concentration")
fig = plot_results(ms_MLD, HTO_var)


# %% [markdown]
# In the first time step, the semi-analytical and numerical solutions do not agree so well. As the time step progresses, they begin to agree more closely.

# %% [markdown]
# **Mass flux calculation**

# %% [markdown]
# Here is a sketch that shows how we calculate the molar flux at the node.

# %%
display(Image(filename="./sketch_molar_flux_calculation.jpg", width=100))


# %% [markdown]
# Additionally, we compute the molar flux profiles at $t$ = 10$^3$, 10$^4$, 10$^5$, and 10$^6$ years. The implementation of molar flux output can be viewed <a href="https://gitlab.opengeosys.org/ogs/ogs/-/merge_requests/4006">at this link</a>.

# %%
HTO_flux_var = ot.variables.Scalar("HTOFlux", "mol/m$^2$/s")
fig = plot_results(ms_MLD, HTO_flux_var)
fig.axes[0].set_yscale("log")


# %% [markdown]
# ## References

# %% [markdown]
# <a href="https://web.archive.org/web/20220314123055/https://backend.nagra.ch/system/files/attachments/ntb/ntb2002/e_ntb02-006.pdf">Nagra, 2002</a>. Project Opalinus Clay: Models, Codes and Data for Safety Assessment. Technical Report NTB 02–06. Nagra, Switzerland.
#
# E. J. Carr and I. W. Turner (2016), A semi-analytical solution for multilayer diffusion in a composite medium consisting of a large number of layers, Applied Mathematical Modelling, 40: pp. 7034–7050. http://dx.doi.org/10.1016/j.apm.2016.02.041

# %% [markdown]
# **Credits:**

# %% [markdown]
# Renchao Lu, Dmitri Naumov, Lars Bilke, Christoph Lehmann, Haibing Shao
