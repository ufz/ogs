# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: .venv
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "McWhorter & Sunada Problem"
# date = "2022-10-19"
# author = "Norbert Grunwald"
# image = "figures/placeholder_mcWhorter.png"
# web_subsection = "th2m"
# coupling = "h2"
# weight = 4
# +++
#

# %% [markdown]
# |<div style="width:330px"><img src="https://www.ufz.de/static/custom/weblayout/DefaultInternetLayout/img/logos/ufz_transparent_de_blue.png" width="300"/></div>|<div style="width:330px"><img src="https://discourse.opengeosys.org/uploads/default/original/1X/a288c27cc8f73e6830ad98b8729637a260ce3490.png" width="300"/></div>|<div style="width:330px"><img src="https://github.com/nagelt/Teaching_Scripts/raw/9d9e29ecca4b04eaf7397938eacbf116d37ddc93/Images/TUBAF_Logo_blau.png" width="300"/></div>|
# |---|---|--:|

# %% [markdown]
# # McWhorter Problem
# <cite>[McWhorter and Sunada][1]</cite> propose an analytical solution to the two-phase flow equation. A one-dimensional problem was considered which describes the flow of two incompressible, immiscible fluids through a porous medium, where the wetting phase (water) displaces the non-wetting fluid (air or oil) in the horizontal direction (without the influence of gravity).
#
#
# ![mcWhorter_concept_schematic.png](figures/mcWhorter_concept.png)
#
# ## Analytical solution
#
#
# A detailed semi-analytical solution and a convenient tool for calculating the solution for different material parameters can be found [here](https://mmg.fjfi.cvut.cz/~fucik/index.php?page=exact).
#
# ### Material Parameters
#
# | Property                                                       | Symbol                 | Value                 | Unit   |
# | -------------------------------------------------------------: | :--------------------- | --------------------- | ------ |
# | Porosity                                                       | $\phi$                 | 0.15                  | 1      |
# | Intrinsic permeability                                         | $K$                    | $$1.0\cdot 10^{-10}$$ | $m^2$  |
# | Residual saturation of the wetting phase                       | $$s_\mathrm{L}^{res}$$ | 0.02                  | 1      |
# | Residual saturation of the non-wetting phase                   | $$s_\mathrm{G}^{res}$$ | 0.001                 | 1      |
# | Dynamic viscosity of the wetting phase                         | $\mu_\mathrm{L}$       | $$1.0\cdot 10^{-3}$$  | Pa s   |
# | Dynamic viscosity of the non-wetting pha                       | $\mu_\mathrm{G}$       | $$5.0\cdot 10^{-3}$$  | Pa s   |
# | Brooks and Corey model parameter: entry pressure               | $p_b$                  | 5000                  | Pa     |
# | Brooks and Corey model parameter: pore size distribution index | $\lambda$              | 3.0                   | 1      |
#
# ### Problem Parameters
#
# | Property                      | Symbol                | Value   | Unit   |
# | ----------------------------: | :-------------------: | ------- | ------ |
# | Initial saturation            | $$s_\mathrm{L}(t=0)$$ | 0.05    | 1      |
# | Injection boundary saturation | $$s_\mathrm{L}(x=0)$$ | 0.8     | 1      |
#
#
# [1]: https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/WR026i003p00399?casa_token=6yzGcmrd7dkAAAAA:E6QsKTxrf12GO-0CY6qgu4XEcX6iFM4O_mnaVV2gWBO8voVnxYXxLOtnAdUnBskEOPZiwaFAggWnmqpg

# %% [markdown]
# ## Exact Solution
#
# The exact solution is not yet calculated in this notebook, instead the [online tool](https://mmg.fjfi.cvut.cz/~fucik/index.php?page=exact) by Radek Fučík is used. This tool calculates the solution and outputs the results with arbitrary accuracy as CSV files, which are plotted below.

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pyvista as pv

# %%
# Import analytical solution from a CSV file
reference = np.loadtxt("data/ref_solution_saturation.csv", delimiter=",")
ref_x, ref_saturation = reference.T

# %% [markdown]
# ## Numerical Solution

# %%
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

# %%
# run OGS
model = ot.Project(input_file="mcWhorter_h2.prj", output_file="mcWhorter_h2.prj")
model.run_model(logfile=f"{out_dir}/out.txt", args=f"-o {out_dir}")

# %%
# Read the OGS results and sample the last timestep on the same points as in
# the reference.
ms = ot.MeshSeries(f"{out_dir}/result_McWhorter_H2.pvd")
pts = pv.PointSet([(x, 0, 0) for x in ref_x])
sL_num = pts.sample(ms.mesh(-1))["saturation"]

# Absolute and relative errors
err_abs = ref_saturation - sL_num
err_rel = err_abs / ref_saturation
assert np.all(np.abs(err_rel) < 0.5)

# %%
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4))
ax3 = ax2.twinx()

ax1.plot(ref_x, sL_num, "b", label=r"$s_\mathrm{L}$ numerical")
ax1.plot(ref_x, ref_saturation, "g", label=r"$s_\mathrm{L}$ exact")
lns_abs = ax2.plot(ref_x, err_abs, "-b", lw=2, label=r"absolute error")
lns_rel = ax3.plot(ref_x, err_rel, "-g", lw=2, label=r"relative error")

lns = lns_abs + lns_rel
ax1.legend()
ax2.legend(lns, [label.get_label() for label in lns])

ax1.set_xlabel("$x$ / m", fontsize=12)
ax1.set_ylabel(r"$s_\mathrm{L}$", fontsize=12)
ax2.set_xlabel("$x$ / m", fontsize=12)
ax2.set_ylabel(r"$\epsilon_\mathrm{abs}$", fontsize=12)
ax3.set_ylabel(r"$\epsilon_\mathrm{rel}$", fontsize=12)
time = ms.timevalues[-1]
fig.suptitle(f"Liquid saturation and errors at t={time} seconds")

# %% [markdown]
# The numerical approximation fits the analytical solution very well in the whole area. Only in the area of the saturation front are deviations recognisable. These deviations are mainly due to the grid resolution and are well known in multiphase simulations. The error can be reduced almost arbitrarily by lowering the size of the grid elements.
