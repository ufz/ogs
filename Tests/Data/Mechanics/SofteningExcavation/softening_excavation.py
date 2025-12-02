# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: 'Python 3.10.8 (''.venv'': venv)'
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "Material Softening Excavation"
# date = "2025-11-25"
# author = "Florian Zill"
# image = "figures/softening_excavation.png"
# web_subsection = "small-deformations"
# +++
#

# %% [markdown]
# ## Problem description
#
# In the planning of tunnel excavations, one possible dimensioning approach involves softening the to be excavated tunnel material until a predefined threshold — such as displacement, stress, or damage — is exceeded.
# To accurately account for stabilizing effects (e.g., shear stiffness) during the relaxation phase, the conventional NodalForceRelease approach is insufficient.
# Instead, a specialized material model can be applied to the tunnel region.
# This model governs the material softening process through the user-defined parameter `strength`.
# In this example the tunnel is softened down to 0% material strength.

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

# %%
prj_file = "softening_excavation.prj"
prj = ot.Project(input_file=prj_file, output_file=prj_file)
prj.run_model(logfile=f"{out_dir}/out.txt", args=f"-o {out_dir}")
results = ot.MeshSeries(f"{out_dir}/SD_softening_excavation.pvd")

# %% [markdown]
# ## Plot results
#
# ### Displacement evolution of the tunnel top point

# %%
fig, ax = plt.subplots(figsize=[8, 4])

tunnel_top = [0.0, -898.0, 0.0]
probe = ot.MeshSeries.extract_probe(results, tunnel_top, "displacement").scale(
    time=("s", "d")
)
ot.plot.line(probe, "time", "displacement", ax=ax, fontsize=10, marker=".")
assert np.isclose(probe[-1]["displacement"][:, 1], -0.10548506)

# %% [markdown]
# ### Contourplots of final timestep (zero strength)

# %%
ot.plot.setup.min_ax_aspect = None
ot.plot.setup.combined_colorbar = False
fig, axs = plt.subplots(1, 3, sharey=True, figsize=(12, 4))
ot.plot.contourf(results[-1], ot.variables.stress.trace, fig=fig, ax=axs[0])
ot.plot.contourf(results[-1], ot.variables.displacement["x"], fig=fig, ax=axs[1])
ot.plot.contourf(results[-1], ot.variables.displacement["y"], fig=fig, ax=axs[2])
ot.plot.utils.update_font_sizes(fig.axes, 12)
_ = axs[0].set(ylim=(-915, -885))

# %%
curve = prj.tree.find("./curves/curve[name='strength_curve']")
curve_x = [float(s) for s in curve.find("coords").text.split(" ")]
curve_y = [float(s) for s in curve.find("values").text.split(" ")]
target_strength = np.interp(np.asarray(results.timevalues), curve_x, curve_y)

center = results[0].threshold([1, 1], scalars="MaterialIDs").center
sig_trace = results.probe(center, ot.variables.stress.trace)
res_strength = sig_trace / sig_trace[0]

np.testing.assert_allclose(res_strength, target_strength)
