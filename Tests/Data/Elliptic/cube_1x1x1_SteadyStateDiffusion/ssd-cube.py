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
# title = "SteadyStateDiffusion Cube Test"
# date = "2021-11-09"
# author = "Lars Bilke"
# web_subsection = "elliptic"
# draft = true
# +++
#

# %% [markdown]
# This notebook is just for testing and does not show up on the web page because `draft = true` in the metadata.

# %%
import os
from pathlib import Path

import pyvista as pv

# %%
# On CI out_dir is set to the notebooks directory inside the build directory
# similar to regular benchmark tests. On local testing it will output to the
# notebooks source directory under a _out-subdirectory.
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

if "CI" in os.environ:
    pv.set_jupyter_backend("static")
else:
    pv.set_jupyter_backend("client")


# %%
resolution = "2e4"
# ! ogs cube_{resolution}.prj -o {out_dir} > {out_dir}/log.txt


# %%
reader = pv.get_reader(f"{out_dir}/cube_{resolution}.pvd")
reader.set_active_time_value(1.0)  # go to 1 s
mesh = reader.read()[0]

plotter = pv.Plotter(notebook=True)
plotter.add_mesh(mesh, scalars="v")  # pressure
plotter.show()
