# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3.9.13 64-bit
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "SimplePETSc"
# date = "2021-11-09"
# author = "Lars Bilke"
# web_subsection = "elliptic"
# +++
#

# %% [markdown]
# The following shows running a simple steady-state diffusion benchmark running on 2 cores.

# %%
import os
from pathlib import Path
from subprocess import run

import numpy as np
import ogstools as ot

# %%
prj_name = "square_1e1_neumann"
data_dir = os.environ.get("OGS_DATA_DIR", "../../../Data")
prj_file = f"{data_dir}/EllipticPETSc/{prj_name}.prj"

# %%
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

command = f"mpirun --bind-to none -np 2 ogs {prj_file} > out.txt"
print(command)
run(command, shell=True, check=True)

# %%
mesh_series = ot.MeshSeries(f"{prj_name}.pvd").scale(time=("s", "a"))
points_coords = np.array([[0.3, 0.5, 0.0], [0.24, 0.21, 0.0]])
labels = [f"{label} linear interpolated" for label in ["pt0", "pt1"]]

ms_pts = ot.MeshSeries.extract_probe(mesh_series, points_coords)
fig = ot.plot.line(
    ms_pts, "time", ot.variables.pressure, labels=labels, colors=["b", "r"]
)
