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
from datetime import datetime
from pathlib import Path
from subprocess import run

import matplotlib.pyplot as plt
import vtuIO

# %%
prj_name = "square_1e1_neumann"
data_dir = os.environ.get("OGS_DATA_DIR", "../../../Data")
prj_file = f"{data_dir}/EllipticPETSc/{prj_name}.prj"

# %%
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
if not out_dir.exists():
    out_dir.mkdir(parents=True)

print(f"mpirun --bind-to none -np 2 ogs {prj_file} > out.txt")
run(f"mpirun --bind-to none -np 2 ogs {prj_file} > out.txt", check=True)

print(datetime.now())


# %%
pvdfile = vtuIO.PVDIO(f"{prj_name}.pvd", dim=2)
time = pvdfile.timesteps
points = {"pt0": (0.3, 0.5, 0.0), "pt1": (0.24, 0.21, 0.0)}
pressure_linear = pvdfile.read_time_series("pressure", points)

# %%
plt.plot(time, pressure_linear["pt0"], "b-", label="pt0 linear interpolated")
plt.plot(time, pressure_linear["pt1"], "r-", label="pt1 linear interpolated")
plt.legend()
plt.xlabel("t")
plt.ylabel("p")
