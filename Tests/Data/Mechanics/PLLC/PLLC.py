# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "Power Law Linear Creep"
# date = "2023-01-02"
# author = "Florian Zill"
# web_subsection = "small-deformations"
# +++
#

# %% [markdown]
# ### Power Law Linear Creep
#
# This benchmark shows the increased creep rate of salt rock at lower deviatoric stress. A two component power law (which we call Power Law Linear Creep, or short PLLC) provides an easy way to capture the power law behaviour (dislocation creep) and the linear behaviour (pressure solution creep). For more details have a look at (Zill et al., 2022).
#

# %%
import contextlib
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot

# %%
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

prj_file = Path("uniax_compression.prj")
ogs_model = ot.Project(input_file=prj_file, output_file=out_dir / prj_file)

# %% [markdown]
# ### Experimental data
#
# A nice overview for the strain rates of salt for different temperatures and differential stresses can be found in (Li et al., 2021).
#
# ### Parameters
#
# This set of parameters gives a good fit with the experimental data. The grain size is a bit larger than the usual grain size of roughly 1 cm.

# %%
A1 = 0.18  # d^-1
Q1 = 54e3  # kJ / mol
A2 = 6.5e-5  # m^3 K d^−1 # noqa: RUF003
Q2 = 24.5e3  # kJ / mol
dGrain = 5e-2  # m
sref = 1.0  # MPa


def BGRa(sig, T):
    return A1 * np.exp(-Q1 / (8.3145 * (273.15 + T))) * np.power(sig / sref, 5.0)


def PLLC(sig, T):
    T_K = T + 273.15
    return (
        A1 * np.exp(-Q1 / (8.3145 * T_K)) * np.power(sig / sref, 5.0)
        + A2 * np.exp(-Q2 / (8.3145 * T_K)) * sig / sref / np.power(dGrain, 3) / T_K
    )


# %% [markdown]
# ### Simulation and plot
#
# The experimental data is compared against the model results (analytically and numerically)

# %%
lo_stresses = np.array([0.2e6, 0.6e6])
hi_stresses = np.array([2e6, 10e6])
Exps = {7.8: "blue", 14.3: "orange", 25: "lime", 60: "red", 100: "gray", 200: "purple"}

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.set_xlabel("$\\sigma_\\mathrm{ax}$ / MPa")
ax.set_ylabel("$\\dot{\\epsilon}_{zz}$ / d$^{-1}$")
ax.set(xlim=(0.15, 30), ylim=(1e-15, 1e1))
ax.minorticks_on()
ax.grid(which="major", color="lightgrey")
ax.grid(which="minor", color="0.95")

sigs = np.logspace(-1, 2, 100)
for temp, col in Exps.items():
    if temp >= 25:  # plot analytical curves
        ax.plot(sigs, BGRa(sigs, temp), color=col, ls="--")
        ax.plot(sigs, PLLC(sigs, temp), color=col, ls="-", label=f"{temp}°C")

    eps_dot = []
    ogs_model.replace_parameter_value("T_ref", str(temp + 273.15))
    stresses = hi_stresses if temp >= 25 else lo_stresses

    for stress in stresses:
        ogs_model.replace_parameter_value("sigma_ax", str(-stress))
        ogs_model.write_input()

        with contextlib.redirect_stdout(None):  # hide output
            ogs_model.run_model(args=f"-m . -o {out_dir}", logfile=out_dir / "out.txt")
            results = ot.MeshSeries(out_dir / f"{prj_file.stem}.pvd")

        eps_zz = results.point_data["epsilon"][:, 0, 2]  # pt0 z-component
        eps_zz_dot = np.abs(np.diff(eps_zz)) / np.diff(results.timevalues)
        eps_dot += [np.mean(eps_zz_dot[1:])]

    np.testing.assert_allclose(eps_dot, PLLC(1e-6 * stresses, temp), rtol=1e-4)
    ax.loglog(1e-6 * stresses, eps_dot, "o", c=col, mec="k")

ax.plot([], [], c="k", label="PLLC")
ax.plot([], [], c="k", ls="--", label="BGRa")
ax.plot([], [], c="w", ls="None", marker="o", mec="k", label="OGS")

markers = {"WIPP": "^", "DeVries": "s", "Berest": "P"}
with Path("reference_data.json").open() as ref_data_file:
    ref_data: dict[str, list[float]] = json.load(ref_data_file)
    for key, data in ref_data.items():
        temp = key.rsplit(" ", 1)[-1]
        stresses, eps_dot = np.array(data).T
        marker = markers[key.split(" ", 1)[0]]
        ax.loglog(stresses, eps_dot, marker, c=Exps[float(temp)])

ax.plot([], [], c="k", ls="None", marker="s", label="DeVries (1988)")
ax.plot([], [], c="k", ls="None", marker="^", label="WIPP CS")
ax.plot([], [], c="b", ls="None", marker="P", label="Bérest (2017) 7.8°C")
ax.plot([], [], c="orange", ls="None", marker="P", label="Bérest (2015) 14.3°C")
ax.legend()
_ = fig.tight_layout()


# %% [markdown]
# ### References
#
# Zill, Florian, Wenqing Wang, and Thomas Nagel. Influence of THM Process Coupling and Constitutive Models on the Simulated Evolution of Deep Salt Formations during Glaciation. The Mechanical Behavior of Salt X. CRC Press, 2022. https://doi.org/10.1201/9781003295808-33.
#
# Li, Shiyuan, and Janos Urai. Numerical Studies of the Deformation of Salt Bodies with Embedded Carbonate Stringers. Online, print. Publikationsserver der RWTH Aachen University, 2012. http://publications.rwth-aachen.de/record/211523/files/4415.pdf
