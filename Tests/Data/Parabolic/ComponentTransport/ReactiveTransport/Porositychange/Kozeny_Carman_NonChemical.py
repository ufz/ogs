# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python (.venv)
#     language: python
#     name: venv
# ---

# %% [raw]
# +++
# title = "Kozeny-Carman for ComponentTransport"
# date = "2026-02-24"
# author = "Mostafa Mollaali"
# web_subsection = "reactive-transport"
# weight = 3
# +++

# %% [markdown]
# This notebook defines a minimal 1D ComponentTransport problem to verify the Kozeny-Carman constitutive update under controlled hydraulic loading.
#
# The pressure boundary conditions are prescribed as
# $$
# p_{\text{inlet}} = 110000\ \text{Pa}, \qquad
# p_{\text{outlet}} = 100000\ \text{Pa}.
# $$
#
# Permeability is evaluated using the Kozeny-Carman relation
# $$
# k_{\text{KozenyCarman}} = k_0
# \left(\frac{1-\phi_0}{1-\phi}\right)^2
# \left(\frac{\phi}{\phi_0}\right)^3,
# $$
# with initial values $k_0$ and $\phi_0$, and time-dependent porosity
# $$
# \phi(t)=
# \begin{cases}
# 0.10, & t < 0.51\ \text{s},\\
# 0.50, & t \ge 0.51\ \text{s}.
# \end{cases}
# $$
#

# %%
import os
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot

# %%
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## Run OpenGeoSys

# %%
prj_in = "Kozeny_Carman_NonChemical.prj"


model = ot.Project(input_file=prj_in, output_file=out_dir / prj_in)
model.write_input()
model.run_model(
    logfile=out_dir / "ogs.log",
    args=f"-o {out_dir} -m .",
)

# %% [markdown]
# ## Post-processing

# %% [markdown]
# ### Read project constants

# %%
tree = ET.parse(prj_in)
root = tree.getroot()


def get_parameter_value(name: str) -> float:
    for p in root.findall("./parameters/parameter"):
        pname = p.findtext("name")
        if pname == name:
            return float(p.findtext("value"))
    msg = f"Parameter '{name}' not found."
    raise KeyError(msg)


phi0 = get_parameter_value("poro0")
k0 = get_parameter_value("kappa0")

print(f"phi0 = {phi0}")
print(f"k0   = {k0:.3e} m^2")

# %% [markdown]
# ### Load OGS output

# %%
pvd_file = out_dir / "simple1d_kozeny.pvd"
ms = ot.MeshSeries(pvd_file)

times = np.asarray(ms.timevalues, dtype=float)
phi_vals = np.nanmean(ms.cell_data["porosity_avg"], axis=1)
k_ogs = np.nanmean(ms.cell_data["permeability_avg"], axis=1)
k_theory = k0 * ((1.0 - phi0) / (1.0 - phi_vals)) ** 2 * (phi_vals / phi0) ** 3
err_rel = np.abs(k_ogs - k_theory) / np.maximum(np.abs(k_theory), 1e-30)

# %% [markdown]
# ### Observation point over time

# %%
obs_point = ms.mesh(0).cell_centers().points[[0]]
phi_obs = ot.MeshSeries.probe(ms, obs_point, data_name="porosity_avg").values(
    "porosity_avg"
)
k_obs = ot.MeshSeries.probe(ms, obs_point, data_name="permeability_avg").values(
    "permeability_avg"
)

fontsize = 14
fig_obs, (ax_phi, ax_k) = plt.subplots(1, 2, figsize=(16, 5.5))
ax_phi.plot(times, np.ravel(phi_obs), marker="o", ms=5)
ax_k.plot(times, np.ravel(k_obs), marker="s", ms=5)

ax_phi.set_ylabel("porosity / 1")
ax_k.set_ylabel("permeability / m$^2$")
ax_k.set_yscale("log")

for ax in (ax_phi, ax_k):
    ax.set_xlabel("time / s")
    ot.plot.utils.update_font_sizes(ax, fontsize=fontsize)

fig_obs.tight_layout()
fig_obs.subplots_adjust(wspace=0.45)
plt.show()

# %% [markdown]
# ## Compare Kozeny-Carman and OGS permeability
#
# At each time step, the theoretical value $k_{\text{KozenyCarman}}$ is compared against the simulated value $k_{\text{OGS}}$ (from `permeability_avg`) through the relative error
# $$
# e_{\text{rel}} =
# \frac{\left|k_{\text{OGS}} - k_{\text{KozenyCarman}}\right|}
# {k_{\text{KozenyCarman}}}.
# $$
#
# The verification objective is constitutive consistency: $k_{\text{OGS}} \approx k_{\text{KozenyCarman}}$ and $e_{\text{rel}}$ should remain close to zero.

# %%
if not (np.all(np.isfinite(k_ogs)) and np.all(np.isfinite(k_theory))):
    msg = (
        "permeability contains non-finite values; this benchmark's mesh has no "
        "inactive elements, so a NaN/inf here indicates a regression."
    )
    raise RuntimeError(msg)

fig_cmp, ax_cmp = plt.subplots(figsize=(8, 5.5))
ax_cmp.plot(
    times,
    k_theory,
    "k-o",
    lw=2,
    ms=5,
    label="Kozeny-Carman model",
)
ax_cmp.plot(times, k_ogs, "s", ms=5, label="OGS")
ax_cmp.set_xlabel("time / s")
ax_cmp.set_ylabel("permeability / m$^2$")
ax_cmp.set_yscale("log")
ax_cmp.grid(True, which="both", alpha=0.3)
ax_cmp.legend()
ot.plot.utils.update_font_sizes(ax_cmp, fontsize=fontsize)
fig_cmp.tight_layout()
plt.show()


# %%
np.testing.assert_array_less(
    err_rel,
    1e-12,
    err_msg="Kozeny-Carman relative error exceeds tolerance 1e-12.",
)

print("max relative error:", f"{np.max(err_rel):.3e}")
