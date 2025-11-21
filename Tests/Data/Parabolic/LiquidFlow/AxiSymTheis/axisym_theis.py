# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "H process: Theis solution (Pumping well)"
# date = "2022-08-24"
# author = "Wenqing Wang, Olaf Kolditz"
# web_subsection = "liquid-flow"
# +++

# %% [markdown]
# <img class="special-img-class" style="width: 25%; float: left" src="./figures/ogs-jupyter-lab.png" />
# <img class="special-img-class" style="width: 10%; float: left" src="./figures/h-tet-1.png" />
# <div style="clear: both;"></div>

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
from scipy.special import exp1

# %% [markdown]
# ## H process: Theis solution

# %% [markdown]
# **Problem description**
#
# Theis' problem examines the transient lowering of the water table induced by a pumping well.
# The assumptions required by the Theis solution are:
#
# The aquifer:
# - is homogeneous, isotropic, confined, infinite in radial extent,
# - has uniform thickness, horizontal piezometric surface.
#
# The well:
# - is fully penetrating the entire aquifer thickness,
# - has a constant pumping rate,
# - well storage effects can be neglected,
# - no other wells or long term changes in regional water levels.

# %% [markdown]
# **Analytical solution**

# %% [markdown]
# The analytical solution of the drawdown as a function of time and distance is expressed by:
# $$
# s(r,t) = h_0 - h(r,t) = \frac{Q}{4\pi T}W(u), \quad \mathrm{where}\quad u = \frac{r^2S}{4Tt}.
# $$
#
# where
# - $s$ [$L$] is the _drawdown_ or change in hydraulic head,
# - $h_0$ is the constant initial hydraulic head,
# - $h$ is the hydrauic head at distance $r$ at time $t$
# - $Q$ [$L^3T^{-1}$] is the constant pumping (discharge) rate
# - $S$ [$-$] is the aquifer storage coefficient (volume of water released per unit decrease in $h$ per unit area)
# - $T$ [$L^2T^{-1}$] is the transmissivity (a measure of how much water is transported horizontally per unit time).
#
# The _Well Function_, $W(u)$ is the exponential integral, $E_1(u).$
# $W(u)$ is defined by an infinite series:
# $$
# W(u) = - \gamma - \ln u + \sum_{k=1}^\infty \frac{(-1)^{k+1} u^k}{k \cdot k!}
# $$
# where $\gamma=0.577215664$ is the Euler-Mascheroni constant
#
# For practical applications an approximation to the exponential integral is used often:
# $$W(u) \approx -\gamma - \ln u$$
#
# This results in an expression for $s(r,t)$ known as the Jacob equation:
# $$
# s(r,t) = -\frac{Q}{4\pi T}\left(\gamma + \ln u \right).
# $$
# For more details we refer to Srivastava and Guzman-Guzman (1998).
#
# The following analytical solution is inspired by [linear and non linear fitting of the theis equation ](https://scipython.com/blog/linear-and-non-linear-fitting-of-the-theis-equation/).


# %%
def calc_u(r, S, T, t):
    """Calculate and return the dimensionless time parameter, u."""
    return r**2 * S / 4 / T / t


def theis_drawdown(t, S, T, Q, r):
    """Calculate and return the drawdown s(r,t) for parameters S, T.

    This version uses the Theis equation, s(r,t) = Q * W(u) / (4.pi.T),
    where W(u) is the Well function for u = Sr^2 / (4Tt).
    S is the aquifer storage coefficient,
    T is the transmissivity (m2/day),
    r is the distance from the well (m), and
    Q is the pumping rate (m3/day).
    """
    u = calc_u(r, S, T, t)
    return Q / 4 / np.pi / T * exp1(u)


# %% [markdown]
# ## Analytical solution setup
#
# Here, we choose some representative parameter values to evaluate the Theis analytical solution.
#
# We set the pumping rate to $Q=2000\, \mathrm{m^3/day}$ and consider a point located at a distance of $r=10\,\mathrm{m}$ from the well.
# The time grid spans from 1 to 100 days.
#
# For the aquifer properties, we use a storage coefficient of $S=3 \times 10^{-4}$ and a transmissivity of  $T=1000\, \mathrm{m^{2}/day}$
#
# Using these values, we compute the analytical drawdown $s(r,t)$ with the `theis_drawdown()` function.

# %%
Q = 2000  # Pumping rate from well (m3/day)
r = 10  # Distance from well (m)

# Time grid, days.
t = np.array([1, 2, 4, 8, 12, 16, 20, 30, 40, 50, 60, 70, 80, 90, 100])

# Calculate some synthetic data to fit.
S, T = 0.0003, 1000
s = theis_drawdown(t, S, T, Q, r)

# Plot the data
fig, ax = plt.subplots(figsize=(16, 10))

ax.plot(t, s, label=f"r = {r} m")
ax.set(
    xlim=(0, 100),
    ylim=(1.6, 2.6),
    title="Theis: Analytical solution",
    xlabel=r"$t\;/\;\mathrm{days}$",
    ylabel=r"$s\;/\;\mathrm{m}$",
)
ax.legend()
ax.grid()
ot.plot.utils.update_font_sizes(ax)


# %% [markdown]
# In this step, we recalculate the drawdown using SI units for consistency with numerical simulations.
#
# The pumping rate is converted to $Q=0.016\,\mathrm{m^{3}/s}$ and the time is set to $t=864000\,\mathrm{s}$ (equivalent to 10 days).
#
# We evaluate the drawdown over radial distances ranging from $r=1\,\mathrm{m}$ to $40\,\mathrm{m}$.

# %%
# Recalculation from days in sec
Q = 0.016  # Pumping rate from well (m3/s)
t = 864000  # Time in s.

# Distance from well (m)
r = np.arange(1, 41, 1)

# Calculate some synthetic data to fit.
S = 0.001
T = 9.2903e-4
u = calc_u(r, S, T, t)
s = theis_drawdown(t, S, T, Q, r)
s = s - 5  # reference head

# Plot the data
fig, ax = plt.subplots(figsize=(16, 10))
ax.plot(r, s, label=f"t = {t} days")
ax.set(
    xlim=(1, 40),
    title="Theis: Analytical solution",
    xlabel=r"$r\;/\mathrm{m}$",
    ylabel=r"$hydraulic head\;/\;\mathrm{m}$",
)
ax.legend()
ax.grid()
ot.plot.utils.update_font_sizes(ax)


# %% [markdown]
# ## Numerical solution

# %%
mesh = ot.MeshSeries("axisym_theis.vtu")


# %%
fig = ot.plot.contourf(mesh, "OGS5_pressure")
fig.get_axes()[0].set_title("Initial pressure")
ot.plot.utils.update_font_sizes(fig.get_axes())


# %% [markdown]
# ## Running OGS

# %%
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

model = ot.Project(input_file="axisym_theis.prj", output_file="axisym_theis.prj")
model.run_model(logfile=out_dir / "out.log", args=f"-o {out_dir}")


# %% [markdown]
# ## Spatial Profiles

# %%
ms = ot.MeshSeries(f"{out_dir}/liquid_pcs.pvd")

# %%
# %%capture --no-display

xaxis = np.column_stack((np.linspace(1.0, 40, 40), np.zeros((40, 2))))
pressure = ot.variables.pressure.replace(
    data_unit="m", output_unit="m", output_name="hydraulic head", symbol=""
)

ms_probe = ot.MeshSeries.extract_probe(ms, xaxis)

# %%
labels = [f"$t={np.round(x, 2)}s$" for x in ms_probe[1:].timevalues]
fig, ax = plt.subplots(figsize=(18, 12))

ot.plot.line(ms_probe[1:], "x", pressure, labels=labels, ax=ax)
ot.plot.line(
    ms_probe[ms_probe.closest_timestep(1728)],
    "x",
    "OGS5_pressure",
    marker="x",
    markersize=8,
    linestyle="",
    label="OGS, $t=1728 s$",
    ax=ax,
)
fig.tight_layout()

# %% [markdown]
# ### Comparing Numerical and Analytical Solutions
#
# In this section, we compare the computed OGS solutions with the [analytical solution](#analytical-solution-setup) we computed earlier.

# %%
time = 864000
timestep = ms_probe.closest_timestep(time)

ana_sol = np.zeros_like(ms_probe.point_data[pressure])
ana_sol[timestep] = s

mask = np.nonzero(ana_sol)
abs_error = np.zeros_like(ms_probe.point_data[pressure])
abs_error[mask] = np.abs(ana_sol[mask] - ms_probe.point_data[pressure][mask])

rel_error = np.zeros_like(ms_probe.point_data[pressure])
rel_error[mask] = np.abs(abs_error[mask] / (ana_sol[mask]))

np.testing.assert_array_less(abs_error, 0.9)
np.testing.assert_array_less(rel_error, 0.07)

# %%
ms_probe.point_data[pressure.anasol.data_name] = ana_sol
ms_probe.point_data[pressure.abs_error.data_name] = abs_error
ms_probe.point_data[pressure.rel_error.data_name] = rel_error

# %%
fig, ax = plt.subplots(ncols=3, figsize=(40, 12))

ot.plot.line(
    ms_probe[timestep],
    "x",
    pressure,
    marker="x",
    markersize=8,
    linestyle="",
    label="numerical solution (ogs6)",
    ax=ax[0],
)
ot.plot.line(
    ms_probe[timestep],
    "x",
    pressure.anasol,
    marker="o",
    label="analytical solution",
    ax=ax[0],
    color="C1",
)
ot.plot.line(ms_probe[timestep], "x", pressure.abs_error, marker="x", ax=ax[1])
ot.plot.line(ms_probe[timestep], "x", pressure.rel_error, marker="x", ax=ax[2])

fig.suptitle("Theis: Comparison analytical and numerical solution", fontsize=38)
fig.tight_layout()

# %%
max_rel_error = (ms_probe[timestep].point_data[pressure.rel_error.data_name]).max()
print(f"Max relative error: {max_rel_error*100:.2f}%")

# %% [markdown]
# The max relative error observed is nearly 6.35%.
# The observed differences between the analytical and numerical results arise from the distinct boundary conditions applied in the two setups.
#
# Overall, the numerical solution closely matches the analytical results.

# %% [markdown]
# ## OGS links
# - Project file: [axisym_theis.prj](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/AxiSymTheis/axisym_theis.prj)

# %% [markdown]
# ## References
# - Rajesh Srivastava and Amado Guzman-Guzman (1998): Practical Approximations of the Well Function. Groundwater, 36(5): 844-848, [doi.org/10.1111/j.1745-6584.1998.tb02203.x](https://doi.org/10.1111/j.1745-6584.1998.tb02203.x)

# %% [markdown]
# ## Credits
# - Christian for the analytical solution in Python: [linear and non linear fitting of the theis-equation](https://scipython.com/blog/linear-and-non-linear-fitting-of-the-theis-equation/)
# - Wenqing Wang for set-up the OGS benchmark: [liquid flow theis problem](https://www.opengeosys.org/stable/docs/benchmarks/liquid-flow/liquid-flow-theis-problem/)
