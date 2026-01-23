# %% [raw]
# +++
# title = "CT process: Axisymmetric Theis solution (Pumping well)"
# date = "2025-12-25"
# author = "Philipp Selzer", adapted and corrected from the LiquidFlow-test and -benchmark
# +++

# %%

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
from scipy.special import exp1

# %% [markdown]
# **Problem description**
#
# Theis problem examines the transient lowering of the water table induced by a pumping well.
# The assumptions required by the Theis solution are:
#
# The aquifer:
# - is homogeneous, isotropic, confined, and infinite in radial extent,
# - has a uniform thickness and a horizontal piezometric surface.
#
# The well:
# - is fully penetrating the entire aquifer thickness,
# - has a constant pumping rate,
# - well storage effects can be neglected,
# - no other wells or long-term changes in regional water levels.

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
# - $h$ is the hydraulic head at distance $r$ at time $t$
# - $Q$ [$L^3T^{-1}$] is the constant pumping (discharge) rate
# - $S$ [$L$] is the depth-integrated aquifer storage coefficient (volume of water released per unit decrease in $h$ per unit area)
# - $T$ [$L^2T^{-1}$] is the transmissivity (a measure of how much water is transported horizontally per unit time).
#
# The _Well Function_, $W(u)$ is the exponential integral, $E_1(u).$
# $W(u)$ is defined by an infinite series:
# $$
# W(u) = - \gamma - \ln u + \sum_{k=1}^\infty \frac{(-1)^{k+1} u^k}{k \cdot k!}
# $$
# where $\gamma=0.577215664$ is the Euler-Mascheroni constant
#
# For practical applications, an approximation to the exponential integral is often used:
# $$W(u) \approx -\gamma - \ln u$$
#
# This results in an expression for $s(r,t)$ known as the Jacob equation:
# $$
# s(r,t) = -\frac{Q}{4\pi T}\left(\gamma + \ln u \right).
# $$
# For more details we refer to Srivastava and Guzman-Guzman (1998).
#
# The following analytical solution is inspired by [linear and non linear fitting of the Theis equation ](https://scipython.com/blog/linear-and-non-linear-fitting-of-the-theis-equation/).


# %%
def calc_u(r, S, T, t):
    """Calculate and return the dimensionless time parameter, u."""
    return r**2 * S / 4 / T / t


def theis_drawdown(t, S, T, Q, r):
    """Calculate and return the drawdown s(r,t) for parameters S, T.

    This version uses the Theis equation, s(r,t) = Q * W(u) / (4.pi.T),
    where W(u) is the well function for u = Sr^2 / (4Tt).
    S is the aquifer storage coefficient: S_0 * thickness,
    T is the transmissivity (m2/s): hydraulic conductivity * thickness,
    r is the distance from the well (m), and
    Q is the pumping rate (m3/s).
    """
    u = calc_u(r, S, T, t)
    return Q / 4 / np.pi / T * exp1(u)


# %%
# Recalculation from days to seconds
Q = 0.016  # Pumping rate from well (m3/s)
t = 1728.0  # Time in s.

# Distance from well (m)
r = np.arange(1, 41, 1)

# Calculate some synthetic data to fit.
S = 0.001
T = 9.2903e-4
u = calc_u(r, S, T, t)
s_ana = theis_drawdown(t, S, T, Q, r)  # reference head as h_0 = 0

# Compute all values

time_vals = [8.64, 86.4, 1728.0, 24192.0, 172800.0, 604800.0, 864000.0]
len_time_vals = len(time_vals)

s_all = np.zeros((40, len(time_vals)))

for ii in range(len_time_vals):
    u = calc_u(r, S, T, time_vals[ii])
    s = theis_drawdown(time_vals[ii], S, T, Q, r)
    s_all[:, ii] = s


# %%
mesh = ot.MeshSeries("axisym_theis.vtu")


# %%
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", ""))
out_dir.mkdir(parents=True, exist_ok=True)

model = ot.Project(input_file="axisym_theis_CT.prj", output_file="axisym_theis_CT.prj")


# %%
ms = ot.MeshSeries(f"{out_dir}/CT_hydraulic_head.pvd")

# %%

xaxis = np.column_stack((np.linspace(1.0, 40, 40), np.zeros((40, 2))))
pressure = ot.variables.pressure.replace(
    data_unit="m", output_unit="m", output_name="hydraulic head", symbol=""
)

ms_probe = ot.MeshSeries.extract_probe(ms, xaxis)

ms_trial = ot.MeshSeries(f"{out_dir}/CT_hydraulic_head_ts_30_t_1728.000000.vtu")

# %%
labels = [f"$t={np.round(x, 2)}s$" for x in ms_probe[1:].timevalues]
fig, ax = plt.subplots(figsize=(18, 12))
#
ot.plot.line(ms_probe[1:], "x", pressure, labels=labels, ax=ax)
ax.plot(r, s_all, "+", markersize=16, label=f"t = {t} s")
ax.set(
    xlim=(1, 40),
    xlabel=r"$r\;/\; \mathrm{m}$",
    ylabel=r"$hydraulic head \; /\; \mathrm{m}$",
)


fig.tight_layout()


# %%
time = 1728.0
timestep = ms_probe.closest_timestep(time)

ana_sol = np.zeros_like(ms_probe.point_data[pressure])
ana_sol[timestep] = s_ana

mask = np.nonzero(ana_sol)
abs_error = np.zeros_like(ms_probe.point_data[pressure])
abs_error[mask] = np.abs(ana_sol[mask] - ms_probe.point_data[pressure][mask])

rel_error = np.zeros_like(ms_probe.point_data[pressure])
rel_error[mask] = np.abs(abs_error[mask] / (ana_sol[mask]))

np.testing.assert_array_less(abs_error, 0.2)
np.testing.assert_array_less(rel_error, 0.07)

# %%
ms_probe.point_data[pressure.anasol.data_name] = ana_sol
ms_probe.point_data[pressure.abs_error.data_name] = abs_error
ms_probe.point_data[pressure.rel_error.data_name] = rel_error

# %% [markdown]

# Solutions at later time steps computed by OGS6 deviate from the analytical solution because
# the solution hits the boundary of the domain leading to a steady-state in the numerical solution
# compared to a transient behaviour in the analytical solution.
