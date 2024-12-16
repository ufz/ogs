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
# author = "Boyan Meng and Yonghui Huang"
# date = "2022-07-01"
# title = "Heat pipe problem"
# web_subsection = "thermal-two-phase-flow"
# +++
#

# %% [markdown]
# ## Introduction
#
# When an unsaturated porous medium is subject to a constant heat flux and the temperature is sufficiently high, water is heated and vaporizes. Vapor flows under its pressure gradient towards the cooler end where it condenses. Vaporization and condensation produce a liquid saturation gradient, creating a capillary pressure gradient inside the porous medium. Condensate flows towards the hot end under the influence of a capillary pressure gradient. This is a heat pipe in an unsaturated porous medium.
#
# A benchmark regarding the heat pipe problem was proposed by Udell and Fitch (1985). A semi-analytical solution was provided for a non-isothermal water-gas system in a porous medium, in which heat convection, conduction, and latent heat transport as well as capillary effects play a key role.

# %% [markdown]
# ## Physical Scenario
#
# As shown in the below figure, the heat pipe was represented by a 2D horizontal column (2.4 m in length and 0.2 m in width) of porous media, which was partially saturated with a liquid phase ($S_w$ = 0.41) at the beginning. A heater is installed on the right-hand-side of the horizontal column generating a constant heat flux of 100 $\mathrm{W/m^2}$ and raises the temperature gradually above the boiling point. At the left-hand boundary, we impose the constant gas phase pressure ($P_g$ = 101934 Pa), constant liquid saturation ($S_w$ = 0.97) and constant temperature ($T$ = 343 K) as Dirichlet boundary conditions.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import vtuIO
from IPython.display import Image, display
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# %%
display(Image(filename="./model_domain.jpg", width=1000))

# %% [markdown]
# ## Model parameters and numerical settings
#
# In this benchmark, the thermal conductivity for an unsaturated medium is given as:
# \begin{equation}
#     \lambda(S_w)=\lambda_{S_w=0}+\sqrt {S_w}\left(\lambda_{S_w = 1}-\lambda_{S_w = 0})\right.
# \end{equation}
# The capillary pressure is dependent on the liquid saturation via the Leverett (Leverett et al. (1941)) function:
# \begin{equation}
#         P_c(S_w)=\sqrt{\frac{\phi}{K}}\gamma\left(1.417(1-S_w)-2.12(1-S_w)^2+1.263(1-Sw)^3)\right.
# \end{equation}
# where $\gamma$ = 0.05878 N/m stands for the surface tension of water. The relative permeabilities are calculated using the Udell (Udell and Fitch (1985)) model:
# \begin{equation}
#     k_{rL}=S_w^3,
# \end{equation}
# \begin{equation}
#     k_{rG}=(1-S_w)^3.
# \end{equation}
# The rest of the parameters used in this benchmark are listed in the following table.
#
# | Parameter | Value | Unit |
# | :-: | :-: | :-: |
# | Intrinsic permeability $K$ | 1e-12 | m$^2$ |
# | Porosity                                        $\phi$          | 0.4                 |   -                     |
# | Thermal conductivity of dry porous medium       $\lambda_{S_w=0}$ | 0.582               | W/m/K   |
# | Thermal conductivity of saturated porous medium $\lambda_{S_w=1}$ | 1.14                | W/m/K  |
# | Specific heat capacity of soil grain            $c_{p,s}$       | 700                 | J/kg/K   |
# | Specific heat capacity of air                   $c_{v,a}$       | 733                 | J/kg/K   |
# | Specific heat capacity of water                 $c_{p,w}$       | 4187                | J/kg/K   |
# | Density of water                                $\rho_w$        |1000                | kg/m$^3$          |
# | Density of soil grain                           $\rho_s$        | 2650                | kg/m$^3$          |
# | Dynamic viscosity of the liquid phase           $\mu_{L}$       | 2.938e-4  | Pa s                |
# | Dynamic viscosity of the gas phase              $\mu_{G}$       | 1.8e-5   | Pa s                |
# | Diffusion coefficient in free gas               $D_{0a}$        | 2.6e-6  | m$^2$/s         |
# | Diffusion coefficient in free water             $D_{0w}$        | 3.0e-9   | m$^2$/s           |
# | Latent heat of water vaporization               $h_{\Delta e}$  | 2.258e6   | J/kg          |
#

# %% [markdown]
# ## Results and analysis
#
# In the CTEST-small, the comparison is made for the time of 10000 seconds. The profiles of saturation and temperature are plotted below.

# %%
plt.rcParams["legend.fontsize"] = 20
plt.rcParams["font.size"] = 20


# %%
filename = "ref_t_10000.000000.vtu"
x = np.array([i * 0.012 for i in range(201)])
r = np.array([[i, 0.1, 0.0] for i in x])

f = vtuIO.VTUIO(filename, nneighbors=100, dim=2)
resp = {}
resp[0] = f.get_set_data("saturation", pointsetarray=r)
resp[1] = f.get_set_data("temperature", pointsetarray=r)

fig, ax = plt.subplots(ncols=2, figsize=(20, 8))
for i in range(2):
    ax[i].plot(x, resp[i], lw=2, label="OGS, $t$ = 10000s")
    ax[i].set_xlim([0, 2.4])
    ax[i].set_xlabel("$x$ / m")
    ax[i].legend()
ax[0].set_ylabel("$S_w$ / -")
ax[1].set_ylabel("$T$ / K")
ax[0].set_title("saturation")
ax[1].set_title("temperature")
fig.tight_layout()


# %% [markdown]
# In the CTEST-large, the comparison is made for the time of 1.4e6 seconds. Around this time, the water is fully evaporated from the heating boundary (right hand side), and single phase zone of gas phase is formulated, while the temperature at this part begins to increase significantly, as shown below.

# %%
filename = "ref_t_1400000.000000.vtu"

f = vtuIO.VTUIO(filename, nneighbors=100, dim=2)
resp = {}
resp[0] = f.get_set_data("saturation", pointsetarray=r)
resp[1] = f.get_set_data("temperature", pointsetarray=r)

fig, ax = plt.subplots(ncols=2, figsize=(20, 8))
for i in range(2):
    ax[i].plot(x, resp[i], lw=2, label="OGS, $t$ = 1.4e6s")
    ax[i].set_xlim([0, 2.4])
    ax[i].set_xlabel("$x$ / m")
    ax[i].legend(fontsize=20)
ax[0].set_ylabel("$S_w$ / -")
ax[1].set_ylabel("$T$ / K")
ax[0].set_ylim([0, 1])
ax[0].set_title("saturation")
ax[1].set_title("temperature")
fig.tight_layout()


# %% [markdown]
# After the gas phase appearance, it is recommended to change to an adaptive time stepping scheme (e.g. Evolutionary PID Controller or Iteration Number Based) to assure the numerical stability. In the case of Iteration Number Based Time Stepping, the time step size is kept around 175 s with 4.5 iterations on average.
#
# For the steady-state solution of this problem, a semi-analytical solution was derived by Udell and Fitch (1985) and extended by Huang et al. (2015). Here we provide the semi-analytical solution as a <a href="https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/ThermalTwoPhaseFlowPP/HeatPipe/SemianalyticalSolutionHeatPipe.zip">MATLAB script</a> which enables us to compute the steady-state gas pressure, saturation and temperature profiles along the $x$-direction (see calculated values in <a href="https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/ThermalTwoPhaseFlowPP/HeatPipe/SemianalyticalSolutionResults.csv">SemianalyticalSolutionResults.csv</a>). In the following, the numerical solution by OpenGeoSys at quasi-steady state ($t$ = 2e7 s) is plotted against the semi-analytical solution for comparison. In addition, the absolute and relative errors are also illustrated.

# %%
result_file = "SemianalyticalSolutionResults.csv"
soln = pd.read_csv(
    result_file,
    sep=",",
    header=None,
    skiprows=0,
    names=["x", "saturation", "temperature", "pressure"],
    index_col=False,
)

filename = "ref_steady_status.vtu"

f = vtuIO.VTUIO(filename, nneighbors=100, dim=2)
resp = {}
resp[0] = f.get_set_data("saturation", pointsetarray=r)
resp[1] = f.get_set_data("temperature", pointsetarray=r)
resp[2] = f.get_set_data("gas_pressure", pointsetarray=r)


# %%
plt.rcParams["legend.fontsize"] = 14
plt.rcParams["font.size"] = 14
fig, ax = plt.subplots(ncols=3, figsize=(18, 6))

ax[0].plot(x, soln["saturation"], lw=1.5, label="semianalytical")
ax[0].plot(
    x,
    resp[0],
    lw=1.5,
    marker="o",
    linestyle="",
    markevery=5,
    color="r",
    label="OGS steady state",
)
ax[1].plot(x, soln["saturation"] - resp[0], lw=1.5)
ax[2].plot(x, (soln["saturation"] - resp[0]) / soln["saturation"], lw=1.5)

for i in range(3):
    ax[i].set_xlim([0, 2.4])
    ax[i].set_xlabel("$x$ / m")
ax[0].set_ylabel("$S_w$ / -")
ax[1].set_ylabel(r"$\Delta S_w$ / -")
ax[2].set_ylabel(r"$\Delta S_w/S_{w, analytical}$")
ax[0].set_ylim([0, 1])
ax[0].set_title("Saturation")
ax[1].set_title("Absolute error")
ax[2].set_title("Relative error")
ax[0].legend()
fig.tight_layout()

ax2 = inset_axes(
    ax[0],
    width="100%",
    height="100%",
    loc="center",
    bbox_to_anchor=[0.45, 0.4, 0.5, 0.4],
    bbox_transform=ax[0].transAxes,
)

ax2.plot(
    x,
    resp[0],
    lw=1.5,
    marker="o",
    linestyle="",
    markevery=1,
    color="r",
    label="OGS steady state",
)
ax2.plot(x, soln["saturation"], lw=1.5, label="semianalytical")
ax2.set_xlim(1.57, 1.63)
ax2.set_ylim(0, 0.1)
ax2.set_yticks(np.arange(0, 0.15, 0.05))

patch, pp1, pp2 = mark_inset(ax[0], ax2, loc1=3, loc2=4, fc="none", ec="0.5")
pp1.loc2 = 2
pp2.loc2 = 1

# %%
fig, ax = plt.subplots(ncols=3, figsize=(18, 6))

ax[0].plot(x, soln["temperature"], lw=1.5, label="semianalytical")
ax[0].plot(
    x,
    resp[1],
    lw=1.5,
    marker="o",
    linestyle="",
    markevery=5,
    color="r",
    label="OGS steady state",
)
ax[1].plot(x, soln["temperature"] - resp[1], lw=1.5)
ax[2].plot(x, (soln["temperature"] - resp[1]) / soln["temperature"], lw=1.5)

for i in range(3):
    ax[i].set_xlim([0, 2.4])
    ax[i].set_xlabel("$x$ / m")
ax[0].set_ylabel("$T$ / K")
ax[1].set_ylabel(r"$\Delta T$ / K")
ax[2].set_ylabel(r"$\Delta T/T_{analytical}$")
ax[0].set_title("Temperature")
ax[1].set_title("Absolute error")
ax[2].set_title("Relative error")
ax[0].legend()
fig.tight_layout()


# %%
fig, ax = plt.subplots(ncols=3, figsize=(18, 6))

ax[0].plot(x, soln["pressure"], lw=1.5, label="semianalytical")
ax[0].plot(
    x,
    resp[2],
    lw=1.5,
    marker="o",
    linestyle="",
    markevery=5,
    color="r",
    label="OGS steady state",
)
ax[1].plot(x, soln["pressure"] - resp[2], lw=1.5)
ax[2].plot(x, (soln["pressure"] - resp[2]) / soln["pressure"], lw=1.5)

for i in range(3):
    ax[i].set_xlim([0, 2.4])
    ax[i].set_xlabel("$x$ / m")
ax[0].set_ylabel("$P_g$ / Pa")
ax[1].set_ylabel(r"$\Delta P_g$ / Pa")
ax[2].set_ylabel(r"$\Delta P_g/P_{g, analytical}$")
ax[0].set_title("Gas pressure")
ax[1].set_title("Absolute error")
ax[2].set_title("Relative error")
ax[0].legend()
fig.tight_layout()


# %% [markdown]
# From the above results, it can be seen that a very good agreement is obtained with respect to the variables saturation, temperature and gas pressure, especially for the latter two. For the saturation, there is only one data point that the divergence between numerical and semi-analytical solutions is obvious, which situates at the end of the two-phase zone (see the embedded subplot above). This might be due to the sharp saturation change around this point which necessitates further mesh refinement locally. Nevertheless, the extent of the heat pipe region at steady state was modeled accurately. The disappearance of the water phase associated with a change of the phase state was carried out well. Note that the OGS solution allows a region near the heated boundary to completely dry out, thus creating increased temperatures (superheated steam) in comparison to the semi-analytical results which assumes coexistence of the liquid and gas phases.
#
# ## References
#
# [1] K. Udell and J. Fitch. Heat and mass transfer in capillary porous media considering evaporation, condensation, and non-condensible
# gas effects. 23rd ASME/AIChE national heat transfer conference, Denver, CO. 1985, pp. 103-110.
#
# [2] Leverett M et al. (1941) Capillary behavior in porous solids. Trans AIME 142(01):152-169
#
# [3] Y. Huang, O. Kolditz, and H. Shao. Extending the persistent primary variable algorithm to simulate non-isothermal two-phase two-component flow with phase change phenomena. Geothermal Energy 3 (1) (2015). http://dx.doi.org/10.1186/s40517-015-0030-8.

# %%
