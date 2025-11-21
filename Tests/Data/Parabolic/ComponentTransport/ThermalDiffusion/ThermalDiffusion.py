# %% [markdown]
# +++
# date = "2019-07-08"
# title = "Tracer diffusion in a thermal gradient"
# weight = 151
# project = ["Parabolic/ComponentTransport/ThermalDiffusion/TemperatureField_transport.prj"]
# author = "Jaime Garibay-Rodriguez, Renchao Lu, Vanessa Montoya, Noor Hasan"
# image = "./figures/thermalDiffusion_ModelSetup.png"
# web_subsection = "reactive-transport"
# +++


# %% [markdown]
# ## Overview
#
# The diffusion of a non-reactive tracer in Opalinus Clay subject to a non-constant temperature field is
# evaluated with OGS-6.
# The temperature dependency over the diffusion process is captured in a two-step model strategy.
# First, the temperature field is obtained by simulating the heating of the domain for a period of one year.
# It is assumed that a quasi-stable temperature field is attained during this period.
# Second, the temperature-dependent diffusion coefficients are calculated with the temperature field of
# step one.
# With this, a new model is set up to obtain the concentration profile of a non-reactive tracer over the
# same domain and period of time.
# Note that, in this case, isotropic diffusion is assumed in this example.
#
# ## Numerical approach and model setup
#
# The benchmark follows the [HT](/docs/processes/thermal-processes/hydro-thermal/) process in a first model to generate
# a temperature field (TemperatureField.prj).
# A 8 x 4 m 2D domain is selected with finer elements on the left side (borehole) for which a Dirichlet
# boundary condition is applied with a temperature of 353.15 K on the upper half (2 m).
# The initial temperature of the porous media is 289.15 K (see Figure 1) and all the other boundaries
# are defined as Neumann (no-flux).
# Opalinus Clay is selected as porous media and full saturation and instantaneous thermal equilibrium
# with its porewater is assumed.
# Properties of the Opalinus clay used in the first-step model are shown in Table 1.
#
# In a second step, a non-reactive tracer is diffused through the domain with a initial concentration
# of $1 * 10^{-8}$ mol/L over 1 m located at the left boundary (see Figure 1).
# For the second-step model, the [Component Transport](/docs/processes/component-transport/hydro-component/)
# process in OGS-6 is used.
#
# **Table 1**: Parameters of the porous media (Opalinus Clay).
#
# | Parameter | Description | Value | Unit |
# | :---: | :---: | :---: | :---: |
# | $\phi$ | Porosity | 0.15 | - |
# | $\rho_{\text {clay }}$ | Dry bulk density | 2720 | $\mathrm{~kg} \mathrm{~m}^{-3}$ |
# | $c_{p\_clay}$  | Specific heat capacity | 800 | $\mathrm{~J} \mathrm{~kg}^{-1} \mathrm{~K}^{-1}$ |
# | $k_{\text {clay }}$ | Thermal conductivity | 0.955 | $\mathrm{~W} \mathrm{~m}{ }^{-1} \mathrm{~K}^{-1}$ |
#
# The $0 \leq x \leq 8$ and $0 \leq y \leq 4$ triangle mesh with finer elements close to the borehole
# lateral
# (left-boundary) is used to i) avoid numerical errors and ii) to gain favorable computational performance.
# An implicit Euler scheme is used for the time discretization with a time horizon of 1 y and fixed time
# steps of 1 d.
#
# ![Domain Mesh](./figures/thermalDiffusion_ModelSetup.png "Figure 1: Representation of the model setup.")
#
# The pore diffusion coefficients (D) used in the model are scalar functions of the temperature in each
# element.
# Arrhenius equation is used as relationship between the effective diffusion coefficient and
# the temperature [1]:
#
# \begin{equation}
# D\left(T_2\right)=D\left(T_1\right) \exp \left[\frac{E_a}{R}\left(\frac{1}{T_1}-\frac{1}{T_2}\right)\right]
# \end{equation}
#
# where $E_a$ = $17000$ J/mol is the activation energy corresponding to bulk water [2] and R is the universal
# gas constant and equal to $8.314$ J/(K mol).
# A value of $D$($298.15$ K) = $2 * 10^{-11}$ $m^2 s^{-1}$ is used for the calculations.
#
# The nodal temperature output is used to compute element-wise diffusion coefficients for the domain.
# This is implemented in ParaView by using a 'Python calculator' filter where the Arrhenius equation
# is used as a scalar function of temperature.
# Afterwards, an additional 'Point data to cell data' ParaView filter is used to carry out the
# element-wise interpolation of the diffusion coefficients.
# The newly generated mesh is then used in the TemperatureField_transport.prj model as input for the
# temperature-dependent molecular diffusion parameters.

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import ogstools as ot
from ogstools.meshlib import MeshSeries

# %%
# creating output path if it doesn't exist already.

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

# %%
# Initiate the OGS

model_temp = ot.Project(
    input_file="TemperatureField.prj", output_file="TemperatureField.prj"
)

# %%
model_temp.run_model(logfile=Path(out_dir) / "log.txt")

# %%
model_conc = ot.Project(
    input_file="TemperatureField_transport.prj",
    output_file="TemperatureField_transport.prj",
)

# %%
model_conc.run_model(logfile=Path(out_dir) / "log2.txt")

# %% [markdown]
# ## Result
#
# Figure 2 shows both, the temperature field (top) and the concentration profile in the Opalinus-Clay
# porewater (bottom) computed after 1 year.
# The results clearly shown that the diffusion of the tracer is confined into the “Borehole-Clay interface”.
# The tracer penetrates at $\approx 0.2$ $m$ in the $x$ direction when measured in the center of the
# diffusion
# range boundary.
# Recall that isotropic diffusion is assumed in this example.

# %%
mesh = MeshSeries("TemperatureField.pvd").mesh(5)
mesh2 = MeshSeries("TemperatureField_transport.pvd").mesh(5)

# %%
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

mesh.plot_contourf("T", fig=fig, ax=ax1, fontsize=15)
ax1.set_ylabel("Borehole heated lateral / m")
cbar1 = fig.axes[-1]
cbar1.set_ylabel("Temperature / K", fontsize=15)

mesh2.plot_contourf("Cs", fig=fig, ax=ax2, fontsize=15, vmin=0)
ax2.set_ylabel("Tracer diffusion boundary / m", fontsize=15)
cbar2 = fig.axes[-1]
cbar2.set_ylabel(r"Tracer concentration / $mol/m^3$", fontsize=15)

fig.suptitle(
    "Figure 2: Temperature (top) and tracer concentration (bottom) profiles at 1 y."
)
plt.tight_layout()
plt.show()

# %%
sample_temp = mesh.sample_over_line([0, 1.75, 0], [0.5, 4, 0])
sample_conc = mesh2.sample_over_line([0, 1.75, 0], [0.5, 3.25, 0])

fig2, ax1 = plt.subplots(figsize=(10, 6))

# First line in blue
ot.plot.line(sample_temp, "T", ax=ax1, color="blue", linewidth=0.5, fontsize=15)
ax1.set_ylabel("Temperature / K", color="blue")
ax1.tick_params(axis="y", labelcolor="blue")

# Second line in red
ax2 = ax1.twinx()
ot.plot.line(sample_conc, "Cs", ax=ax2, color="red", linewidth=0.5, fontsize=15)
ax2.set_ylabel(r"Tracer concentration / $mol/m^3$", color="red")
ax2.tick_params(axis="y", labelcolor="red")

ax2.grid(False, which="both")

plt.tight_layout()
plt.show()

# %% [markdown]
# ## Domain size sensitivity analysis
#
# A sufficient domain size was used in this exercise to neglect the effect of the no-flux boundary
# conditions on the temperature profile.
# However, the top, bottom and right boundaries do have a small effect on the computed temperature
# profile (i.e., temperature at the right boundary is different than the initial temperature).
# Nonetheless, the size of the domain in this benchmark was selected due to its favorable
# computational properties.
#
# For comparison purposes, Figure 3 shows the temperature profiles (along a line in the center of the
# diffusion boundary) for the 8 x 4 m domain used in this description and a bigger domain size
# (10 x 15 m), which is of sufficient dimensions to neglect the Neumann boundaries on the
# temperature profile.
# Using the bigger domain size results in a computational time more than 5 times greater than
# the current benchmark.
# In Figure 3, however, note that the temperature difference between the two 3 cases is minimal
# at the diffusion zone of the tracer (0 - 0.2 m).
#
# ![Different domain](./figures/differentDomain.png "Figure 3: Temperature profiles for different domain sizes at the center of the diffusion boundary.")


# %% [markdown]
# ## Reference
# [1] LR Van Loon, W Müller, and K Iijima. Activation energies of the self-diffusion of hto, 22na+ and
# 36cl- in a highly compacted argillaceous rock (opalinus clay). _Applied Geochemistry_, 20(5):961–
# 972, 2005.
#
# [2] Fátima González Sánchez, Luc R Van Loon, Thomas Gimmi, Andreas Jakob, Martin A Glaus, and
# Larryn W Diamond. Self-diffusion of water and its dependence on temperature and ionic strength
# in highly compacted montmorillonite, illite and kaolinite. _Applied Geochemistry_, 23(12):3840–
# 3851, 2008.

# %%
