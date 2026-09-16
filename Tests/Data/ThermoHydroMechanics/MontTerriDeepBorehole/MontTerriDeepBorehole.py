# %% [raw]
# +++
# title = "Thermo osmosis under geothermal gradient - Mont Terri Deep Borehole experiment"
# date = "2026-09-18"
# author = "Feliks Kiszkurno, Fabien Magri and Thomas Nagel"
# web_subsection = "thermo-hydro-mechanics"
# weight = 3
# +++


# %% [markdown]
# # Introduction
# This benchmark is based on the Mont Terri Deep Borehole experiment (MTDB) and
# the work presented by Gonçalvès et al., 2023, and extended by
# Kiszkurno et al., 2026.
# While presenting the same experiment, this benchmark differs from the paper.
# Changes were made to how the conceptual model was expanded.
# They were made to allow presenting the conclusions in a more concise format
# without distorting them.
# For full details, the reader is referred to the original papers.
#
# ## What is thermo-osmosis?
# Thermo-osmosis (TO) is a fluid flux that occurs in porous media (in this
# benchmark, it is clay rock) in the presence of a thermal gradient (e.g.:
# geothermal gradient, thermal disturbance introduced by a nuclear waste
# container) and under steady-state conditions
# (Soler, 2001; Gonçalvès et al., 2018).
#
# ## Experiment geometry
# The MTDB consists of a single borehole 301 meters long.
# It crosses the Opalinus Clay formation and the surrounding geological layers.
#
# ## Material properties
# For an overview of the material properties, readers are referred to the
# original papers.
# The following table will present only the parameters crucial to
# thermo-osmosis: thermo-osmotic coefficient (given in two forms: permeability
# $\epsilon_T$ and coefficient $k_T$), thermal gradient and permeability ($k$).
#
# | Unit            | Subunit   | $\epsilon_T$ / PaK$^{-1}$ | $k_T$ / m$^{2}$K$^{-1}$s$^{-1}$ | Permeability k / m$^2$ | grad T / Km$^{-1}$ |
# |-----------------|-----------|-------------------------|-----------------------|-----------------------|--------------|
# | Hauptrogenstein | -         | 0                       | 0                     | 1.00e-10              | -0.0556      |
# | Passwang        | -         | 0                       | 0                     | 1.00e-19              | -0.0556      |
# | Opalinus clay   | Sandy I   | 1.31e5                  | 6.93e-13              | 3.74e-21              | -0.0556      |
# | Opalinus clay   | Shaly I   | 1.08e5                  | 7.27e-12              | 4.76e-20              | -0.0556      |
# | Opalinus clay   | Sandy II  | 2.70e5                  | 1.43e-12              | 3.74e-21              | -0.0556      |
# | Opalinus clay   | Carbonate | 0.55e5                  | 7.05e-13              | 9.06e-21              | -0.0556      |
# | Opalinus clay   | Shaly II  | 1.60e5                  | 1.76e-11              | 7.8e-20               | -0.0556      |
# | Staffelegg      | -         | 1.44e5                  | 2.00e-9               | 9.8e-18               | -0.0556      |
#
# ## Boundary conditions
#
# ### Displacement
#
# | Component | Boundary     | Type      | Value |
# |-----------|--------------|-----------|-------|
# | 0         | Whole domain | Initial   | 0     |
# | 1         | Whole domain | Initial   | 0     |
# | 0         | Whole domain | Dirichlet | 0     |
# | 1         | Whole domain | Dirichlet | 0     |
#
#
# ## Pressure
#
# | Boundary        | Type      | Value / Expression        |
# |-----------------|-----------|---------------------------|
# | Whole domain    | Initial   | `-0.006320e6*y + 0.562e6` |
# | Hauptrogenstein | Dirichlet | `-0.006320e6*y + 0.562e6` |
# | Upper boundary  | Dirichlet | `0.562e6 Pa`              |
# | Bottom boundary | Dirichlet | `2.463e6 Pa`              |
#
# ### Temperature
#
# | Boundary        | Type      | Value                 |
# |-----------------|-----------|-----------------------|
# | Whole domain    | Initial   | "Temperature profile" |
# | Whole domain    | Dirichlet | "Temperature profile" |
# | Bottom boundary | Dirichlet | `302.32 K`            |
# | Upper boundary  | Dirichlet | `285.58 K`            |
#
# The "temperature profile" applied as initial and boundary conditions is based
# on the observation points with one value of geothermal gradient per layer.
# The value from observation point boundary between Opalinus carbonate and
# Opalinus shaly II units has been skipped, as including it would result in
# negative gradient.
#
# ## Analytical solution
#
# The analytical solution is based on a system of $n$ equations, where $n$
# indicates the number of model units.
#
# \begin{aligned}
# & 0=v_z+k_{\mathrm{f}}^i \frac{\Delta h_i}{\Delta z_i}+\epsilon_T^i k_{\mathrm{f}}^i \partial_z T \quad \forall i \in[1, n] \\
# & \Delta h=\sum_{i=1}^n \Delta h_i
# \end{aligned}
#
# Thermo-osmosis is represented by the thermo-osmotic permeability ($\epsilon_T$).
# The relation between thermo-osmotic permeability and coefficient ($k_T$)
# is described by following equation:
#
# \begin{aligned}
# k_T=\frac{k_{\mathrm{f}}*\epsilon_{T}}{\mu}
# \end{aligned}
#
# where $k_{\mathrm{f}}$ is the permeability and $\mu$ is water viscosity.
# Further details can be found in the paper (Kiszkurno et al., 2026).
#

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
from _helper_functions import (
    DOMAIN_LENGTH,
    add_thermo_osmosis,
    analytical_solution_pressure,
    analytical_solution_temperature,
    format_figure,
    load_observation_data_pressure,
    load_observation_data_temperature,
    test_against_reference_data,
)
from _mesh import generate_mesh

# %%
# Creating output directory if it doesn't exist already.
# It will be used to store results and auxiliary files required to run
# the simulations.
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# In the following cell, the observation data for pressure and temperature
# is loaded.

# %%
p_ref_p, p_ref_z = load_observation_data_pressure()
t_ref_t, t_ref_z = load_observation_data_temperature()

# %% [markdown]
# Now, observation points are defined.
# They will be used for all model variants presented in this benchmark.

# %%
points_p = np.array([5 * np.ones_like(p_ref_z), p_ref_z, np.zeros_like(p_ref_z)]).T
points_T = np.array(
    [5 * np.ones_like(t_ref_z), t_ref_z - DOMAIN_LENGTH, np.zeros_like(t_ref_z)]
).T

# %% [markdown]
# The analytical solution and numerical model are compared against observational data
# to verify that:
# - OpenGeoSys is capable of modelling the considered system accurately
# - the numerical model is correctly implemented and reflects the investigated
# system
# - the analytical solution has been implemented correctly.
#
# # Model extensions
#
# The original study has utilised one global value to characterise the
# geothermal profile.
# While this approach is generally valid, the thermal gradient is crucial to
# thermo-osmosis; therefore it is modelled and investigated in more detail.
# The geothermal gradient is added to the model with values based on the
# observed data and fitting process (more details can be found in the original
# paper).
# According to data presented by Yu et al., 2017, the permeability in the Passwang
# model unit varies with depth.
# This is in contrast to other units in this experiment.
# To include this characteristic, the Passwang unit is divided into 4
# sub-units designed to be representative.
# Additionally, TO will be enabled in it.
# The following table presents the geometry and parameterisation of the introduced
# division of the Passwang Unit.
#
# | Unit            | Subunit | t / m | Upper boundary - z / m | Lower boundary - z / m | k / m |
# |-----------------|---------|-------|------------------------|------------------------|-------|
# | Hauptrogenstein | -       | 37    | 301                    | 264                    | 1e-10 |
# | Passwang        | I       | 18    | 264                    | 246                    | 1e-20 |
# | Passwang        | II      | 17    | 246                    | 229                    | 1e-21 |
# | Passwang        | III     | 17    | 229                    | 212                    | 1e-20 |
# | Passwang        | IV      | 17    | 212                    | 195                    | 1e-19 |
# | Opalinus clay   | Sandy I | 32    | 195                    | 163                    | 3.74e-21 |
#
# Further modification of the model involves fitting the values of TO permeability.
# The exact parameterisation of TO is uncertain; therefore, the potential
# range of valid values is explored with a fitting process (more details can be
# found in the original paper).
# In the following table, the fitted values of geothermal gradient and
# TO permeability are presented.
#
# | Unit            | Subunit   | $\nabla T$ / Km$^{-1}$ | Fitted $k_T$ / m$^2$K$^{-1}$s$^{-1}$ |
# |-----------------|-----------|---------------------|---------------------------|
# | Hauptrogenstein | -         | -0.0544             | 0                         |
# | Passwang        | I         | -0.0544             | 2.78e-12                  |
# | Passwang        | II        | -0.0544             | 2.78e-13                  |
# | Passwang        | III       | -0.0544             | 2.78e-12                  |
# | Passwang        | IV        | -0.0544             | 2.78e-11                  |
# | Opalinus clay   | Sandy I   | -0.0489             | 6.32e-13                  |
# | Opalinus clay   | Shaly I   | -0.0934             | 3.64e-12                  |
# | Opalinus clay   | Sandy II  | -0.0328             | 7.14e-13                  |
# | Opalinus clay   | Carbonate | -0.0328             | 3.52e-13                  |
# | Opalinus clay   | Shaly II  | -0.0328             | 1.39e-11                  |
# | Staffelegg      | -         | -0.0619             | 2.77e-9                   |
#
# Now, a mesh with subdivided Passwang can be generated, and the simulations
# with fitted $k_T$ can be executed.
#

# %%
mesh_dir = Path(out_dir, "Modified", "Mesh")
mesh_dir.mkdir(parents=True, exist_ok=True)
generate_mesh(mesh_dir)

# %% [markdown]
# With the mesh ready, the numerical simulation can be executed in two variants:
# - pure thermo-hydro-mechanical process (THM)
# - THM process extended with TO.

# %%
# Run THM model and read results
out_dir = Path(out_dir, "Modified")
model_THM_mod = ot.Project(
    input_file="MTDB_THM.prj",
    output_file=out_dir / "MTDB_THM_Mod_modified.prj",
)
model_THM_mod.save(out_dir, overwrite=True)
model_THM_mod.run_model(
    logfile=Path(out_dir) / "log_THM_Modified.txt",
    args=f"-o {out_dir} -m {mesh_dir!s}",
)
ms_THM_mod = ot.MeshSeries(out_dir / "MTDB_THM_Modified.pvd").scale(time="d")

# %%
# Run THM+TO model and read results
model_THM_TO_mod = ot.Project(
    input_file="MTDB_THM.prj",
    output_file=out_dir / "MTDB_THM_TO_Mod_modified.prj",
)
model_THM_TO_mod = add_thermo_osmosis(model_THM_TO_mod)
model_THM_TO_mod.save(out_dir, overwrite=True)
model_THM_TO_mod.run_model(
    logfile=Path(out_dir) / "log_THM_TO_Modified.txt",
    args=f"-o {out_dir} -m {mesh_dir!s}",
)
ms_THM_TO_mod = ot.MeshSeries(out_dir / "MTDB_THM_TO_Modified.pvd").scale(time="d")

# %%
# Probe data - pressure
THM_p_mod = ms_THM_mod.probe(points_p)[-1]["pressure_interpolated"]
THM_TO_p_mod = ms_THM_TO_mod.probe(points_p)[-1]["pressure_interpolated"]

# %%
# Probe data - temperature
THM_T_mod = ms_THM_mod.probe(points_T)[-1]["temperature_interpolated"] - 273.15
THM_TO_T_mod = ms_THM_TO_mod.probe(points_T)[-1]["temperature_interpolated"] - 273.15

# %% [markdown]
# The modifications of the model are also reflected in the analytical solution.

# %%
# Prepare analytical solutions
THM_analytical_p, THM_analytical_z = analytical_solution_pressure(
    base=False, TO=False, gradT=True
)
THM_TO_analytical_p, THM_TO_analytical_z = analytical_solution_pressure(
    base=False, TO=True, gradT=True
)

THM_analytical_T, THM_analytical_T_z = analytical_solution_temperature(gradT=True)

# %% [markdown]
# Overall, the model extensions have led to a better
# match between the numerical model/analytical solution and the observational
# data.
# The observation data (red crosses) presented in the following figure was
# published in the paper by Gonçalvès et al., 2023 (see README.txt note in
# Observation_data folder for more details).

# %%
# Plot comparison between numerical model, analytical solution and data
plt.rcParams.update({"font.size": 16})
fig_mod, axs_mod = plt.subplots(
    1, 3, sharey=True, figsize=(12, 5), width_ratios=[1.25, 1.25, 0.5]
)

### Subplot 0 - Pressure
## Numerical models
axs_mod[0].plot(
    THM_p_mod / 1e6, DOMAIN_LENGTH + p_ref_z, label="NUM: THM", marker="x", c="blue"
)
axs_mod[0].plot(
    THM_TO_p_mod / 1e6,
    DOMAIN_LENGTH + p_ref_z,
    label="NUM: THM+TO",
    marker="x",
    c="black",
)
## Analytical solution
axs_mod[0].plot(
    THM_analytical_p,
    THM_analytical_z,
    label="ANA: THM",
    marker="x",
    c="orange",
    ls="--",
)
axs_mod[0].plot(
    THM_TO_analytical_p,
    THM_TO_analytical_z,
    label="ANA: THM+TO",
    marker="x",
    c="grey",
    ls="--",
)
## Obs data
axs_mod[0].scatter(
    p_ref_p, DOMAIN_LENGTH + p_ref_z, label="Obs data", marker="x", c="red", zorder=99
)

### Subplot 1 - Temperature
## Numerical models
axs_mod[1].plot(THM_T_mod, t_ref_z, label="NUM: THM", marker="x", c="blue")
axs_mod[1].plot(THM_TO_T_mod, t_ref_z, label="NUM: THM+TO", marker="x", c="black")
# ## Analytical
axs_mod[1].plot(
    THM_analytical_T,
    THM_analytical_T_z,
    label="ANA: THM",
    marker="x",
    c="orange",
    ls="--",
)
## Observation data
axs_mod[1].scatter(t_ref_t, t_ref_z, label="Obs data", marker="x", c="red", zorder=99)

format_figure(axs_mod)

# %% [markdown]
# # Summary
# This benchmark shows an excellent agreement between the analytical model and
# numerical simulations performed with OpenGeoSys, proving its suitability
# for researching thermo-osmosis.
# An overall improvement in agreement between the model and observation data has
# been achieved by the introduction of a per-unit geothermal gradient, including
# permeability gradient and sub-units in the Passwang model unit and exploring
# parameter uncertainty using fitting methods.


# %% [markdown]
# # Bibliography:
# Gonçalvès, J., Matray, J.-M., & Yu, C. J. (2023). Assessing relevant transport processes in Opalinus Clay at the Mont Terri rock laboratory using excess-pressure, concentration and temperature profiles. Applied Clay Science, 242, 107016. https://doi.org/10.1016/j.clay.2023.107016
#
# Kiszkurno, F. K., Magri, F., & Nagel, T. (2026). Learning from data—Calibration and improvement of modeling thermo-osmosis effects in THM simulations based on the Mont Terri Deep Borehole experiment. International Journal of Rock Mechanics and Mining Sciences, 202, 106513. https://doi.org/10.1016/j.ijrmms.2026.106513
#
# Soler, J. M. (2001). The effect of coupled transport phenomena in the Opalinus Clay and implications for radionuclide transport. Journal of Contaminant Hydrology, 53(1), 63–84. https://doi.org/10.1016/S0169-7722(01)00140-1
#
# Gonçalvès, J., Ji Yu, C., Matray, J.-M., & Tremosa, J. (2018). Analytical Expressions for Thermo-Osmotic Permeability of Clays. Geophysical Research Letters, 45(2), 691–698. https://doi.org/10.1002/2017GL075904
#
# Yu, C., Matray, J.-M., Gonçalvès, J., Jaeggi, D., Gräsle, W., Wieczorek, K., Vogt, T., & Sykes, E. (2017). Comparative study of methods to estimate hydraulic parameters in the hydraulically undisturbed Opalinus Clay (Switzerland). Swiss Journal of Geosciences, 110(1), 85–104. https://doi.org/10.1007/s00015-016-0257-9


# %% [markdown]
# # Acknowledgments
# We gratefully acknowledge Prof. Julio Gonçalves for the data from the Mont Terri Deep Borehole experiment, which was crucial for this work.
#
# The ThORN project was supported by the German Federal Office for the Safety of Nuclear Waste Management (BASE) through the research grant “Experimental investigations on thermo-osmotic flow in argillaceous materials relevant to deep geological repositories for radioactive waste” (Grant No. 4723F00104) and by the French National Agency for Radioactive Waste Management (Andra).


# %% [markdown]
# # Testing
#
# Check pressure and temperature against reference results.

# %%
test_against_reference_data(THM_p_mod, THM_T_mod, THM_TO_p_mod, THM_TO_T_mod)
