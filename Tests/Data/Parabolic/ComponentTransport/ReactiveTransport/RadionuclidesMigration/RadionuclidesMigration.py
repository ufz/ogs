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
# title = "Radionuclides migration in Opalinus clay"
# date = "2022-05-13"
# author = "Jaime Garibay-Rodriguez, Chaofan Chen, Haibing Shao, Lars Bilke, Olaf Kolditz, Vanessa Montoya, Renchao Lu"
# web_subsection = "reactive-transport"
# +++
#

# %% [markdown]
# # Radionuclides migration in Opalinus clay
#
# The diffusion of Cs-135 and U-238 through a 1-D column of fully water-saturated Opalinus clay for 1 million years is simulated using a mechanistic approach via the coupling of OGS to PHREEQC (OGS#iPHREEQC). Radioactive decay has not been considered due to the long half-lives (Nagra, 2002) of Cs-135 (2.3e6 years) and U-238 (4.5e9 years). The diffusion of radionuclides into the clay porewater is modeled with a Dirichlet boundary condition applied at the left side of the domain, whereas the right-side is modeled as a closed boundary (Neumann). The 1-D domain length is 40 m, which is sufficient to ensure that Cs does not reach the right-side boundary.
#
# The following Opalinus clay properties are used:
#
# |                   Parameter                    | Value  |  Unit   |
# | :--------------------------------------------: | :----: | :-----: |
# |                Porosity $\phi$                 |  0.15  |    -    |
# |   Pore diffusion coefficient $\mathrm{D_p}$    | 1e-10  | m$^2$/s |
# | Effective diffusion coefficient $\mathrm{D_e}$ | 1.5-11 | m$^2$/s |
# |                  Pressure $P$                  |   1    |   atm   |
# |                Temperature $T$                 |   25   |    C    |
#
# The domain is discretized in 200 linear elements (avg. element size = 0.2 m) with a non-uniform mesh refinement at the inlet boundary (to account for the high concentration gradient produced due to the constant concentration). The minimum and maximum element sizes are 0.04 m and 0.35 m, respectively (element length multiplier = 1.00625). Mesh-independence has been verified by running the simulations with a mesh with 400 elements (0.1 m element size). The results did not show a significant difference. The temporal discretization follows the maximum time-step size given by the von Neumann stability critetion for diffusive transport:
#
# $$
# \Delta t < \frac{{\Delta x}^2}{2 D_p},
# $$
#
# therefore, a uniform time-stepping scheme is adopted with $\Delta t = \text{2e8}~s~(\approx 6.34~\text{years})$. This means that the simulation needs a total of 157680 time-steps to estimate the migration of 1 million years. This causes that large CPU times are needed when using OGS#iPHREEQC.
#
# Finally, the total concentrations in Opalinus clay of the main species, including background concentrations of Cs and U, are as follows:
#
# | Species (mol/kgw) |           Value           |
# | :---------------: | :-----------------------: |
# |     HCO$_3^-$     | 4.659 $\times$ 10$^{-4}$  |
# |     Ca$^{2+}$     | 2.651 $\times$ 10$^{-2}$  |
# |      Cl$^-$       | 3.058 $\times$ 10$^{-1}$  |
# |       K$^+$       | 1.631 $\times$ 10$^{-3}$  |
# |     Mg$^{2+}$     | 1.733 $\times$ 10$^{-2}$  |
# |      Na$^+$       | 2.447 $\times$ 10$^{-1}$  |
# |    SO$_4^{2-}$    | 1.417 $\times$ 10$^{-2}$  |
# |     Sr$^{2+}$     | 4.577 $\times$ 10$^{-4}$  |
# |      Cs$^+$       | 1.019 $\times$ 10$^{-10}$ |
# | U$_\text{total}$  | 2.549 $\times$ 10$^{-9}$  |
# |        pH         |            7.6            |
# |        pE         |            -2             |

# %% [markdown]
# ## Cs-135 migration
#
# A Cs concentration at the left boundary of 1.019e-7 mol/kgw is set constant throughout the simulation. The sorption model is taken from Águila et al., 2021. Cs-135 migration is simulated by using OGS#iPHREEQC via the `<surface>` keyword with the following parameters:
#
# | Site name |              Description               | Capacity (mol/kgw) |
# | :-------: | :------------------------------------: | :----------------: |
# |   Su\_    |       Double layer planar sites        |       1.425        |
# |  Su\_ii   | Type II sites of intermediate strength |      1.38e-1       |
# |  Su\_fes  |           Frayed edge sites            |       1.8e-3       |
#
# The corresponding reactions of the sorption model are added to the latest version (2020) of the PSI/Nagra thermodynamic database (`PSINagra2020v1-1_davies.dat`).
#
# In the `Cs.prj` file, the surface input under the `<chemical_system>` keyword is:
# ```XML
# <chemical_system chemical_solver="Phreeqc">
#     ...
#         <surface site_unit="mole">
#             <site>
#                 <name>Su_</name>
#             </site>
#             <site>
#                 <name>Su_fes</name>
#             </site>
#             <site>
#                 <name>Su_ii</name>
#             </site>
#         </surface>
#     ...
# ```
#
# The amounts of each site are given through the `<media>` keyword as follows:
# ```XML
#                 <phase>
#                     <type>Solid</type>
#                     <components>
#                         <component>
#                             <name>Su_</name>
#                             <properties>
#                                 <property>
#                                     <name>molality</name>
#                                     <type>Constant</type>
#                                     <value>1.4250e+00</value>
#                                 </property>
#                             </properties>
#                         </component>
#                         <component>
#                             <name>Su_fes</name>
#                             <properties>
#                                 <property>
#                                     <name>molality</name>
#                                     <type>Constant</type>
#                                     <value>1.8000e-03</value>
#                                 </property>
#                             </properties>
#                         </component>
#                         <component>
#                             <name>Su_ii</name>
#                             <properties>
#                                 <property>
#                                     <name>molality</name>
#                                     <type>Constant</type>
#                                     <value>1.3800e-01</value>
#                                 </property>
#                             </properties>
#                         </component>
#                     </components>
#                 </phase>
# ```
#
# Where the site names should be an exact match in both `<media>` and `<chemical_system>` keywords. This means that it is possible to setup complex simulations comprised of multiple physical domains with distinct chemical systems (e.g., different sorption sites, diffusion, porosity) such as the multiple barriers in a nuclear waste repository.

# %% [markdown]
# ### Solve the numerical model
#
# The simulation takes approximately 23 hrs. to complete when t = 1e6 years. Therefore, only the first 50 time-steps are simulated in this notebook. To run the full simulation, the time loop parameters can be easily adapted.

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)


def plot_results(var: ot.variables.Scalar, x_max: float) -> plt.Figure:
    nucl = var.data_name
    this_ms = ot.MeshSeries(out_dir / f"{nucl}.pvd")
    full_ms = ot.MeshSeries(f"{nucl}_full_simulation/out.pvd")[1:]
    mesh_list = list(this_ms) + list(full_ms)
    timevalues = np.append(this_ms.timevalues, full_ms.timevalues)
    ms = ot.MeshSeries.from_data(mesh_list, timevalues).scale(time=("s", "a"))
    difference = ms[nucl][:-1] <= ms[nucl][1:]
    assert np.all((difference >= 0.0) | np.isclose(difference, 0.0))
    labels = [f"{tv:.1e} years" for tv in ms.timevalues]
    fig = ot.plot.line(ms, var, marker="o", labels=labels, figsize=[8, 6], fontsize=10)
    fig.axes[0].set_xlim((0.0, x_max))
    fig.axes[0].set_yscale("log")
    return fig


# %% Run model
model = ot.Project(input_file="Cs.prj", output_file="Cs.prj")
model.write_input()
model.run_model(logfile=out_dir / "out.txt", args=f"-o {out_dir} -m .")


# %% [markdown]
# ### Plot the results
# Here we plot the last timestep of the above simulation and the timesteps of the full simulation.

# %%
cs_var = ot.variables.Scalar("Cs", output_name="Cs concentration / mol/kgw")
cs_fig = plot_results(cs_var, x_max=5)

# %% [markdown]
# The results show that Cs-135 does not migrate for more than 3 m even after 1 million years.

# %% [markdown]
# ## U-238 migration
#
# Similarly, a U concentration at the left boundary of 1.060e-6 mol/kgw is applied to the domain. In this case, the sorption model is taken from Hennig et al., 2020. U-238 migration is simulated by using OGS#iPHREEQC via the `<surface>` keyword definition with the following parameters:
#
# | Site name |         Description         | Capacity (mol/kgw) |
# | :-------: | :-------------------------: | :----------------: |
# | Ill\_wOH  |      Illite weak site       |      6.254e-1      |
# | Ill\_sOH  |     Illite strong site      |      1.384e-2      |
# | Kln\_aOH  |   Kaolinite Aluminol site   |      7.351e-2      |
# | Kln\_siOH |   Kaolinite Silanol site    |      7.351e-2      |
# | Mll\_wOH  |  Montmorillonite weak site  |      1.137e-1      |
# | Mll\_sOH  | Montmorillonite strong site |      2.516e-3      |
# |   ChlOH   |          Chlorite           |      9.414e-3      |
#
# The Opalinus clay parameters remain the same as in the Cs case, as well as the space-time discretization. In the same way, the reactions of the sorption model are added to the themodynamic database `PSINagra2020v1-1_davies.dat`.

# %% [markdown]
# ### Solve the numerical model
#
# The model is solved for the first 50 time steps only to minimize the CPU time of this notebook. The time loop parameters in `U.prj` can be adapted accordingly to cover the full simulation time of 1 million years.

# %%
model = ot.Project(input_file="U.prj", output_file="U.prj")
model.write_input()
model.run_model(logfile=out_dir / "out.txt", args=f"-o {out_dir} -m .")


# %% [markdown]
# ### Plot the results
#
# Again, the full simulation results are provided in `/U_full_simulation`. They are plotted to show the migration of U after 1 million years.

# %%
U_var = ot.variables.Scalar("U", output_name="U concentration / mol/kgw")
U_fig = plot_results(U_var, x_max=30)


# %% [markdown]
# In this case, it is clear that the migration of U-238 is much more pronounced than Cs-135. After 1 million years, U-238 shows a breakthrough to around 20 m.

# %% [markdown]
# ## References

# %% [markdown]
# Águila, J. F., Montoya, V., Samper, J., Montenegro, L., Kosakowski, G., Krejci, P., & Pfingsten, W. (2021). Modeling cesium migration through Opalinus clay: a benchmark for single-and multi-species sorption-diffusion models. Computational Geosciences, 25(4), 1405-1436.
#
# Garibay-Rodriguez, J., Chen, C., Shao, H., Bilke, L., Kolditz, O., Montoya, V., & Lu, R. (2022). Computational framework for radionuclide migration assessment in clay rocks. Frontiers in Nuclear Engineering (under review).
#
# Hennig, T., Stockmann, M., & Kühn, M. (2020). Simulation of diffusive uranium transport and sorption processes in the Opalinus Clay. Applied Geochemistry, 123, 104777.
#
# Hummel W. Thoenen T. (2021): The PSI Chemical Thermodynamic Database 2020, Nagra NTB 21-03.
#
# <a href="https://web.archive.org/web/20220314123055/https://backend.nagra.ch/system/files/attachments/ntb/ntb2002/e_ntb02-006.pdf">Nagra, 2002.</a> Project Opalinus Clay: Models, Codes and Data for Safety Assessment. Technical Report NTB 02–06. Nagra, Switzerland.

# %% [markdown]
# **Credits:**

# %% [markdown]
# Jaime Garibay-Rodriguez, Chaofan Chen, Haibing Shao, Lars Bilke, Olaf Kolditz,
# Vanessa Montoya, Renchao Lu
