# %% [markdown]
# +++
# date = "2017-02-16"
# title = "Richards Flow"
# weight = 141
# project = ["Parabolic/Richards/RichardsFlow_2d_small.prj"]
# author = "Yonghui Huang"
# image = ""
# web_subsection = "richards-flow"
# +++

# %% [markdown]
# ## Introduction
#
# The Richards equation is often used to mathematically describe
# water movement in the unsaturated zone.
# It has been introduced by Richards (1931) who suggested that
# Darcy's law under consideration of the mass conservation principle,
# is also appropriate for unsaturated flow conditions in porous media.
# The pressure based formulation of this governing equation (Eq. 1), which
# elects the unknown primary variable as $p$, can be written as:
#
# \begin{equation}
# \phi \rho_w \frac{\partial S}{\partial p_c} \frac{\partial p_c}{\partial t}+\nabla \cdot\left(\rho_w \frac{k_{r e l} \mathbf{k}}{\mu_w}\left(\nabla p_w-\rho_w \mathbf{g}\right)\right)=Q_w
# \end{equation}
#
# where $\phi$ is porosity, $t$ is time, $\rho_w$ is the liquid density,
# $\mu_w$ is the liquid viscosity, $p_c$ is the capillary pressure with
# $p_c = -p_w$ , $p_w$ is the water pressure, $S$ is the water saturation,
# $g$ is gravity acceleration vector, $Q_w$ is the source term, $k_{rel}$ is
# the relative permeability and $k$ is the intrinsic permeability which is
# related to the hydraulic conductivity $K$ with
#
# \begin{equation}
# \mathbf{k}=\frac{\mu_w}{\rho_w \mathrm{~g}} \mathbf{K}
# \end{equation}
#
# In an unsaturated porous media, the capillary pressure is fundamentally
# related to the saturation of the gas and liquid phase.
# If e.g. the water saturation decreases and hence the saturation of air
# increases, then the water retreats to smaller pores and the capillary
# pressure increases.
# The capillary pressure can be seen as a function of the effective
# saturation $S_{eff}$ .
# This relationship is primarily determined by the nature of the pore
# space geometry and interconnectivity and is highly non-linear.
# Brooks and Corey (1964) and Van Genuchten (1980), among many other
# scientists, derived functional correlations which contain empiric shape
# parameters that characterize pore-specific properties.
# With the Van Genuchten parameterization the capillary pressure
# can be described as
#
# \begin{equation}
# p_c=\frac{\rho_w \mathrm{~g}}{\alpha}\left[S_{\mathrm{eff}}^{-1 / m}-1\right]^{1 / n}
# \end{equation}
#
# where $\alpha$ $[1/m]$ is a conceptualized parameter related to the
# air entry pressure, $n$ is a dimensionless pore size distribution
# index and $m = 1-(1/n)$.
# These parameters are usually used to fit the saturation dependent
# curves of capillary pressure and hydraulic conductivity to
# experimental data.
# The relative permeability can be given as
#
# \begin{equation}
# k_{r e l}=S_{\mathrm{eff}}^{1 / 2}\left[1-\left(1-S_{\mathrm{eff}}^{1 / m}\right)^m\right]^2
# \end{equation}
#
# The effective saturation is
#
# \begin{equation}
# S_{\mathrm{eff}}=\frac{S-S_r}{S_{\mathrm{max}}-S_r}
# \end{equation}
#
# with $S_{max}$ and $S_r$ as the maximum and residual saturation.
#
# ## Infiltration in homogenous soil
# ### Definition
#
# This inﬁltration problem refers to a classical ﬁeld experiment
# described by Warrick et al. (1971), who examined simultaneous
# solute and water transfer in unsaturated soil within the Panoche
# clay loam, an alluvial soil of the Central Valley of California.
# A quadratic 6.10 m plot, which had an average initial saturation
# of 0.455, was wetted for 2.8 h with 0.076 m of 0.2 N $CaCl_2$,
# followed by 14.7 h inﬁltration of 0.229 m solute-free water.
# The soil-water pressure was monitored by duplicate tensiometer
# installations at 0.3, 0.6, 0.9, 1.2, 1.5 and 1.8 m below surface.
# Two ﬁxed pressure boundary conditions are used in the ﬂow equation
# with a uniform initial saturation in the whole domain of 45.5%.
# At the top, the 2 m high soil column is open to the atmosphere,
# i.e. the capillary pressure is 0 Pa.
# The bottom of the column has a capillary pressure of 21,500 Pa.
# Homogeneous material properties are assumed within the whole domain.
# The average saturated moisture content, which is equal to the
# porosity of the soil, is 0.38.
# The saturated permeability is $9.35^{-12}$ $m^2$.
# The relative permeability and capillary pressure vs. saturation
# data are ﬁtted by the soil characteristic functions respectively.
#
# ### Results
#
# The simulated and experimental saturation data at various time
# steps are plotted in Fig.1.
# The OGS6 simulated infiltration front propagates through the soil
# column and resembles well the saturation results of OGS5.
#
# To ensure the consistency of different interpolation functions,
# the soil column has been discretized by two dimensional geometrical
# models, which contain accordingly consistent finite element types
# such as triangles or quadrilaterals.
# Fig. 1 shows the saturation contours after 2, 9 and 17 hours for
# structured meshes.


# %%
# importing libraries
import os
from pathlib import Path

import matplotlib.pyplot as plt
import ogstools as ot

# %%
# Creating output path if it doesn't exist already.
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

# %%
# Initiating OGS object
model = ot.Project(
    input_file="RichardsFlow_2d_small.prj",
    output_file="RichardsFlow_2d_small_modified.prj",
)
# %%
# Modifying convergence criterion
model.replace_text(
    "Residual", xpath="./time_loop/processes/process/convergence_criterion/type"
)

model.replace_text(
    1e-8, xpath="./time_loop/processes/process/convergence_criterion/abstol"
)

# Modifying simulation end time and time stepping
model.replace_text(61200, xpath="./time_loop/processes/process/time_stepping/t_end")

model.remove_element(xpath="./time_loop/processes/process/time_stepping/timesteps/pair")

model.time_loop.add_time_stepping_pair(process="GW23", repeat="612", delta_t="100")

# Modifying output time stepping
model.remove_element(xpath="./time_loop/output/timesteps")

model.add_element(
    parent_xpath="./time_loop/output", tag="fixed_output_times", text="7200 32400 61200"
)

# Write the prj file
model.write_input()

# %%
# Run OGS
model.run_model(logfile=Path(out_dir) / "log.txt")

# %%
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 7), sharey=True)

ms = ot.MeshSeries("richards.pvd").scale(time=("s", "h"))
ms_ogs5 = ot.MeshSeries("h_us_line_Warrick/h_us_line_Warrick_RICHARDS_FLOW.pvd").scale(
    time=("s", "h")
)

# Define shared resolution and colors
res = 200
color = ["C0", "C1", "C2"]
mesh_pairs = [(1, 7), (2, 8), (3, 9)]

for i, (ogs6, ogs5) in enumerate(mesh_pairs):
    line6 = ms.mesh(ogs6).sample_over_line([0.02, 0, 0], [0.02, 2, 0], resolution=res)
    sat6 = line6["saturation"]
    y_coords = line6.points[:, 1]

    line5 = ms_ogs5.mesh(ogs5).sample_over_line([0, 0, 0], [0, 0, 2], resolution=res)
    sat5 = line5["SATURATION1"]

    ax1.plot(sat6, y_coords, c=color[i], label=f"ogs6_{ms.timevalues[ogs6]} h")
    ax1.scatter(
        sat5,
        y_coords,
        c=color[i],
        marker="x",
        label=f"ogs5_{ms_ogs5.timevalues[ogs5]} h",
    )

    difference = sat6 - sat5
    ax2.plot(difference, y_coords, c=color[i], label=f"{ms.timevalues[ogs6]} h")

ax1.set_xlabel("Saturation")
ax1.set_ylabel("y / m")
ax1.set_title("Infiltration (OGS6 - line; OGS5 - scatter)")
ax1.legend(loc="lower center", bbox_to_anchor=(0.5, -0.3), ncol=3)

ax2.set_xlabel(r"Difference ($\Delta S$)")
ax2.set_title("Difference (OGS6 - OGS5)")
ax2.axvline(0, color="black", linestyle="--", linewidth=1, alpha=0.7)  # Zero Reference
ax2.legend(loc="lower center", bbox_to_anchor=(0.5, -0.3), ncol=3)

fig.suptitle("Figure 1: Comparison of OGS6 and OGS5")
plt.tight_layout()
plt.show()
