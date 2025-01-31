# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: 'Python 3.10.8 (''.venv'': venv)'
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "Phase Appearance/Disappearance"
# date = "2022-10-19"
# author = "Norbert Grunwald"
# image = "figures/placeholder_bourgeat.png"
# web_subsection = "th2m"
# coupling = "h2"
# weight = 3
# +++
#

# %% [markdown]
# |<div style="width:300px"><img src="https://www.ufz.de/static/custom/weblayout/DefaultInternetLayout/img/logos/ufz_transparent_de_blue.png" width="300"/></div>|<div style="width:300px"><img src="https://discourse.opengeosys.org/uploads/default/original/1X/a288c27cc8f73e6830ad98b8729637a260ce3490.png" width="300"/></div>|<div style="width:330px"><img src="https://github.com/nagelt/Teaching_Scripts/raw/9d9e29ecca4b04eaf7397938eacbf116d37ddc93/Images/TUBAF_Logo_blau.png" width="300"/></div>|
# |---|---|--:|

# %% [markdown]
# In <cite>Bourgeat et al. (2013)</cite> a numerical experiment is presented which investigates the transport behaviour of hydrogen in clay rock. This test example was investigated by several institutions using different simulation tools and the results were then compared in detail.
#
# To verify the phase transition, the diffusive mass transport and the two-phase behaviour, this numerical experiment was recalculated with the TH2M process class. The material parameters and boundary conditions used essentially correspond to those of the experiment and are summarised as follows:
#
#
# | Parameter                      | Symbol           | Value    | Unit                   |
# |--------------------------------|:----------------:|----------|------------------------|
# | binary diffusion coefficient   | $D$              | 3.0e-9   | Pa                     |
# | viscosity liquid               | $\mu_\text{LR}$  | 1.0e-3   | Pa s                   |
# | viscosity gas                  | $\mu_\text{GR}$  | 9.0e-6   | Pa s                   |
# | Henry-coefficient              | $H$              | 7.65e-6  | mol Pa$^{-1}$m$^{-3}$  |
# | molar mass hydrogen            | $M_\mathrm{H_2}$ | 2.0e-3   | kg mol$^{-1}$          |
# | molar mass water               | $M_\mathrm{H_2O}$| 1.0e-2   | kg mol$^{-1}$          |
# | density liquid                 | $\rho_\text{LR}$ | eq. (1)  | kg m$^{-3}$            |
# | intrinsic permeability         | $\mathbf{k}$     | 5.0e-20  | m$^2$                  |
# | porosity                       | $\phi$           | 0.15     | 1                      |
#
# The van Genuchten model was used as the saturation relation with $m=$0.329 and $p_b$=2.0e6 Pa, $s_\mathrm{L}^\mathrm{res}=$0.4 and $s_\mathrm{G}^\mathrm{res}=$0.0. The relative permeabilities $k^\mathrm{rel}_\mathrm{L}$ and $k^\mathrm{rel}_\mathrm{G}$ were determined using the van Genuchten and van Genuchten-Mualem models, respectively.
#
# Unlike in the cited article, the density of the liquid phase was not considered constant, but was represented by a linear equation of state of the form
#
# $$
# \rho_\text{LR}=\rho^\text{ref}_\text{LR}\left(1+\chi_{c,\text{L}}c^\text{C}_\text{L}\right) \;\;\;(1)
# $$
#
# Here, $\rho^\text{ref}_\text{LR}$=1000 kg$^{-3}$ is a reference density and $\chi_{c,\text{L}}=M_{H_2}{\rho^\text{ref}_\text{LR}}^{-1}$=2.0e-6 m$^{3}$mol$^{-1}$ is the slope of the function with respect to the hydrogen concentration.
#
# This deviation from the specifications in <cite>Bourgeat et al. (2013)</cite> is necessary because the phase transition model in TH2M evaluates the ratio of the densities of solvent and solution to calculate the dissolved gas fraction.
# Therefore, the equation of state must take into account the concentration of dissolved gases, even if the actual influence on the density seems negligible.
#
# ![](figures/Bourgeat_concept.png#one-half)
#
# The experiment considers a 200 m long, porous quasi-1D region, at the left edge of which a hydrogen mass flow is injected. This injection is kept constant for 500,000 years and then stopped.
# The chosen mass flow is low enough that all the injected hydrogen is initially dissolved in the liquid phase and is transported there primarily by diffusion.
#
# The validity of the test is carried out by comparing the variables gas saturation, gas pressure and water pressure with the results of various codes from <cite>Bourgeat et al. (2013)</cite>.
#
# ---
#
# Bourgeat, A. P., Granet, S., & Smaï, F. (2013). Compositional two-phase flow in saturated–unsaturated porous media: benchmarks for phase appearance/disappearance. Simulation of Flow in Porous Media, 1–29. https://doi.org/10.1515/9783110282245.81

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pandas as pd

# %%
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
if not out_dir.exists():
    out_dir.mkdir(parents=True)

# %%
model = ot.Project(input_file="bourgeat.prj", output_file=f"{out_dir}/modified.prj")
# This Jupyter notebook version of this test runs not as far as its cTest-counterpart,
# it'll stop after approx. 800 ka
timestepping = "./time_loop/processes/process/time_stepping"
model.replace_text(2.5e13, xpath=f"{timestepping}/t_end")


# The cTest version shows only a few output-timesteps while this version is supposed to
# output much higher resolution time steps to be able to compare the results. Thus, every
# timestep will be written and the maximum timestep-size of the adaptive time stepping
# method will be reduced drastically
time_end = 1e11
model.replace_text(time_end, xpath=f"{timestepping}/maximum_dt")

# The following for loop generates a text with output times,which is then replaced by
# the project file API (ogstools.Project) in the project file.
new_line = "\u000A"
timesteps = str(0.4 * 1e6 * 86400 * 365.25) + new_line
for t in np.arange(0.6, 1.1, 0.1):
    timesteps += str(t * 1e6 * 86400 * 365.25) + new_line

model.replace_text(1, xpath="./time_loop/output/timesteps/pair/each_steps")
model.replace_text(timesteps, xpath="./time_loop/output/fixed_output_times")
model.replace_text("XDMF", xpath="./time_loop/output/type")
model.write_input()

# Run OGS
model.run_model(logfile=f"{out_dir}/out.txt", args=f"-o {out_dir} -m .")

# %%
# Read the results and compare them with the reference data
ms = ot.MeshSeries(f"{out_dir}/result_bourgeat_domain.xdmf")

saturation = ot.variables.Scalar("saturation", "", "%", symbol="s_{G}")
gas_pressure = ot.variables.Scalar("gas_pressure", "Pa", "MPa", symbol="p_{GR}")
liquid_pressure = ot.variables.Scalar(
    "liquid_pressure_interpolated", "Pa", "MPa", "liquid pressure", "p_{LR}"
)


def plot_results(var: ot.variables.Scalar, ref: str, max_err: float) -> None:
    "Plot evolution at [0, 0, 0], a contourplot at t=10^5 and a timeslice."
    # computes gas saturation from liquid saturation in numerical results
    var_OGS = var.replace(func=lambda x: 1.0 - x) if var == saturation else var

    # === Test for valid results ========================================
    df_refs = pd.read_csv(f"references/bourgeat_{ref}.csv")
    num_vals = var_OGS.transform(ms.probe([0, 0, 0], var.data_name)[:, 0])
    mean_ref_vals = var.transform(df_refs.drop(columns="time").aggregate("mean", 1))
    num_vals_interp = np.interp(df_refs["time"], ms.timevalues("a"), num_vals)
    mean_rel_err = mean_ref_vals - num_vals_interp
    assert np.all(np.abs(mean_rel_err) < max_err)

    # === Line plot =====================================================
    fig1, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(12, 4))
    for ax in [ax1, ax2]:
        for institute in sorted(df_refs.drop(columns="time")):
            ref_vals = var.transform(df_refs[institute])
            ax.plot(df_refs["time"], ref_vals, "-", label=institute)
        ax.plot(ms.timevalues("a"), num_vals, "--k", label="OGS-TH$^2$M")
        ax.set_xlabel("time / a")
    ax1.legend()
    fig1.suptitle(var.output_name + r" vs. time at $x=0$m")
    ax1.set_xscale("log")
    ax1.set_xlim(left=10)
    ax1.set_ylabel(var.get_label())

    # === Contourplot ==================================================
    mesh = ms.mesh(ms.closest_timestep(1e5 * 86400 * 365.25))
    fig2 = ot.plot.contourf(mesh, var_OGS, fontsize=20)
    fig2.suptitle(rf"{var.output_name} at t=$10^5$ a", fontsize=20)

    # === Timeslice=====================================================
    pts = np.linspace([0, 0, 0], [200, 0, 0], num=500)
    fig3 = ms.plot_time_slice(
        var_OGS, pts, time_unit="ka", figsize=[20, 7], interpolate=False, fontsize=20
    )
    fig3.suptitle(f"{var.output_name} over time and x", fontsize=20)
    fig3.tight_layout()


# %% [markdown]
# # Saturation

# %%
plot_results(var=saturation, ref="sG", max_err=0.15)

# %% [markdown]
# # Gas pressure

# %%
plot_results(var=gas_pressure, ref="pGR", max_err=0.35)


# %% [markdown]
# # Liquid pressure

# %%
plot_results(var=liquid_pressure, ref="pLR", max_err=0.1)


# %% [markdown]
# After about 10,000 years, the dissolution capacity of the liquid phase is exhausted and a separate gas phase is formed. From this point on, an increase in water pressure can be observed. This is shown in the liquid pressure plot
# by comparison with the results shown in <cite>Bourgeat et al. (2013)</cite>.
#
# Regarding the accuracy (especially during times t<10,000 years), it should be noted that the reference graphs from the original paper were digitised manually; in the paper, the graphs were only linearly plotted (like those shown here in the right column), so no real statement can be made about the accuracy of the results for times t<10,000 years.
#
# The comparison of the time evolution of the water pressure, gas pressure und gas saturation shows that the numerical solution calculated by OGS6-TH2M is very close to the results shown in <cite>Bourgeat et al. (2013)</cite>.
# Both the time of formation of the gas phase (approximately at t=15,000 a and the shape and magnitude of the resulting pressure rise and pressure drop due to switching off the source term are within the differences of the reference solutions.
#
# It can be concluded that both the transport and the phase transition behaviour of the TH2M process have been successfully derived and implemented. The self-controlled appearance and disappearance of a free gas phase should be emphasised.
