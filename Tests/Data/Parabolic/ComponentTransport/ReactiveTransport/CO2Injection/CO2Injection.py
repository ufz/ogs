# %% [markdown]
# +++
# title = "OGS-PHREEQC-Benchmark: CO2 injection into Opalinus Clay"
# date = "2023-09-07"
# author = "Shuang Chen, Vinay Kumar, Jobst Maßmann"
# web_subsection = "reactive-transport"
# image = "figures/results_comparision_sce1.png"
# +++

# %% [markdown]
# $$\require{mhchem}$$
#
# ## Overview
#
# The present numerical work is motivated by the CO2LPIE project (shortened as CL-experiment) [1], which is an in-situ experiment that is being conducted at the Mont Terri
# rock laboratory.
# In the experiment, carbon dioxide ($\ce{CO2}$) is injected into the Opalinus Clay leading to changes in its hydraulic, mechanical and chemical properties.
# In general, these processes are of great interest in the evaluation of barrier integrity.
#
# Two scenarios are considered in this benchmark.
# In the first scenario, the pure $\ce{CO2}$ gas induced Calcite dissolution is simulated by OGS-6-IPhreeqc and the results are verified with the related PHREEQC example presented by Appelo and Postma [2].
# A comprehensive information regarding the computational procedure of OGS-6-IPhreeqc can be found in Lu et al. [3].
# In the second scenario, the simulation considers a more accurate representation of chemical environments based on the CL Project.
# This includes the incorporation of primary minerals typically found in Opalinus Clay and the relevant composition of species present in the pore water.
#
# ## Scenario #1: Kinetics of $\ce{CO2}$ induced calcite dissolution
#
# ### Problem description
#
# In this case, carbon dioxide ($\ce{CO2}$) gas with a partial pressure $\mathrm{10^{-1.5}}$ $\mathrm{bar}$ is injected into a fluid containing calcite.
# The temperature is maintained at 10 $\mathrm{°C}$, the fluid`s pH value is set at 6 and the initial concentration of calcite is 1 $\mathrm{mol\cdot kg^{-1}\cdot water}$ (shortened as $\mathrm{mol/kgw}$).
# The dissolution of the $\ce{CO2}$ gas in water together with the calcite dissolution pathways can be described with the following reactions (Plummer et al. [4], Appelo and Postma [2]).
#
# $$ \ce{CO2_{(g)} + H2O -> H2CO3^{\ast}} $$
# $$ \ce{CaCO3 + H+ <=> Ca^2+ + HCO3-} $$
# $$ \ce{CaCO3 + H2CO3^{\ast} <=> Ca^2+ + 2HCO3-} $$
# $$ \ce{CaCO3 + H2O <=> Ca^2+ + HCO3- + OH-} $$
#
# Two kinetic rates are adopted to describe the calcite dissolution.
# The first approach involves a simplified rate calculation derived directly from the current concentration of calcium,
# employing the following formula: $\text{rate} = 10^{-6.91} - 10^{-1.52} * (c_{\ce{Ca}})^2$, with $c_{\ce{Ca}}$ in $\mathrm{mol/kgw}$.
# And another one is the well known PWP rate, which is proposed by Plummer, Wigley and Parkhurst (denoted further as “PWP”) in 1978 for calcite dissolution based on three dissolution reactions [4].
# The formulae of the PWP approach is available to be found in the PHREEQC database e.g. `phreeqc.dat`.
# In the numerical experiment, the total simulation time is 30 000 $\mathrm{s}$.
# Details about the case study is described in the example 5.9 from Appelo and Postma [2].
# The related PHREEQC script is available online to be found [5].
#
# ### Model and results
#
# A simple 1D line-element model with one element is constructed in this work.
# The coupled Hydraulic-Chemical (HC) Process is adopted for the OGS-6-IPhreeqc simulation.
# To match the PHREEQC example, the advection and diffusion have been set to zero in the modelling.
# The initial conditions and chemical parameters are provided as outlined in the associated example.
# The Fig. 1 depicts the comparison between the computed results obtained from the OGS-6-IPhreeqc model and the results derived from PHREEQC, which shows a very
# good agreement with each other.
# Preparation of the simulation and results evaluation are presented in what follows.
#
# ### 1) Solve the numerical model

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pandas as pd

# %%
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

_ = os.system(
    f"ogs calcite_simple.prj -o {out_dir} 2>&1 > {out_dir / 'out_simple.txt'}"
)
_ = os.system(f"ogs calcite_pwp.prj -o {out_dir} 2>&1 > {out_dir / 'out_pwp.txt'}")

# %% [markdown]
### Comparsion of OGS-6 simulation and PHREEQC results


# %%
ot.variables.u_reg.define("kgw = kg")
ca = ot.variables.Scalar("Ca", "mol/kgw", "mmol/kgw", "Ca+")
h = ot.variables.Scalar("H", output_name="H+")
ph = ot.variables.Scalar("H", output_name="pH", func=lambda x: -np.log10(x))

ogs_results: dict[str, ot.MeshSeries] = {}
phreeqc_results: dict[str, pd.DataFrame] = {}
for kind in ["simple", "pwp"]:
    ogs_results[kind] = ot.MeshSeries(f"{out_dir}/calcite_{kind}.pvd").extract(0)
    phreeqc_results[kind] = pd.read_csv(
        f"./PHREEQC_results_{kind}.txt", sep=r"\s+", header=0, skipinitialspace=True
    )

fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
axs[0].plot([], [], "-k", label="OGS6")
axs[0].plot([], [], "xk", label="PHREEQC")
for i, kind in enumerate(["simple", "pwp"]):
    ms_pt = ogs_results[kind]
    phreeqc = phreeqc_results[kind]
    c = f"C{i}"
    for j, var in enumerate([ca, ph]):
        ot.plot.line(ms_pt, "time", var, ax=axs[j], fontsize=10, c=c, label=kind)
        axs[j].scatter(
            phreeqc["x"] * 1e3, phreeqc[var.output_name[:2]], marker="x", c=c
        )
    assert np.all(ms_pt.timevalues >= 0.0)
    assert np.all(ca.transform(ms_pt) >= 0.0)
axs[0].set(xlim=(0, None), ylim=(0, None))
_ = axs[1].legend(*fig.axes[0].get_legend_handles_labels())

# %% [markdown]
# Fig. 1: pH and calcium increase with kinetic dissolution of Calcite.

## L2 difference

# %%
phreeqc_results["pwp"]["H"] = 10 ** (-phreeqc_results["pwp"]["pH"])
t = phreeqc_results["pwp"]["x"][1:] * 1e3

fig, ax = plt.subplots(figsize=(6, 4), sharex=True)
for color, var in zip(["g", "k"], [ca, h], strict=True):
    OGS_data = var.transform(ogs_results["pwp"][24:])[:, 0]
    phreeqc_data = phreeqc_results["pwp"][var.output_name[:-1]][1:]
    error = np.log10(((OGS_data - phreeqc_data) ** 2) ** 0.5 / phreeqc_data)
    np.testing.assert_array_less(error, -1.8)
    ax.plot(t, error, c=color, marker="o", markersize=4, label=var.output_name)

    assert np.all((t >= 1500) & (t <= 3.0e4))
    assert np.all((error >= -5) & (error <= 0.0))

ax.set(xlim=(1500.0, 3.0e4), ylim=(-5, 0.0), xticks=np.arange(0, 30001, 5000))
ax.set_xlabel("time / s")
ax.set_ylabel(
    r"Log($\frac{|| \mathbf{c}^\mathrm{OGS-6} - \mathbf{c}^{\mathrm{PHREEQC}}||_{2}}{\mathbf{c}^{\mathrm{PHREEQC}}}$)"
)
ax.grid(color="gray", linestyle="dashed")
_ = ax.legend(loc="best", fontsize=11)

# %% [markdown]
# Fig. 2:  L2 relative difference norm between the obtained results from OGS-6 and PHREEQC for the scenario with the PWP rate.
#
# ## Scenario #2: Modelling of the $\ce{CO2}$ injection into Opalinus Clay
#
# ### Model description
#
# The identical 1D line-element model as in scenario #1 is used in the OGS-6-IPhreeqc simulation.
# In the model, the initial porosity is set to 0.15. Similarly to the scenario #1, hydraulic advection and diffusion are excluded from the simulation.
# In this case, chemical environments based on the CL-experiment are considered.
# For the porous medium in the OGS-6 model, in terms of the solid component, the following mineral composition is assumed in this work: 36% illite, 24% kaolinite, 7.5% calcite and 2.5% dolomite-dis (namely sedimentary (disordered) dolomite).
# The remaining 30% of the mineral is mostly quartz and is not considered in this simulation scenario.
# A more detailed description of the mineral composition of Opalinus Clay can be found in the work of Thury [6].
# In OGS-6-IPhreeqc, the molar amount of each reactive solid component per kilogram of water ${b_m}$ [$\mathrm{mol/kgw}$] is calculated by
#
# $$ {b_m} = \frac{\phi_m}{\rho^{l}\phi V_m},$$
#
# where $\mathrm{\phi_m}$ is the volume fraction of the corresponding mineral ${m}$.
# The $\mathrm{\rho^{l}}$ is the density of the fluid in $\mathrm{kg/m^3}$.
# The $\mathrm{\phi}$ denotes the porosity of the porous medium.
# And the ${V_m}$ is the molar volume of the corresponding mineral ${m}$ in $\mathrm{m^3/mol}$.
# Details of the equation description can be found in Lu et al. [3].
# In the simulation, the used parameters of the Opalinus Clay minerals are listed in table 1.
#
# During the simulation, a constant carbon dioxide ($\ce{CO2}$) gas partial pressure of $10^{1.5}$ $\mathrm{bar}$ is applied to the model domain, same as in the equilibrium phase.
# Notice in this scenario, the $\ce{CO2}$ partial pressure is much higher than it in scenario #1.
# Consequently, the dissolved $\ce{CO2}$ results in an alteration in the acidity of the pore water and leads to the different chemical reactions of each mineral.
# The related reaction formulae are described in the PHREEQC database `llnl.dat`.
# The adopted kinetic dissolution rate of each mineral are referenced from the transition state theory-based reaction mechanism following the work by Palandri and Kharaka [7]. The general equation formula reads
#
# $${{\text{Rate}_\mathrm{mineral}}} = [ k_{\mathrm{acid_{mineral}}}^\mathrm{298.15K}{a_\mathrm{H^{+}}} + k_{\mathrm{neutral_{mineral}}}^\mathrm{298.15K} + k_{\mathrm{base_{mineral}}}^\mathrm{298.15K}{a_\mathrm{{H_{2}CO_{3}}^{*}}} ](1 - \mathrm{SR_{mineral}}),
# \label{eq:transition_state}$$ with
# $$\mathrm{SR_{mineral}} = \frac{\mathrm{IAP_{mineral}}}{K_\mathrm{eq,mineral}},$$
#
# where $a$ denotes the activity of the ion, and $\mathrm{SR}$ is the abbreviation for the saturation ratio of a phase, which describes the ion activity product $\mathrm{IAP}$ divided by equilibrium constant $K_\mathrm{eq}$.
# In the simulation, the main species composition within the Opalinus Clay pore water are considered and the corresponding values are listed in table 2, following the
# work of Wersin et al. [8].
#
# | Minerals     |    Volume fraction [-]  | Molar volume in [$\mathrm{m^3/mol}$] | Molality in [$\mathrm{mol/kgw}$]      |
# |:--------------:|:-------------------------:|:--------------------------------------:|:---------------------------------------:|
# | Calcite      | 0.06375                 | 3.6934e-5                            | 11.507                                |
# | Dolomite-dis | 0.02125                 | 6.439e-5                             | 2.2                                   |
# | Illite       | 0.306                   | 1.4148e-4                            | 14.419                                |
# | Kaolinite    | 0.204                   | 9.935e-5                             | 13.689                                |
#
# Table 1: Parameter of the Opalinus Clay minerals used in the modeling
#
# | Species | Value in [mol/kgw]      |
# |:---------:|:-------------------------:|
# | C(4)    | $3.89 \times 10^{-3}$   |
# | Ca      | $1.89 \times 10^{-2}$   |
# | Mg      | $2.197 \times 10^{-2}$  |
# | Cl      | $3.2667 \times 10^{-1}$ |
# | K       | $1.92 \times 10^{-3}$   |
# | Na      | $2.8067 \times 10^{-1}$ |
# | S(6)    | $1.679 \times 10^{-2}$  |
# | Al      | $2.295 \times 10^{-10}$ |
# | Si      | $1.7 \times 10^{-5}$    |
# | Sr      | $4.6 \times 10^{-5}$    |
#
# Table 2: concentration of the species in the Opalinus Clay pore water
#
# For verification purposes, a corresponding PHREEQC model was constructed with identical parameter settings.
#
# ### Results
#
# On the Fig. 3 a) the evolution of calcium and magnesium molality during the $\ce{CO2}$ injection computed from the OGS-6 and PHREEQC model is shown.
# And the Fig. 3 b) displays the L2 relative difference norm between the results computed by the two software.
# The results of the two software programs match perfectly.
# On the Fig. 4 the evolution of the dissoluted molality of each minerals during the simulation is illustrated.
# Only calcite and dolomite have been partially dissolved due to $\ce{CO2}$ injection.
# In contrast, clay minerals, the illite and kaolinite were hardly affected.
# The kinetic dissolution rate of calcite and dolomite is mostly controlled by their saturation ratio states.
# When the SR value of calcite and dolomite reaches 1, the dissolution process of the correspondences is stopped.
#
# **Note:**
# The scenario #2 is not directly calculated in the notebook due to the lack of the large input file `llnl.dat` database in the commit. Also the post-processing of the scenario #2 is similar to which is showed in the scenario #1.
#
# ![Figure 3](figures/results_cl_comparision_sce2.png)
# Fig. 3: a) evolution of the calcium and magnesium molality during the $\ce{CO2}$ injection; b) L2 relative difference norm between the results computed by OGS-6 and PHREEQC.
#
# ![Figure 4](figures/results_disso_SR.png)
# Fig. 4: Evolution of the dissoluted minerals molality and the related saturation ratio of calcite and dolomite over the time.
#
# ## Literature
#
# <!-- vale off -->
#
# [1] BGR, CO2LPIE project, 2023. URL: [https://www.bgr.bund.de/DE/Themen/Endlagerung/Projekte/Wirtsgesteine_geotechnische_Barrieren/laufend/Nur-Deutsch/mont_terri_experimente.html?nn=1542156](https://www.bgr.bund.de/DE/Themen/Endlagerung/Projekte/Wirtsgesteine_geotechnische_Barrieren/laufend/Nur-Deutsch/mont_terri_experimente.html?nn=1542156).
#
# [2] Appelo, C. A. J. and Postma, D. (2004). Geochemistry, groundwater and pollution. CRC press.
#
# [3] Lu, R., Nagel, T., Poonoosamy, J., Naumov, D., Fischer, T., Montoya, V., Kolditz, O., and Shao, H.
# (2022). A new operator-splitting finite element scheme for reactive transport modeling in saturated porous
# media. Computers and Geosciences, 163(April 2021):105106. URL: https://doi.org/10.1016/j.cageo.2022.105106.
#
# [4] Plummer, L. N., Wigley, T. M., and Parkhurst, D. L. (1978). KINETICS OF CALCITE DISSOLUTION
# IN CO2-WATER SYSTEMS AT 5 degree TO 60 degree C AND 0. 0 TO 1. 0 ATM CO2.
#
# [5] Appelo, C. A. J. and Postma, D. (2023). Example 5.9: Kinetic dissolution of calcite. URL: https://www.hydrochemistry.eu/a&p/5/ex_5_9.phr.
#
# [6] Thury, M. (2002). The characteristics of the Opalinus Clay investigated in the Mont Terri underground
# rock laboratory in Switzerland. doi:10.1016/S1631-0705(02)01372-5.
#
# [7] Palandri, J. L. and Kharaka, Y. K. (2004). A compilation of rate parameters of water-mineral interaction
# kinetics for application to geochemical modeling.
#
# [8] Wersin, P., Mazurek, M., and Gimmi, T. (2022). Porewater chemistry of Opalinus Clay revisited:
# Findings from 25 years of data collection at the Mont Terri Rock Laboratory. Applied Geochemistry,
# 138(November 2021):105234. URL: https://doi.org/10.1016/j.apgeochem.2022.105234.
