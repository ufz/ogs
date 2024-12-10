# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = Original Modified Cam Clay model versus a basic version with constant elastic parameters under isotropic compression
# author = Prof. Dr. Thomas Nagel and Dr. Christian Silbermann
# note = Dummy frontmatter, not shown on website, only tested via wheels
# +++
#

# %% [markdown]
# # Interactive comparison:
# ## Original Modified Cam Clay model versus a basic version with constant elastic parameters under isotropic compression
#
# Comments to:
#
# *Prof. Dr. Thomas Nagel, Dr. Christian Silbermann
# Chair of Soil Mechanics and Foundation Engineering
# Geotechnical Institute
# Technische Universität Bergakademie Freiberg.*

# %%
# HIDDEN
import site

import matplotlib.pyplot as plt
import mtest
import numpy as np

# %%
# Path to MFront library including the material models
lib_path = site.getsitepackages()[0] + "/lib64/libOgsMFrontBehaviour.so"


# %%
# Some plot settings
plt.rcParams["lines.linewidth"] = 2.0
plt.rcParams["lines.color"] = "black"
plt.rcParams["legend.frameon"] = True
plt.rcParams["figure.figsize"] = (14, 14)
plt.rcParams["font.family"] = "serif"
plt.rcParams["legend.fontsize"] = 16
plt.rcParams["font.size"] = 16


# %%
# Set stress unit to Pa
kPa = 1e3


# %%
# Material parameters taken from Triax test
nu = 0.30  # Poisson ratio
la = 7.7e-2  # slope of the virgin consolidation line
ka = 6.6e-3  # slope of the swelling line
M = 1.2  # slope of the critical state line (CSL)
v0 = 1.785714286
phi0 = 1 - 1 / v0  # Initial porosity
pc0 = 200.0e3  # Initial pre-consolidation pressure in Pa


# %% [markdown]
# ## Analytical solution
#
# The essential MCC evolution equations for the pressure and the pre-consolidation pressure, i.e.
# $$
#   \frac{\dot{p}}{p} = -\left(\frac{v_0}{\kappa}\right) \dot{\varepsilon}_\text{e}^\text{V}
#   \quad,\quad
#   \frac{\dot{p}_\text{c}}{p_\text{c}} = -\left(\frac{v_0}{\lambda - \kappa}\right)\dot{\varepsilon}_\text{p}^\text{V}
#   \ . \\
# $$
# can be integrated easily by separation of variables. Assuming zero initial volume strains and the initial values $p_0$ and $p_\text{c0}$ we get the expressions
# $$
#   {\varepsilon}_\text{e}^\text{V} = -\left(\frac{\kappa}{v_0}\right)\ln\left(\frac{p}{p_0}\right)
#   \quad,\quad
#   {\varepsilon}_\text{p}^\text{V} = -\left(\frac{\lambda-\kappa}{v_0}\right) \ln\left(\frac{p_\text{c}}{p_\text{c0}}\right)
#   \ . \\
# $$
# With the additive composition of the total volumetric strain, i.e. ${\varepsilon}^\text{V} = {\varepsilon}_\text{e}^\text{V} + {\varepsilon}_\text{p}^\text{V}$ and with the linear kinematic relation $v - v_0 = v_0 \varepsilon^\text{V}$ we finally get
# $$
#   v = v_0 - \kappa\ln\left(\frac{p}{p_0}\right) - (\lambda-\kappa) \ln\left(\frac{p_\text{c}}{p_\text{c0}}\right)
#   \ . \\
# $$
# This is an analytical solution holding true for small changes of the volume ratio, i.e. $v\approx v_0$.
# Now, monotonic loading is considered for a stress-controlled isotropic compression test. Starting from the state $(p_0, v_0)$ with
# $$
# p=p_0, p_\text{c}=p_\text{c0} \quad\rightarrow\quad v = v_0
# $$
# the pressure is increased monotonically and the behavior is elastic until reaching the pre-consolidation pressure (the yield stress under isotropic compression) when
# $$
# p=p_\text{c0}, p_\text{c}=p_\text{c0} \quad\rightarrow\quad v = v_0 - \kappa\ln\left(\frac{p_\text{c0}}{p_0}\right) \ .
# $$
# Further increasing the pressure in the elastic-plastic region up to some end value $p=p_\text{E}>p_\text{c0}$, the final state is
# $$
# p=p_\text{E}, p_\text{c}=p_\text{E} \quad\rightarrow\quad v = v_0 - \kappa\ln\left(\frac{p_\text{E}}{p_0}\right) - (\lambda-\kappa) \ln\left(\frac{p_\text{E}}{p_\text{c0}}\right) \ .
# $$
# Between these corner points, the solution is a straight line in the semi-log space. To be precise, we choose $p_0 = p_\text{c0} / 4$, $p_\text{E} = 2 p_\text{c0}$. Thus, we obtain three states defining the analytical solution:
# $$
# \text{State 0:} \quad p=\frac{p_\text{c0}}{4} ,\quad v = v_0 \\
# \text{State 1:} \quad p=p_\text{c0} ,\quad v = v_0 - \kappa\ln(4) \\
# \text{State 2:} \quad p=2p_\text{c0} ,\quad v = v_0 - \kappa\ln(8) - (\lambda-\kappa) \ln(2) \\
# $$

# %%
# Loading
p_0 = pc0 / 4  # initial pressure
p_E = 2 * pc0  #  final pressure
OCR = pc0 / p_0  # highest experienced pressure over current/initial pressure

t_discrete = np.linspace(0, 1, 500)


# %%
# Analytical Solution: two straight lines in semi-log-space
def MCC_isotropic_compression_analytical(la, ka, v0, pc0, p_0, p_E):
    # swelling line (elastic)
    Pp_swl = [p_0, pc0]
    Pv_swl = [v0, v0 - ka * np.log(pc0 / p_0)]
    # normal compression line (elastoplastic)
    Pp_ncl = [pc0, p_E]
    Pv_ncl = [
        v0 - ka * np.log(pc0 / p_0),
        v0 - ka * np.log(p_E / p_0) - (la - ka) * np.log(p_E / pc0),
    ]

    return np.array(Pp_swl), np.array(Pv_swl), np.array(Pp_ncl), np.array(Pv_ncl)


# %% [markdown]
# ## Numerical solution


# %%
# MFront test definition
def MCC_isotropic_compression_numerical(
    mcc_model, lib_path, nu, la, ka, M, v0, phi0, pc0, p_0, p_E, t_discrete
):
    m = mtest.MTest()
    mtest.setVerboseMode(mtest.VerboseLevel.VERBOSE_QUIET)
    m.setMaximumNumberOfSubSteps(10)
    m.setBehaviour("generic", lib_path, mcc_model)
    m.setExternalStateVariable("Temperature", 293.15)

    m.setStress([-p_0, -p_0, -p_0, 0.0, 0.0, 0.0])
    m.setStrain([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    m.setImposedStress("SXX", {0: -p_0, 1: -p_E})
    m.setImposedStress("SYY", {0: -p_0, 1: -p_E})
    m.setImposedStress("SZZ", {0: -p_0, 1: -p_E})

    if mcc_model == "ModCamClay_semiExpl_constE":
        E = 3 * (1 - 2 * nu) / (1 - phi0) * p_0 / ka
        m.setMaterialProperty("YoungModulus", E)
        m.setParameter("AmbientPressure", 0.0)
        print("Young Modulus set to E =", E / 1e6, " MPa")
    if mcc_model in ("ModCamClay_semiExpl", "ModCamClay_semiExpl_absP"):
        m.setMaterialProperty("InitialVolumeRatio", v0)

    m.setMaterialProperty("PoissonRatio", nu)
    m.setMaterialProperty("CriticalStateLineSlope", M)
    m.setMaterialProperty("SwellingLineSlope", ka)
    m.setMaterialProperty("VirginConsolidationLineSlope", la)
    m.setMaterialProperty("CharacteristicPreConsolidationPressure", pc0)

    m.setInternalStateVariableInitialValue("PreConsolidationPressure", pc0)
    m.setInternalStateVariableInitialValue("VolumeRatio", v0)

    s = mtest.MTestCurrentState()
    wk = mtest.MTestWorkSpace()
    m.completeInitialisation()
    m.initializeCurrentState(s)
    m.initializeWorkSpace(wk)

    tau_xy = np.array([s.s0[3] / np.sqrt(2.0)])
    eps_xy = np.array([s.e0[3] / np.sqrt(2.0)])
    tau_yz = np.array([s.s0[5] / np.sqrt(2.0)])
    eps_yz = np.array([s.e0[5] / np.sqrt(2.0)])
    tau_xz = np.array([s.s0[4] / np.sqrt(2.0)])
    eps_xz = np.array([s.e0[4] / np.sqrt(2.0)])
    eps_xx = np.array([s.e0[0]])
    sig_xx = np.array([s.s0[0]])
    eps_yy = np.array([s.e0[1]])
    sig_yy = np.array([s.s0[1]])
    eps_zz = np.array([s.e0[2]])
    sig_zz = np.array([s.s0[2]])
    pc = np.array([pc0])
    vr = np.array([v0])
    for i in range(len(t_discrete) - 1):
        m.execute(s, wk, t_discrete[i], t_discrete[i + 1])
        eps_xy = np.append(eps_xy, s.e1[3] / np.sqrt(2.0))  # Kelvin mapping backwards!
        tau_xy = np.append(tau_xy, s.s1[3] / np.sqrt(2.0))  # Kelvin mapping backwards!
        eps_xz = np.append(eps_xz, s.e1[4] / np.sqrt(2.0))  # Kelvin mapping backwards!
        tau_xz = np.append(tau_xz, s.s1[4] / np.sqrt(2.0))  # Kelvin mapping backwards!
        eps_yz = np.append(eps_yz, s.e1[5] / np.sqrt(2.0))  # Kelvin mapping backwards!
        tau_yz = np.append(tau_yz, s.s1[5] / np.sqrt(2.0))  # Kelvin mapping backwards!
        sig_xx = np.append(sig_xx, s.s1[0])
        eps_xx = np.append(eps_xx, s.e1[0])
        sig_yy = np.append(sig_yy, s.s1[1])
        eps_yy = np.append(eps_yy, s.e1[1])
        sig_zz = np.append(sig_zz, s.s1[2])
        eps_zz = np.append(eps_zz, s.e1[2])
        pc = np.append(pc, s.getInternalStateVariableValue("PreConsolidationPressure"))
        vr = np.append(vr, s.getInternalStateVariableValue("VolumeRatio"))

    p = -(sig_xx + sig_yy + sig_zz) / 3
    return p, pc, vr


# %%
# Analytical solution:
[p_swl, v_swl, p_ncl, v_ncl] = MCC_isotropic_compression_analytical(
    la, ka, v0, pc0, p_0, p_E
)

# Numerical solution:
mcc_model = "ModCamClay_semiExpl_constE"
[p_data0, pc_data0, vr_data0] = MCC_isotropic_compression_numerical(
    mcc_model, lib_path, nu, la, ka, M, v0, phi0, pc0, p_0, p_E, t_discrete
)

mcc_model = "ModCamClay_semiExpl"
[p_data1, pc_data1, vr_data1] = MCC_isotropic_compression_numerical(
    mcc_model, lib_path, nu, la, ka, M, v0, phi0, pc0, p_0, p_E, t_discrete
)

mcc_model = "ModCamClay_semiExpl_absP"
[p_data2, pc_data2, vr_data2] = MCC_isotropic_compression_numerical(
    mcc_model, lib_path, nu, la, ka, M, v0, phi0, pc0, p_0, p_E, t_discrete
)


# %% [markdown]
# ### Consistency of initial elastic constants
#
# Given the swelling line slope $\kappa$ and the initial hydrostatic pressure $p_0$, the compression modulus is determined according to
# $$
# K = v_0 \frac{p_0}{\kappa} \ .
# $$
#
# With the Poisson's ratio $\nu$, Young's modulus can be obtained as
# $$
# E = 3(1-2\nu) K = v_0\cdot 3(1-2\nu) \frac{p_0}{\kappa} = 3 \frac{(1-2\nu)}{(1-\phi_0)} \frac{p_0}{\kappa} \ .
# $$
# In order to make the MCC models with pressure-dependent elasticity (the original MCC) and with constant elastic parameters *initially* consistent, $E$ cannot be chosen arbitrarily, but according to the formula above. Nevertheless, the hardening behaviour will be different, which is seen in the next figure.

# %%
fig, ax = plt.subplots(figsize=(12, 10))

ax.plot(
    p_swl / kPa,
    v_swl,
    color="black",
    lw=4,
    label="analytical solution: normal consolidation line",
)
ax.plot(p_ncl / kPa, v_ncl, color="black", lw=4)

ax.plot(
    p_data0 / kPa,
    vr_data0,
    label="numerical solution: MCC constant elastic parameters",
    color="red",
    ls="--",
)
ax.plot(
    p_data1 / kPa,
    vr_data1,
    label="numerical solution: MCC pressure-dependent elasticity incremental",
    color="orange",
    ls="-.",
)
# ax.plot(p_data2 / kPa, vr_data2, label='numerical solution: MCC pressure-dependent elasticity absolute', color='green' ,ls=':')
ax.legend()
ax.grid()
ax.set_title(f"Isotropic compression test with OCR={OCR}")
ax.set_xlabel("$p$ / kPa")
ax.set_ylabel("$v$")
ax.set_xscale("log")  # ,base=np.e)
fig.tight_layout()

fig.savefig(f"Comparison_NCL_MCC_new_OCR={OCR}.pdf")
plt.show()


# %% [markdown]
# ## Difference between the two implementation with pressure-dependent elasticity
#
# Note that for extrem high initial overconsolidation, the absolute formulation is superior to the fully incremental one. The latter one produces an offset right in the first time step, which persists throughout the simulation.
# To demonstrate this we use `p_0 = pc0/400` and plot only `data1` and `data2`.
#
# However, in practical applications with normal consolidation this seems hardly relevant.

# %%
# Loading
p_0 = pc0 / 400
p_E = 2 * pc0
OCR = pc0 / p_0
t_discrete = np.linspace(0, 1, 500)

# Analytical solution:
[p_swl, v_swl, p_ncl, v_ncl] = MCC_isotropic_compression_analytical(
    la, ka, v0, pc0, p_0, p_E
)

mcc_model = "ModCamClay_semiExpl"
[p_data3, pc_data3, vr_data3] = MCC_isotropic_compression_numerical(
    mcc_model, lib_path, nu, la, ka, M, v0, phi0, pc0, p_0, p_E, t_discrete
)

mcc_model = "ModCamClay_semiExpl_absP"
[p_data4, pc_data4, vr_data4] = MCC_isotropic_compression_numerical(
    mcc_model, lib_path, nu, la, ka, M, v0, phi0, pc0, p_0, p_E, t_discrete
)


# %%
fig, ax = plt.subplots(figsize=(12, 10))

ax.plot(
    p_swl / kPa,
    v_swl,
    color="black",
    lw=4,
    label="analytical solution: normal consolidation line",
)
ax.plot(p_ncl / kPa, v_ncl, color="black", lw=4)

ax.plot(
    p_data3 / kPa,
    vr_data3,
    label="numerical solution: MCC pressure-dependent elasticity incremental",
    color="orange",
    ls="-.",
)
ax.plot(
    p_data4 / kPa,
    vr_data4,
    label="numerical solution: MCC pressure-dependent elasticity absolute",
    color="green",
    ls=":",
)
ax.legend()
ax.grid()
ax.set_title(f"Isotropic compression test with OCR={OCR}")
ax.set_xlabel("$p$ / kPa")
ax.set_ylabel("$v$")
ax.set_xscale("log")
plt.show()
