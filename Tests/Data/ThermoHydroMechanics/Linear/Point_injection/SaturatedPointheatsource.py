# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: venv-with-ogs
#     language: python
#     name: venv-with-ogs
# ---

# %% [raw]
# +++
# author = "Jörg Buchwald and Kata Kurgyis"
# date = "2022-11-02"
# title = "Point-Heatsource Problem"
# weight = 8
# image = "figures/placeholder_pointheatsource.png"
# web_subsection = "th2m"
# coupling = "thm"
# +++
#

# %% [markdown]
# ## Problem description
#
# The problem describes a heat source embedded in a fully fluid-saturated porous
# medium. The spherical symmetry is modeled using a 10 m x 10 m disc with a
# point heat source ($Q=150\;\mathrm{W}$) placed at one corner ($r=0$) and a
# curved boundary at $r=10\;\mathrm{m}$. Applying rotational axial symmetry at
# one of the linear boundaries, the model region transforms into a half-space
# configuration of the spherical symmetrical problemcorresponding to the
# analytical solution.
# The initial temperature and the excess pore pressure are 273.15 K and 0 Pa,
# respectively. The axis-normal displacements along the symmetry (inner)
# boundaries were set to zero, whereas the excess pore pressure, as well as the
# temperature are set to their initial values along the outer (curved) boundary.
# The heat coming from the point source is propagated through the medium,
# causing the fluid and the solid to expand at different rates. The resulting
# pore pressure (gradient) is triggering a thermally driven consolidation
# process caused by the fluid flowing away from the heat source until
# equilibrium is reached.
#
# ![PointHeatSourceSchematic.png](figures/PointHeatSourceSchematic.png)

# %% [markdown]
# # Governing equations
#
# For this problem we consider the following assumptions:
#
# * No thermal adverction is considered: $\rho_\text{w}c_\text{w}T_{,i} v_i = 0$.
#
# * Gravitational forces are neglected: $\rho g = 0$.
#
# * Both fluid and solid phases are intrinsically incompressible: $\alpha_B = 1$; $\beta = 0$.
#
# * No external fluid sink or source term: $q_H = 0$.
#
# * The porous medium is isotropic and homogeneous.
#
# These assumptions lead to the following set of governing equation describing
# the system behavior:
#
# **Energy balance**
#
# $$
# \begin{gather}
#     m \dot T - (K T_{,i})_{,i} = q_T
# %
# \\
# %
#     \text{where}\nonumber
# %
# \\
# %
#     m = \phi \rho_w c_w + (1-\phi) \rho_s c_s
# %
# \\
# %
#     K = \phi K_w + (1 - \phi) K_s
# %
# \\
# %
#     v_i = -\dfrac{k_s}{\eta} (p_{,i})
# \end{gather}
# $$
#
# **Mass balance**
#
# $$
# \begin{gather}
#     - a_u \dot T+ \dot u_{i,i} + v_{i,i} = 0
# %
# \\
# %
#     \text{where}\nonumber
# %
# \\
# %
#     a_u = \phi a_w + (1-\phi) a_s
# \end{gather}
# $$
#
# **Momentum balance**
#
# $$
# \begin{equation}
#     \sigma_{ij} = \sigma^\prime_{ij} - p \delta_{ij} = 0
# \end{equation}
# $$
#
# A detailed description about the problem formulation and equation derivation
# can be found in the original work of Booker and Savvidou (1985) or Chaudhry
# et al. (2019).
#
# ## Input parameters
#
# We considered the following set of values as input parameters:
#
# ![PointHeatSourceInput.png](figures/PointHeatSourceInput.png)
#

# %% [markdown]
# # The analytical solution
#
#
# The analytical solution of the coupled THM consolidation problem is derived
# in the original work of Booker and Savvidou (1985). In Chaudhry et al. (2019),
# a corrected solution is given for the effective stress term.
#
# For clarification, the equations below are based on the solid mechanics sign
# convention (tensile stress is positive). Furthermore, temporal partial
# derivative is indicated by the dot convention, while spatial partial
# derivatives are expressed by the comma convention, i.e.
# $(\cdot)_{,i}=\partial (\cdot)/\partial x_i$.
#
# The analytical solution for the three primary variables are expressed as:
#
# **Temperature**
#
# $$
# \begin{equation}
#     \Delta T = \dfrac{Q}{4 \pi K r} f^{\kappa}
# \end{equation}
# $$
#
# **Pore pressure**
#
# $$
# \begin{equation}
#     p = \dfrac{X Q}{(1 - \frac{c}{\kappa}) 4 \pi K r} (f^{\kappa}-f^{c})
# \end{equation}
# $$
#
# **Displacement of the solid skeleton**
#
# $$
# \begin{equation}
#     u_{i} = \dfrac{Q a_u x_i}{4 \pi K r} g^{\ast}
# \end{equation}
# $$
#
# In the above equations, the following derived parameters are used:
#
# $$
# \begin{align}
#     \kappa &= \dfrac{K}{m}
# %
# \\
# %
#     c &= \dfrac{k_s}{\eta}(\lambda + 2G)
# %
# \\
# %
#     r &= \sqrt{x_{1}^{2}+x_{2}^{2}+x_{3}^{2}}
# %
# \\
# %
#     X &= a_\text{u}\left(\lambda+2G\right)-b^{\prime}
# %
# \\
# %
#     Y &= \dfrac{1}{\lambda+2G}\left(\dfrac{X}{\left(1-\dfrac{c}{\kappa}\right)a_\text{u}}+\dfrac{b^{\prime}}{a_\text{u}}\right)
# %
# \\
# %
#     Z &= \dfrac{1}{\lambda+2G}\left(\dfrac{X}{\left(1-\dfrac{c}{\kappa}\right)a_\text{u}}\right)
# %
# \\
# %
#     f^{A} &= \text{erfc}\left(\dfrac{r}{2\sqrt{At}}\right),\quad A=\kappa,c
# %
# \\
# %
#     g^{A} &= \dfrac{At}{r^{2}}+\left(\frac{1}{2}-\dfrac{At}{r^{2}}\right)f^{A}-\sqrt{\dfrac{At}{\pi r^{2}}} \exp\left(-\dfrac{r^{2}}{4At}\right)
# %
# \\
# %
#     g^{\ast} &= Yg^{\kappa}-Zg^{c}
# %
# \\
# %
#     g^{A}_{,i} &= \frac{2x_{i}At}{r^{4}}\left(f^{A}-1+\frac{r}{\sqrt{\pi At}}\exp\left(-\frac{r^{2}}{4At}\right)\right),\quad i=1,2,3
# %
# \\
# %
#     g^{\ast}_{,i} &= Yg^{\kappa}_{,i}-Zg^{c}_{,i}
# \end{align}
# $$
#
# The corrected form of the effective stress:
#
# $$
# \begin{align}
#     \sigma^{\prime}_{ij|j=i} &= \frac{Q a_\text{u}}{4\pi Kr}\left( 2G\left[g^{\ast}\left(1-\frac{x^{2}_{i}}{r^{2}}\right)+x_{i}g^{\ast}_{,i}\right]+\lambda \left[x_{i}g^{\ast}_{,i}+2g^{\ast}\right]\right)-b^{\prime}\Delta T
# %
# \\
# %
#     \sigma^\prime_{ij|j \neq i} &= \frac{Q a_\text{u}}{4\pi Kr}\left( G\left[x_{i}g^{\ast}_{,j}+x_{j}g^{\ast}_{,i}-2g^{\ast}\dfrac{x_{i}x_{j}}{r^{2}}\right]\right)
# \end{align}
# $$

# %%
import concurrent.futures
import os
import platform
from pathlib import Path
from subprocess import run
from timeit import default_timer as timer

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pyvista as pv
from scipy.special import erfc


# %%
class ANASOL:
    def __init__(self):
        # material parameters
        self.phi = 0.16  # porosity of soil
        self.k = 2e-20  # coefficient of permeability
        self.eta = 1e-3  # viscosity water at 20 deg
        self.E = 5.0e9  # Youngs modulus
        self.nu = 0.3  # Poisson ratio
        self.rho_w = 999.1  # density of pore water
        self.c_w = 4280  # specific heat of pore water
        self.K_w = 0.6  # thermal conductivity of pore water
        self.rho_s = 2290.0  # density of solid matrix
        self.c_s = 917.654  # specific heat capacity of solid matrix
        self.K_s = 1.838  # themal conductivity of solid matrix
        # volumetric expansivity of matrix - conversion from linear volumetric
        self.a_s = 3 * 1.5e-5
        # coefficient of volume expansion of pore water (beta_w)
        self.a_w = 4.0e-4

        # initial and boundary condition
        # strength of heat source - value corrected to account for domain size
        self.Q = 2 * 150  # [Q] = W
        self.T0 = 273.15  # initial temperature

        self.Init()

    # derived parameters
    def f(self, ka, R, t):
        return erfc(R / (2 * np.sqrt(ka * t)))

    def g(self, ka, R, t):
        return (
            ka * t / R**2
            + (1 / 2 - ka * t / R**2) * erfc(R / (2 * np.sqrt(ka * t)))
            - np.sqrt(ka * t / (np.pi * R**2)) * np.exp(-(R**2) / (4 * ka * t))
        )

    def gstar(self, R, t):
        return self.Y * self.g(self.kappa, R, t) - self.Z * self.g(self.c, R, t)

    def R(self, x, y, z):
        return np.sqrt(x**2 + y**2 + z**2)

    def dg_dR(self, ka, i, R, t):
        return (2 * i / R**3) * np.sqrt(ka * t / np.pi) * np.exp(
            -R * R / (4 * ka * t)
        ) + (2 * i * ka * t / R**4) * (self.f(ka, R, t) - 1)

    def dgstar_dR(self, i, R, t):  # Subscript R means derivative w.r.t R
        return self.Y * self.dg_dR(self.kappa, i, R, t) - self.Z * self.dg_dR(
            self.c, i, R, t
        )

    # corrected form of effective stress
    def sigma_ii(self, x, y, z, t, ii):  # for normal components
        R = self.R(x, y, z)
        index = {"xx": x, "yy": y, "zz": z}
        return (self.Q * self.a_u / (4 * np.pi * self.K * R)) * (
            2
            * self.G
            * (
                self.gstar(R, t) * (1 - index[ii] ** 2 / R**2)
                + index[ii] * self.dgstar_dR(index[ii], R, t)
            )
            + self.lambd
            * (
                x * self.dgstar_dR(x, R, t)
                + y * self.dgstar_dR(y, R, t)
                + z * self.dgstar_dR(z, R, t)
                + 2 * self.gstar(R, t)
            )
        ) - self.bprime * (self.temperature(x, y, z, t) - self.T0)

    def sigma_ij(self, x, y, z, t, i, j):  # for shear components
        R = self.R(x, y, z)
        index = {"x": x, "y": y, "z": z}
        return (self.Q * self.a_u / (4 * np.pi * self.K * R)) * (
            2
            * self.G
            * (
                index[i] * self.dgstar_dR(index[j], R, t) / 2
                + index[j] * self.dgstar_dR(index[i], R, t) / 2
                - index[i] * index[j] * self.gstar(R, t) / R**2
            )
        )

    # primary variables
    def temperature(self, x, y, z, t):
        R = self.R(x, y, z)
        return self.Q / (4 * np.pi * self.K * R) * self.f(self.kappa, R, t) + self.T0

    def porepressure(self, x, y, z, t):
        R = self.R(x, y, z)
        return (
            self.X
            / (1 - self.c / self.kappa)
            * self.Q
            / (4 * np.pi * self.K * R)
            * (self.f(self.kappa, R, t) - self.f(self.c, R, t))
        )

    def u_i(self, x, y, z, t, i):
        R = self.R(x, y, z)
        index = {"x": x, "y": y, "z": z}
        return (
            self.a_u * index[i] * self.Q / (4 * np.pi * self.K * R) * self.gstar(R, t)
        )

    def Init(self):
        # derived constants
        self.lambd = (
            self.E * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))
        )  # Lame constant
        self.G = self.E / (2 * (1 + self.nu))  # shear constant

        self.K = (
            self.phi * self.K_w + (1 - self.phi) * self.K_s
        )  # average thermal conductivity
        self.m = (
            self.phi * self.rho_w * self.c_w + (1 - self.phi) * self.rho_s * self.c_s
        )
        self.kappa = self.K / self.m  # scaled heat conductivity
        self.c = self.k / self.eta * (self.lambd + 2 * self.G)

        self.aprime = self.a_s
        self.a_u = self.a_s * (1 - self.phi) + self.a_w * self.phi
        self.bprime = (self.lambd + 2 * self.G / 3) * self.aprime

        self.X = self.a_u * (self.lambd + 2 * self.G) - self.bprime
        self.Y = (
            1
            / (self.lambd + 2 * self.G)
            * (self.X / ((1 - self.c / self.kappa) * self.a_u) + self.bprime / self.a_u)
        )
        self.Z = (
            1
            / (self.lambd + 2 * self.G)
            * (self.X / ((1 - self.c / self.kappa) * self.a_u))
        )


ana_model = ANASOL()

# %% [markdown]
#
# ## The numerical solutions
#
# For the numerical solution we compare the Thermal-Hydro-Mechanical (THM -
# linear and quadratic mesh), Thermal-2-Phase-Hydro-Mechanical (TH2M) and
# Thermal-Richard-Mechanical (TRM - quadratic mesh) formulation of OGS.
#
# The TH2M and TRM formulation methods have essential differences when applied
# to an unsaturated media where a gas phase is also present along side the
# aqueous phase. The difference originates from the way how the two mobile phases are treated specifically in the equation system: in the TH2M formulation, both the gas phase and the liquid phase is explicitely present and each phase is comprised of the two distinct component of aqueous component and non-aqueous component. In this case, the gas phase has a variable pressure solved explicitely in the governing equations. On the other hand, the TRM model assumes that the gas phase mobility is high and fast enough that gas drainage can occur significantly faster than the other processes in the system and hence, gas pressure doesn't build up. This leads to the simplification, that no gas pressure is calculated in the TRM model explicitely.
#
# The THM model is a simplified form of the general TH2M model, where there is
# no gas phase, only the aqueous phase is present in the equation system.
#
# In addition to the different formulation, we also compare the performance of
# the THM formulation with a linear and a quadratic mesh as well.

# %% [markdown]
# Preparing the different models:

# %%
data_dir = os.environ.get("OGS_DATA_DIR", "../../..")
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
if not out_dir.exists():
    out_dir.mkdir(parents=True)

dir_th2m = f"{data_dir}/TH2M/THM/sphere/"
dir_trm = f"{data_dir}/ThermoRichardsMechanics/PointHeatSource/"
dirs = {"thm_lin": "./", "thm_quad": "./", "th2m": dir_th2m, "trm": dir_trm}

prj_files = {
    "thm_lin": "pointheatsource_linear-mesh.prj",
    "thm_quad": "pointheatsource_quadratic-mesh.prj",
    "th2m": f"{dir_th2m}point_heatsource.prj",
    "trm": f"{dir_trm}point_heat_source_2D.prj",
}
models = {
    process_key: ot.Project(
        input_file=prj,
        output_file=f"{out_dir}/point_heatsource_{process_key}.prj",
    )
    for process_key, prj in prj_files.items()
}
models["th2m"].replace_text(
    "150", xpath="./parameters/parameter[name='temperature_source_term']/value"
)

# %%
# Simulation time
t_end = 2e6  # <= was originally 5e6
for process_key, model in models.items():
    model.set(t_end=t_end)
    model.set(output_prefix=f"point_heatsource_{process_key}")
    model.write_input()

# %%
# Run models in parallel via concurrent.futures
ogs_models = [
    {
        "prj": model.prjfile,
        "logfile": f"{out_dir}/{process_key}-out.txt",
        "args": f"-o {out_dir} -m {dirs[process_key]} -s {dirs[process_key]}",
    }
    for process_key, model in models.items()
]


def run_ogs(model):
    print(f"Starting {model['prj']} ...\n")
    start_sim = timer()
    # Starting via ogs6py does not work ("cannot pickle lxml"), at least on mac.
    command = f"ogs {model['prj']} {model['args']} > {model['logfile']}"
    run(command, shell=True, check=True)
    runtime = timer() - start_sim
    return [f"Finished {model['prj']} in {runtime} s", runtime]


if platform.system() == "Darwin":
    import multiprocessing as mp

    mp.set_start_method("fork")

runtimes = []
start = timer()
with concurrent.futures.ProcessPoolExecutor() as executor:
    outputs = executor.map(run_ogs, ogs_models)
    for output in outputs:
        print(output[0])
        runtimes.append(output[1])
print(f"Elapsed time for all simulations: {timer() - start} s")

# %% [markdown]
# ## Evaluation and Results
#
# The analytical expressions together with the numerical model can now be
# evaluated at different points as a function of time (time series) or for a
# given time as a function of their spatial coordinates (along radial axis).

# %%
results = {
    process_key: ot.MeshSeries(f"{out_dir}/point_heatsource_{process_key}.pvd")
    for process_key in models
}


# %% [markdown]
# ### Time series plots for temperature, pressure and displacement
#
# Comparison between the analytical solution and the numerical solution shows
# very good agreement, as displayed below in the figures.

# %%
plt.rcParams["lines.linewidth"] = 2.0
plt.rcParams["font.family"] = "serif"
plt.rcParams["legend.fontsize"] = 14
plt.rcParams["font.size"] = 14

t_i, p_i = ("temperature_interpolated", "pressure_interpolated")
data_names = {
    "thm_lin": {"T": "temperature", "p": "pressure", "u": "displacement"},
    "thm_quad": {"T": t_i, "p": p_i, "u": "displacement"},
    "th2m": {"T": t_i, "p": f"gas_{p_i}", "u": "displacement"},
    "trm": {"T": t_i, "p": p_i, "u": "displacement"},
}
variables: dict[str, ot.variables.Variable] = {
    "T": ot.variables.temperature,
    "p": ot.variables.pressure,
    "u": ot.variables.displacement["x"].replace(output_unit="mm"),
}
style = {"thm_lin": "r+", "thm_quad": "rx", "th2m": "b+", "trm": "g+"}
pt = [0.5, 0.5, 0.0]
t = np.linspace(1, 50000 * 200, num=201, endpoint=True)


def plot_over_time(
    var_key: str, ana_sol: np.ndarray, err_bounds: list[float]
) -> plt.Figure:
    fig, axs = plt.subplots(ncols=2, figsize=(16, 6))
    var = variables[var_key]
    # only do unit conversion here
    ana_sol = var.replace(func=lambda _: _).transform(ana_sol)
    axs[0].plot(t, ana_sol, "k", label="analytical")

    for key, ms in results.items():
        values = var.transform(ms.probe(pt, data_names[key][var_key])).ravel()
        error = values - np.interp(ms.timevalues, t, ana_sol)
        assert np.all((error > err_bounds[0]) & (error < err_bounds[1]))

        axs[0].plot(ms.timevalues, values, style[key], label=f"ogs6_{key}")
        axs[1].plot(ms.timevalues, error, style[key], label=f"ogs6_{key}")

    axs[0].set_ylabel(var.get_label())
    axs[1].set_ylabel(f"error (numerical - analytical) / {var.difference.output_unit}")
    for ax in axs:
        ax.set_xlim(4.0e4, 3.0e6)
        ax.set_xscale("log")
        ax.set_xlabel("t / s")
        ax.legend()
    return fig


# %%
T_fig = plot_over_time("T", ana_model.temperature(*pt, t), [-0.06, 0.2])
T_fig.tight_layout()

# %%
p_fig = plot_over_time("p", ana_model.porepressure(*pt, t), [-0.06, 0.1])
p_fig.tight_layout()

# %%
u_fig = plot_over_time("u", ana_model.u_i(*pt, t, "x"), [-0.0035, 0.0005])
u_fig.tight_layout()

# %% [markdown]
# ### Plots for temperature, pressure and displacement along the radial axis
#
# The comparison between the analytical and the numerical results along the
# radial axis generally shows good agreement. The differences observed can be
# primarily explained by mesh discretization and finite size effects. This is
# particularly the case for the th2m simulation results, where the differences
# are slightly more emphasized which is the results of larger time steps.


# %%
# Time stamp for the results along the radial axis
# Radial coordinates for plotting
x = np.linspace(start=0.0001, stop=2.0, num=50)
t_i = 1.0e5


def plot_over_line(
    var_key: str, ana_sol: np.ndarray, err_bounds: list[float]
) -> plt.Figure:
    fig, axs = plt.subplots(ncols=2, figsize=(16, 6))
    var = variables[var_key]
    # only do unit conversion here
    ana_sol = var.replace(func=lambda _: _).transform(ana_sol)
    axs[0].plot(x, ana_sol, "k", label="analytical")

    for key, ms in results.items():
        mesh = ms.mesh(ms.closest_timestep(t_i))
        obs_pts = pv.PolyData([(x_i, 0, 0) for x_i in x[1:]]).sample(mesh)
        values = var.transform(obs_pts[data_names[key][var_key]]).ravel()
        pts_x = obs_pts.points[:, 0]
        error = values - np.interp(pts_x, x, ana_sol)
        # do not check first entry, which corresponds to the origin
        assert np.all((error[2:] > err_bounds[0]) & (error[2:] < err_bounds[1]))

        axs[0].plot(pts_x, values, style[key], label=f"ogs6_{key}")
        axs[1].plot(pts_x, error, style[key], label=f"ogs6_{key}")

    axs[0].set_ylim([-0.2, min(130, max(ana_sol) * 1.25)])
    axs[0].set_ylabel(var.get_label())
    axs[1].set_ylabel(f"error (numerical - analytical) / {var.difference.output_unit}")
    for ax in axs:
        ax.set_xlim([0.0, 2.0])
        ax.set_xlabel("r / m")
        ax.legend()
    return fig


# %%
T_fig = plot_over_line("T", ana_model.temperature(x, 0, 0, t_i), [-2.5, 0.5])
T_fig.tight_layout()

# %%
p_fig = plot_over_line("p", ana_model.porepressure(x, 0, 0, t_i), [-1.0, 2.5])
p_fig.tight_layout()

# %%
u_fig = plot_over_line("u", ana_model.u_i(x, 0, 0, t_i, "x"), [-0.015, 0.01])
u_fig.tight_layout()

# %% [markdown]
# ## Execution times
#
# To compare the performance of the different numerical solutions implemented in OGS6, we compare the execution time of the simulations. The linear thm and trm solutions perform best, while the quadratic thm and th2m solutions take significantly longer time to run. It is also important to mention here, that the time step size selected for the th2m solution are twice as big as the other 3 implementation, yet simulation time still takes longer than any of the other solution.

# %%
fig, ax = plt.subplots(figsize=(16, 6))
ax.bar(["thm linear", "thm quadratic", "th2m", "trm"], runtimes)
ax.set_ylabel("exec. time / s")
fig.tight_layout()

# %% [markdown]
# ## References
#
# [1] Booker, J. R.; Savvidou, C. (1985), Consolidation around a point heat
# source. International Journal for Numerical and Analytical Methods in
# Geomechanics, 1985, 9. Jg., Nr. 2, S. 173-184.
#
# [2] Chaudhry, A. A.; Buchwald, J.; Kolditz, O. and Nagel, T. (2019),
# Consolidation around a point heatsource (correction & verification).
# International Journal for Numerical and Analytical Methods in Geomechanics,
# 2019, <https://doi.org/10.1002/nag.2998>.

# %%
