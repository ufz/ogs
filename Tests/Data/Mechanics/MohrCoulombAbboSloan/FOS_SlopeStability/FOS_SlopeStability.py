# ---
# jupyter:
#   jupytext:
#     notebook_metadata_filter: note
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python (.venv)
#     language: python
#     name: venv
# ---

# %% [raw]
# +++
# title = "Shear strength reduction - Slope stability"
# date = "2026-01-28"
# author = "Mostafa Mollaali, Thomas Nagel"
# image = "figures/schematic.png"
# web_subsection = "small-deformations"
# weight = 3
# +++

# %%
import os
from pathlib import Path

import numpy as np
import ogstools as ot
import pyvista as pv
import slope_postprocessing as sp
from slope_mesh_stage import CaseSpec, MeshStage
from slope_runner import MultiCaseRunner

# %%
ot.plot.setup.show_region_bounds = False

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## Shear strength reduction ($\varphi$–$c$ reduction)
#
# We compute the slope factor of safety using the **shear strength reduction (SSR)** method, in which the Mohr–Coulomb shear strength is reduced by a trial factor ($F$) until static equilibrium can no longer be achieved (typically observed as non-convergence and/or formation of a continuous plastic shear band). In the classical formulation, reduction is applied as
# \begin{equation}
# c_{\text{trial}}=\frac{c_0}{F}, \qquad
# \tan\varphi_{\text{trial}}=\frac{\tan\varphi_0}{F},
# \qquad
# \left(\Rightarrow\ \varphi_{\text{trial}}=\arctan\left(\frac{\tan\varphi_0}{F}\right)\right),
# \end{equation}
# and the limiting value of ($F$) at failure is interpreted as the factor of safety.
#
# ### Staged procedure in OpenGeoSys
#
# To avoid numerical artefacts from an unrealistic initial stress state, loading and reduction are applied in stages using pseudo-time. This staged procedure follows the OGS slope-stability workflow where gravity and surface load are ramped first and strength reduction starts later.
#
# **Stage 1 — Gravity ramp ($t=0\rightarrow1$)**
# Gravity is activated by ramping the **solid density** from 0 to its target value using a time-dependent function,
# \begin{equation}
# \rho(t)=\rho_0\cdot r_\rho(t), \qquad
# r_\rho(t)=\max\left(0,\ \min\left(1,\ \frac{t}{t_\rho}\right)\right),
# \qquad
# \rho_0=1600~\mathrm{kg/m^3},
# \end{equation}
# with $t_\rho=1$ and fixed body force $\mathbf{b}=(0,-9.81)\ \mathrm{m/s^2}$.
#
# **Stage 2 — Surcharge (top load) ramp ($t=1\rightarrow2$)**
# A distributed top load is applied using a time-dependent ramp,
# \begin{equation}
# q(t)=q_0\cdot r_q(t), \qquad
# r_q(t)=\max\left(0,\ \min\left(1,\ \frac{t-t_{q0}}{t_{q1}-t_{q0}}\right)\right),
# \qquad
# q_0=-30~\mathrm{kPa},
# \end{equation}
# with $t_{q0}=1$ and $t_{q1}=2$. Thus $q=0$ for $t\le 1$ and reaches full magnitude at $t=2$.
#
# **Stage 3 — Hold / equilibration ($t=2\rightarrow 3$)**
# After gravity and surcharge are fully applied, the model is held for **1** unit of pseudo-time with no strength reduction before weakening starts.
#
# **Stage 4 — Strength reduction (SSR) ($t\ge 3$)**
# Strength reduction is prescribed as a time-dependent multiplier $s(t)\in[s_{\min},1]$, and the trial factor is defined as
# \begin{equation}
# F(t)=\frac{1}{s(t)}.
# \end{equation}
# We use a linear decrease (clipped to $[s_{\min},1]$) during the SSR interval $t\in[t_0,t_1]$,
# \begin{equation}
# s(t)=
# \left(
# s_{\min},
# \quad
# \left(
# 1,
# \quad
# 1-(1-s_{\min})\frac{t-t_0}{t_1-t_0}
# \right)\right),
# \qquad
# t_0=3,\quad t_1=9,\quad s_{\min}=0.4.
# \end{equation}
#
# * **Cohesion reduction**
#   Cohesion is reduced linearly as
#   \begin{equation}
#   c(t)=s(t)\cdot c_0, \qquad c_0=5000~\mathrm{Pa}.
#   \end{equation}
#
# * **Friction and dilatancy reduction (SSR-consistent)**
#   In SSR the reduction is defined on $\tan\varphi$, i.e.
#   \begin{equation}
#   \tan\varphi(t)=\frac{\tan\varphi_0}{F(t)}=\tan\varphi_0\cdot s(t),
#   \qquad
#   \Rightarrow\quad
#   \varphi(t)=\arctan\big(\tan\varphi_0\cdot s(t)\big).
#   \end{equation}
#   An analogous relation is used for the dilatancy angle ($\psi$):
#   \begin{equation}
#   \psi(t)=\arctan\big(\tan\psi_0\cdot s(t)\big).
#   \end{equation}
#
# Failure is identified by loss of numerical convergence, corroborated by the formation of a continuous plastic shear band.
#

# %% [markdown]
# # Input data

# %%
# ----------------------------
# Project basename and template
# ----------------------------
BASENAME = "slope"
TEMPLATE_PRJ = Path(f"{BASENAME}.prj")
SRC_MESH_DIR = Path("mesh")

# ----------------------------
# MASTER: global configuration
# ----------------------------
MASTER = {
    "mesh_name": "slope_case",
    "geometry": {
        # Domain size / slope geometry [m]
        "L_bottom": 18.0,
        "H_top": 6.5,
        "H_bench": 2.0,
        "x_toe": 3.0,
        "x_crest": 12.0,
        # x-coordinates of top polyline "break points" [m]
        # (order matters)
        "top_breaks": [18.0, 15.5, 12.5, 12.0],
        # Target mesh size [m]
        "h": 0.25,
    },
    "mesh": {
        "algo_2d": 5,
        "recombine": False,  # True -> quads where possible
        "generate_mesh": True,
    },
    "physics": {
        # Body force [m/s^2] (gravity in -y)
        "specific_body_force": [0.0, -9.81]
    },
    "material": {
        "parameters": {
            "scalar": {
                # Elasticity (constants)
                "YoungModulus": 26e5,
                "PoissonRatio": 0.3,
                # Baseline strength constants (the only strength numbers user edits)
                "cohesion_rs0": 5000,  # [Pa]
                "FrictionAngle_rs0": 20,  # [deg]
                "DilatancyAngle_rs0": 10,  # [deg]
                "TransitionAngle": 27,  # [deg]
                "TensionCutOffParameter_rs0": 300,  # [Pa]
                # Baseline loads/constants
                "rho_sr0": 1600,
                "top_load0": -30000,
                "dirichlet0": 0,
            },
            "vector": {
                "displacement0": [0.0, 0.0],
                "T_ref": [293.15],
            },
        },
    },
    # staging schedule (user edits times only)
    "stages": {
        "rho_t1": 1.0,  # density ramp 0->1 on [0, rho_t1]
        "load_t0": 1.0,  # load starts ramping at t=load_t0
        "load_t1": 2.0,  # load reaches full at t=load_t1
    },
    # SSR schedule (user edits only these 3)
    "ssr": {
        "t0": 3.0,  # SSR starts
        "t1": 9.0,  # SSR reaches minimum
        "smin": 0.4,  # minimum strength multiplier s(t)=1/F(t)
    },
    "time_stepping": {
        # Time window
        "t_initial": 0.0,
        "t_end": 9.0,
        # dt controls
        "initial_dt": 0.1,
        "minimum_dt": 1e-3,
        "maximum_dt": 0.5,
        "number_iterations": "1 8 13 20",
        "multiplier": "1.2 1.0 0.9 0.5",
    },
    "output": {"prefix": "slope", "suffix": "_ts_{:timestep}_t_{:time}"},
}

# ----------------------------
# Load configurations
# ----------------------------
# These dictionaries are used by your geometry/BC generation logic.

LOAD_ALL = {
    "mode": "all",  # apply load on the default/top region (your code decides which part of the top is “eligible”)
    "L_load": 3.0,  # [m] patch length along the top boundary (used by edge/crest modes; may be ignored in "all")
    "clamp": True,
    "L_top_1": 2.5,  # [m] length of top segment 1 (only used when mode="top_segments")
    "L_top_2": 3.0,  # [m] length of top segment 2
    "L_top_3": 0.5,  # [m] length of top segment 3
    "loaded_zones": [
        2
    ],  # segment indices to load (only used when mode="top_segments"), e.g. [2]=only segment 2
}

LOAD_RIGHT = {
    **LOAD_ALL,
    "mode": "from_right_edge",  # load a patch of length L_load anchored at the right end of the top boundary
    "L_load": 3.0,
}

LOAD_CREST = {
    **LOAD_ALL,
    "mode": "from_slope_crest",  # load a patch of length L_load anchored at the slope crest (near x_crest)
    "L_load": 3.0,
}

LOAD_SEG2 = {
    **LOAD_ALL,
    "mode": "top_segments",  # split top into 3 segments (L_top_1/2/3) and load only 'loaded_zones'
    "L_top_1": 2.5,
    "L_top_2": 3.0,
    "L_top_3": 0.5,
    "loaded_zones": [2],
}


# Reference mesh size convenience
h_ref = float(MASTER["geometry"]["h"])

# ----------------------------
# Case list
# ----------------------------
# Notes on mesh handling:
# - Case 1 ("reference") uses a pre-existing mesh from the local ./mesh folder.
#   This avoids CI/test failures caused by small, OS-dependent differences in Gmsh mesh generation
#   (e.g., Linux vs macOS).
# - Case 2 ("load_right") generates a new mesh into _out/<case>/mesh to test the full workflow
#   (mesh generation + simulation) for a different loading configuration.
cases = [
    CaseSpec(
        label="reference",
        h=h_ref,
        element_order=2,
        nonlinear_reltol=1e-14,
        minimum_dt=1e-4,
        load=LOAD_ALL,
        load_tag="all",
        mesh_mode="use",
        src_mesh_dir=SRC_MESH_DIR,
        mesh_dir=SRC_MESH_DIR,
    ),
    CaseSpec(
        label="load_right",
        h=h_ref,
        element_order=2,
        nonlinear_reltol=1e-14,
        minimum_dt=1e-4,
        load=LOAD_RIGHT,
        load_tag="right",
        mesh_mode="generate",  # generate a new mesh into _out/<case>/mesh
        regenerate=True,  # force rebuild
        plot_mesh=True,
    ),
]


# %% [markdown]
# # Mesh generation and pre-processing

# %%
mesh_stage = MeshStage(out_root=out_dir)
mesh_results = mesh_stage.build_all(
    basename="slope",
    template_prj=TEMPLATE_PRJ,
    master=MASTER,
    cases=cases,
)


# %% [markdown]
# # Prepare the prj file and run the simulations

# %%
runner = MultiCaseRunner(out_root=out_dir)

ogs_results = runner.run_all(
    basename="slope",
    template_prj=TEMPLATE_PRJ,
    master=MASTER,
    cases=cases,
)


# %% [markdown]
# # Post-processing
# Plots: factor of safety vs pseudo-time, slope crown displacement vs time, and field maps.
#
# The first row  of profiles shows baseline fields at the at end of Stage 3 (hold/equilibration, ($t=3$)): mean stress ($\pi$),
# minimum principal stress ($\sigma_{\min}$), and trace strain ($\mathrm{tr}(\boldsymbol{\varepsilon})$).
#
# The second row shows end-of-SSR fields ($t=t_{\mathrm{end}}$): relative displacement magnitude
# ($\|\mathbf{u}-\mathbf{u}_0\|$), equivalent plastic strain ($\bar{\varepsilon}^{\,p}$),
# and trace strain ($\mathrm{tr}(\boldsymbol{\varepsilon})$).

# %%
sp.configure_plots(ogs_fontsize=16)
user_ranges = {}

sp.run_ssr_postprocessing(
    ogs_results,
    MASTER,
    ranges=user_ranges,
    auto_compute_missing_ranges=True,
    colorbar_mode="case",  # case or golbal
)

# %% [markdown]
# ## Reference check
# %% [markdown]
# We compare displacement fields at two key stages against stored VTUs
# (elastic at $t=3$, plastic during SSR at $t=4.7$). We run the current case
# and compare it with the reference results to ensure the comparison is consistent.
# %%
case_dir = (
    out_dir / "slope__ord-quad__h-0p25__dtmin-1e-04__tol-1e-14__load-all__reference"
)
results_dir = case_dir / "results"
expected_dir = Path("expected")

tests = [
    # Elastic stage (pre-SSR / before failure)
    ("3.000000", "slope_ts_17_t_3.000000.vtu", 1e-5, 1e-7, "elastic"),
    # Plastic stage (during SSR / closer to failure)
    ("4.700000", "slope_ts_41_t_4.700000.vtu", 5e-3, 1e-5, "plastic"),
]

for t, ref_name, rtol, atol_factor, stage in tests:
    new_vtu = next(iter(sorted(results_dir.glob(f"slope_ts_*_t_{t}.vtu"))), None)
    if new_vtu is None:
        msg = f"no vtu at t={t} in {results_dir.resolve()}"
        raise FileNotFoundError(msg)

    ref_vtu = expected_dir / ref_name
    if not ref_vtu.exists():
        msg = f"missing reference vtu: {ref_vtu.resolve()}"
        raise FileNotFoundError(msg)

    n = pv.read(new_vtu)
    r = pv.read(ref_vtu)
    if "displacement" not in n.point_data or "displacement" not in r.point_data:
        msg = f"missing point_data['displacement'] in new or ref at t={t}"
        raise KeyError(msg)

    s = r.sample(pv.PolyData(n.points))
    m = s.point_data.get("vtkValidPointMask")
    if m is not None and not bool(np.all(m)):
        msg = f"some points could not be sampled for t={t}"
        raise ValueError(msg)

    a = np.asarray(n.point_data["displacement"], float)
    b = np.asarray(s.point_data["displacement"], float)

    scale = max(float(np.nanmax(np.abs(b))), float(np.nanmax(b) - np.nanmin(b)), 1.0)
    np.testing.assert_allclose(a, b, rtol=rtol, atol=atol_factor * scale)

    print(f"OK ({stage}) t={t}")


# %% [markdown]
# # References
#
# 1. Abbo, A. J., & Sloan, S. W. (1995).
#    **A smooth hyperbolic approximation to the Mohr–Coulomb yield criterion.**
#    *Computers & Structures*, 54(3), 427–441.
#
# 2. Nagel, T., Minkley, W., Böttcher, N., Naumov, D., Görke, U.-J., & Kolditz, O. (2017).
#    **Implicit numerical integration and consistent linearization of inelastic constitutive models of rock salt.**
#    *Computers & Structures*, 182, 87–103.
#
# 3. Marois, G., Nagel, T., Naumov, D., & Helfer, T. (2020).
#    **Invariant-based implementation of the Mohr-Coulomb elasto-plastic model in OpenGeoSys using MFront.**
#    https://doi.org/10.13140/RG.2.2.34335.10403
#
# 4. Dawson, E. M., Roth, W. H., & Drescher, A. (1999).
#    **Slope stability analysis by strength reduction.**
#    *Géotechnique*, 49(6), 835–840.
#
# 5. Matsui, T., & San, K.-C. (1992).
#    **Finite element slope stability analysis by shear strength reduction technique.**
#    *Soils and Foundations*, 32(1), 59–70.
#
# 6. Zienkiewicz, O. C., Humpheson, C., & Lewis, R. W. (1975).
#    **Associated and non-associated visco-plasticity and plasticity in soil mechanics.**
#    *Géotechnique*, 25(4), 671–689.
#
# 7. Hergl, C., Kern, D., & Nagel, T. (2025).
#    **From Problem to Failure – Insights from the Eigenvalue Problem of the Stiffness Matrix in Non-linear Structural Analysis.**
#    *GAMM Archive for Students*, 7(1).
#
# 8. Kafle, L., Xu, W.-J., Zeng, S.-Y., & Nagel, T. (2022).
#    **A numerical investigation of slope stability influenced by the combined effects of reservoir water level fluctuations and precipitation: A case study of the Bianjiazhai landslide in China.**
#    *Engineering Geology*, 297(May 2021), 106508.
