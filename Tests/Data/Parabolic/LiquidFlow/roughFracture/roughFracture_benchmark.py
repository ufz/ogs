# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python (.venv)
#     language: python
#     name: venv
# ---

# %% [raw]
# +++
# title = "Liquid flow in rough fracture"
# author = "Mostafa Mollaali, Thomas Nagel"
# date = "2026-05-29"
# web_subsection = "liquid-flow"
# +++

# %%
from __future__ import annotations

import contextlib
import importlib
import os
from pathlib import Path

import numpy as np
import ogstools as ot
import pyvista as pv
import roughfracture_postprocessing as rp
import roughfracture_preprocessing as rpre
import roughfracture_runner as rr

# %% [markdown]
# ## Introduction
#
# A rough fracture does not behave like a slot with a single width. The gap
# between the two walls (the aperture) varies from zero (where the walls
# touch) to several hundred micrometres (where they don't). Under normal
# stress the contact patches grow and the conductive network shrinks. Flow
# concentrates in a handful of preferential paths, and the rest of the
# fracture barely participates.
#
# We quantify this for three roughness classes (JRC 4, 7, 10) and three
# normal stress levels (0.2, 2, 10 MPa). The surfaces are synthetic fracture
# profiles generated to match a target JRC value using the method of
# Stigsson & Mas Ivars (2025). Steady-state Darcy flow is computed via
# transient time-marching to convergence using the OGS `LiquidFlow` process.

# %% [markdown]
# ## Problem setup
#
# ### Domain and mesh
#
# The fracture is a 100 mm x 100 mm square. It is represented by its
# midplane (the geometric average of the two rough walls), triangulated
# at 1 mm node spacing. That gives **10 201 nodes** and **20 000 triangular
# elements**. Because the walls are rough, the midplane is not flat: node
# elevations vary by a few millimetres across the domain.
#
# ### Boundary conditions
#
# | Face | Coordinate | Condition | Value |
# |------|-----------|-----------|-------|
# | Inlet | x = x_min = 0.0 m | Dirichlet (pressure) | $p_\mathrm{in} = 5 \times 10^5$ Pa |
# | Outlet | x = x_max = 0.1 m | Dirichlet (pressure) | $p_\mathrm{out} = 10^5$ Pa |
# | Bottom | y = y_min = 0.0 m | no-flow |  |
# | Top | y = y_max = 0.1 m | no-flow | |
#
# Prescribing pressure on both inlet and outlet guarantees a physically bounded
# pressure field regardless of fracture permeability. The resulting flow rate
# is an output — analogous to a laboratory constant-head permeability test.
# A fixed-flux (Neumann) inlet would require arbitrarily large pressures in
# nearly-closed fractures (high JRC, high stress), causing numerical blowup.
#
# ### Body force
#
# Gravity acts in the -z direction: **0, 0, -9.81 m s$^{-2}$**. The fracture
# plane is roughly horizontal and flow is driven laterally, so gravity
# affects the solution only through the fluid density. Over a 100 mm domain
# the hydrostatic contribution to the lateral pressure gradient is
# negligible.
#
# ### Fluid
#
# Water, with a linearly pressure-dependent density (slightly compressible):
#
# $$
# \rho(p) = \rho_0\bigl[1 + \beta(p - p_0)\bigr],
# \quad \rho_0 = 1000\;\text{kg m}^{-3},
# \quad \beta = 4.5 \times 10^{-7}\;\text{Pa}^{-1},
# \quad p_0 = 10^5\;\text{Pa}
# $$
#
# The mass-balance form of Darcy's law is used (`equation_balance_type =
# mass`), consistent with a slightly compressible fluid. Dynamic viscosity
# is 1e-3 Pa s (water at ~20 C).
#
# ### How the aperture is calculated
#
# The aperture field is produced using the
# **midplane method** of Stigsson & Mas Ivars (2025), combined with the
# **Barton–Bandis hyperbolic closure law** (Bandis et al. 1983) for
# normal closure. Local permeability is then assigned cell by cell via
# the **cubic law** (Witherspoon et al. 1980):
#
# $$k = \frac{w^2}{12}$$
#
# Contact cells (where the two walls touch) have zero mechanical aperture.
# A minimum aperture $w_{\min} = 10^{-8}$ m is applied before computing
# permeability, representing the residual hydraulic void that persists even
# under full contact. This regularises the permeability field and avoids
# a singular system.

# %% [markdown]
# ## Setup

# %%
with contextlib.suppress(Exception):
    ot.plot.setup.show_region_bounds = False

if pv.global_theme.trame.server_proxy_enabled:
    pv.set_jupyter_backend("client")
else:
    pv.set_jupyter_backend("static")

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

try:
    _HERE = Path(__file__).parent.resolve()
except NameError:
    _HERE = Path.cwd()

BASE_PRJ = _HERE / "roughFracture_synthesis_LF.prj"
SRC_MESH_DIR = _HERE / "meshes"

JRC_LIST = [4, 7, 10]
sigmas_mpa = [0.2, 2.0, 10.0]

# %% [markdown]
# ## Input parameters

# %%
user_parameters_rough = {
    "prefix": "roughFracture_synthesis_LF",
    "t_end": "1000",
    "initial_dt": ".1",
    "minimum_dt": "1e-3",
    "maximum_dt": "100",
    "specific_body_force": "0 0 -9.81",
    "initial_pressure": 100000,
    "outlet_pressure": 100000,  # Pa -- p_out, Dirichlet at outlet (x = x_max)
    "inlet_pressure": 500000,  # Pa -- p_in,  Dirichlet at inlet  (x = x_min); Delta_P = 4e5 Pa
    "porosity_value": "1.0",
    "fluid_density": {
        "type": "Linear",
        "reference_value": "1000",
        "variable_name": "liquid_phase_pressure",
        "reference_condition": "1e5",
        "slope": "4.5e-7",
    },
    "equation_balance_type": "mass",
    "permeability_mesh_field": "permeability_cubic",
    "fracture_thickness_mesh_field": "aperture_closed",
    "user_w_min_aperture": "1e-8",  # minimum aperture [m] applied to contact cells
}

# %% [markdown]
# ## Meshes
#
# Nine input meshes in `meshes/`, one per (JRC, s_n) pair. Each carries
# `aperture_closed` and `contact_closed` as cell fields from the synthesis.

# %%
rr.print_mesh_inventory(SRC_MESH_DIR)

# %% [markdown]
# JRC 10 at 0.2 MPa, the roughest, most open case. Vertical scale
# exaggerated 2x to make the surface relief visible.

# %%
importlib.reload(rpre)

rpre.plot_fracture_surface_3d(
    SRC_MESH_DIR / "joint_JRC10_sigma_0p2MPa.vtu",
    z_scale=2,
)

# %% [markdown]
# Aperture fields for all nine cases (rows = JRC, columns = normal stress).
# Increasing stress (left to right) causes contact patches to grow and
# open channels to narrow. Increasing roughness (top to bottom) changes
# the spatial pattern of open and closed regions.

# %%
rpre.plot_aperture_grid(
    jrc_list=JRC_LIST,
    sigmas_mpa=sigmas_mpa,
    mesh_dir=SRC_MESH_DIR,
)

# %% [markdown]
# ## Fracture geometry
#
# Aperture decreases with stress and contact area grows, as expected
# from the Barton-Bandis closure law. Rougher surfaces (higher JRC) show
# larger spatial variability at any given stress.

# %%
importlib.reload(rpre)

rpre.plot_fracture_geometry_summary(
    jrc_list=JRC_LIST,
    sigmas_mpa=sigmas_mpa,
    mesh_dir=SRC_MESH_DIR,
)

# %% [markdown]
# ## Simulation

# %%
PREPARE_CASES = True
RUN_OGS = True

importlib.reload(rr)

prepared_cases = []

if PREPARE_CASES:
    prepared_cases = rr.prepare_all_cases(
        out_dir=out_dir,
        base_prj=BASE_PRJ,
        src_mesh_dir=SRC_MESH_DIR,
        user_parameters=user_parameters_rough,
        jrc_list=JRC_LIST,
        sigmas_mpa=sigmas_mpa,
        w_min=float(user_parameters_rough["user_w_min_aperture"]),
    )

if RUN_OGS:
    if not prepared_cases:
        prepared_cases = rr.prepare_all_cases(
            out_dir=out_dir,
            base_prj=BASE_PRJ,
            src_mesh_dir=SRC_MESH_DIR,
            user_parameters=user_parameters_rough,
            jrc_list=JRC_LIST,
            sigmas_mpa=sigmas_mpa,
            w_min=float(user_parameters_rough["user_w_min_aperture"]),
        )
    rr.run_all_cases(prepared_cases)

# %% [markdown]
# ## Results
#
# Three fields are shown for all nine cases: pressure with streamlines,
# aperture, and contact indicator. All cases are driven by the same pressure
# difference ($\Delta P = 4 \times 10^5$ Pa); the resulting flow rate varies with the
# fracture conductivity. The streamlines trace where the flow actually goes,
# channelling around contact zones through the wider aperture paths. At
# high normal stress the open network shrinks to a few narrow corridors.

# %%
RUN_POSTPROCESSING = True

importlib.reload(rp)

if RUN_POSTPROCESSING:
    rp.run_roughfracture_postprocessing(
        out_dir=out_dir,
        jrc_list=JRC_LIST,
        sigmas_mpa=sigmas_mpa,
    )

# %% [markdown]
# ## Assertion
#
# The pressure field at the last timestep is checked against stored reference
# results for all nine (JRC, sigma) combinations. References were produced with
# Dirichlet pressure BCs on both boundaries (p_in = 5e5 Pa, p_out = 1e5 Pa)
# and w_min = 1e-8 m.

# %%

REF_DIR = _HERE / "expected"

for jrc in JRC_LIST:
    for sigma in sigmas_mpa:
        s_tag = rr.sigma_tag(sigma)
        prefix = f"roughFracture_synthesis_LF_JRC{jrc}_sigma_{s_tag}"
        pvd = (
            out_dir / "runs" / f"JRC_{jrc}" / f"sigma_{s_tag}" / "out" / f"{prefix}.pvd"
        )

        ref_matches = sorted(REF_DIR.glob(f"{prefix}_ts_*_t_1000.000000.vtu"))
        if not ref_matches:
            msg = f"No reference file found in {REF_DIR} for {prefix}"
            raise FileNotFoundError(msg)
        ref = pv.read(ref_matches[0])

        p_comp = np.asarray(ot.MeshSeries(str(pvd))[-1].point_data["pressure"])
        p_ref = np.asarray(ref.point_data["pressure"])

        np.testing.assert_allclose(p_comp, p_ref, rtol=2e-3, atol=1e-3)
        print(
            f"JRC={jrc}, sigma={sigma:g} MPa: max|dp| = {np.abs(p_comp - p_ref).max():.3g} Pa -- OK"
        )

# %% [markdown]
# ## References
#
# Stigsson, M. & Mas Ivars, D. (2025). Synthetic rough fracture surfaces
# from JRC and fractal parameters. *Rock Mechanics and Rock Engineering*
# (in preparation).
#
# Bandis, S., Lumsden, A. C. & Barton, N. R. (1983). Fundamentals of rock
# joint deformation. *International Journal of Rock Mechanics and Mining
# Sciences*, 20(6), 249–268.
#
# Barton, N., Bandis, S. & Bakhtar, K. (1985). Strength, deformation and
# conductivity coupling of rock joints. *International Journal of Rock
# Mechanics and Mining Sciences*, 22(3), 121–140.
#
# Witherspoon, P. A. et al. (1980). Validity of cubic law for fluid flow
# in a deformable rock fracture. *Water Resources Research*, 16(6),
# 1016–1024.
