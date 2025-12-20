#!/usr/bin/env python
# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause


# This script computes the reference solution of the "only_grad_c" test case.
# Only a concentration gradient is applied over a quasi-one-dimensional domain.
# The pressure gradient is zero. The system is in its steady state.

# The left half of the domain and the right half of the domain have different
# porosities. However, the mass flux through the domain must be the same
# everywhere (due to the steady state).

# In this example only diffusive mass flux plays a role.

import numpy as np
import pyvista as pv

c1 = 1
c0 = 0
phi1 = 0.65
phi0 = 0.15
L = 0.8
D = 1e-9

a0 = (c1 - c0) / (1 + phi0 / phi1) * 2 / L
a1 = (c1 - c0) / (1 + phi1 / phi0) * 2 / L
b0 = c0
b1 = (2 * c0 + c1 * (phi1 / phi0 - 1)) / (1 + phi1 / phi0)


# @np.vectorize
def c_ana(x):
    # return  a0 * x + b0 if x < L/2 else a1 * x + b1
    return np.where(x < L / 2, a0 * x + b0, a1 * x + b1)


# @np.vectorize
def flux_ana(x):
    # grad_c = a0 if x < L/2 else a1
    # phi = phi0 if x < L/2 else phi1

    # return (-D * phi * grad_c, 0)
    flux_x = -D * np.where(x < L / 2, a0 * phi0, a1 * phi1)
    return np.c_[flux_x, np.zeros_like(x)]


mesh = pv.read("mesh_2D.vtu")

xs = mesh.points[:, 0]
cs = c_ana(xs)
fluxes = flux_ana(xs)
min_flux_x = np.min(fluxes[:, 0])
max_flux_x = np.max(fluxes[:, 0])

# mass flux must be constant on the entire mesh
assert (max_flux_x - min_flux_x) / abs(max_flux_x) < 1e-15

print(f"concentration range: [{np.min(cs)}, {np.max(cs)}]")
print(f"flux_x range:        [{min_flux_x}, {max_flux_x}]")

mesh.point_data["pressure"] = np.zeros_like(xs)
mesh.point_data["C"] = cs
mesh.point_data["CFlux"] = fluxes

# spatially homogenous flux can be extrapolated exaclty
mesh.cell_data["CFlux_residual"] = np.zeros((mesh.n_cells, fluxes.shape[1]))

mesh.point_data["darcy_velocity"] = np.zeros((mesh.n_points, fluxes.shape[1]))
mesh.cell_data["darcy_velocity_residual"] = np.zeros((mesh.n_cells, fluxes.shape[1]))

mesh.save("only_grad_c_ts_4_t_400000000.000000.vtu")
