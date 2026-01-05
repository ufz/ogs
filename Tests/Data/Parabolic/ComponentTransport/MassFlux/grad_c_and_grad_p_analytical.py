#!/usr/bin/env python
# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause


import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


class Solution:
    """
    Computes the analytical solution to a diffusion-advection-reaction equation
    in steady state on a 1D domain subject to a pressure gradient from p0 (left)
    to p1 (right) and a concentration gradient from c0 (left) to c1 (right).

    The porosities in the left and right halves of the domain differ.

    The analytical solution works for non-zero Darcy velocities only.
    """

    _phi1 = 0.65
    _phi0 = 0.15
    L = 0.8
    _D = 1e-9

    def __init__(self, *, c0, c1, p0, p1, r):
        phi0 = Solution._phi0
        phi1 = Solution._phi1
        D = Solution._D
        L = Solution.L

        mu = 1  # viscosity
        k = 1e-9  # permeability
        v_x = -k * (p1 - p0) / (mu * L)  # Darcy velocity

        print(f"Darcy velocity is {v_x} m/s")

        a0 = v_x / (phi0 * D)
        a1 = v_x / (phi1 * D)

        eps0 = np.exp(a0 * L / 2)
        eps1 = np.exp(a1 * L / 2)

        phit = phi0 / phi1
        rho = -r / v_x

        d0 = (rho / a1 * (1 - phit) * (1 - eps1) + c1 - c0 + rho * L) / (
            (eps0 - 1) / a0 + phit * eps0 * (eps1 - 1) / a1
        )

        d1 = (phit * d0 * eps0 + rho * (1 - phit)) / eps1

        f0 = c0 - d0 / a0
        f1 = c1 - d1 / a1 * eps1**2 + rho * L

        self._a0 = a0
        self._a1 = a1
        self._d0 = d0
        self._d1 = d1
        self._f0 = f0
        self._f1 = f1
        self.v_x = v_x
        self.p0 = p0
        self.p1 = p1
        self._r = r

        # check BCs
        assert self.c_ana(0) == c0
        # print(abs(self.c_ana(L) - c1))
        assert abs(self.c_ana(L) - c1) < 1e-14
        # check concenctration continuity
        eps = 1e-8
        # print(abs(self.c_ana(L/2 - eps) - self.c_ana(L/2 + eps)))
        assert abs(self.c_ana(L / 2 - eps) - self.c_ana(L / 2 + eps)) < 1e-7
        # check flux continuity
        # print(abs(self.flux_ana(L/2 - eps)[0] - self.flux_ana(L/2 + eps))[0])
        assert np.all(
            np.abs(self.flux_ana(L / 2 - eps) - self.flux_ana(L / 2 + eps)) < 1e-15
        )

    # general solution
    def _gs(self, d, a, f, x):
        r = self._r
        v = self.v_x
        return d / a * np.exp(a * x) + f + r * x / v

    # general gradient
    def _gg(self, d, a, x):
        r = self._r
        v = self.v_x
        return d * np.exp(a * x) + r / v

    # analytical solution (concentration profile)
    def c_ana(self, x):
        a0 = self._a0
        a1 = self._a1
        d0 = self._d0
        d1 = self._d1
        f0 = self._f0
        f1 = self._f1
        L = Solution.L

        return np.where(x < L / 2, self._gs(d0, a0, f0, x), self._gs(d1, a1, f1, x))

    # analytical solution (concentration gradient)
    def flux_ana(self, x):
        a0 = self._a0
        a1 = self._a1
        d0 = self._d0
        d1 = self._d1
        v_x = self.v_x
        phi0 = Solution._phi0
        phi1 = Solution._phi1
        D = Solution._D
        L = Solution.L

        phi = np.where(x < L / 2, phi0, phi1)
        grad_c = np.where(x < L / 2, self._gg(d0, a0, x), self._gg(d1, a1, x))
        c = self.c_ana(x)
        flux_x = v_x * c - phi * D * grad_c

        return np.c_[flux_x, np.zeros_like(x), np.zeros_like(x)]


def create_quasi_1D_mesh(L, n):
    mesh = pv.ImageData()
    mesh.dimensions = (n + 1, 2, 2)
    mesh.origin = (0, 0, 0)
    mesh.spacing = np.array([L, L, L]) / n
    mesh = mesh.cast_to_unstructured_grid()
    mesh.points_to_double()
    return mesh


def plot_solution_to_files(sol):
    xs = np.linspace(0, Solution.L, 1000)
    cs = sol.c_ana(xs)
    fluxes = sol.flux_ana(xs)

    fig, ax = plt.subplots()

    ax.plot(xs, cs)
    ax.set_xlabel("$x$ / m")
    ax.set_ylabel("$c$")
    ax.axvline(Solution.L / 2, ls=":", color="0.5")
    fig.subplots_adjust(left=0.18)

    fig.savefig("cs.pdf")

    fig, ax = plt.subplots()

    ax.plot(xs, fluxes[:, 0])
    ax.set_xlabel("$x$ / m")
    ax.set_ylabel("total mass flux")
    ax.axvline(Solution.L / 2, ls=":", color="0.5")
    fig.subplots_adjust(left=0.18)

    fig.savefig("fluxes.pdf")


def generate_meshes(sol, file_prefix, num_cells):
    L = Solution.L

    # mesh with initial data
    mesh = create_quasi_1D_mesh(L, num_cells)
    xs = mesh.points[:, 0]
    mesh.point_data["C_ini"] = sol.c_ana(xs)
    p0 = sol.p0
    p1 = sol.p1
    mesh.point_data["p_ini"] = (p0 * (L - xs) + p1 * xs) / L
    mesh.cell_data["MaterialIDs"] = (mesh.cell_centers().points[:, 0] > L / 2).astype(
        np.int32
    )
    mesh.point_data["bulk_node_ids"] = np.arange(len(xs), dtype=np.uint64)
    mesh.save(file_prefix + ".vtu")

    # extract boundary meshes
    eps = L * 1e-6
    surf = mesh.extract_surface()
    surf_ctr_xs = surf.cell_centers().points[:, 0]
    surf_left = surf.extract_cells(np.nonzero(surf_ctr_xs < eps)[0])
    surf_left.save(file_prefix + "_left.vtu")
    surf_right = surf.extract_cells(np.nonzero(surf_ctr_xs > L - eps)[0])
    surf_right.save(file_prefix + "_right.vtu")

    # mesh with reference data
    mesh = create_quasi_1D_mesh(L, num_cells)

    xs = mesh.points[:, 0]
    mesh.point_data["pressure"] = (p0 * (L - xs) + p1 * xs) / L
    mesh.point_data["C"] = sol.c_ana(xs)

    fluxes = sol.flux_ana(xs)
    mesh.point_data["CFlux"] = fluxes

    # The constant and linear flux profiles from these tests can be extrapolated
    # exactly on a FEM mesh.
    dim = fluxes.shape[1]
    mesh.cell_data["CFlux_residual"] = np.zeros((mesh.n_cells, dim))

    v = np.pad((sol.v_x,), (0, dim - 1))
    mesh.point_data["darcy_velocity"] = np.tile(v, (mesh.n_points, 1))
    mesh.cell_data["darcy_velocity_residual"] = np.zeros((mesh.n_cells, dim))

    mesh.save(file_prefix + "_ts_4_t_400000000.000000.vtu")


#### Advection only

# This test case has spatially homogeneous concentration and mass flux profiles.
# Therefore, it can be represented exactly on a low resolution FEM mesh.
sol = Solution(c0=1, c1=1, p0=0.3, p1=0, r=0)
generate_meshes(sol, "only_grad_p", 10)
# plot_solution_to_files(sol)


#### Advection and diffusion

# Spatially homogeneous flux profile, but varying concentrations.
sol = Solution(c0=0, c1=1, p0=0.3, p1=0, r=0)
generate_meshes(sol, "grad_c_and_grad_p", 1000)
# plot_solution_to_files(sol)


#### Advection and diffusion and reaction

# The spatially homogeneous reaction rate term leads to a constant slope flux
# profile, which will be computed rather exactly by OGS.
sol = Solution(c0=0.5, c1=1, p0=0.3, p1=0, r=-5e-10)
generate_meshes(sol, "grad_c_and_grad_p_and_r", 1000)
# plot_solution_to_files(sol)
