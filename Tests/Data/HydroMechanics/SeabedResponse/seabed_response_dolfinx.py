from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import basix
import dolfinx as dolx
import numpy as np
import ufl
from dolfinx import fem, mesh
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.io import XDMFFile
from mpi4py import MPI

dolx.log.set_log_level(dolx.log.LogLevel.INFO)
print(
    f"DOLFINx version: {dolx.__version__} based on GIT commit: {dolx.git_commit_hash} of https://github.com/FEniCS/dolfinx/"
)

subprocess.run([sys.executable, "-m", "pip", "install", "pyvista"], check=True)


def plot_dolfinx_mesh(domain, result_folder: Path) -> None:
    import pyvista as pv

    pv.set_plot_theme("document")
    pv.set_jupyter_backend("static")

    topology, cell_types, geometry = dolx.plot.vtk_mesh(domain, domain.geometry.dim)
    grid = pv.UnstructuredGrid(topology, cell_types, geometry)

    plotter = pv.Plotter()
    plotter.add_mesh(grid, show_edges=True)
    plotter.view_xy()
    if pv.OFF_SCREEN:
        plotter.screenshot(str(result_folder / "mesh.png"))
    else:
        plotter.show()


def run_dolfinx_simulation(result_folder: Path) -> None:
    height = 100
    width = 200
    ny_cells = 25
    nx_cells = 50

    domain = mesh.create_rectangle(
        MPI.COMM_WORLD,
        [np.array([0, 0]), np.array([width, height])],
        [nx_cells, ny_cells],
        mesh.CellType.triangle,
    )
    cell_name = domain.topology.cell_name()
    gdim = domain.geometry.dim

    fdim = domain.topology.dim - 1

    vue = basix.ufl.element("Lagrange", cell_name, degree=2, shape=(gdim,))
    vpe = basix.ufl.element("Lagrange", cell_name, degree=1)
    mixed = basix.ufl.mixed_element([vue, vpe])
    v = fem.functionspace(domain, mixed)

    n_reference_fa3 = ufl.as_vector((0, 1))
    dimension = "2D"

    plot_dolfinx_mesh(domain, result_folder)

    def bottom(x):
        return np.isclose(x[1], 0)

    def top(x):
        return np.isclose(x[1], height)

    def left(x):
        return np.isclose(x[0], 0)

    def right(x):
        return np.isclose(x[0], width)

    bottom_facets = mesh.locate_entities_boundary(domain, fdim, bottom)
    top_facets = mesh.locate_entities_boundary(domain, fdim, top)
    left_facets = mesh.locate_entities_boundary(domain, fdim, left)
    right_facets = mesh.locate_entities_boundary(domain, fdim, right)

    marked_facets = np.hstack([bottom_facets, right_facets, top_facets, left_facets])
    marked_values = np.hstack(
        [
            np.full_like(bottom_facets, 1),
            np.full_like(right_facets, 2),
            np.full_like(top_facets, 3),
            np.full_like(left_facets, 4),
        ]
    )
    sorted_facets = np.argsort(marked_facets)
    facet_tag = mesh.meshtags(
        domain, fdim, marked_facets[sorted_facets], marked_values[sorted_facets]
    )

    def lame(e_modulus, poisson):
        mu_ = e_modulus / (2 * (1 + poisson))
        lambda_ = e_modulus * poisson / ((1 + poisson) * (1 - 2 * poisson))
        return mu_, lambda_

    e = fem.Constant(domain, 1e7)
    nu = fem.Constant(domain, 0.0)
    mu, lambda_ = lame(e, nu)

    rho_sr = fem.Constant(domain, 2650.0)
    n_f_ini = fem.Constant(domain, 0.4)

    k_s = ufl.as_tensor([[1.33e-11, 0], [0, 1.33e-11]])
    mu_fr = fem.Constant(domain, 1.3e-3)

    beta_m = fem.Constant(domain, 4.81e-10)
    p_ref = fem.Constant(domain, 101325.0)
    rho_fr_ref = fem.Constant(domain, 1000.0)

    b = fem.Constant(domain, (0.0, 0.0))

    t_end = 10
    n_steps = 100
    dt_value = t_end / n_steps
    dt = fem.Constant(domain, dt_value)

    u_all = fem.Function(v)
    u, p = ufl.split(u_all)

    test = ufl.TestFunction(v)
    vu, vp = ufl.split(test)

    trial = ufl.TrialFunction(v)
    _du, _dp = ufl.split(trial)

    ramp_steps = 5

    amplitude_wave = 0.1e5
    freq = 0.1
    omega = 2 * np.pi * freq
    wave_length = 100

    a_t = fem.Constant(domain, 0.0)

    v_p, _ = v.sub(1).collapse()
    p_wave = fem.Function(v_p)

    top_dofs_p = fem.locate_dofs_topological(
        (v.sub(1), v_p), facet_tag.dim, facet_tag.find(3)
    )

    v_ux, _ = v.sub(0).sub(0).collapse()
    u_x = fem.Function(v_ux)

    v_uy, _ = v.sub(0).sub(1).collapse()
    u_y = fem.Function(v_uy)

    bottom_dofs_uy = fem.locate_dofs_topological(
        (v.sub(0).sub(1), v_uy), facet_tag.dim, facet_tag.find(1)
    )
    right_dofs_ux = fem.locate_dofs_topological(
        (v.sub(0).sub(0), v_ux), facet_tag.dim, facet_tag.find(2)
    )
    left_dofs_ux = fem.locate_dofs_topological(
        (v.sub(0).sub(0), v_ux), facet_tag.dim, facet_tag.find(4)
    )

    bc_h1 = fem.dirichletbc(p_wave, top_dofs_p, v.sub(1))
    bc_m1 = fem.dirichletbc(u_y, bottom_dofs_uy, v.sub(0).sub(1))
    bc_m2 = fem.dirichletbc(u_x, right_dofs_ux, v.sub(0).sub(0))
    bc_m3 = fem.dirichletbc(u_x, left_dofs_ux, v.sub(0).sub(0))

    def kinematics(displacement):
        dim = len(displacement)
        identity = ufl.Identity(dim)
        f_def = identity + ufl.grad(displacement)
        jacobian = ufl.det(f_def)
        cof_f = jacobian * ufl.inv(f_def.T)
        b_s = ufl.dot(f_def, f_def.T)
        b_inv = ufl.inv(b_s)
        c = ufl.dot(f_def.T, f_def)
        c_inv = ufl.inv(c)
        return dim, identity, f_def, jacobian, cof_f, b_s, b_inv, c, c_inv

    def neumann_condition(dim, nc, n_reference_fa, p_fa):
        n_reference = n_reference_fa
        if nc == "constant":
            traction = 1 * p_fa * n_reference
        else:
            msg = "nc not defined"
            raise ValueError(msg)

        return ufl.as_vector((traction[0], traction[1])) if dim == "2D" else traction

    nc = "constant"
    x_ufl = ufl.SpatialCoordinate(domain)
    p_fa3 = -a_t * ufl.cos(((2 * ufl.pi) / wave_length) * x_ufl[0])

    q = fem.Constant(domain, 0.0)

    u_old = fem.Function(v, name="Initial")
    uold, pold = ufl.split(u_old)

    u_derivative = (u - uold) / dt
    p_derivative = (p - pold) / dt

    def porosity_development(displacement, initial_porosity):
        _d, _i, _f_s, j_s, _cof_f_s, _b_s, _b_inv, _c_s, _c_s_inv = kinematics(
            displacement
        )
        n_f = 1 - ((1 - initial_porosity) / j_s)
        e_local = n_f / (1 - n_f)
        n_s = 1 - n_f
        return n_f, e_local, n_s

    def density_development(displacement, initial_porosity, rho_solid, rho_fluid):
        n_f, _e, n_s = porosity_development(displacement, initial_porosity)
        rho_s = n_s * rho_solid
        rho_f = n_f * rho_fluid
        return rho_s + rho_f

    def second_piola_kirchhoff(displacement, mu_, lambda_local):
        _d, identity, _f, j_local, _cof_f, _b_s, _b_inv, _c, c_inv = kinematics(
            displacement
        )
        return mu_ * (identity - c_inv) + lambda_local * ufl.ln(j_local) * c_inv

    def cauchy(displacement, mu_, lambda_local):
        _d, _i, f_def, j_local, _cof_f, _b_s, _b_inv, _c, _c_inv = kinematics(
            displacement
        )
        s_eff = second_piola_kirchhoff(displacement, mu_, lambda_local)
        return ufl.inv(j_local) * ufl.dot(ufl.dot(f_def, s_eff), f_def.T)

    def d_green_lagrange(displacement, test_u):
        grad_v = ufl.grad(test_u)
        grad_u = ufl.grad(displacement)
        return (1 / 2) * (
            grad_v + grad_v.T + ufl.dot(grad_u.T, grad_v) + ufl.dot(grad_v.T, grad_u)
        )

    def variational_form_hm_mech(
        displacement,
        pressure,
        test_u,
        mu_,
        lambda_local,
        initial_porosity,
        rho_solid,
        rho_fluid_ref,
        body_force,
        p_fa,
        nc_local,
        n_reference,
        dim_local,
        dx,
        ds,
    ):
        _d, _i, _f_s, j_s, _cof_f_s, _b_s, _b_inv, _c_s, c_s_inv = kinematics(
            displacement
        )
        rho_fluid = rho_fluid_ref * ufl.exp(beta_m * (pressure - p_ref))
        rho_bulk = density_development(
            displacement, initial_porosity, rho_solid, rho_fluid
        )
        rho_0 = j_s * rho_bulk
        s_eff = second_piola_kirchhoff(displacement, mu_, lambda_local)
        traction = neumann_condition(dim_local, nc_local, n_reference, p_fa)
        d_egl_term = d_green_lagrange(displacement, test_u)
        vol_body_force = rho_0 * body_force
        return (
            ufl.inner(s_eff, d_egl_term) * dx
            - pressure * ufl.inner(j_s * c_s_inv, d_egl_term) * dx
            - ufl.inner(vol_body_force, test_u) * dx
            - ufl.inner(traction, test_u) * ds(3)
        )

    def variational_form_hm_hydro(
        displacement,
        displacement_dt,
        pressure,
        pressure_dt,
        test_p,
        beta_m_local,
        rho_fluid_ref,
        pressure_ref,
        permeability,
        viscosity,
        initial_porosity,
        body_force,
        fluid_flux,
        dx,
        ds,
    ):
        _d, _i, f_s, j_s, _cof_f_s, _b_s, _b_inv, _c_s, c_s_inv = kinematics(
            displacement
        )
        n_f, _e, _n_s = porosity_development(displacement, initial_porosity)
        rho_fluid = rho_fluid_ref * ufl.exp(beta_m_local * (pressure - pressure_ref))
        k_perm = permeability / viscosity
        k_spatial = j_s * ufl.dot(ufl.dot(ufl.inv(f_s), k_perm), ufl.inv(f_s.T))
        grad_p = ufl.grad(pressure)
        grad_vp = ufl.grad(test_p)
        d_egl = d_green_lagrange(displacement, displacement_dt)
        return (
            rho_fluid * ufl.inner(j_s * c_s_inv, d_egl) * test_p * dx
            + rho_fluid * ufl.dot(ufl.dot(grad_vp, k_spatial), grad_p) * dx
            - (rho_fluid**2)
            * ufl.dot(grad_vp, ufl.dot(k_spatial, (ufl.dot(body_force, f_s))))
            * dx
            + rho_fluid * j_s * n_f * beta_m_local * pressure_dt * test_p * dx
            - rho_fluid * ufl.inner(fluid_flux, test_p) * ds(3)
        )

    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tag)
    dx = ufl.Measure("dx", domain=domain, metadata={"quadrature_degree": 4})

    mech_form = variational_form_hm_mech(
        u,
        p,
        vu,
        mu,
        lambda_,
        n_f_ini,
        rho_sr,
        rho_fr_ref,
        b,
        p_fa3,
        nc,
        n_reference_fa3,
        dimension,
        dx,
        ds,
    )
    hydro_form = variational_form_hm_hydro(
        u,
        u_derivative,
        p,
        p_derivative,
        vp,
        beta_m,
        rho_fr_ref,
        p_ref,
        k_s,
        mu_fr,
        n_f_ini,
        b,
        q,
        dx,
        ds,
    )
    ff = mech_form + hydro_form

    petsc_options = {
        "snes_type": "newtonls",
        "snes_linesearch_type": "none",
        "snes_rtol": 1e-9,
        "snes_monitor": None,
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "ksp_type": "preonly",
        "ksp_rtol": 1.0e-18,
        "ksp_atol": 2e-7,
        "ksp_max_it": 2000,
        "ksp_monitor": None,
    }

    jc = ufl.derivative(ff, u_all, trial)

    problem = NonlinearProblem(
        ff,
        u_all,
        bcs=[bc_h1, bc_m1, bc_m2, bc_m3],
        petsc_options=petsc_options,
        petsc_options_prefix="seabed_",
        J=jc,
    )

    file_results_u = XDMFFile(
        domain.comm, str(result_folder / f"seabed_{nc}_u.xdmf"), "w"
    )
    file_results_p = XDMFFile(
        domain.comm, str(result_folder / f"seabed_{nc}_p.xdmf"), "w"
    )
    file_results_stress = XDMFFile(
        domain.comm, str(result_folder / f"seabed_{nc}_stress.xdmf"), "w"
    )
    file_results_u.write_mesh(domain)
    file_results_p.write_mesh(domain)
    file_results_stress.write_mesh(domain)

    k_degree = v.ufl_element().degree
    v_stress = fem.functionspace(domain, ("CG", k_degree - 1, (gdim, gdim)))
    stress_expr = cauchy(u, mu, lambda_)
    expr = fem.Expression(stress_expr, v_stress.element.interpolation_points)
    stress_h = fem.Function(v_stress, name="Cauchy")
    stress_h.interpolate(expr)

    v_u_out = fem.functionspace(
        domain, basix.ufl.element("Lagrange", cell_name, degree=1, shape=(gdim,))
    )
    v_p_out = fem.functionspace(
        domain, basix.ufl.element("Lagrange", cell_name, degree=1)
    )

    u_out = fem.Function(v_u_out, name="Displacement")
    p_out = fem.Function(v_p_out, name="Porewater pressure")

    time_values = np.linspace(0, t_end, n_steps + 1)
    for n, time_ in enumerate(time_values):
        print("Steps:", n)
        print("Time:", time_)

        ramp_factor = min(n / ramp_steps, 1.0)
        a_t.value = ramp_factor * amplitude_wave * np.sin(omega * time_)

        p_wave.interpolate(
            lambda x: a_t.value * np.cos(((2 * np.pi) / wave_length) * x[0])
        )

        problem.solve()
        converged = problem.solver.getConvergedReason()
        num_iter = problem.solver.getIterationNumber()
        assert converged > 0, f"Solver did not converge, got {converged}."
        print(
            f"Solver converged after {num_iter} iterations with converged reason {converged}."
        )

        u_out.interpolate(u_all.sub(0))
        p_out.interpolate(u_all.sub(1))

        stress_h.interpolate(expr)

        file_results_u.write_function(u_out, time_)
        file_results_p.write_function(p_out, time_)
        file_results_stress.write_function(stress_h, time_)

        u_old.x.array[:] = u_all.x.array
        u_old.x.scatter_forward()

    file_results_u.close()
    file_results_p.close()
    file_results_stress.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--result-folder",
        default="_out",
        help="Directory for simulation output files.",
    )
    args = parser.parse_args()

    output_dir = Path(args.result_folder)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"DOCKER: {output_dir}")
    run_dolfinx_simulation(output_dir)
