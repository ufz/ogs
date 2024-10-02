import pyvista as pv

mesh = pv.read("axisymmetric_bbar_ts_1_t_1.000000.vtu")
E = 1.0e10
nv = 0.499
p = 1e7

sigma_analytical = mesh.point_data["sigma"]
sigma_analytical[:, 0] = 0.0
sigma_analytical[:, 1] = -p
sigma_analytical[:, 2] = 0.0
sigma_analytical[:, 3] = 0.0
mesh.point_data["analytic_sigma"] = sigma_analytical

eps_analytical = mesh.point_data["epsilon"]
eps_1 = -p / E
eps_0 = -eps_1 * nv
eps_2 = eps_0
eps_3 = 0.0
eps_analytical[:, 0] = eps_0
eps_analytical[:, 1] = eps_1
eps_analytical[:, 2] = eps_2
eps_analytical[:, 3] = eps_3
mesh.point_data["analytic_eps"] = eps_analytical


u_analytical = mesh.point_data["displacement"]
u_analytical[:, 0] = eps_0 * mesh.points[:, 0]
u_analytical[:, 1] = eps_1 * mesh.points[:, 1]
mesh.point_data["analytic_displacement"] = u_analytical

mesh.save("axisymmetric_bbar_ts_1_t_1.000000.vtu")
# mesh.plot(scalars="analytic_displacement")
