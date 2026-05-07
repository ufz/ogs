# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "DFN generating using Porepy"
# date = "2025-04-28"
# author = "Mostafa Mollaali, Thomas Nagel"
# web_subsection = "reactive-transport"
# +++

# %%
import json
import os
from pathlib import Path

import numpy as np
import porepy as pp
import pyvista as pv
from numpy.random import default_rng

# %%

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
if not out_dir.exists():
    out_dir.mkdir(parents=True)
orig_dir = Path.cwd()


# %% [markdown]
# # DFN generating using Porepy
#
# > **Note:** This script generates the mesh for
# > [DFNbyPorePy_to_OGS.py](DFNbyPorePy_to_OGS.py)

# %% [markdown]
# **PorePy** is an open-source Python library developed by the University of Bergen, designed for simulating multiphysics processes in fractured and porous media. It emphasizes discrete fracture network (DFN) modeling through a mixed-dimensional approach.
#
# - Official site: [PorePy at University of Bergen](https://www.uib.no/en/rg/pmg/143656/porepy)
# - Source code: [PorePy GitHub Repository](https://github.com/pmgbergen/porepy)
#
# ### Installation
#
# Install PorePy within the same virtual environment as OGS:
#
# ```bash
# pip install porepy
# ```
#
# For full installation instructions, see the [PorePy Installation Guide](https://github.com/pmgbergen/porepy/blob/develop/Install.md).
#
#

# %% [markdown]
# ### Input data

# %%
domain_size = 100.0
mesh_size_boundary = 0.1 * domain_size
mesh_size_fracture = 0.05 * domain_size
mesh_size_min = 0.01 * domain_size

# %% [markdown]
# ### Domain seteup
# Defines the minimum and maximum coordinates of a 3D cubic domain.

# %%
mins = np.array([0.0, 0.0, 0.0])
maxs = np.array([domain_size, domain_size, domain_size])
bounding_box = {
    "xmin": mins[0],
    "xmax": maxs[0],
    "ymin": mins[1],
    "ymax": maxs[1],
    "zmin": mins[2],
    "zmax": maxs[2],
}
domain = pp.Domain(bounding_box=bounding_box)

# %% [markdown]
#
# ### Fracture generation loop
#
# - The loop iterates `n_fractures` times, with each iteration generating a single fracture.
# - For every fracture:
#   - A `random center` is generated within the domain by scaling a uniformly distributed random vector.
#   - A `random radius` is selected within the `radius_\text{R}ange`, setting its size between the allowed minimum and maximum values.
#   - Both the `major axis` and `minor axis` are set equal to the radius—while PorePy supports ellipses, we're using circular fractures here.
#   - `Strike and dip angles` are randomly drawn from the interval $[-\frac{\pi}{2}, \frac{\pi}{2}]$.
#     - These angles define the fracture plane's orientation (see [link](https://github.com/pmgbergen/porepy/blob/main/src/porepy/fracs/plane_fracture.py)):
#       - **Strike angle**: Sets the direction of the rotation axis in the horizontal plane, measured from the x-axis.
#       - **Dip angle**: Describes the tilt of the fracture plane from the horizontal.

# %%
use_saved_fractures = True
fracture_params_file = Path("fracture_params.json")
radius_range = np.array([50, 80])
n_fractures = 10

# borehole
borehole_height = 60
borehole_radius = maxs[0] * 0.005  # radius in x/y
z_center = maxs[2]
y_center = 0.5 * (mins[1] + maxs[1])
x_left = 0.2 * maxs[0]
x_right = 0.8 * maxs[0]

# %%
fractures = []
if use_saved_fractures and fracture_params_file.exists():
    fracture_params = json.loads(fracture_params_file.read_text())
    print("Loaded existing fracture parameters.")
else:
    rng = default_rng(12345)
    fracture_params = []
    for _ in range(n_fractures):
        center = rng.random(3) * (maxs - mins) + mins

        radius = rng.random() * (radius_range[1] - radius_range[0]) + radius_range[0]
        major_axis = minor_axis = radius

        major_axis_angle = rng.uniform(0, 2 * np.pi)
        strike_angle = rng.uniform(-0.5 * np.pi, 0.5 * np.pi)
        dip_angle = rng.uniform(-0.5 * np.pi, 0.5 * np.pi)

        fracture_params.append(
            {
                "center": center.tolist(),
                "major_axis": major_axis,
                "minor_axis": minor_axis,
                "major_axis_angle": major_axis_angle,
                "strike_angle": strike_angle,
                "dip_angle": dip_angle,
            }
        )

    fracture_params_file.write_text(json.dumps(fracture_params, indent=4))
    print("New fracture parameters generated and saved explicitly.")

for p in fracture_params:
    frac = pp.create_elliptic_fracture(
        center=np.array(p["center"]),
        major_axis=p["major_axis"],
        minor_axis=p["minor_axis"],
        major_axis_angle=p["major_axis_angle"],
        strike_angle=p["strike_angle"],
        dip_angle=p["dip_angle"],
    )
    fractures.append(frac)


#  vertical boreholes as long, thin fractures
center_left = np.array([x_left, y_center, z_center])
borehole_left = pp.create_elliptic_fracture(
    center=center_left,
    major_axis=borehole_height,
    minor_axis=borehole_radius,
    major_axis_angle=np.pi / 2,
    strike_angle=0.0,
    dip_angle=-0.5 * np.pi,
)
fractures.append(borehole_left)

center_right = np.array([x_right, y_center, z_center])
borehole_right = pp.create_elliptic_fracture(
    center=center_right,
    major_axis=borehole_height,
    minor_axis=borehole_radius,
    major_axis_angle=np.pi / 2,
    strike_angle=0.0,
    dip_angle=-0.5 * np.pi,
)
fractures.append(borehole_right)

# %% [markdown]
# ### Fracture network generation
#
# Once all fractures are generated, they're passed as a list to `create_fracture_network()`, along with the domain's bounding box. This creates a 3D Fracture Network object, representing the full 3D fracture system. The object takes care of intersection detection, applies the orientation from strike/dip angles, and prepares the geometry for meshing and simulation (see [FractureNetwork3d documentation](https://pmgbergen.github.io/porepy/html/docsrc/porepy/porepy.fracs.fracture_network_3d.html)).
#
# The system consists of two vertical fractures representing inlet and outlet boreholes, modeled as long, thin elliptic fractures. The minor axis is set to a small value (`borehole_Radius = maxs[0] * 0.005`) to reflect the narrow borehole dimensions.
# 1. **Left Borehole**: Located at `[x_Left, y_Center, zmax]`, with a major axis of `borehole_height` and a minor axis of `borehole_Radius`, oriented vertically (`major_axis_angle =`$\frac{\pi}{2}$, `dip_angle =`$\frac{\pi}{2}$).
# 2. **Right Borehole**: Located at `[x_Right, y_Center, zmax]`, with identical dimensions and orientation as the left borehole.
#
# Both fractures are added to the `fractures` list.
#

# %%
network = pp.create_fracture_network(fractures=fractures, domain=domain)

# %% [markdown]
# ### Meshing the fracture network
#
# We generate the computational mesh using PorePy's mixed-dimensional approach:
#
# - `cell_size_boundary`: maximum cell size near the domain boundaries.
# - `cell_size_fracture`: target cell size along the fracture surfaces.
# - `cell_size_min`: minimum allowed cell size anywhere in the mesh.
# - `simplex`: enables triangular (2D) or tetrahedral (3D) elements.
#
# This produces a **mixed-dimensional grid (`mdg`)** with 3D rock matrix cells, 2D fracture surfaces, and lower-dimensional intersections.
#
# For fracture-only simulations, the 3D matrix is removed, leaving just the 2D fractures and their intersections.

# %%
mesh_args = {
    "cell_size_boundary": mesh_size_boundary,
    "cell_size_fracture": mesh_size_fracture,
    "cell_size_min": mesh_size_min,
    "export": True,
    "filename": "mdg2d",
}

os.chdir(out_dir)
mdg = pp.create_mdg("simplex", mesh_args, network)
os.chdir(orig_dir)

mdg2d = mdg.copy()
for sd in list(mdg2d.subdomains()):
    if sd.dim == 3:
        mdg2d.remove_subdomain(sd)

for g in list(mdg2d.subdomains(dim=1)):
    cn = g.cell_nodes()
    if g.num_cells == 0 or cn.indices.size == 0 or cn.indptr.size <= 1:
        mdg2d.remove_subdomain(g)

# %% [markdown]
# ## Export the mesh to VTU format
#
# Creates a VTK file of the 2D mixed-dimensional grid (mdg2d) using `Exporter` (see [link](https://pmgbergen.github.io/porepy/html/docsrc/porepy/porepy.viz.exporter.html#module-porepy.viz.exporter)). The output is compatible with visualization tools like ParaView, and serves as input for OpenGeoSys.

# %%
exporter = pp.Exporter(mdg2d, "mixed_dimensional_grid", folder_name=out_dir).write_vtu()

# %%
os.chdir(orig_dir)
try:
    pv.set_jupyter_backend("static")
except Exception as e:
    print("PyVista backend not set:", e)

plot_mesh = pv.read(out_dir / "mixed_dimensional_grid_2.vtu")
scalar_name = "MaterialIDs" if "MaterialIDs" in plot_mesh.cell_data else "subdomain_id"

plotter = pv.Plotter(off_screen=True)
plotter.add_mesh(
    plot_mesh,
    scalars=scalar_name,
    cmap="tab20",
    opacity=0.4,
    scalar_bar_args={"title": scalar_name, "vertical": True},
)
plotter.add_mesh(plot_mesh.outline(), color="black", line_width=2)
plotter.show_axes()
plotter.enable_parallel_projection()
plotter.view_isometric()

output_path = out_dir / "mesh_overview.png"
plotter.screenshot(str(output_path))
plotter.show()
plotter.close()


# %% [markdown]
# ## Mesh-based regression checks
#
# The following assertions validate fracture geometry using mesh-derived quantities
# (orientation and in-plane metrics), not node-to-node equality.
#
# With PorePy v1.12, different operating systems can generate meshes with
# different numbers of nodes/cells for the same DFN input (e.g., Linux vs macOS).
# Using geometry-invariant checks keeps the validation robust across platforms.


# %%
def _expected_normal(strike_angle, dip_angle):
    return np.array(
        [
            -np.sin(strike_angle) * np.sin(dip_angle),
            np.cos(strike_angle) * np.sin(dip_angle),
            -np.cos(dip_angle),
        ]
    )


def _angle_diff_mod_pi(a, b):
    """Smallest absolute difference of two line orientations (period pi)."""
    return np.abs((a - b + 0.5 * np.pi) % np.pi - 0.5 * np.pi)


mesh_candidates = [
    out_dir / "mixed_dimensional_grid_2.vtu",
    out_dir / "mixed_dimensional_grid.vtu",
]
mesh_path = next((p for p in mesh_candidates if p.exists()), None)
assert mesh_path is not None, "Exported mesh file not found for fracture checks."

mesh = pv.read(mesh_path)
if "MaterialIDs" in mesh.cell_data:
    material_ids = np.asarray(mesh.cell_data["MaterialIDs"], dtype=int)
elif "subdomain_id" in mesh.cell_data:
    subdomain_id = np.asarray(mesh.cell_data["subdomain_id"], dtype=int)
    material_ids = subdomain_id - subdomain_id.min()
else:
    msg = "Neither 'MaterialIDs' nor 'subdomain_id' found in exported mesh."
    raise KeyError(msg)

fracture_params_from_json = json.loads(fracture_params_file.read_text())

normal_alignment_tol = 1e-8
strike_err_tol_deg = 1e-3
dip_abs_err_tol_deg = 1e-3
plane_dist_tol = 1e-8
center_plane_dist_tol = 1e-5
radius_upper_tol = 1e-8
radius_touch_tol = 1e-6

detail_rows = []
for mat_id, p in enumerate(fracture_params_from_json):
    cell_idx = np.flatnonzero(material_ids == mat_id)
    assert cell_idx.size > 0, f"No mesh cells found for fracture MaterialID {mat_id}."

    frac_mesh = (
        mesh.extract_cells(cell_idx)
        .extract_surface(algorithm=None)
        .triangulate()
        .compute_normals(
            cell_normals=True,
            point_normals=False,
            consistent_normals=False,
            auto_orient_normals=False,
            inplace=False,
        )
        .compute_cell_sizes(area=True)
    )

    cell_normals = np.asarray(frac_mesh.cell_data["Normals"], dtype=float)
    cell_areas = np.asarray(frac_mesh.cell_data["Area"], dtype=float)

    normal_expected = _expected_normal(
        float(p["strike_angle"]),
        float(p["dip_angle"]),
    )
    dots = cell_normals @ normal_expected

    # Normal direction on triangles can flip; align signs before averaging.
    cell_normals[dots < 0.0] *= -1.0
    normal_avg = (cell_normals * cell_areas[:, None]).sum(axis=0)
    normal_avg /= np.linalg.norm(normal_avg)

    alignment = np.abs(
        np.dot(normal_avg, normal_expected)
        / (np.linalg.norm(normal_avg) * np.linalg.norm(normal_expected))
    )

    angle_err_deg = float(np.degrees(np.arccos(np.clip(alignment, -1.0, 1.0))))

    strike_ref = np.mod(float(p["strike_angle"]), np.pi)
    strike_mesh = np.mod(np.arctan2(-normal_avg[0], normal_avg[1]), np.pi)
    strike_err_deg = float(np.degrees(_angle_diff_mod_pi(strike_mesh, strike_ref)))

    dip_ref_abs_deg = float(np.degrees(np.abs(float(p["dip_angle"]))))
    dip_mesh_abs_deg = float(np.degrees(np.arccos(np.clip(-normal_avg[2], -1.0, 1.0))))
    dip_abs_err_deg = abs(dip_mesh_abs_deg - dip_ref_abs_deg)

    pts = np.asarray(frac_mesh.points, dtype=float)
    center_ref = np.asarray(p["center"], dtype=float)
    signed_plane_dist = (pts - center_ref) @ normal_expected
    max_plane_dist = float(np.max(np.abs(signed_plane_dist)))
    center_to_plane_dist = float(np.abs((center_ref - pts[0]) @ normal_avg))

    # Radius checks are done in-plane around the JSON centre.
    vec = pts - center_ref
    vec_plane = vec - np.outer(vec @ normal_expected, normal_expected)
    radius_max = float(np.linalg.norm(vec_plane, axis=1).max())
    axis_major = float(p["major_axis"])
    axis_minor = float(p["minor_axis"])

    assert alignment > 1.0 - normal_alignment_tol, (
        f"Fracture {mat_id}: |dot(n_mesh, n_json)|={alignment:.12f} "
        f"is below {1.0 - normal_alignment_tol:.12f}."
    )
    assert strike_err_deg < strike_err_tol_deg, (
        f"Fracture {mat_id}: strike mismatch {strike_err_deg:.3e} deg "
        f"(tol {strike_err_tol_deg:.3e} deg)."
    )
    assert dip_abs_err_deg < dip_abs_err_tol_deg, (
        f"Fracture {mat_id}: |dip| mismatch {dip_abs_err_deg:.3e} deg "
        f"(tol {dip_abs_err_tol_deg:.3e} deg)."
    )
    assert max_plane_dist < plane_dist_tol, (
        f"Fracture {mat_id}: points are off JSON plane by {max_plane_dist:.3e} m "
        f"(tol {plane_dist_tol:.3e} m)."
    )
    assert center_to_plane_dist < center_plane_dist_tol, (
        f"Fracture {mat_id}: JSON centre is off mesh plane by "
        f"{center_to_plane_dist:.3e} m (tol {center_plane_dist_tol:.3e} m)."
    )
    assert radius_max <= axis_major + radius_upper_tol, (
        f"Fracture {mat_id}: mesh radius {radius_max:.6e} exceeds JSON major axis "
        f"{axis_major:.6e} (tol {radius_upper_tol:.1e})."
    )
    if np.isclose(axis_major, axis_minor, rtol=0.0, atol=1e-12):
        assert abs(radius_max - axis_major) < radius_touch_tol, (
            f"Fracture {mat_id}: circular fracture radius from mesh ({radius_max:.6e}) "
            f"does not match JSON radius ({axis_major:.6e}) within {radius_touch_tol:.1e}."
        )

    detail_rows.append(
        {
            "strike_err_deg": strike_err_deg,
            "dip_err_deg": dip_abs_err_deg,
            "axis_err_m": abs(radius_max - axis_major),
            "normal_err_deg": angle_err_deg,
            "max_plane_dist_m": max_plane_dist,
            "center_plane_dist_m": center_to_plane_dist,
        }
    )

max_strike_err = max(r["strike_err_deg"] for r in detail_rows)
max_dip_err = max(r["dip_err_deg"] for r in detail_rows)
max_axis_err = max(r["axis_err_m"] for r in detail_rows)
max_normal_err = max(r["normal_err_deg"] for r in detail_rows)
max_plane_dist = max(r["max_plane_dist_m"] for r in detail_rows)
max_center_dist = max(r["center_plane_dist_m"] for r in detail_rows)
print(
    f"Fracture validation passed ({len(detail_rows)}/{len(fracture_params_from_json)}). "
    f"Max errors: strike={max_strike_err:.3e} deg, dip={max_dip_err:.3e} deg, "
    f"axis={max_axis_err:.3e} m, normal={max_normal_err:.3e} deg, "
    f"plane={max_plane_dist:.3e} m, centre={max_center_dist:.3e} m."
)
