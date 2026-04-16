# %%
import json
import os
from pathlib import Path

import numpy as np
import porepy as pp
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
fracture_params_path = Path(out_dir, "fracture_params.json")
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

# --- Remove the 3D matrix cells to keep only fractures + intersections ---
mdg2d = mdg.copy()
for sd in list(mdg2d.subdomains()):
    if sd.dim == 3:
        mdg2d.remove_subdomain(sd)

# --- Clean up empty 1D subdomains & interfaces (avoid IndexError on export) ---
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
