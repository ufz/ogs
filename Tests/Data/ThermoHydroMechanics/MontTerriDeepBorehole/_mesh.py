from pathlib import Path

import gmsh
import numpy as np
import ogstools as ot
from _helper_functions import DOMAIN_LENGTH, LAYER_DZ, temperature_profile

width = 10.0


def _add_next_layer(z, h, p_right, p_left, l_top):
    p3_new = gmsh.model.geo.addPoint(width, z + h, 0)
    p4_new = gmsh.model.geo.addPoint(0, z + h, 0)

    l2 = gmsh.model.geo.addLine(p_right, p3_new)
    l3_new = gmsh.model.geo.addLine(p3_new, p4_new)
    l4 = gmsh.model.geo.addLine(p4_new, p_left)

    cl = gmsh.model.geo.addCurveLoop([-l_top, l2, l3_new, l4])
    s = gmsh.model.geo.addPlaneSurface([cl])

    gmsh.model.geo.mesh.setTransfiniteCurve(l2, -h)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3_new, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(l4, -h)

    gmsh.model.geo.mesh.setTransfiniteSurface(s)
    gmsh.model.geo.mesh.setRecombine(2, s)

    z += h

    return z, p3_new, p4_new, l3_new


def _insert_T_init(meshes):
    THM_analytical_T, THM_analytical_T_z = temperature_profile(gradT=True)

    THM_analytical_T_z = THM_analytical_T_z - DOMAIN_LENGTH
    THM_analytical_T = THM_analytical_T + 273.15

    def T_init(points):
        return np.flip(
            np.interp(
                np.flip(points[:, 1]),
                np.flip(THM_analytical_T_z),
                np.flip(THM_analytical_T),
            )
        )

    meshes["domain"].point_data["T_init"] = T_init(meshes["domain"].points)


def generate_mesh(savepath):

    gmsh.initialize()
    gmsh.model.add("MTDB_Modified")
    gmsh.option.setNumber("General.Terminal", 0)

    # each rectangle height; LAYER_DZ is the same layer thicknesses used by
    # _helper_functions.temperature_profile for the mesh's initial condition
    heights = (-LAYER_DZ).tolist()
    z = 0.0

    # Make rectangle for first layer
    h = heights[0]

    p1 = gmsh.model.geo.addPoint(0, z, 0)
    p2 = gmsh.model.geo.addPoint(width, z, 0)
    p3 = gmsh.model.geo.addPoint(width, z + h, 0)
    p4 = gmsh.model.geo.addPoint(0, z + h, 0)

    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    s = gmsh.model.geo.addPlaneSurface([cl])

    gmsh.model.geo.mesh.setTransfiniteCurve(l1, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, -h)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(l4, -h)

    gmsh.model.geo.mesh.setTransfiniteSurface(s)
    gmsh.model.geo.mesh.setRecombine(2, s)

    z += h

    p_right = p3
    p_left = p4
    l_top = l3

    for h in heights[1:]:
        z, p_right, p_left, l_top = _add_next_layer(z, h, p_right, p_left, l_top)

    gmsh.model.geo.synchronize()

    ### Physical groups
    ## Surfaces
    gmsh.model.addPhysicalGroup(2, [1], name="Hauptrogenstein")
    gmsh.model.addPhysicalGroup(2, [2], name="PasswangI")
    gmsh.model.addPhysicalGroup(2, [3], name="PasswangII")
    gmsh.model.addPhysicalGroup(2, [4], name="PasswangIII")
    gmsh.model.addPhysicalGroup(2, [5], name="PasswangIV")
    gmsh.model.addPhysicalGroup(2, [6], name="Opa_sandy1")
    gmsh.model.addPhysicalGroup(2, [7], name="Opa_shaly1")
    gmsh.model.addPhysicalGroup(2, [8], name="Opa_sandy2")
    gmsh.model.addPhysicalGroup(2, [9], name="Opa_carbo")
    gmsh.model.addPhysicalGroup(2, [10], name="Opa_shaly2")
    gmsh.model.addPhysicalGroup(2, [11], name="Staffelegg")

    ## Boundaries
    gmsh.model.addPhysicalGroup(1, [1], name="Upper_boundary")
    gmsh.model.addPhysicalGroup(1, [l_top], name="Bottom_boundary")

    gmsh.model.geo.synchronize()

    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.model.mesh.generate(2)

    msh_path = str(Path(savepath, "MTDB_Modified.msh"))

    gmsh.write(msh_path)
    gmsh.finalize()

    meshes = ot.Meshes.from_gmsh(filename=msh_path, reindex=True, log=False)
    meshes.identify_subdomain(include_domain=True)

    _insert_T_init(meshes)

    meshes.save(savepath, overwrite=True)
