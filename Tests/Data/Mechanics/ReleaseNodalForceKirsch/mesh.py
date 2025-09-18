# %%
from pathlib import Path

import gmsh
import ogstools as ot
import pyvista as pv


# %%
class MeshGenerator:
    """Mesh generator for the Kirsch problem"""

    def __init__(self, gmsh_model):
        self.gmsh_model = gmsh_model

        if not gmsh.isInitialized():
            gmsh.initialize()

        self.out_dir = None

        # Parameters
        L = 70  # Square size
        r = 6.5  # Radius of quarter circle
        y0 = -857
        lc = 1.0  # Mesh size (used for transfinite divisions)
        lc_o = 5.0

        # --- Define points ---
        # Quarter circle points
        p0 = gmsh.model.geo.addPoint(0, y0, 0, lc)
        p1 = gmsh.model.geo.addPoint(r, y0, 0, lc)
        p2 = gmsh.model.geo.addPoint(0, r + y0, 0, lc)

        # Square points
        p3 = gmsh.model.geo.addPoint(L, y0, 0, lc_o)
        p3a = gmsh.model.geo.addPoint(0.5 * L, y0, 0, lc)

        p4 = gmsh.model.geo.addPoint(L, L + y0, 0, lc_o)
        p5 = gmsh.model.geo.addPoint(0, L + y0, 0, lc_o)
        p5a = gmsh.model.geo.addPoint(0, 0.5 * L + y0, 0, lc)

        # --- Define curves ---
        arc1 = gmsh.model.geo.addCircleArc(p1, p0, p2)
        arc2 = gmsh.model.geo.addCircleArc(p3a, p0, p5a)

        l1 = gmsh.model.geo.addLine(p2, p5a)
        l1a = gmsh.model.geo.addLine(p5a, p5)
        l2 = gmsh.model.geo.addLine(p5, p4)
        l3 = gmsh.model.geo.addLine(p4, p3)
        l4a = gmsh.model.geo.addLine(p3, p3a)
        l4 = gmsh.model.geo.addLine(p3a, p1)

        self.tunnel_boundary = [arc1]

        self.left_boundary = [l1, l1a]
        self.right_boundary = [l3]

        self.bottom_boundary = [l4a, l4]
        self.top_boundary = [l2]

        # Line loop and surface
        loop1 = gmsh.model.geo.addCurveLoop([arc1, l1, -arc2, l4])
        surface1 = gmsh.model.geo.addPlaneSurface([loop1])

        loop2 = gmsh.model.geo.addCurveLoop([arc2, l1a, l2, l3, l4a])
        surface2 = gmsh.model.geo.addPlaneSurface([loop2])

        l5 = gmsh.model.geo.addLine(p0, p2)
        l6 = gmsh.model.geo.addLine(p1, p0)
        loop3 = gmsh.model.geo.addCurveLoop([l5, -arc1, l6])
        self.surface3 = gmsh.model.geo.addPlaneSurface([loop3])
        self.left_hole = l5
        self.bottom_hole = l6

        # --- Transfinite setup ---
        # Define transfinite curves (divisions based on geometry and smoothness)
        n_arc = 20  # number of divisions on arc
        n_side = 50  # for square edges

        gmsh.model.geo.mesh.setTransfiniteCurve(arc1, n_arc)
        gmsh.model.geo.mesh.setTransfiniteCurve(arc2, n_arc)
        gmsh.model.geo.mesh.setTransfiniteCurve(l1, n_side)
        gmsh.model.geo.mesh.setTransfiniteCurve(l4, n_side)

        ## Automatically set transfinite surface
        gmsh.model.geo.mesh.setTransfiniteSurface(surface1)
        gmsh.model.geo.mesh.setRecombine(2, surface1)  # Turn triangles into quads

        self.surfaces = [surface1, surface2]

        # --- Synchronize and add physical groups ---
        self.gmsh_model.geo.synchronize()

    def generate_meshes(self, out_dir="", order=1, with_cavern=False):
        self.out_dir = Path(out_dir)
        if not self.out_dir.exists():
            self.out_dir.mkdir(parents=True)

        gmsh.option.setNumber("Mesh.Algorithm", 5)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.5)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 10)
        gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 1)
        gmsh.option.setNumber("Mesh.Format", 1)

        # ot.cli().NodeReordering does not work for higher order elements generated from Gmsh
        # gmsh.option.setNumber("Mesh.ElementOrder", order)
        # gmsh.option.setNumber("Mesh.HighOrderOptimize", order)

        # --- Mesh generation ---
        self.gmsh_model.mesh.generate(2)

        # Optional GUI
        # if "-nopopup" not in sys.argv:
        #    gmsh.fltk.run()

        self.gmsh_model.removePhysicalGroups()
        bulk_mesh_name = "domain"
        self.gmsh_model.addPhysicalGroup(2, self.surfaces, 1, name=bulk_mesh_name)
        bulk_mesh_name = Path(self.out_dir, bulk_mesh_name + ".vtu")

        mesh_dim = 1
        self.gmsh_model.addPhysicalGroup(mesh_dim, self.tunnel_boundary, name="arc")
        left_boundary = self.left_boundary
        bottom_boundary = self.bottom_boundary
        if with_cavern:
            self.gmsh_model.addPhysicalGroup(2, [self.surface3], 2, name="cavern")
            left_boundary.append(self.left_hole)
            bottom_boundary.append(self.bottom_hole)

        self.gmsh_model.addPhysicalGroup(mesh_dim, left_boundary, name="left")
        self.gmsh_model.addPhysicalGroup(mesh_dim, bottom_boundary, name="bottom")
        self.gmsh_model.addPhysicalGroup(mesh_dim, self.top_boundary, name="top")
        self.gmsh_model.addPhysicalGroup(mesh_dim, self.right_boundary, name="right")

        # Set output format to MSH2
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        mesh_name = "domain.msh"
        mesh_file_name = Path(self.out_dir, mesh_name)
        gmsh.write(str(mesh_file_name))

        # `dim` is ommitted for converting all meshes
        meshes = ot.meshes_from_gmsh(
            filename=str(mesh_file_name), reindex=True, log=False
        )

        vtu_names = []
        for name, mesh in meshes.items():
            print(f"{name}: {mesh.n_cells} cells")
            if "physical_group_" in name:
                vtu_name = Path(
                    self.out_dir, f"{name.replace('physical_group_', '')}.vtu"
                )
            else:
                vtu_name = Path(self.out_dir, f"{name}.vtu")
            vtu_names.append(vtu_name)
            pv.save_meshio(vtu_name, mesh)

        ot.cli().NodeReordering("-i", bulk_mesh_name, "-o", bulk_mesh_name)
        if order == 2:
            files_to_process = {bulk_mesh_name, *vtu_names}
            for vtu_file_name in files_to_process:
                print(f"Create quadratic mesh for {vtu_file_name}")
                ot.cli().createQuadraticMesh("-i", vtu_file_name, "-o", vtu_file_name)

        ot.cli().identifySubdomains(
            "-m",
            bulk_mesh_name,
            f"-o {self.out_dir}/",
            "-f",
            "-s 1e-6",
            "--",
            *vtu_names,
        )
