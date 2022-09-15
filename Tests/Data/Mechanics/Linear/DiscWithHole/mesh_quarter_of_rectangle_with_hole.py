# mesh quarter of rectangle with circular hole
# meshing strategy as FEM example in "Hoehere TM" by Kreissig and Benedix
# Dominik Kern
import numpy as np
import gmsh


def run(output_file, a=3.0, b=4.0, r=1.0, R=2.0, lc=0.5, Nx=4, Ny=5, NR=5, Nr=5, P=1.3):
    """
    lc ... characteristic length for meshing
    P ... Progression towards hole
    """

    # Before using any functions in the Python API, Gmsh must be initialized:
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("rectangle")

    # Dimensions
    dim1 = 1
    dim2 = 2

    # Auxiliary
    s45 = np.sin(np.pi / 4)
    c45 = np.cos(np.pi / 4)

    # Outer points (ccw)
    gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
    gmsh.model.geo.addPoint(r, 0, 0, lc, 2)
    gmsh.model.geo.addPoint(R, 0, 0, lc, 3)
    gmsh.model.geo.addPoint(a, 0, 0, lc, 4)
    gmsh.model.geo.addPoint(a, R, 0, lc, 5)
    gmsh.model.geo.addPoint(a, b, 0, lc, 6)
    gmsh.model.geo.addPoint(R, b, 0, lc, 7)
    gmsh.model.geo.addPoint(0, b, 0, lc, 8)
    gmsh.model.geo.addPoint(0, R, 0, lc, 9)
    gmsh.model.geo.addPoint(0, r, 0, lc, 10)
    gmsh.model.geo.addPoint(r * c45, r * s45, 0, lc, 11)

    # Inner point (only one)
    gmsh.model.geo.addPoint(R, R, 0, lc, 12)

    # Outer lines (ccw)
    gmsh.model.geo.addLine(2, 3, 1)
    gmsh.model.geo.addLine(3, 4, 2)
    gmsh.model.geo.addLine(4, 5, 3)
    gmsh.model.geo.addLine(5, 6, 4)
    gmsh.model.geo.addLine(6, 7, 5)
    gmsh.model.geo.addLine(7, 8, 6)
    gmsh.model.geo.addLine(8, 9, 7)
    gmsh.model.geo.addLine(9, 10, 8)
    gmsh.model.geo.addCircleArc(10, 1, 11, 9)
    gmsh.model.geo.addCircleArc(11, 1, 2, 10)

    # Inner lines (pointing inwards)
    gmsh.model.geo.addLine(3, 12, 11)
    gmsh.model.geo.addLine(5, 12, 12)
    gmsh.model.geo.addLine(7, 12, 13)
    gmsh.model.geo.addLine(9, 12, 14)
    gmsh.model.geo.addLine(11, 12, 15)

    # The third elementary entity is the surface. In order to define a surface
    # from the curves defined above, a curve loop has first to be defined (ccw).
    gmsh.model.geo.addCurveLoop([2, 3, 12, -11], 1)
    gmsh.model.geo.addCurveLoop([4, 5, 13, -12], 2)
    gmsh.model.geo.addCurveLoop([6, 7, 14, -13], 3)
    gmsh.model.geo.addCurveLoop([8, 9, 15, -14], 4)
    gmsh.model.geo.addCurveLoop([10, 1, 11, -15], 5)

    # Add plane surfaces defined by one or more curve loops.
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addPlaneSurface([2], 2)
    gmsh.model.geo.addPlaneSurface([3], 3)
    gmsh.model.geo.addPlaneSurface([4], 4)
    gmsh.model.geo.addPlaneSurface([5], 5)

    gmsh.model.geo.synchronize()

    #  Prepare structured grid
    gmsh.model.geo.mesh.setTransfiniteCurve(1, Nr, "Progression", P)
    gmsh.model.geo.mesh.setTransfiniteCurve(2, Nx)
    gmsh.model.geo.mesh.setTransfiniteCurve(3, NR)
    gmsh.model.geo.mesh.setTransfiniteCurve(4, Ny)
    gmsh.model.geo.mesh.setTransfiniteCurve(5, Nx)
    gmsh.model.geo.mesh.setTransfiniteCurve(6, NR)
    gmsh.model.geo.mesh.setTransfiniteCurve(7, Ny)
    gmsh.model.geo.mesh.setTransfiniteCurve(8, Nr, "Progression", -P)  # note direction
    gmsh.model.geo.mesh.setTransfiniteCurve(9, NR)
    gmsh.model.geo.mesh.setTransfiniteCurve(10, NR)
    gmsh.model.geo.mesh.setTransfiniteCurve(11, NR)
    gmsh.model.geo.mesh.setTransfiniteCurve(12, Nx)
    gmsh.model.geo.mesh.setTransfiniteCurve(13, Ny)
    gmsh.model.geo.mesh.setTransfiniteCurve(14, NR)
    gmsh.model.geo.mesh.setTransfiniteCurve(15, Nr, "Progression", P)

    gmsh.model.geo.mesh.setTransfiniteSurface(1, "Alternate")
    gmsh.model.geo.mesh.setTransfiniteSurface(2, "Alternate")
    gmsh.model.geo.mesh.setTransfiniteSurface(3, "Alternate")
    gmsh.model.geo.mesh.setTransfiniteSurface(4, "Alternate")
    gmsh.model.geo.mesh.setTransfiniteSurface(5, "Alternate")
    gmsh.model.geo.mesh.setRecombine(dim2, 1)
    gmsh.model.geo.mesh.setRecombine(dim2, 2)
    gmsh.model.geo.mesh.setRecombine(dim2, 3)
    gmsh.model.geo.mesh.setRecombine(dim2, 4)
    gmsh.model.geo.mesh.setRecombine(dim2, 5)

    gmsh.model.geo.synchronize()

    # Physical groups (only this gets saved to file per default)
    Bottom = gmsh.model.addPhysicalGroup(dim1, [1, 2])
    gmsh.model.setPhysicalName(dim1, Bottom, "Bottom")

    Right = gmsh.model.addPhysicalGroup(dim1, [3, 4])
    gmsh.model.setPhysicalName(dim1, Right, "Right")

    Top = gmsh.model.addPhysicalGroup(dim1, [5, 6])
    gmsh.model.setPhysicalName(dim1, Top, "Top")

    Left = gmsh.model.addPhysicalGroup(dim1, [7, 8])
    gmsh.model.setPhysicalName(dim1, Left, "Left")

    Hole = gmsh.model.addPhysicalGroup(dim1, [9, 10])
    gmsh.model.setPhysicalName(dim1, Hole, "Hole")

    Plate = gmsh.model.addPhysicalGroup(dim2, [1, 2, 3, 4, 5])
    gmsh.model.setPhysicalName(dim2, Plate, "Plate")

    gmsh.model.geo.synchronize()

    gmsh.model.mesh.generate(dim2)
    # gmsh.option.setNumber('Mesh.SecondOrderIncomplete', 1)   # serendipity elements
    # gmsh.model.mesh.setOrder(2)   # higher order elements (quadratic)
    gmsh.write(output_file)

    gmsh.finalize()
