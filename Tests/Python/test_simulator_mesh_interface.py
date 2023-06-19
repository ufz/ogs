import tempfile
import os
import sys

import numpy as np

import pytest

import ogs.mesh as mesh
from ogs import simulator


def crossProduct(v, w):
    u = np.array(
        [
            v[1] * w[2] - v[2] * w[1],
            v[2] * w[0] - v[0] * w[2],
            v[0] * w[1] - v[1] * w[0],
        ]
    )
    return u


def computeVectorFromPoints(a, b):
    v = np.array([b[0] - a[0], b[1] - a[1], b[2] - a[2]])
    return v


def computeTriArea(a, b, c):
    v = computeVectorFromPoints(b, a)
    w = computeVectorFromPoints(c, a)
    return 0.5 * np.linalg.norm(crossProduct(v, w))


def computeQuadArea(point0, point1, point2, point3):
    return computeTriArea(point0, point1, point2) + computeTriArea(
        point1, point2, point3
    )


def comparePointCoordinates(points):
    cnt = 0
    for y in (0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0):
        for x in (0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0):
            if (
                float(x) != points[cnt, 0]
                or float(y) != points[cnt, 1]
                or 1.0 != points[cnt, 2]
            ):
                print(
                    "Python: error: expected point [["
                    + str(x)
                    + " "
                    + str(y)
                    + " 1.0]], point: "
                    + str(points[cnt])
                )
            cnt += 1


# check cell types - should be quads
def checkCells(cells, celltypes, points):
    for celltype in celltypes:
        if celltype != 9:
            print("Python: error: cell isn't a quad")

    # structured mesh with equal size cells
    for c in range(0, len(celltypes)):
        area = computeQuadArea(
            points[cells[5 * c + 1]],
            points[cells[5 * c + 2]],
            points[cells[5 * c + 3]],
            points[cells[5 * c + 4]],
        )
        if abs(area - 1.0 / 9.0) > sys.float_info.epsilon:
            print(
                "Python: error: area for cell "
                + str(c)
                + " is: "
                + str(area)
                + ", expected value is "
                + str(1.0 / 1.9)
            )


def test_simulator():
    arguments = [
        "",
        f"{os.path.abspath(os.path.dirname(__file__))}/../Data/Parabolic/LiquidFlow/Flux/3D/Hex/cuboid_1x1x1_hex_27_Dirichlet_Dirichlet_Python.prj",
        "-o " + tempfile.mkdtemp(),
    ]

    try:
        print("Python: OpenGeoSys.init ...")
        assert 0 == simulator.initialize(arguments)

        top_boundary_grid = simulator.getMesh("cuboid_1x1x1_hex_27_top_boundary")
        # compare grid point coordinates with expected point coordinates
        points = np.array(top_boundary_grid.getPointCoordinates())
        number_of_points = int(len(points) / 3)
        points.shape = (number_of_points, 3)
        comparePointCoordinates(points)
        # set top boundary conditions values for first time step
        bc_values_for_first_time_step = np.ones(number_of_points) * 5e6
        top_boundary_grid.setPointDataArray(
            "values_set_from_python", bc_values_for_first_time_step, 1
        )

        print("Python: OpenGeoSys.executeSimulation ...")
        assert 0 == simulator.executeTimeStep()
        print("Python: simulator.executeTimeStep() done")
        top_boundary_grid = simulator.getMesh("cuboid_1x1x1_hex_27_top_boundary")

        (cells, celltypes) = top_boundary_grid.getCells()
        checkCells(cells, celltypes, points)

        # set values of cell data array and get it back
        bc_values_for_second_time_step = np.ones(number_of_points) * 1e7
        top_boundary_grid.setPointDataArray(
            "values_set_from_python", bc_values_for_second_time_step, 1
        )
        read_back_bc_values = top_boundary_grid.getPointDataArray(
            "values_set_from_python", 1
        )
        # check lengths
        if len(read_back_bc_values) != len(bc_values_for_second_time_step):
            print(
                "Python: error: data array size mismatch: got "
                + str(len(new_bc_values))
                + ", expected "
                + str(len(bc_values))
            )
        comparison = read_back_bc_values == bc_values_for_second_time_step
        if not comparison.all():
            print("Python: error: data arrays contain different values")

        print("Python: simulator.executeTimeStep() ...")
        simulator.executeTimeStep()
        print("Python: simulator.executeTimeStep() done")
        print("Python: update OGS done")

    finally:
        print("Python: OpenGeoSys.finalize() ...")
        simulator.finalize()
