import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest
from ogs import OGSSimulator, mesh


def crossProduct(v, w):
    return np.array(
        [
            v[1] * w[2] - v[2] * w[1],
            v[2] * w[0] - v[0] * w[2],
            v[0] * w[1] - v[1] * w[0],
        ]
    )


def computeVectorFromPoints(a, b):
    return np.array([b[0] - a[0], b[1] - a[1], b[2] - a[2]])


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
                or points[cnt, 2] != 1.0
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
    for c in range(len(celltypes)):
        area = computeQuadArea(
            points[cells[4 * c + 0]],
            points[cells[4 * c + 1]],
            points[cells[4 * c + 2]],
            points[cells[4 * c + 3]],
        )
        if abs(area - 1.0 / 9.0) > sys.float_info.epsilon:
            print(
                "Python: error: area for cell "
                + str(c)
                + " is: "
                + str(area)
                + ", expected value is "
                + str(1.0 / 9.0)
            )


@pytest.mark.skipif("OGS_USE_PATH" in os.environ, reason="Works in wheel only.")
def test_simulator():
    import ogs.OGSMesh as OGSMesh  # noqa: F401, PLC0415
    from ogs import OGSSimulation  # noqa: PLC0415

    current_dir = Path(__file__).parent.resolve()
    arguments = [
        "",
        f"{current_dir}/../Data/Parabolic/LiquidFlow/Flux/3D/Hex/cuboid_1x1x1_hex_27_Dirichlet_Dirichlet_Python.prj",
        "-o " + tempfile.mkdtemp(),
    ]

    print("Python OGSSimulation() ...")
    sim = OGSSimulator.OGSSimulation(arguments)
    top_boundary_grid: mesh = sim.mesh("cuboid_1x1x1_hex_27_top_boundary")
    # compare grid point coordinates with expected point coordinates
    points = np.array(top_boundary_grid.getPointCoordinates())
    number_of_points = int(len(points) / 3)
    points.shape = (number_of_points, 3)
    comparePointCoordinates(points)
    # set top boundary conditions values for first time step
    bc_values_for_first_time_step = top_boundary_grid.dataArray(
        "values_set_from_python", "double"
    )
    bc_values_for_first_time_step[:] = np.ones(number_of_points) * 5e6

    assert sim.execute_time_step() == 0
    print("Python: sim.execute_time_step() done")
    top_boundary_grid: mesh = sim.mesh("cuboid_1x1x1_hex_27_top_boundary")

    (cells, celltypes) = top_boundary_grid.getCells()
    checkCells(cells, celltypes, points)

    # reset values of cell data array and get it back
    bc_values_for_second_time_step = top_boundary_grid.dataArray(
        "values_set_from_python", "double"
    )
    bc_values_for_second_time_step = np.ones(number_of_points) * 1e7
    read_back_bc_values = top_boundary_grid.dataArray(
        "values_set_from_python", "double"
    )

    # check lengths
    if len(read_back_bc_values) != len(bc_values_for_second_time_step):
        print(
            "Python: error: data array size mismatch: got "
            + str(len(read_back_bc_values))
            + ", expected "
            + str(len(bc_values_for_second_time_step))
        )
    comparison = read_back_bc_values == bc_values_for_second_time_step
    if not comparison.all():
        print("Python: error: data arrays contain different values")

    print("Python: sim.execute_time_step() ...")
    sim.execute_time_step()
    print("Python: sim.execute_time_step() done")
    print("Python: sim.close() ...")
    sim.close()
