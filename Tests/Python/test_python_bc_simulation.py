import tempfile
import os

import pytest
import ogs.simulator as sim


def test_HM_ground_equil_TaylorHood_Python():
    srcdir = os.path.join(os.path.dirname(__file__), "..", "..")
    prjpath = os.path.join(
        srcdir,
        "Tests/Data/HydroMechanics/GroundEquilibrium/simHM_ground_quadBCu_python.prj",
    )

    with tempfile.TemporaryDirectory() as tmpdirname:
        arguments = ["ogs", prjpath, "-o", tmpdirname]

        try:
            print("Python OpenGeoSys.init ...")
            assert 0 == sim.initialize(arguments)
            print("Python OpenGeoSys.executeSimulation ...")
            assert 0 == sim.executeSimulation()
        finally:
            print("Python OpenGeoSys.finalize() ...")
            sim.finalize()

        assert os.path.exists(
            os.path.join(
                tmpdirname, "simHM_ground_quadBCu_python_ts_10_t_1000000.000000.vtu"
            )
        )
