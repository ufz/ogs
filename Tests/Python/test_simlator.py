import tempfile
import os

import pytest
import ogs.simulator as sim


def test_simulator():
    arguments = [
        "",
        f"{os.path.abspath(os.path.dirname(__file__))}/../Data/Parabolic/LiquidFlow/Flux/cube_1e3_calculatesurfaceflux.prj",
        "-o " + tempfile.mkdtemp(),
    ]

    print("Python OpenGeoSys.init ...")
    sim.initialize(arguments)
    print("Python OpenGeoSys.executeSimulation ...")
    sim.executeSimulation()
    print("Python OpenGeoSys.finalize() ...")
    sim.finalize()
