import os
import tempfile

import ogs.simulator as sim


def test_simulator():
    arguments = [
        "",
        f"{os.path.abspath(os.path.dirname(__file__))}/../Data/Parabolic/LiquidFlow/Flux/cube_1e3_calculatesurfaceflux.prj",
        "-o " + tempfile.mkdtemp(),
    ]

    try:
        print("Python OpenGeoSys.init ...")
        assert sim.initialize(arguments) == 0
        print("Python OpenGeoSys.executeSimulation ...")
        assert sim.executeSimulation() == 0
    finally:
        print("Python OpenGeoSys.finalize() ...")
        sim.finalize()
