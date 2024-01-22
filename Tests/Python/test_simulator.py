import tempfile
from pathlib import Path

import ogs.simulator as sim


def test_simulator():
    current_dir = Path(__file__).parent.resolve()
    arguments = [
        "",
        f"{current_dir}/../Data/Parabolic/LiquidFlow/Flux/cube_1e3_calculatesurfaceflux.prj",
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
