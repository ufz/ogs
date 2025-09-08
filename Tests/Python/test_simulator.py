import os
import tempfile
from pathlib import Path

import pytest


@pytest.mark.skipif("OGS_USE_PATH" in os.environ, reason="Works in wheel only.")
def test_simulator():
    from ogs import OGSSimulator  # noqa: PLC0415

    current_dir = Path(__file__).parent.resolve()
    arguments = [
        "",
        f"{current_dir}/../Data/Parabolic/LiquidFlow/Flux/cube_1e3_calculatesurfaceflux.prj",
        "-o " + tempfile.mkdtemp(),
    ]

    try:
        print("Python OpenGeoSys.init ...")
        sim = OGSSimulator.OGSSimulation(arguments)
        print("Python OpenGeoSys.executeSimulation ...")
        assert sim.executeSimulation() == 0
    finally:
        print("Python OpenGeoSys.finalize() ...")
        sim.finalize()
