import os
import tempfile
from pathlib import Path

import pytest


@pytest.mark.skipif("OGS_USE_PATH" in os.environ, reason="Works in wheel only.")
def test_simulator():
    from ogs import OGSSimulation  # noqa: PLC0415

    current_dir = Path(__file__).parent.resolve()
    arguments = [
        "",
        f"{current_dir}/../Data/Parabolic/LiquidFlow/Flux/cube_1e3_calculatesurfaceflux.prj",
        "-o " + tempfile.mkdtemp(),
    ]

    print("Python OGSSimulation() ...")
    sim = OGSSimulation(arguments)
    print("Python sim.execute_simulation() ...")
    assert sim.execute_simulation() == 0
    print("Python sim.close() ...")
    sim.close()
