import tempfile
from pathlib import Path

import ogs.simulator as sim


def test_HM_ground_equil_TaylorHood_Python():
    srcdir = Path(__file__).parent.parent.parent
    prjpath = (
        srcdir
        / "Tests/Data/HydroMechanics/GroundEquilibrium/simHM_ground_quadBCu_python.prj"
    )

    with tempfile.TemporaryDirectory() as tmpdirname:
        arguments = ["ogs", str(prjpath), "-o", tmpdirname]

        try:
            print("Python OpenGeoSys.init ...")
            assert sim.initialize(arguments) == 0
            print("Python OpenGeoSys.executeSimulation ...")
            assert sim.executeSimulation() == 0
        finally:
            print("Python OpenGeoSys.finalize() ...")
            sim.finalize()

        assert (
            Path(tmpdirname) / "simHM_ground_quadBCu_python_ts_10_t_1000000.000000.vtu"
        ).exists()
