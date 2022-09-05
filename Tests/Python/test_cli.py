import OpenGeoSys
import tempfile
import os


def test_init():
    arguments = [
        "",
        f"{os.path.abspath(os.path.dirname(__file__))}/../Data/Parabolic/LiquidFlow/Flux/cube_1e3_calculatesurfaceflux.prj",
        "-o " + tempfile.mkdtemp(),
    ]

    print("Python OpenGeoSys.init ...")
    OpenGeoSys.initialize(arguments)
    print("Python OpenGeoSys.executeSimulation ...")
    OpenGeoSys.executeSimulation()
    print("Python OpenGeoSys.finalize() ...")
    OpenGeoSys.finalize()
