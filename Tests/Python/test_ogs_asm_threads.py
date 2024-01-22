import platform
import tempfile
from pathlib import Path

import ogs.simulator as sim
import pytest


def run(prjpath, outdir, expect_successful):
    arguments = ["ogs", prjpath, "-o", outdir]

    try:
        while True:
            print("Python OpenGeoSys.init ...")
            status = sim.initialize(arguments)

            if status != 0:
                break

            print("Python OpenGeoSys.executeSimulation ...")
            status = sim.executeSimulation()

            break
    finally:
        print("Python OpenGeoSys.finalize() ...")
        sim.finalize()

    status_expected = 0 if expect_successful else 1
    assert status == status_expected


def check_simulation_results_exist(outdir):
    assert (outdir / "anisotropic_thermal_expansion_ts_1_t_1000000.000000.vtu").exists()


@pytest.mark.parametrize(
    "asm_threads_parameter",
    [
        # first entry: how to set OGS_ASM_THREADS
        # second entry: should OGS run successfully?
        (False, True),  # do not set env var
        ("1", True),  # positive integer
        ("1  ", True),  # positive integer with trailing spaces
        ("   1", True),  # positive integer with leading spaces
        ("    1   ", True),  # positive integer surrounded by spaces
        #
        ("", False),  # set but empty
        (" ", False),  # blank string
        ("1.1", False),  # floating point number
        ("0", False),  # zero
        ("-1", False),  # minus one (has a special meaning sometimes)
        ("-13", False),  # another negative number
        ("4 5", False),  # list of integers
        ("4x", False),  # integer and garbage
        ("zxyf", False),  # garbage
        ("!23", False),  # garbage and integer
    ],
)
def test_ogs_asm_threads_env_var(monkeypatch, asm_threads_parameter):
    srcdir = Path(__file__).parent.parent.parent
    # fast running model with TRM process (OpenMP parallelized)
    prjpath = (
        srcdir
        / "Tests/Data/ThermoRichardsMechanics/anisotropic_thermal_expansion/aniso_expansion.prj"
    )

    asm_threads_setting, expect_ogs_success = asm_threads_parameter

    if platform.system() == "Windows" and asm_threads_setting == "":
        # Empty env var not supported on Windows!
        return

    # https://docs.pytest.org/en/6.2.x/reference.html#pytest.MonkeyPatch
    with tempfile.TemporaryDirectory() as tmpdirname, monkeypatch.context() as ctx:
        # prepare environment
        if asm_threads_setting is not False:
            ctx.setenv("OGS_ASM_THREADS", asm_threads_setting)

        ctx.chdir(tmpdirname)

        # run and test
        run(str(prjpath), tmpdirname, expect_ogs_success)

        if expect_ogs_success:
            check_simulation_results_exist(Path(tmpdirname))
