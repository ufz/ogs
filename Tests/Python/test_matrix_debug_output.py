import os
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
    assert (outdir / "square_1e2_ts_0_t_0.000000.vtu").exists()
    assert (outdir / "square_1e2_ts_4_t_1.000000.vtu").exists()


def check_local_matrix_output_files_exist(outdir, file_should_exist=True):
    logfile = outdir / "ogs_local_matrix.log"

    assert file_should_exist == logfile.exists()

    if file_should_exist:
        num_lines = sum(1 for line in logfile.open())

        # the number is a "reference value", it should not change unless the test is changed
        assert num_lines == 230


def check_global_matrix_output_files_exist(outdir, file_should_exist=True):
    # the numbers are reference values, they should not change unless the test is changed
    map_filename_to_num_lines = {
        "ogs_global_matrix_cnt_001_t_0.25_pcs_0_b.vec": 243,
        "ogs_global_matrix_cnt_001_t_0.25_pcs_0_Jac.mat": 3846,
        "ogs_global_matrix_cnt_001_t_0.25_pcs_0_K.mat": 2,  # global K and M matrices are empty
        "ogs_global_matrix_cnt_001_t_0.25_pcs_0_M.mat": 2,
        "ogs_global_matrix_cnt_002_t_0.25_pcs_0_b.vec": 243,
        "ogs_global_matrix_cnt_002_t_0.25_pcs_0_Jac.mat": 3846,
        "ogs_global_matrix_cnt_002_t_0.25_pcs_0_K.mat": 2,
        "ogs_global_matrix_cnt_002_t_0.25_pcs_0_M.mat": 2,
        "ogs_global_matrix_cnt_003_t_0.5_pcs_0_b.vec": 243,
        "ogs_global_matrix_cnt_003_t_0.5_pcs_0_Jac.mat": 3846,
        "ogs_global_matrix_cnt_003_t_0.5_pcs_0_K.mat": 2,
        "ogs_global_matrix_cnt_003_t_0.5_pcs_0_M.mat": 2,
        "ogs_global_matrix_cnt_004_t_0.75_pcs_0_b.vec": 243,
        "ogs_global_matrix_cnt_004_t_0.75_pcs_0_Jac.mat": 3846,
        "ogs_global_matrix_cnt_004_t_0.75_pcs_0_K.mat": 2,
        "ogs_global_matrix_cnt_004_t_0.75_pcs_0_M.mat": 2,
        "ogs_global_matrix_cnt_005_t_1_pcs_0_b.vec": 243,
        "ogs_global_matrix_cnt_005_t_1_pcs_0_Jac.mat": 3846,
        "ogs_global_matrix_cnt_005_t_1_pcs_0_K.mat": 2,
        "ogs_global_matrix_cnt_005_t_1_pcs_0_M.mat": 2,
    }

    for f, num_lines_expected in map_filename_to_num_lines.items():
        logfile = outdir / f

        print(logfile)
        assert file_should_exist == logfile.exists()

        if file_should_exist:
            num_lines_actual = sum(1 for line in logfile.open())

            assert num_lines_actual == num_lines_expected


@pytest.mark.parametrize(
    "prefix_parameter",
    [
        # first entry: how to set OGS_LOCAL_MAT_OUT_PREFIX
        # second entry: is the matrix logfile writable?
        # third entry: is the environment variable set properly?
        (False, True, False),  # do not set env var
        (True, True, True),  # set env var to an absolute path
        (f".{os.sep}", True, True),  # set env var to a relative path
        (
            f".{os.sep}dir_does_not_exist{os.sep}",
            False,
            True,
        ),  # set env var to a directory that does not exist
        ("", True, True),  # set env var to an empty string
    ],
)
@pytest.mark.parametrize(
    "elements_parameter",
    [
        # first entry: how to set OGS_LOCAL_MAT_OUT_ELEMENTS
        # second entry: is the environment variable set properly?
        (False, False),  # do not set env var
        ("2 3", True),  # set env var properly
        ("x", False),  # set env var to a wrong value
    ],
)
def test_local_matrix_debug_output(monkeypatch, prefix_parameter, elements_parameter):
    srcdir = Path(__file__).parent.parent.parent
    prjpath = srcdir / "Tests/Data/Mechanics/Linear/square_1e2.prj"

    assert "OGS_LOCAL_MAT_OUT_PREFIX" not in os.environ
    assert "OGS_LOCAL_MAT_OUT_ELEMENTS" not in os.environ

    prefix_setting, matrix_logfile_writable, prefix_expect_output = prefix_parameter
    elements_setting, elements_expect_output = elements_parameter

    # Empty env var not supported on Windows!
    if platform.system() == "Windows" and (
        prefix_setting == "" or elements_setting == ""
    ):
        return

    # whether a logfile is expected based on the environment variable
    # settings
    expect_matrix_logfile_output = prefix_expect_output and elements_expect_output

    # ... but the entire OGS simulation might be cancelled because the
    # log file is not writable
    expect_ogs_success = (not expect_matrix_logfile_output) or matrix_logfile_writable

    # https://docs.pytest.org/en/6.2.x/reference.html#pytest.MonkeyPatch
    with tempfile.TemporaryDirectory() as tmpdirname, monkeypatch.context() as ctx:
        tmpdirpath = Path(tmpdirname)
        # prepare environment
        if prefix_setting is False:
            pass
        elif prefix_setting is True:
            ctx.setenv("OGS_LOCAL_MAT_OUT_PREFIX", tmpdirname + os.sep)
        elif prefix_setting == "" or prefix_setting.startswith("."):
            ctx.setenv("OGS_LOCAL_MAT_OUT_PREFIX", prefix_setting)
            # change to the temporary directory such that log files will be written there.
            ctx.chdir(tmpdirname)

        if elements_setting is not False:
            ctx.setenv("OGS_LOCAL_MAT_OUT_ELEMENTS", elements_setting)

        # run and test
        check_local_matrix_output_files_exist(tmpdirpath, False)

        run(str(prjpath), str(tmpdirpath), expect_ogs_success)

        if expect_ogs_success:
            check_simulation_results_exist(tmpdirpath)

        check_local_matrix_output_files_exist(
            tmpdirpath, expect_matrix_logfile_output and expect_ogs_success
        )


@pytest.mark.parametrize(
    "prefix_parameter",
    [
        # first entry: how to set OGS_GLOBAL_MAT_OUT_PREFIX
        # second entry: is the matrix logfile writable?
        # third entry: is the environment variable set properly?
        (False, True, False),  # do not set env var
        (True, True, True),  # set env var to an absolute path
        (f".{os.sep}", True, True),  # set env var to a relative path
        (
            f".{os.sep}dir_does_not_exist{os.sep}",
            False,
            False,
        ),  # set env var to a directory that does not exist
        ("", True, True),  # set env var to an empty string
    ],
)
def test_global_matrix_debug_output(monkeypatch, prefix_parameter):
    srcdir = Path(__file__).parent.parent.parent
    prjpath = srcdir / "Tests/Data/Mechanics/Linear/square_1e2.prj"

    assert "OGS_GLOBAL_MAT_OUT_PREFIX" not in os.environ

    prefix_setting, expect_ogs_success, prefix_expect_output = prefix_parameter

    if platform.system() == "Windows" and prefix_setting == "":
        # Empty env var not supported on Windows!
        return

    # https://docs.pytest.org/en/6.2.x/reference.html#pytest.MonkeyPatch
    with tempfile.TemporaryDirectory() as tmpdirname, monkeypatch.context() as ctx:
        tmpdirpath = Path(tmpdirname)
        # prepare environment
        if prefix_setting is False:
            pass
        elif prefix_setting is True:
            ctx.setenv("OGS_GLOBAL_MAT_OUT_PREFIX", tmpdirname + os.sep)
        elif prefix_setting == "" or prefix_setting.startswith("."):
            ctx.setenv("OGS_GLOBAL_MAT_OUT_PREFIX", prefix_setting)
            # change to the temporary directory such that log files will be written there.
            ctx.chdir(tmpdirname)

        # run and test
        check_global_matrix_output_files_exist(tmpdirpath, False)

        run(str(prjpath), str(tmpdirpath), expect_ogs_success)

        if expect_ogs_success:
            check_simulation_results_exist(tmpdirpath)

        check_global_matrix_output_files_exist(tmpdirpath, prefix_expect_output)
