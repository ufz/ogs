import filecmp
import os
import shutil
import subprocess
from pathlib import Path

import pytest

if os.name != "posix":
    pytest.skip("scripts run need a unix shell and grep", allow_module_level=True)


@pytest.fixture(scope="module")
def doc_scripts_dir(ogs_src_dir):
    return ogs_src_dir / "scripts" / "doc"


@pytest.fixture(scope="module")
def res_dir():
    return Path(__file__).parent / "res" / "test_docu_scripts"


# since fixtures are module-scoped, tmpdir cannot be used and temporary directories have to be created "manually"
@pytest.fixture(scope="module")
def tmp_path(tmp_path_factory):
    return tmp_path_factory.mktemp("test_docu_scripts")


@pytest.fixture(scope="module")
def get_prj_params_output(doc_scripts_dir, res_dir, tmp_path):
    outfile = tmp_path / "get-project-params-output.txt"

    shutil.copyfile(
        res_dir / "CreateTwoPhaseComponentialFlowProcess.cpp.test_input.txt",
        tmp_path / "CreateTwoPhaseComponentialFlowProcess.cpp",
    )

    with outfile.open("w") as fh:
        # the functionality of this script is tested in this module
        subprocess.run(
            [doc_scripts_dir / "get-project-params.sh", tmp_path], check=True, stdout=fh
        )

    return outfile


@pytest.fixture(scope="module")
def normalize_param_cache_output(doc_scripts_dir, tmp_path, get_prj_params_output):
    outfile = tmp_path / "normalize-param-cache-output.txt"

    with get_prj_params_output.open() as fh_in, outfile.open("w") as fh_out:
        # the functionality of this script is tested in this module
        subprocess.run(
            [doc_scripts_dir / "normalize-param-cache.py"],
            check=True,
            stdin=fh_in,
            stdout=fh_out,
        )

    return outfile


def grep_count_fixed_string(string, path):
    res = subprocess.run(
        ["grep", "-c", "-F", string, path], check=True, capture_output=True, text=True
    )

    return int(res.stdout)


@pytest.mark.parametrize(
    ("string", "occurrences"),
    [
        ("getConfigSubtree(", 1),
        ("getConfigParameter<", 9),
        ("checkConfigParameter(", 1),
        (r"\ogs_file_param{", 11),
        (r"\ogs_file_param_special{", 5),
        # check that we didn't accidentally remove the double white space from the test input file
        (r"//!  \ogs", 2),
        (r"//!\ogs", 1),
        # check that we didn't accidentally remove the trailing white space from the test input file
        ("h2o_in_gas} ", 1),
    ],
)
def test_get_prj_params_output_contains_string(
    get_prj_params_output, string, occurrences
):
    assert grep_count_fixed_string(string, get_prj_params_output) == occurrences


def test_get_prj_params_output_full(get_prj_params_output, res_dir):
    assert filecmp.cmp(
        get_prj_params_output,
        res_dir / "get-project-params-ref-output.txt",
        shallow=False,
    )


@pytest.mark.parametrize(
    ("string", "occurrences"),
    [("getConfigSubtree", 1), ("getConfigParameter", 9), ("checkConfigParameter", 1)],
)
def test_normalize_param_cache_output_contains_string(
    normalize_param_cache_output, string, occurrences
):
    assert grep_count_fixed_string(string, normalize_param_cache_output) == occurrences


def test_normalize_param_cache_output_full(normalize_param_cache_output, res_dir):
    assert len(normalize_param_cache_output.open().readlines()) == 16

    assert filecmp.cmp(
        normalize_param_cache_output,
        res_dir / "normalize-param-cache-ref-output.txt",
        shallow=False,
    )
