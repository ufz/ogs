# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import multiprocessing as mp

import ogs
import pytest
from ogs.OGSSimulator import OGSSimulation


def init_ogs(*arguments, expect_initialized):
    str_args = [""] + [str(arg) for arg in arguments]
    sim = OGSSimulation(str_args)
    assert sim.initialized == expect_initialized
    assert sim.status == 0
    return sim


def run_ogs(sim):
    assert sim.execute_simulation() == 0
    assert sim.initialized
    assert sim.status == 0


def close_ogs(sim):
    sim.close()
    assert not sim.initialized
    assert sim.status == 0


def init_and_close_ogs(*arguments, expect_initialized):
    sim = init_ogs(*arguments, expect_initialized=expect_initialized)
    close_ogs(sim)
    close_ogs(sim)  # closing a second time should not be a problem


def init_and_run_and_close_ogs(*arguments, expect_initialized=True):
    sim = init_ogs(*arguments, expect_initialized=expect_initialized)
    run_ogs(sim)
    close_ogs(sim)
    close_ogs(sim)  # closing a second time should not be a problem


def run_ogs_after_close(*arguments):
    sim = init_ogs(*arguments, expect_initialized=True)
    run_ogs(sim)
    close_ogs(sim)
    run_ogs(sim)


def maybe_run_in_separate_process(fct, *args, **kwargs):
    if ogs.OGS_USE_PETSC == "ON":
        # We need separate processes, otherwise we might run into this error:
        # The MPI_Init() function was called after MPI_FINALIZE was invoked.
        ctx = mp.get_context("spawn")
        p = ctx.Process(target=fct, args=args, kwargs=kwargs)
        p.start()
        p.join()
        if p.exitcode != 0:
            # Using RuntimeError here, because this error is triggered in these
            # unit tests by code that in the non-petsc builds would raise a
            # RuntimeError.
            # Which error condition? -> prj file does not exist, see
            # test_prj_does_not_exist()
            msg = f"Exit code of the spawned process was {p.exitcode} (expected 0)."
            raise RuntimeError(msg)
    else:
        fct(*args, **kwargs)


def test_simulation(tmp_path, ogs_src_dir):
    arguments = [
        ogs_src_dir
        / "Tests/Data/Parabolic/LiquidFlow/Flux/cube_1e3_calculatesurfaceflux.prj",
        "-o",
        tmp_path,
    ]

    maybe_run_in_separate_process(init_and_run_and_close_ogs, *arguments)

    assert (
        tmp_path
        / "LF_cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3_complete_surface.pvd"
    ).exists()
    assert (
        tmp_path / "LF_cube_1e3_calculatesurfaceflux_cube_1x1x1_hex_1e3.pvd"
    ).exists()


@pytest.mark.parametrize("arg", ["--help", "-h", "--version"])
def test_help_version(arg):
    maybe_run_in_separate_process(init_and_close_ogs, arg, expect_initialized=False)


@pytest.mark.xfail(reason="Prj file does not exist", raises=RuntimeError, strict=True)
def test_prj_does_not_exist(tmp_path):
    maybe_run_in_separate_process(init_and_run_and_close_ogs, [tmp_path / "some.prj"])


@pytest.mark.xfail(
    reason="OGS run although not initialized", raises=RuntimeError, strict=True
)
@pytest.mark.parametrize("arg", ["--help", "-h", "--version"])
def test_run_uninitialized(arg):
    maybe_run_in_separate_process(
        init_and_run_and_close_ogs, arg, expect_initialized=False
    )


@pytest.mark.xfail(reason="OGS run after close()", raises=RuntimeError, strict=True)
def test_run_after_close(tmp_path, ogs_src_dir):
    arguments = [
        ogs_src_dir
        / "Tests/Data/Parabolic/LiquidFlow/Flux/cube_1e3_calculatesurfaceflux.prj",
        "-o",
        tmp_path,
    ]

    maybe_run_in_separate_process(run_ogs_after_close, *arguments)
