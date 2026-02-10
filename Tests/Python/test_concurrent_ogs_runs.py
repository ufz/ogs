# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import random
import threading
import time
from pathlib import Path

import ogs
import pytest
from ogs.OGSSimulator import OGSSimulation


def create(prj_path: Path, output_dir: Path):
    output_dir.mkdir(exist_ok=True)

    arguments = [
        "",
        str(prj_path),
        "-o",
        str(output_dir),
    ]

    return OGSSimulation(arguments)


def create_and_wait_and_close(prj_path: Path, output_dir: Path):
    sim = create(prj_path, output_dir)
    assert sim.initialized
    assert sim.status == 0
    time.sleep(random.uniform(0.05, 0.1))
    print(f"== closing simulation {output_dir.name}")
    sim.close()
    assert sim.status == 0
    assert not sim.initialized


@pytest.mark.xfail(
    ogs.OGS_USE_PETSC == "ON",
    reason="requires multiple MPI setups (not possible in OGS)",
    # code will trigger OGS_FATAL, which throws a std::runtime_error, which in
    # turn is converted to a RuntimeError by pybind11.
    raises=RuntimeError,
    strict=True,
)
def test_subsequent_runs(tmp_path, ogs_src_dir):
    prj_path = (
        ogs_src_dir
        / "Tests/Data/Parabolic/LiquidFlow/LineDirichletNeumannBC/line_dirichlet_neumannBC.prj"
    )

    sim1 = create(prj_path, tmp_path / "1")
    sim1.close()
    assert sim1.status == 0

    sim2 = create(prj_path, tmp_path / "2")
    sim2.close()
    assert sim2.status == 0


@pytest.mark.xfail(
    ogs.OGS_USE_PETSC == "ON",
    reason="requires multiple MPI setups (not possible in OGS)",
    # code will trigger OGS_FATAL, which throws a std::runtime_error, which in
    # turn is converted to a RuntimeError by pybind11.
    raises=RuntimeError,
    strict=True,
)
@pytest.mark.parametrize("num_simulations", [2, 3, 4, 8])
def test_parallel_runs(tmp_path, ogs_src_dir, num_simulations):
    prj_path = (
        ogs_src_dir
        / "Tests/Data/Elliptic/square_1x1_SteadyStateDiffusion/square_1e0.prj"
    )

    # many simulations are open at the same time
    sims = [create(prj_path, tmp_path / str(i)) for i in range(num_simulations)]
    for s in sims:
        s.close()

    for s in sims:
        assert s.status == 0


# Necessary to communicate exceptions raised in the thread outside
class ThreadStoringException(threading.Thread):
    def __init__(self, *, target, args):
        super().__init__(target=self._my_wrapped_target, args=args)
        self._my_target = target
        self.exception = None

    def _my_wrapped_target(self, *args):
        try:
            self._my_target(*args)
        except Exception as e:
            self.exception = e


@pytest.mark.xfail(
    ogs.OGS_USE_PETSC == "ON",
    reason="requires multiple MPI setups (not possible in OGS)",
    raises=AssertionError,
    strict=True,
)
@pytest.mark.parametrize("num_simulations", [2, 3, 4, 8])
def test_parallel_runs_threads(tmp_path, ogs_src_dir, num_simulations):
    prj_path = (
        ogs_src_dir
        / "Tests/Data/Elliptic/square_1x1_SteadyStateDiffusion/square_1e0.prj"
    )

    threads = [
        ThreadStoringException(
            target=create_and_wait_and_close, args=(prj_path, tmp_path / str(i))
        )
        for i in range(num_simulations)
    ]

    for t in threads:
        t.start()

    for t in threads:
        t.join()

    for t in threads:
        assert t.exception is None
