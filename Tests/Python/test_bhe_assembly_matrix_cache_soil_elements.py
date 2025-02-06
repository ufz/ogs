# ruff: noqa: E402

import pytest

ogstools = pytest.importorskip("ogstools")

from pathlib import Path

import pandas as pd
from ogs import cli
from ogstools.logparser.common_ogs_analyses import (
    analysis_time_step,
    fill_ogs_context,
    time_step_vs_iterations,
)
from ogstools.logparser.log_parser import parse_file


def parse_ogs_log(logfile):
    records = parse_file(logfile)

    dfa = pd.DataFrame(records)
    dfb = fill_ogs_context(dfa)
    dfc = time_step_vs_iterations(dfb)
    dfd = analysis_time_step(dfa)
    # Removing MPI_process (index=0) from result (all are 0) for serial log.
    dfd = dfd.loc[0]
    return dfd.join(dfc)


@pytest.mark.performance_test
def test_assembly_optimization(tmp_path, capfd, monkeypatch):
    srcdir = Path(__file__).parent.parent.parent
    testsrcdir = srcdir / "Tests/Data/Parabolic/T/3D_Beier_sandbox"
    cases = [("base", "beier_sandbox.prj"), ("linear", "beier_sandbox_linear.xml")]

    # run OGS and save logs
    for name, prj_or_patch in cases:
        outdir = tmp_path / name
        outdir.mkdir()
        with monkeypatch.context() as ctx:
            ctx.setenv("OMP_NUM_THREADS", "1")  # enforce a single thread
            status = cli.ogs(testsrcdir / prj_or_patch, o=outdir)

        captured = capfd.readouterr()
        with (outdir / "ogs.log").open("w") as fh:
            fh.write(captured.out)

        assert status == 0  # OGS run successful

    # parse logs
    mean_timings = {}
    for name, _ in cases:
        outdir = tmp_path / name
        df_timing = parse_ogs_log(outdir / "ogs.log")

        # remove timestep 0 (only output) and 1 (assembly is run in <linear> case, too)
        df_timing = df_timing.drop([0, 1])

        # make sure we will average over several measurements in the next step
        assert df_timing.shape[0] >= 5

        mean_timings[name] = df_timing.mean()

    # negative values in rel_timings mean <linear> is faster than the base case
    rel_timings = mean_timings["linear"] / mean_timings["base"] - 1
    rt = rel_timings
    mtb = mean_timings["base"]
    mtb_sol = mtb["time_step_solution_time"]

    ## actual checks
    # the assembly must be significantly faster in the <linear> case while not
    # sacrificing performance elsewhere

    # assembly must be at least twice as fast for the <linear> case
    assert rt["assembly_time"] < -0.5
    # assembly must have a significant runtime share
    assert mtb["assembly_time"] / mtb_sol > 0.05

    # no regression allowed, any possible regression in other parts of OGS must
    # be at least compensated by the accelerated assembly
    assert rt["time_step_solution_time"] < 0.0

    # some regression allowed due to measurement uncertainty
    # note: the linear solver time check is important since a messed up assembly
    # could in hypothetically lead to worse linear solver performance; the
    # output time limit is not that important
    assert rt["linear_solver_time"] < 0.03
    assert rt["output_time"] < 0.03
    assert rt["dirichlet_time"] < 0.03

    # just assert that the Dirichlet BC setting has a small runtime share
    assert mtb["dirichlet_time"] / mtb_sol < 0.01
