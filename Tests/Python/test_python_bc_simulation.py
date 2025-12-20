# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import tempfile
from pathlib import Path

from ogs import cli


def test_HM_ground_equil_TaylorHood_Python():
    srcdir = Path(__file__).parent.parent.parent
    prjpath = (
        srcdir
        / "Tests/Data/HydroMechanics/GroundEquilibrium/simHM_ground_quadBCu_python.prj"
    )

    with tempfile.TemporaryDirectory() as tmpdirname:
        status = cli.ogs(str(prjpath), o=tmpdirname)
        assert status == 0
        assert (
            Path(tmpdirname) / "simHM_ground_quadBCu_python_ts_10_t_1000000.000000.vtu"
        ).exists()
