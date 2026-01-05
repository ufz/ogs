# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import tempfile
from pathlib import Path

import ogs


def test_pvtu2vtu_has_version():
    cli = ogs.cli
    cli.pvtu2vtu("--version")


def test_pvtu2vtu_based_on_output_of_parallel_LiquidFlow_on_cube_1e3():
    current_dir = Path(__file__).parent.resolve()
    path = current_dir / Path("../Data/EllipticPETSc/")
    original_vtu = path / Path("cube_1x1x1_hex_1e3.vtu")

    pvtu = path / Path("cube_1e3_ts_1_t_1_000000.pvtu")
    vtu = tempfile.mkdtemp() / Path("cube_1e3_ts_1_t_1_000000.vtu")
    cli = ogs.cli
    cli.pvtu2vtu(i=pvtu, m=original_vtu, o=vtu)

    cli.vtkdiff(vtu, a="pressure", b="Linear_1_to_minus1", abs=1e-12, rel=1e-12)


def test_pvtu2vtu_based_on_output_of_parallel_LiquidFlow_on_cube_1e3_neumann():
    current_dir = Path(__file__).parent.resolve()
    path = current_dir / Path("../Data/EllipticPETSc/")
    original_vtu = path / Path("cube_1x1x1_hex_1e3.vtu")

    pvtu = path / Path("cube_1e3_neumann_ts_1_t_1_000000.pvtu")
    vtu = tempfile.mkdtemp() / Path("cube_1e3_neumann_ts_1_t_1_000000.vtu")
    cli = ogs.cli
    cli.pvtu2vtu(i=pvtu, m=original_vtu, o=vtu)

    cli.vtkdiff(vtu, a="pressure", b="D1_left_front_N1_right", abs=1e-2, rel=1e-2)
