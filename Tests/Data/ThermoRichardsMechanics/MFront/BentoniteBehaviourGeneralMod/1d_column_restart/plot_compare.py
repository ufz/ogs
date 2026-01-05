#!/usr/bin/env python
# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause


from pathlib import Path

import numpy as np
import postprocessing_utils as pu

outdir = Path("out")
outdir.mkdir(exist_ok=True)
with (outdir / ".gitignore").open("w") as fh:
    fh.write("*\n")


def run():
    ts, line_mesh, _dfs_aggregation, profiles = pu.extract_data(
        pu.IterPvd("out/ogs/bentonite_column.pvd")
    )

    ts_ref, _line_mesh_ref, _dfs_aggregation_ref, profiles_ref = pu.extract_data(
        pu.IterVtus(
            [
                "bentonite_column_ts_0_t_41242.9_sec.vtu",
                "bentonite_column_ts_1_t_44842.9_sec.vtu",
            ],
            [41242.9, 44842.9],
        ),
        line_mesh=line_mesh,
    )

    assert np.all(ts == ts_ref)

    profile_quantities = [
        "sigma_total_avg[2]",
        "epsilon_avg[2]",
        "porosity_avg",
        "saturation_avg",
        "pressure_interpolated",
        "displacement[2]",
        "eM",
        "re",
    ]

    pu.plot_profiles(
        ts,
        line_mesh,
        profiles,
        profile_quantities,
        "out/plot_profiles_compare",
        profiles_ref=profiles_ref,
    )


if __name__ == "__main__":
    run()
