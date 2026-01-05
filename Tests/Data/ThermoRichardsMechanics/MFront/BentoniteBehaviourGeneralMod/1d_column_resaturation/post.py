#!/usr/bin/env python
# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause


from pathlib import Path

import numpy as np
from postprocessing_utils import (
    IterPvd,
    check_zero_and_constant_quantities,
    extract_data,
    is_checked_via_another_quantity,
    plot_profiles,
    plot_time_series,
    select_profiles,
)

outdir = Path("out")
outdir.mkdir(exist_ok=True)
with (outdir / ".gitignore").open("w") as fh:
    fh.write("*\n")


def run():
    ts, line_mesh, dfs_aggregation, profiles = extract_data(
        IterPvd("out/ogs/bentonite_column.pvd")
    )

    success, successfully_checked_quantities = check_zero_and_constant_quantities(
        dfs_aggregation
    )

    ignored_quantities = {
        "MaterialIDs",
        "bulk_node_ids",
        "bulk_element_ids",
        "cell_ids",
        "IntegrationPointMetaData",
        "OGS_VERSION",
    }

    all_quantities = set(dfs_aggregation["min"].columns) - ignored_quantities
    print("all quantities in the result meshes: ", sorted(all_quantities))

    quantities_not_checked_via_another_quantity = {
        q
        for q in sorted(all_quantities)
        if not is_checked_via_another_quantity(all_quantities, q)
    }
    assert quantities_not_checked_via_another_quantity

    interesting_quantities = quantities_not_checked_via_another_quantity - set(
        successfully_checked_quantities
    )

    plot_time_series(dfs_aggregation, interesting_quantities)

    ts_filtered, profiles_filtered = select_profiles(
        ts, profiles, ts[-1] * np.array([0, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0])
    )

    profile_quantities = [
        "sigma_total_avg[2]",
        "epsilon_avg[2]",
        "porosity_avg",
        "saturation_avg",
        "pressure_interpolated",
        "displacement[2]",
    ]

    plot_profiles(
        ts_filtered,
        line_mesh,
        profiles_filtered,
        profile_quantities,
        "out/plot_profiles",
    )

    assert success


if __name__ == "__main__":
    run()
