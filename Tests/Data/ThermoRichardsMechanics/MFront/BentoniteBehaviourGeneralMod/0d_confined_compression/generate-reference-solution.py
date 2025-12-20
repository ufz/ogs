#!/usr/bin/env python
# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause


# invoke with
# LD_LIBRARY_PATH="<OGS_BUILD_DIR>/lib:$LD_LIBRARY_PATH" ./generate-reference-solution.py


from pathlib import Path

import mtest
import numpy as np
import tfel.tests  # noqa: F401 necessary for Boost.Python type converters
from ogs_and_mfront import (
    aggregate_ogs_results,
    generate_reference_meshes,
    plot,
    read_mtest_results,
    read_single_element_ogs_results,
    report_and_plot_summary,
    to_dataframes_with_mfront_names,
    write_test_definition_snippet,
)

out_dir = Path("out")
out_dir.mkdir(exist_ok=True)
with Path("out/.gitignore").open("w") as fh:
    fh.write("*\n")
(out_dir / "plot_mtest").mkdir(exist_ok=True)
(out_dir / "plot_mtest_vs_ogs").mkdir(exist_ok=True)

dt = 0.0125
mtest_output_file = "out/confined_compression_mtest.res"
ogs_output_file = "out/ogs/confined_compression.pvd"  # OGS must be run manually, e.g., run ctest and link ctest output dir
mesh_filename = "meshes/cube_1e0_1x1x1_hex8.vtu"

check_times = np.arange(0.5, 6 + 0.5 / 2, 0.5)

ref_mesh_filename_template = "confined_compression_ts_{ts}_t_{t:.6f}.vtu"


def exec_mtest():
    m = mtest.MTest()
    m.setBehaviour(
        "generic",
        "libOgsMFrontBehaviourBentoniteGeneralModForCTestsOnly",
        "BentoniteBehaviour",
    )
    m.setScalarInternalStateVariableInitialValue("e", 0.9)
    m.setScalarInternalStateVariableInitialValue("em", 0)
    m.setScalarInternalStateVariableInitialValue("eM", 0)
    m.setScalarInternalStateVariableInitialValue("SrM", 0)
    m.setScalarInternalStateVariableInitialValue("a_scan", 0)
    m.setScalarInternalStateVariableInitialValue("re", 0)

    # initial value of the gradients (strain+suction)
    m.setGradientsInitialValues([0, 0, 0, 0, 0, 0, -110e6])
    # initial values of the thermodynamic forces (total stress+Saturation)
    m.setThermodynamicForcesInitialValues([-2e6, -2e6, -2e6, 0, 0, 0, 0.345])

    m.setImposedThermodynamicForce("StressXX", {0: -2e6, 6: -4e6})
    # m.setImposedGradient('StrainXX', 0)
    m.setImposedGradient("StrainYY", 0)
    m.setImposedGradient("StrainZZ", 0)
    m.setImposedGradient("StrainXY", 0)
    m.setImposedGradient("StrainXZ", 0)
    m.setImposedGradient("StrainYZ", 0)
    m.setImposedGradient("LiquidPressure", -110e6)
    #
    m.setExternalStateVariable("AirPressure", 0)
    m.setExternalStateVariable("Temperature", 273.15)

    # output file
    m.setOutputFileName(mtest_output_file)
    m.setOutputFilePrecision(15)

    # // Imposing the time
    time = np.arange(
        0, 6 + dt / 2, dt
    ).tolist()  # conversion to list necessary to satisfy mtest Python bindings
    m.setTimes(time)

    print("## Executing MTest...")
    m.execute()


exec_mtest()

df_ref = read_mtest_results(mtest_output_file)

report_and_plot_summary(df_ref, file_prefix="out/plot_mtest/")

generate_reference_meshes(
    mesh_filename, df_ref, check_times, ref_mesh_filename_template
)

ts_ogs, recs_ogs = read_single_element_ogs_results(ogs_output_file)

recs_ogs_agg = aggregate_ogs_results(recs_ogs)

dfs_ogs_agg_mfront = to_dataframes_with_mfront_names(recs_ogs_agg, ts_ogs)

print("\n## Plotting comparison between OGS and MTest...")
plot(
    df_ref,
    dfs_ogs_agg_mfront["min"],
    dfs_ogs_agg_mfront["max"],
    labels=["mtest", "ogs min", "ogs max"],
    file_prefix="out/plot_mtest_vs_ogs/",
    quantities="intersection",
)

write_test_definition_snippet(
    df_ref.columns.values,
    dfs_ogs_agg_mfront["min"].columns.values,
    "confined_compression",
    "out/test_definition.xml",
)
