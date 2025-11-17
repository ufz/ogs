+++
date = "2024-09-10T10:57:31+02:00"
title = "Debugging Options"
author = "Christoph Lehmann"
weight = 1
+++

If you encounter an error or unexpected behaviour during an OGS simulation there
are a couple of things you can do to investigate the error. These are, e.g...

## Increasing the Log Level

Running OGS with the [command line option](/docs/userguide/basics/cli-arguments/)
`-l debug` will make OGS print a lot of debug output, which might help locating
the cause of the error.

## Writing local and/or global Vectors and Matrices

When the environment variables `OGS_LOCAL_MAT_OUT_PREFIX="<some prefix>"` and
`OGS_LOCAL_MAT_OUT_ELEMENTS="0 1 5 16 123"` (or `"*"` to select all elements)
are set, OGS will write the local matrices and vectors it has assembled to
`<some prefix>ogs_local_matrix.log`. That file can be checked, if the assembled
Jacobian is singular, among others. `<some prefix>` can contain slashes, making
output to an arbitrary directory possible.

Likewise, output of global matrices can be requested via
`OGS_GLOBAL_MAT_OUT_PREFIX="<some prefix>"`. That will make OGS write a series
of files, one for each global matrix/vector for each non-linear iteration.

Local matrix output is available for all processes in OGS, as of now
(2024-09-10) global matrix output is available only for the `LARGE_DEFORMATION`,
`SMALL_DEFORMATION`, `ThermalTwoPhaseFlowWithPP`, `ThermoRichardsMechanics`,
`TH2M` and `ThermoHydroMechanics` processes, and only for serial (non-PETSc) runs.

## Generating more output

When the [`<time_loop/output/variables>`](https://doxygen.opengeosys.org/stable/d2/d5b/ogs_file_param__prj__time_loop__output__variables)
list in the OGS project file is empty, OGS will write all available data to its
output files.

Furthermore, you can set
[`<time_loop/output/output_iteration_results>`](https://doxygen.opengeosys.org/stable/dc/d50/ogs_file_param__prj__time_loop__output__output_iteration_results.html)
to true. Then OGS will generate output files after each non-linear iteration,
which might help to find the cause of a non-linear solver divergence.
This option is available for VTK output only (as of Sep. 2024).

When interpreting secondary variables it makes sense to enable
[`<time_loop/output/output_extrapolation_residuals>`](https://doxygen.opengeosys.org/stable/d1/d30/ogs_file_param__prj__time_loop__output__output_extrapolation_residuals.html).
The generated extrapolation residuals provide a measure of the extrapolation
errors of secondary variables.

## Enabling Floating Point Exceptions

If you are running OGS on Linux you can pass it the [command line option](/docs/userguide/basics/cli-arguments/)
`--enable-fpe` while running OGS in a debugger. That command line option will
raise a signal (i.e., make OGS crash) if there is a divide by zero or other
illegal floating point operation. The debugger will notice that situation and
stop at the respective source location. That, of course, requires an OGS build
with debug information (this is different from debug output above) enabled. See
[here](/docs/devguide/getting-started/build-configuration/#available-cmake-presets)
for more information.
