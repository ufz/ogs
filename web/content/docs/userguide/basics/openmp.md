+++
date = "2024-09-09T13:52:17+02:00"
title = "OpenMP parallelization"
author = "Christoph Lehmann"
weight = 91


[menu]
  [menu.userguide]
    parent = "basics"
+++

OpenGeoSys supports shared memory parallelization (multi threading) to some
extent: by default OpenGeoSys is built with OpenMP support and OpenMP-enabled
linear solvers will automatically make use of it.

The number of threads a linear solver uses can be controlled with the
environment variable
[`OMP_NUM_THREADS`](https://www.openmp.org/spec-html/5.0/openmpse50.html).

Note:

* Some linear solvers and some preconditioners from the Eigen library run with a
  single thread only, see
  [here](https://eigen.tuxfamily.org/dox/TopicMultiThreading.html).
* For Eigen's CG solver, multi-threading can be exploited using both triangular
  matrices (see
  [here](https://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient.html).
  This can be set in OGS using `<triangular_matrix>LowerUpper</triangular_matrix>`
  as linear solver option.
* When using distributed memory parallelization (MPI) it might not make sense to
  use OpenMP at the same time: the different threads might compete for the same
  CPU resources.
* There is an interference between the Intel MKL and OpenMP parallelization: If
  you use OGS together with the former, you might break the latter. In that case
  you might want to change the
  [`OGS_EIGEN_PARALLEL_BACKEND`](/docs/devguide/advanced/configuration-options/).

In addition, the assembly of some processes in OGS can be run
OpenMP-parallelized when setting the environment variable `OGS_ASM_THREADS` to
the number of threads to be used. At the moment (2024-09-09) that's supported by
the `HydroMechanics`, `TH2M`, `ThermoRichardsMechanics`, and
`ThermoRichardsFlow` processes.
