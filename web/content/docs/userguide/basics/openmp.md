+++
date = "2024-09-09T13:52:17+02:00"
title = "Shared-memory parallelism"
author = "Christoph Lehmann"
weight = 91

[menu]
  [menu.userguide]
    parent = "basics"
+++

OGS can use multiple CPU threads for linear solvers, element assembly, and parts
of VTK mesh I/O. These stages have separate controls: setting the solver thread
count does not automatically set the assembly or VTK thread count.

## Thread settings overview

Set these environment variables before starting OGS:

| Variable | Controls | Default when unset |
| --- | --- | --- |
| `OMP_NUM_THREADS` | OpenMP-enabled Eigen solver operations | OpenMP runtime default |
| `MKL_NUM_THREADS` | MKL operations, when enabled | MKL/runtime default |
| `OGS_ASM_THREADS` | Parallel element assembly in supported processes | 1 |
| `VTK_SMP_MAX_THREADS` | [VTK's multiprocessing](https://docs.vtk.org/en/latest/design_documents/smptools.html) | 1 for `TBB`- / `STDThread`-backend; OpenMP backend also uses `OMP_NUM_THREADS` |

The available parallelism depends on the libraries and backends included in
your OGS build. A thread count is a limit or request, not a promise that every
operation will use that many threads.

## Linear solvers

OGS normally builds with OpenMP support. For Eigen operations that support
parallel execution, set
[`OMP_NUM_THREADS`](https://www.openmp.org/spec-html/5.0/openmpse50.html):

```bash
OMP_NUM_THREADS=4 ogs model.prj
```

If the variable is unset, Eigen uses the OpenMP runtime's default thread count,
which can depend on the machine and execution environment. Set it explicitly
when running several simulations on the same machine.

Some Eigen solvers and preconditioners only run serially. Increasing the thread
count does not make them parallel. See
[Eigen's list of parallel operations](https://libeigen.gitlab.io/eigen/docs-5.0.1/TopicMultiThreading.html).
The start-up log's `Eigen use ... threads` message reports the configured count,
not whether every solver operation uses it.
For Eigen's conjugate gradient solver, using both triangular parts enables its
parallel implementation. Select this with the following option in the project's
Eigen linear solver configuration:

```xml
<triangular_matrix>LowerUpper</triangular_matrix>
```

See the [conjugate gradient documentation](https://libeigen.gitlab.io/eigen/docs-5.0.1/classEigen_1_1ConjugateGradient.html)
for the matrix requirements.

For builds using Intel MKL with its OpenMP threading layer, set
`MKL_NUM_THREADS` to control MKL solver operations:

```bash
MKL_NUM_THREADS=4 ogs model.prj
```

MKL thread settings take precedence over `OMP_NUM_THREADS` for MKL operations;
they do not configure element assembly or VTK. See
[Intel's threading controls](https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-windows/2024-1/techniques-to-set-the-number-of-threads.html).

The build option `OGS_EIGEN_PARALLEL_BACKEND` selects the Eigen
parallel backend; it is separate from the runtime thread settings above. See
[build configuration options]({{< relref "/docs/devguide/advanced/configuration-options/index.md#advanced-options" >}})
for MKL and OpenMP compatibility details.

## Element assembly

Set `OGS_ASM_THREADS` to a positive integer to request parallel assembly:

```bash
OGS_ASM_THREADS=4 ogs model.prj
```

Assembly defaults to one thread, even if `OMP_NUM_THREADS` is greater than one.
Conversely, setting `OGS_ASM_THREADS` does not change the solver's OpenMP default.

Parallel assembly requires OpenMP support and is available in these processes:

- `ComponentTransport`
- `HeatConduction`
- `HeatTransportBHE`
- `HydroMechanics`
- `RichardsMechanics`
- `SmallDeformation`
- `TH2M`
- `ThermoHydroMechanics`
- `ThermoRichardsFlow`
- `ThermoRichardsMechanics`

The start-up log reports the requested count as
`Threads used for ParallelVectorMatrixAssembler: ...`. During assembly, the
`Number of threads: ...` message reports the actual OpenMP thread count, which may
be constrained by the runtime.

## VTK mesh I/O

VTK uses its [SMP backend](https://docs.vtk.org/en/latest/design_documents/smptools.html) for some work during VTU output, such as computing
array ranges for XML metadata and the maximum cell size. This does not make
all mesh preparation, compression, and file writing parallel.

For the `TBB` and `STDThread` backends, OGS defaults to one VTK thread.
Override this with a positive integer:

```bash
VTK_SMP_MAX_THREADS=4 ogs model.prj
```

`OMP_NUM_THREADS` does not configure `TBB` or `STDThread`: there is no automatic
fallback from `VTK_SMP_MAX_THREADS` to `OMP_NUM_THREADS`.

For the **Sequential backend**, VTK SMP operations remain serial regardless of the
thread limit. Setting a larger limit does not enable a parallel backend.

For the **OpenMP backend**, OGS leaves the VTK thread configuration unchanged to
preserve the OpenMP settings also used by the solver. VTK normally follows the
OpenMP runtime settings, including `OMP_NUM_THREADS`. An explicit
`VTK_SMP_MAX_THREADS` can also change the OpenMP default used by the solver, so
these are not independent limits with this backend.

The backend is determined by your VTK build. The environment variable `VTK_SMP_BACKEND_IN_USE` can select
another backend only if it was compiled into VTK. A build or preset name alone
does not guarantee which VTK backend is available. See
[VTK SMP configuration documentation](https://vtk.org/doc/nightly/html/classvtkSMPTools.html).

## Combining settings and avoiding CPU contention

For four solver threads, four assembly threads, and serial VTK SMP work with
`TBB`/`STDThread`, use:

```bash
OMP_NUM_THREADS=4 OGS_ASM_THREADS=4 VTK_SMP_MAX_THREADS=1 ogs model.prj
```

For an OpenMP-backed VTK, omit `VTK_SMP_MAX_THREADS` to keep the same OpenMP
thread setting for the solver and VTK:

```bash
OMP_NUM_THREADS=4 OGS_ASM_THREADS=4 ogs model.prj
```

Assembly, solving, and output normally run in separate stages. When they use
the same OpenMP runtime, that runtime can reuse workers between stages. `TBB` and
`STDThread` use separate workers; assigning the same limit to both libraries
does not create a shared thread pool. The number of operating-system threads
can therefore exceed the number actively computing at a given time.

These settings do not impose a process-wide thread cap. When running several
OGS processes concurrently, choose thread counts that fit the available CPU
resources across all runs. Start with one thread per stage and compare run times
with representative models; small problems can be slower with more threads.

For notebooks run through CTest, the test's processor requirement must also
reflect its thread usage. See
[registering notebook tests]({{< relref "/docs/devguide/documentation/jupyter-docs/index.md#register-with-ctest" >}}).
