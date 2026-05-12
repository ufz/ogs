+++
title = "Parallel Chemistry with OpenMP"
date = "2026-02-06"
author = "OpenGeoSys Community"
weight = 5
+++

## Overview

In reactive transport simulations, geochemical speciation calculations
performed by PHREEQC typically consume >95% of the total simulation time.
Because these calculations are embarrassingly parallel --- each element's
chemistry can be computed independently --- OGS can distribute them across
multiple OpenMP threads for a significant speedup.

Two usage modes are available:

| | OpenMP-only | Hybrid MPI+OpenMP |
|---|---|---|
| Build preset | `release` / `debug` | `release-petsc` / `debug-petsc` |
| Transport | single process | distributed across MPI ranks |
| Chemistry | multi-threaded | multi-threaded within each rank |
| Use case | desktop / workstation | clusters, large problems |

## Mode 1: OpenMP-only (no MPI)

This is the simplest way to speed up reactive transport. Transport runs on a
single CPU while chemistry is parallelized across OpenMP threads. No MPI
installation is required.

### Build

See the [developer guide — build configuration](/docs/devguide/getting-started/build-configuration/).

### Run

Set the number of chemistry threads via the environment variable or the
project file (see [Configuration](#configuration) below), then run OGS
normally:

```bash
OGS_CHEM_THREADS=8 ogs project.prj
```

Or with the project-file option:

```xml
<chemical_system chemical_solver="Phreeqc">
    <mesh>ReactiveDomain</mesh>
    <database>phreeqc.dat</database>
    <chemistry_threads>8</chemistry_threads>
    <!-- ... other settings ... -->
</chemical_system>
```

```bash
ogs project.prj
```

## Mode 2: Hybrid MPI+OpenMP

For large problems or cluster environments, combine MPI domain decomposition
for transport with OpenMP threads for chemistry within each rank:

- MPI handles transport with optimal communication patterns (typically 2-4 ranks)
- OpenMP threads parallelize chemistry within each rank (8-16+ threads)
- Total chemistry parallelism = MPI ranks x threads per rank

### Build

See the [developer guide — build configuration for MPI and PETSc](/docs/devguide/getting-started/build-configuration-for-mpi-and-petsc/).

### Run

```bash
OGS_CHEM_THREADS=8 mpirun -np 2 ogs project.prj
```

**Total cores used:** 2 ranks x 8 threads = 16 cores during the chemistry phase.

Or equivalently, with `<chemistry_threads>8</chemistry_threads>` in the
project file:

```bash
mpirun -np 2 ogs project.prj
```

### Recommended settings

> **Note:** OpenMP threads use shared memory and therefore cannot span
> multiple cluster nodes. Each MPI rank and all its chemistry threads
> must fit within a **single node**.

**Single-node examples** (all ranks on one node):

| Cores per node | MPI ranks | Chemistry threads/rank | Total cores used |
|----------------|-----------|------------------------|------------------|
| 16             | 2         | 8                      | 16               |
| 32             | 2         | 16                     | 32               |
| 64             | 4         | 16                     | 64               |

**Multi-node jobs:** Use one (or a small number of) MPI rank(s) per node
and fill the remaining cores on each node with chemistry threads. The
key constraint is:

```text
ranks_per_node × threads_per_rank ≤ cores_per_node
```

For example, on a cluster with 32-core nodes and a 4-node job:

| Nodes | Ranks/node | Threads/rank | Total MPI ranks | Total chemistry parallelism |
|-------|-----------|--------------|-----------------|----------------------------|
| 4     | 2          | 16           | 8               | 8 × 16 = 128               |
| 4     | 4          | 8            | 16              | 16 × 8 = 128               |

Chemistry throughput scales with the total parallelism (ranks × threads),
so both knobs contribute equally to the chemistry phase. The two rows
above therefore reach the same chemistry throughput. The trade-off is on
the transport side and on resource usage: more MPI ranks improve
transport scalability but increase communication overhead and per-node
memory (each rank loads its own PHREEQC database), while more threads
per rank reduce both at fixed total parallelism but cannot span nodes.
The optimal balance depends on the ratio of transport to chemistry cost
for your specific problem.

**General guidelines:**

- Keep `ranks_per_node × threads_per_rank = cores_per_node` to fully utilise each node
- Use your MPI launcher's binding options (e.g. `--ntasks-per-node`,
  `--cpus-per-task` in SLURM) to ensure each rank's threads stay on the same node
- Benchmark with a few time steps to find the optimal rank/thread ratio

### Example SLURM job scripts

Single node (2 ranks × 8 threads = 16 cores):

```bash
#!/bin/bash
#SBATCH --job-name=reactive_transport
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2     # 2 MPI ranks on this node
#SBATCH --cpus-per-task=8       # 8 OpenMP threads per rank
#SBATCH --mem=16G
#SBATCH --time=24:00:00

export OGS_CHEM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores

srun ogs your_project.prj -o output
```

Multi-node (4 nodes × 2 ranks/node × 16 threads/rank = 128 chemistry threads):

```bash
#!/bin/bash
#SBATCH --job-name=reactive_transport
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2     # 2 MPI ranks per node
#SBATCH --cpus-per-task=16      # 16 OpenMP threads per rank
#SBATCH --mem=64G
#SBATCH --time=24:00:00

export OGS_CHEM_THREADS=16
export OMP_PROC_BIND=close
export OMP_PLACES=cores

srun ogs your_project.prj -o output
```

Using `--ntasks-per-node` and `--cpus-per-task` together tells SLURM to
pin each rank to a fixed set of cores on one node, so OpenMP threads
never migrate across node boundaries.

## Configuration

There are two ways to set the number of chemistry threads. Both work
identically in OpenMP-only and hybrid MPI+OpenMP mode.

### Environment variable

Set `OGS_CHEM_THREADS` before running OGS:

```bash
export OGS_CHEM_THREADS=8
```

If `OGS_CHEM_THREADS` is not set, the number of assembly threads
(`OGS_ASM_THREADS`) is used as the default. If neither variable is set,
chemistry runs sequentially (single-threaded).

### Project file parameter

Add `<chemistry_threads>` inside `<chemical_system>`:

```xml
<chemical_system chemical_solver="Phreeqc">
    <chemistry_threads>8</chemistry_threads>
    <!-- ... -->
</chemical_system>
```

**Priority:** The project file setting overrides the environment variable.

## Memory Considerations

Each chemistry thread creates an independent PHREEQC instance, which
loads the thermodynamic database into memory. Typical memory usage:

- ~10-50 MB per PHREEQC instance (depends on database size)
- 8 threads = 80-400 MB additional memory per OGS process (in hybrid
  MPI+OpenMP mode, this cost is incurred on each MPI rank)

For memory-constrained systems, balance thread count against available RAM.

## Performance Tips

1. **Profile first**: Run with different thread counts to find the optimal
   configuration for your specific problem
2. **Avoid over-committing CPUs**: Ensure `threads <= physical cores`
   (OpenMP-only) or `ranks x threads <= physical cores` (hybrid)
3. **Consider problem size**: Small problems may not benefit from many threads
   due to overhead
4. **Monitor memory**: Use `htop` or similar to check memory usage when
   increasing threads

## Troubleshooting

### Chemistry threads not working

- Verify OpenMP is enabled: check CMake output for `OpenMP: YES`
- Check environment variable: `echo $OGS_CHEM_THREADS`
- Ensure `use_stream_for_data_exchange` is `true` (default) in your project file

### Slower than expected

- Reduce threads if memory-bound (check swap usage with `free -h`)
- Ensure no CPU over-commitment
- Check if problem is too small for parallel overhead

### Different results with threads

Results should be identical within floating-point precision (typically < 1e-10
relative difference). If results differ significantly, please report this as a
bug.

### Debugging parallel issues

Run with a single thread first to verify the sequential case works:

```bash
# OpenMP-only
OGS_CHEM_THREADS=1 ogs project.prj

# Hybrid
OGS_CHEM_THREADS=1 mpirun -np 2 ogs project.prj
```

Then gradually increase threads to identify where issues occur.

## Technical Details

The parallel chemistry implementation uses:

- **Multiple PHREEQC instances**: Each thread has its own isolated PHREEQC instance
- **Dynamic scheduling**: `schedule(dynamic)` handles load imbalance between chemical systems
- **Thread-safe design**: Each chemical system is processed by exactly one thread
- **Exception handling**: Failed chemical systems are collected into a
  mutex-guarded list inside the parallel region; after it ends, any failure is
  raised collectively across MPI ranks via `BaseLib::MPI::allRanksThrowOrNone`

The parallel execution follows this pattern:

1. Generate input strings for all chemical systems (sequential)
2. Execute PHREEQC calculations in parallel with OpenMP
3. Parse output strings and update shared state (sequential)
