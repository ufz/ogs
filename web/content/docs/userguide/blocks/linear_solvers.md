+++
date = "2018-02-27T11:00:13+01:00"
title = "Linear solvers"
author = "Feliks Kiszkurno"
weight = 9
+++

<div class="note">

### Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

</div>

An important parameter defined for each iterative linear solver is the "error tolerance" (direct solvers don't use the tolerance settings).
Combined with `abstols` from section [Time loop](/docs/userguide/blocks/time_loop/) it defines the acceptable level of error in the obtained result.
Those two parameters are interconnected with each other. Setting either of them too tightly will result in an error message referring to problems, implying that a smaller time step needs to be set up. If the settings are too loose, results may be erroneous.
As the variables are interconnected, it is possible to avoid error when one of the tolerances is set too low, by setting the other one a bit too high.

Generally, the error tolerance of a linear solver should always be tighter than `abstols`.
Unlike `abstol`, "error tolerance" is not specific to "process variables" and is defined as one value.
For most cases value below $10^{-10}$ is recommended.

<!-- TODO: The text above is not completely clear. It needs re-wording and additional clarifications. -->

## Eigen

<!-- TODO: Add description of Eigen -->

## PETSc

See [Running OGS with MPI]({{< ref "parallel_computing_mpi">}})-page for details.
