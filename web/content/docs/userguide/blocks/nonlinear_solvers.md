+++
date = "2018-02-27T11:00:13+01:00"
title = "Nonlinear solvers"
author = "Feliks Kiszkurno"
weight = 8
+++

<div class="note">

### Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

</div>

Following non-linear solvers are available in OpenGeoSys:

- [Newton](/docs/userguide/blocks/nonlinear_solvers/#newton)
- [Picard](/docs/userguide/blocks/nonlinear_solvers/#picard)

## Newton

<!-- TODO: add content -->

The nonlinear solver of "Newton" type is an implementation of the Newton-Raphson method.
The basic definition of the non-linear solver with "Newton" follows this template:

```xml
<nonlinear_solver>
    <name>basic_newton</name>
    <type>Newton</type>
    <max_iter>10</max_iter>
    <linear_solver>linear_solver</linear_solver>
</nonlinear_solver>
```

### Configuration parameters

- `name`: Unique identifier for the nonlinear solver.
- `type`: Must be `Newton` for this solver type.
- `max_iter`: Maximum number of nonlinear iterations.
- `linear_solver`: Reference to a defined linear solver.
- `recompute_jacobian`: (Optional, default: `1`) Frequency of Jacobian recalculation. `1` means recalculate every iteration.
- `damping`: (Optional, default: `1.0`) Damping factor for the Newton update.
- `damping_reduction`: (Optional) Factor by which the damping is reduced when convergence is slow.
- `tikhonov`: (Optional) Tikhonov regularisation configuration (see below).

### Tikhonov regularisation

Tikhonov regularisation can be used to improve convergence for ill-conditioned systems by adding a value to the diagonal of the Jacobian matrix:

```xml
<nonlinear_solver>
    <name>newton_with_tikhonov</name>
    <type>Newton</type>
    <max_iter>100</max_iter>
    <linear_solver>linear_solver</linear_solver>
    <tikhonov>
        <lambda>1e-14</lambda>
        <starting_iteration>11</starting_iteration>
    </tikhonov>
</nonlinear_solver>
```

The `tikhonov` element contains:

- `lambda`: The regularisation parameter (typically 1e-20 to 1e-10).
- `starting_iteration`: The iteration from which regularisation is applied (default: `0`).

## Picard

<!-- TODO: add content -->

The nonlinear solver of "Picard" type is an implementation of the Picard-Iteration method.
The basic definition of the non-linear solver with "Picard" follows this template:

```xml
<nonlinear_solver>
    <name>basic_picard</name>
    <type>Picard</type>
    <max_iter>100</max_iter>
    <linear_solver>linear_solver</linear_solver>
</nonlinear_solver>
```
