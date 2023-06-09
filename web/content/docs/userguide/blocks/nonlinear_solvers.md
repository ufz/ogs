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
The basic definition of the non-linear solver follows this template:

```xml
<nonlinear_solver>
    <name>basic_newton</name>
    <type>Newton</type>
    <max_iter>100</max_iter>
    <linear_solver>linear_solver</linear_solver>
</nonlinear_solver>
```

## Picard

<!-- TODO: add content -->
