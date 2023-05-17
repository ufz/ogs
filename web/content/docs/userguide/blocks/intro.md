+++
title = "Overview"
author = "Feliks Kiszkurno"
weight = 1
+++

## Required blocks

Following blocks with valid content have to be provided in each project file, their content will be discussed in the following pages:

- [Meshes](/docs/userguide/blocks/meshes/)
- [Processes](/docs/userguide/blocks/processes/)
- [Media](/docs/userguide/blocks/media/)
- [Parameters](/docs/userguide/blocks/parameters/)
- [Time loop](/docs/userguide/blocks/time_loop/)
- [Process variables](/docs/userguide/blocks/process_variables/)
- [Nonlinear solvers](/docs/userguide/blocks/nonlinear_solvers/)
- [Linear solvers](/docs/userguide/blocks/linear_solvers/)

## Optional blocks

Blocks in this list are optional but may be convenient to use and are used in many benchmarks:

- [Curves](/docs/userguide/blocks/curves/)
- [Test definitions](/docs/userguide/blocks/test_definitions/)
- [Geometry](/docs/userguide/blocks/geometry/)

If you already have a project file but you are wondering how to modify it automatically or how it can be interacted with from external tools, please have a look at the [Working with project files](/docs/userguide/basics/working_with_project_files/) section.

## Minimum requirements for a valid project file

Generally, there is a certain degree of flexibility, when it comes to the internal structure of the project file.
Still, there are some exceptions.
For instance, in the "Time loop" section, the "set stepping" subsection has to be defined after "add process".
