+++
date = "2018-02-23T15:28:13+01:00"
title = "Missing content in documentation"
author = "Feliks Kiszkurno"
weight = 1027

+++

This list was obtained using grep and added here to provide an  overview of what sections are missing the documentation. In the future it should be generated automatically.

It can be obtained by running script `todo-check.sh` from `web/content/docs/userguide` folder.

## TODOs in `userguide/basics`

## TODOs in `userguide/blocks`

```bash
blocks/curves.md-19-</div>
blocks/curves.md-20-
blocks/curves.md:21:TODO: Add general description
--
blocks/time_loop.md-210-### Time discretizationn
blocks/time_loop.md-211-
blocks/time_loop.md:212:TODO: describe different options for time discretization
--
blocks/time_loop.md-276-
blocks/time_loop.md-277-The list of available variables differs between processes.
blocks/time_loop.md:278:TODO: Create list of variables available in different processes
--
blocks/linear_solvers.md-31-## Eigen
blocks/linear_solvers.md-32-
blocks/linear_solvers.md:33:TODO: Add description of Eigen
--
blocks/media.md-60-The available phases depend on the type of the process.
blocks/media.md-61-
blocks/media.md:62:TODO: Add overview of what mediums are available in different processes
--
blocks/media.md-273-
blocks/media.md-274-There are limitations of what variables can be used inside of the `<expression> </expression>` tags. Only the ones related to the process variables can be called. The values defined for example in parameter block are out of reach.
blocks/media.md:275:TODO: Check if this is correct
--
blocks/media.md-348-### Supported types of properties and special cases
blocks/media.md-349-
blocks/media.md:350:TODO: Add content
--
blocks/misc/constitutive_relations.md-22-
blocks/misc/constitutive_relations.md-23-## CreepBGRa
blocks/misc/constitutive_relations.md:24:TODO: Add content
blocks/misc/constitutive_relations.md-25-
blocks/misc/constitutive_relations.md-26-## Ehlers
blocks/misc/constitutive_relations.md:27:TODO: Add content
blocks/misc/constitutive_relations.md-28-
blocks/misc/constitutive_relations.md-29-## Linear Elastic Isotropic
blocks/misc/constitutive_relations.md:30:TODO: Add content
blocks/misc/constitutive_relations.md-31-
blocks/misc/constitutive_relations.md-32-## Linear Elastic Orthotropic
blocks/misc/constitutive_relations.md:33:TODO: Add content
blocks/misc/constitutive_relations.md-34-
blocks/misc/constitutive_relations.md-35-## Lubby2
blocks/misc/constitutive_relations.md:36:TODO: Add content
--
blocks/boundary_conditions.md-24-## Dirichlet
blocks/boundary_conditions.md-25-
blocks/boundary_conditions.md:26:TODO: add description of Dirichlet boundary condition (including different subvariants) and how it can be used in prj file
--
blocks/boundary_conditions.md-36-## Neumann
blocks/boundary_conditions.md-37-
blocks/boundary_conditions.md:38:TODO: add description of Dirichlet boundary condition (including different subvariants) and how it can be used in prj file
--
blocks/parameters.md-33-## Parameters vs properties
blocks/parameters.md-34-
blocks/parameters.md:35:TODO: describe differences in access to parameters and properties
--
blocks/geometry.md-6-+++
blocks/geometry.md-7-
blocks/geometry.md:8:TODO: describe what geometry block does
blocks/geometry.md-9-
blocks/geometry.md:10:TODO: describe how geometry files are structured
--
blocks/test_definitions.md-6-+++
blocks/test_definitions.md-7-
blocks/test_definitions.md:8:TODO: describe usage of test definitions, reference the `-o` and `-r` command line arguments, `vtkdiff`.
--
blocks/nonlinear_solvers.md-27-## Newton
blocks/nonlinear_solvers.md-28-
blocks/nonlinear_solvers.md:29:TODO: add content
--
blocks/nonlinear_solvers.md-43-## Picard
blocks/nonlinear_solvers.md-44-
blocks/nonlinear_solvers.md:45:TODO: add content
--
blocks/processes.md-155-## Jacobian Assembler
blocks/processes.md-156-
blocks/processes.md:157:TODO: Explanations for each type.
blocks/processes.md-158-
blocks/processes.md:159:The global non-linear equation system can be solved either with Picard fix-point iterations or a Newton scheme. (TODO: Reference NLS scheme)
```

## TODOs in `userguide/features`

```bash
features/python_bc.md-31-## Using python boundary condition in project file
features/python_bc.md-32-
features/python_bc.md:33:TODO: add description of how to call python bc from the boundary condition tag
--
features/mfront.md-29-## Preparing MFront file
features/mfront.md-30-
features/mfront.md:31:TODO: add content
--
features/mfront.md-96-## References
features/mfront.md-97-
features/mfront.md:98:TODO: add content
--
features/mfront.md-100-### Benchmarks using MFront
features/mfront.md-101-
features/mfront.md:102:TODO: add content
--
features/mfront.md-104-### Available MFront models
features/mfront.md-105-
features/mfront.md:106:TODO: add content
--
features/mfront.md-108-## Testing MFront model with MTest
features/mfront.md-109-
features/mfront.md:110:TODO: add content
```

## TODOs in `userguide/troubleshooting`
