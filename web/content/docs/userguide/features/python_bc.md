+++
date = "2018-02-27T11:00:13+01:00"
title = "Defining boundary conditions with Python"
author = "Feliks Kiszkurno"
weight = 3
+++
<div class="note">

### Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

</div>

## Goal

Defining a boundary condition with a Python script allows for a more flexible approach of assigning boundary conditions.
For example it can be used to account for measurements in space and time like accounting for a transient boundary condition via
time series data.

## Prerequisites

This feature requires a basic understanding of classes in Python. Information about them and their syntax can be found in the
[official Python documentation](https://docs.python.org/3/tutorial/classes.html).

## Using a Python boundary condition in a project file

<!-- TODO: add a description of how to call a Python bc from the boundary condition tag -->

In order to define a boundary condition with Python, the path to the `*.py` file, containing the boundary condition script, has
to be provided.
It can be done with the tag `<python_script> </python_script>`.
It is commonly placed in the beginning of the project file between the sections "Meshes" and "Processes".
The path to the file can be defined in relative or absolute terms. For example:

```xml
<python_script>Other_files/Boundary_conditions/BC_North.py</python_script>
```

## Prepare the definition of a boundary condition in a Python script

First, the OpenGeoSys module has to be imported:

```python
try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys
```

This module doesn't need to be installed, if OpenGeoSys is compiled with Python support, then it is part of the Python
environment.
For testing or debugging of the Python script outside of the OpenGeoSys run environment, it is necessary to comment this line
out.

OpenGeoSys will expect a child class of "OpenGeoSys.BoundaryCondition" in the Python script:

```python
class SomeBoundaryCondition(OpenGeoSys.BoundaryCondition)
```

The two types of boundary conditions that can be defined in a Python script are Dirichlet and Neumann.

### Dirichlet boundary condition

OpenGeoSys will obtain values for a Dirichlet boundary condition by calling the method "getDirichletBCValue".

#### Input

The following variables are always passed as an input to the Dirichlet boundary condition method (in brackets the default
variable name used in examples) is given:

- time (`t`),
- spatial coordinates (`coords`),
- node id (`node_id`),
- primary variables (`primary_vars`).

#### Output

The Dirichlet boundary condition method has to return two values:

- Boolean (True)
- Value at the point for which it was called

The order in which those variables are provided is important.

#### Example

### Neumann boundary condition

OpenGeoSys will obtain values for a Neumann boundary condition by calling the method `getFlux`.

#### Input

The following variables are always passed as an input to the Dirichlet boundary condition method (in brackets the default variable name used in examples is given):

- time (`t`),
- spatial coordinates (`coords`),
- primary variables (`primary_vars`).

The order in which those variables are provided is important.

#### Output

The expected output of the Neumann boundary condition method is similar to the one for Dirichlet boundary conditions but it is
extended by values of the Jacobian at a specific point at the boundary:

- Boolean (True)
- Value at point for which it was called
- Jacobian

## Last step: Initialize the boundary condition object

Within the script, an object has to be created by which OpenGeoSys will access the values provided by the Python script.
There can be an arbitrary number of boundary condition objects defined in the Python script.
Each of them can be used for one or multiple instances of boundary conditions defined in the project file.

The name of this object will be included in the project file in the tag `<boundary_condition> <boundary_condition>`

## Benchmarks

Examples of the application of Python boundary conditions can be found in [this group of benchmarks](/docs/benchmarks/python-bc/).

## Full example of a Python boundary condition

```python
try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys

class BoundaryCondition(OpenGeoSys):

    def __init__(self):
        super(OpenGeoSys, BoundaryCondition)
        # Do here all steps, that should be done only once, e.g.: read and preprocess the data from csv file

    def getFlux(self, coords, t, primary_vars):
        # For Neumann BC

        return True, BC_value, BC_jacobian

    def getDirichletBCValue(self, t, coords, node_id, primary_vars):

        return True, BC_value

# Initialize boundary condition object
BC_class = BoundaryCondition()
```
