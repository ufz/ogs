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

Defining boundary condition with Python script allows more flexible approach.
For example it can be used to insert time series data or measurements as boundary conditions.

## Prerequisites

OpenGeoSys has to be compiled with `-DOGS_USE_PYTHON` enabled. For more details about compiling OpenGeoSys, see [developer guide - build configuration](/docs/devguide/getting-started/build-configuration/).
This feature requires basic understanding of classes in Python. Information about them and their syntax can be found in [official Python documentation](https://docs.python.org/3/tutorial/classes.html).

## Using python boundary condition in project file

TODO: add description of how to call python bc from the boundary condition tag

In order to define a boundary condition with python the path `*.py` file, containing the boundary condition script, has to be provided.
It can be done with tag `<python_script> </python_script>`.
It is commonly placed in the beginning of the project file between sections "Meshes" and "Processes".
The path to the file can be defined relatively or absolutely. For example:

```xml
<python_script>Other_files/Boundary_conditions/BC_North.py</python_script>
```

## Prepare definition of boundary condition in python script

First, OpenGeoSys module has to be imported:

```python
import OpenGeoSys
```

This module doesn't have to be installed, it is part of the python environment created when OpenGeoSys is compiled with Python support.
For testing or debugging of the Python script outside of OpenGeoSys run environment, it is necessary to comment this line out.

OpenGeoSys will expect a child class of "OpenGeoSys.BoundaryCondition" in the python script:

```python
class SomeBoundaryCondition(OpenGeoSys.BoundaryCondition)
```

The two types of boundary conditions that can be defined in Python script are Dirichlet and Neumann.

### Dirichlet boundary condition

OpenGeoSys will obtain values for Dirichlet boundary condition by calling method "getDirichletBCValue".

#### Input

Following variables are always passed as an input to Dirichlet boundary condition method (in bracket default variable name used in examples):

- time (t),
- spatial coordinates (coords),
- node id (node_id),
- primary variables (primary_vars).

#### Output

Dirichlet boundary condition method has to return two values:

- Boolean (True)
- Value at point for which it was called

The order in which those variables are provided is important.

#### Example

### Neumann boundary condition

OpenGeoSys will obtain values for Neumann boundary condition by calling method "getFlux".

#### Input

Following variables are always passed as an input to Dirichlet boundary condition method (in bracket default variable name used in examples):

- time (t),
- spatial coordinates (coords),
- primary variables (primary_vars).

The order in which those variables are provided is important.

#### Output

The expected output of Neumann boundary condition method is similar to the one for Dirichlet but is expanded by value of the Jacobian at specific point at the boundary:

- Boolean (True)
- Value at point for which it was called
- Jacobian

## Last step initialize boundary condition object

Within the script, an object has to be created by which OpenGeoSys will access the values provided by the Python script.
There can be an arbitrary number of boundary condition objects defined in the Python script.
Each of them can be used for one or multiple instance of boundary condition defined in the project file.

The name of this object will be included in project file in tag `<boundary_condition> <boundary_condition>`

## Benchmarks

For examples of the application of Python boundary conditions can be found in [this group of benchmarks](/docs/benchmarks/python-bc/).

## Full example of python boundary condition

```python
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
