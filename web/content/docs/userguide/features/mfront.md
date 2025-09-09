+++
title = "Using MFront with OpenGeoSys"
author = "Feliks Kiszkurno"
weight = 2
models = [ "material" ]
+++
<div class="note">

### Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

</div>

## Goal

[MFront](https://tfel.sourceforge.net/) can be used to run simulation with constitutive model not available in OpenGeoSys or
with a modification of one of available ones.

## Prerequisites

OpenGeoSys has to be compiled with the flag `-DOGS_USE_MFRONT` in that case. It will enable MFront support and download
necessary libraries.
For more details about compiling OpenGeoSys, see [developer guide - build configuration](/docs/devguide/getting-started/build-configuration/) and [developer guide - MFront installation](/docs/devguide/packages/mfront/).

## Preparing MFront file

<!-- TODO: add content -->

The details of preparing MFront file will not be discussed here. Only the necessary details will be stated.

### Name of the model

At the beginning of the file, name of the model has to be provided as follows:

```c++
@Behaviour ModelName
```

The MFront file has to be the same as name provided in the line above:

```c++
ModelName.mfront
```

### Material properties

The material properties, that should be read from project file and be used in the MFront model have to be defined using
following syntax:

```c++
@MaterialProperty type MatPropName
```

The project file has to contain references pointing from one of the parameters to material properties used by the MFront model.
For details, see next section.

## Preparing project file

Details on how to use a constitutive relation defined with MFront is described [here](/docs/userguide/blocks/misc/constitutive_relations/#mfront).

## Preparing OpenGeoSys to use a new MFront model

After the MFront file has been prepared, it has to be placed in the folder containing OpenGeoSys source code at the following
path:

```c++
ogs-source-code/MaterialLib/SolidModels/MFront
```

in the same folder there is the `CMakeLists.txt`.
In that file, there is a list containing names of all MFront models stored in the folder.
Name of the newly added model has to be added to that list:

```c++
mfront_behaviours_check_library(
    ...
    ModelName
    ...
)
```

To make the new "ModelName" available rerun the "configure" and "generate" CMake-steps, and recompile OpenGeoSys.
(This process should take less time than the first time, as only new code will be compiled.)

<div class='note'>

### Warning

If a compilation error is encountered, it can be caused by inconsistent use of name of the model.
It is case sensitive and it has to be written in exactly the same way everywhere, where it is used.

</div>

## References

<!-- TODO: add content -->

### Benchmarks using MFront

<!-- TODO: add content -->

### Available MFront models

<!-- TODO: add content -->

## Testing MFront model with MTest

<!-- TODO: add content -->
