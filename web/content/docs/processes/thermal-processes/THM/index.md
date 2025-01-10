+++
author = "Feliks Kiszkurno"
date = "2023-01-10"
title = "Thermo Hydro Mechanics Process"
weight = 1
+++

This page describes Thermo-Hydro-Mechanics Process (THM).

<div class="note">

### Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to the contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

</div>

## Introduction

## Theoretical background

## Implementation

### Supported phases

- Aqueous liquid
- Frozen liquid
- Solid

## Input variables and parameters

List of medium properties required by THM process.

### Medium properties

- bulk modulus
- density
- specific heat capacity

## Input parameters in the project file

THM process has to be declared in project file in the processes block. For example in following way:

```xml
<processes>
    <process>
        <type>THERMO_HYDRO_MECHANICS</type>
    </process>
</processes>
```

### Process variables

Following process variables are available in THM process:

- `temperature`
- `pressure`
- `displacement`

For more details, see [Process variables]({{< ref "process_variables" >}}).

### Example of full section defining THM process

## Features

### Specific body force

### Thermal porosity mixing

THM can automatically obtain thermal conductivity for the medium based on thermal conductivities of phases and porosity.

See [Thermal conductivity: effective porosity mixing]({{% ref "effective-porosity-mixing" %}}) for more information.

#### Examples

## Available benchmarks

To gain more insight into THM process, you can investigate [THM benchmarks]({{< ref "thermo-hydro-mechanics" >}}).

## References
