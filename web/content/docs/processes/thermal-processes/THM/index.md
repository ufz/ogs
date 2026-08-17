+++
author = "Feliks Kiszkurno"
date = "2023-01-10"
title = "Thermo Hydro Mechanics Process"
weight = 1
+++

This page describes Thermo-Hydro-Mechanics Process (THM).

<div class="note">

## Work in progress

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

Those properties are defined on medium level. See [medium properties]({{< ref "media#properties" >}}) for more details on defining them.

| Property name | Mandatory | Constant | Function | Curve | Parameter | Other |
| --- | --- | --- | --- | --- | --- | --- |
| Thermo-osmosis coefficient | No | Yes | No | No | No | - |
| Thermo-osmosis permeability | No | Yes | No | No | No | - |

Thermo-osmosis can be parametrised in two, mutually exclusive, ways:

- `thermal_osmosis_coefficient` sets the thermo-osmotic coefficient $k_T$
  ($[k_T] = m^2/(K \cdot s)$) directly. Use this when a tabulated $k_T$ is
  available.
- `thermal_osmosis_permeability` sets the scalar thermo-osmotic permeability
  $\epsilon_T$ ($[\epsilon_T] = Pa/K$), from which $k_T = \epsilon_T k / \mu$
  is formed using the medium's intrinsic permeability $k$ and the
  AqueousLiquid phase's dynamic viscosity $\mu$. Use this when
  $\epsilon_T$ is the measured quantity and $k_T$ should track a spatially
  varying $k$ and a temperature-dependent $\mu$.

Note that `thermal_osmosis_permeability` is not a permeability and does not
have a permeability's dimension; the name was kept for continuity with the
existing property naming.

Older project files defining `thermal_osmosis_coefficient` on the solid phase
are moved to the medium level by the
[`scripts/dev/move_thermal_osmosis.py`](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/scripts/dev/move_thermal_osmosis.py)
helper script, which with `--convert` also rewrites the property as
`thermal_osmosis_permeability`. A `thermal_osmosis_coefficient` already given
on the medium level is a valid parametrisation and is left unchanged by the
script.

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
