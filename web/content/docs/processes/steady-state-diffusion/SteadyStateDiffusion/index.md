+++
author = "Dmitri Naumov"
date = "2025-11-01"
title = "Steady-State Diffusion Process"
weight = 1
+++

## Introduction

The Steady-State Diffusion Process solves the steady-state diffusion equation for pressure (or concentration) transport in porous media.
This is a linear elliptic partial differential equation that models phenomena such as:

- Steady-state groundwater flow in porous media
- Heat conduction in solids
- Mass diffusion in continuous media
- Electrical potential distribution

### Key Features

The Steady-State Diffusion process is linear, making it suitable for efficient solution using direct or iterative linear solvers.
The diffusion coefficient can be specified as a tensor to model anisotropic diffusion behavior.
The process supports calculation of surface fluxes through the `calculatesurfaceflux` feature, which can be used for mass balance verification.
The surface flux is computed as `specific_flux` on the specified mesh.
The process can output secondary variables such as `darcy_velocity`, the velocity field computed as $-\mathbf{k} \nabla p$.

## Theoretical background

The Steady-State Diffusion Process solves the steady-state diffusion equation in strong form:

$$
\nabla \cdot (-\mathbf{k} \nabla p) = Q
$$

where:

- $p$ is the liquid phase pressure (primary variable)
- $\mathbf{k}$ is the diffusion coefficient tensor (material property)
- $Q$ is the source term (can be zero for no source)
- $\nabla$ is the gradient operator

### Finite Element Discretization

The element stiffness matrix and load vector are computed as:

$$
\mathbf{K}_e = \int_{\Omega^e} (\nabla \mathbf{N})^T \mathbf{k} (\nabla \mathbf{N})  d\Omega
$$

$$
\mathbf{b}_e = \int_{\Omega^e} \mathbf{N}^T Q  d\Omega
$$

where $\mathbf{N}$ contains the shape functions and their gradients.

## Definition in the project file

Steady-State Diffusion process has to be declared in the project file in the processes block.
For example in following way:

```xml
<process>
    <name>SteadyStateDiffusion</name>
    <type>STEADY_STATE_DIFFUSION</type>
    <integration_order>2</integration_order>
    <process_variables>
        <process_variable>pressure</process_variable>
    </process_variables>
    <secondary_variables>
        <secondary_variable name="darcy_velocity"/>
    </secondary_variables>
    <calculatesurfaceflux> <!-- optional -->
        <mesh>boundary_mesh</mesh>
        <property_name>specific_flux</property_name>
    </calculatesurfaceflux>
</process>
```

For more detailed description of tags used in this snippet, please see [Processes]({{< ref "/docs/userguide/blocks/processes" >}}).

### Process variables

The Steady-State Diffusion process requires single scalar process variable, *e.g.* `pressure` in the configuration above.
For more details, see [Process variables]({{< ref "process_variables" >}}).

### Medium properties

Required medium properties for each medium. See [medium properties]({{< ref "media#properties" >}}) for more details on defining them.

| Property name | Mandatory | Constant | Function | Linear | Parameter | Other |
|---|---|---|---|---|---|---|
| `reference_temperature` | Yes | Yes | No | No | No | - |
| `diffusion` | Yes | Yes | Yes | Yes | Yes | - |

## Available benchmarks

To gain more insight into Steady-State Diffusion process, you can investigate [Steady-State Diffusion benchmarks]({{< ref "/docs/benchmarks/elliptic" >}}).
