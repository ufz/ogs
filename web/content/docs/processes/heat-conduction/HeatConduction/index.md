+++
author = "Dmitri Naumov"
date = "2025-11-06"
title = "HeatConduction Process"
weight = 1
+++

## Introduction

This process simulates transient heat conduction to model temperature distribution and heat flow.

It supports various thermal material models and can handle both steady-state and transient heat conduction problems with temperature-dependent properties.

### Key Features

The major features of this process are:

- Temperature-dependent material properties
- Mass lumping option for improved numerical stability
- Heat flux computation and output

The process handles both Dirichlet and Neumann boundary conditions and can account for various heat sources and material interfaces.

## Theoretical background

The heat conduction process solves the transient heat equation in strong form:

$$
\rho c_p \frac{\partial T}{\partial t} = \nabla \cdot (\boldsymbol{\lambda} \nabla T) + Q,
$$

where:

- $\rho$ is the density [M·L⁻³]
- $c_p$ is the specific heat capacity [L²·T⁻²·Θ⁻¹]
- $T$ is the temperature [Θ]
- $\boldsymbol{\lambda}$ is the thermal conductivity tensor [M·L·T⁻³·Θ⁻¹]
- $Q$ is the heat source term [M·L⁻¹·T⁻³]; in SI-units [W·m⁻³].

The heat flux vector is defined by Fourier's law:

$$
\mathbf{q} = -\boldsymbol{\lambda} \nabla T,
$$

where $\mathbf{q}$ is the heat flux vector [M·T⁻³].

### Finite Element Discretization

The heat equation is discretized using the finite element method, resulting in:

$$
\mathbf{M} \frac{d T}{dt} + \mathbf{K} T = \mathbf{f},
$$

where the element matrices are defined as follows.

#### Element Storage Matrix

The element storage matrix $\mathbf{M}_e$ [M·L⁻¹·Θ⁻¹] represents the heat capacity:

$$
\mathbf{M}_e = \int_{\Omega^e} \mathbf{N}^T \rho c_p \mathbf{N}  d\Omega,
$$

where $\mathbf{N}$ is the shape function matrix.

#### Element Conductivity Matrix

The element conductivity matrix $\mathbf{K}_e$ [M·L·T⁻³·Θ⁻¹] represents the thermal diffusion:

$$
\mathbf{K}_e = \int_{\Omega^e} (\nabla \mathbf{N})^T \boldsymbol{\lambda} \nabla \mathbf{N}  d\Omega,
$$

where $\nabla \mathbf{N}$ is the gradient of shape functions.

#### Element Load Vector

The load vector $\mathbf{f}_e$ [M·L·T⁻³] includes heat sources and boundary fluxes:

$$
\mathbf{f}_e = \int_{\Omega^e} \mathbf{N}^T Q  d\Omega + \int_{\Gamma^e} \mathbf{N}^T q_n  d\Gamma,
$$

where:

- $Q$ is the volumetric heat source [M·L⁻¹·T⁻³]
- $q_n = \mathbf{q} \cdot \mathbf{n}$ is the surface heat flux [M·T⁻³] with $\mathbf{n}$ the outward normal to the surface.

## Definition in the project file

The heat conduction process has to be declared in the project file in the processes block.
For example in following way:

```xml
<process>
    <name>HeatConduction</name>
    <type>HEAT_CONDUCTION</type>
    <integration_order>2</integration_order>
    <linear>true</linear>
    <mass_lumping>false</mass_lumping>
    <process_variables>
        <process_variable>temperature</process_variable>
    </process_variables>
    <secondary_variables>
        <secondary_variable name="heat_flux"/>
    </secondary_variables>
</process>
```

For more detailed description of tags used in this snippet, please see [Processes]({{< ref "/docs/userguide/blocks/processes" >}}).

### Process variables

The heat conduction process requires a single scalar temperature process variable.
For more details, see [Process variables]({{< ref "process_variables" >}}).

## Media

The heat conduction process requires properties for the medium only.

### Medium properties

Required medium property

| Property name | Units | SI | Notes |
|---|---|---|---|
| `thermal_conductivity` | M·L·T⁻³·Θ⁻¹ | W·m⁻¹·K⁻¹ | Thermal conductivity of the material |
| `specific_heat_capacity` | L²·T⁻²·Θ⁻¹ | J·kg⁻¹·K⁻¹ | Specific heat capacity at constant pressure |
| `density` | M·L⁻³ | kg·m⁻³ | Mass density of the material |

See [medium properties]({{< ref "media" >}}) for more details on defining them.

## Features

### Mass lumping

Mass lumping can be enabled to improve numerical stability for transient problems:

```xml
<mass_lumping>true</mass_lumping>
```

### Linear solver optimization

The process supports linear solver optimizations for improved performance when solving linear problems:

```xml
<linear>true</linear>
<linear_solver_compute_only_upon_timestep_change>false</linear_solver_compute_only_upon_timestep_change>
```

### Source terms

The heat conduction process supports various source term types including:

- **Volumetric heat sources**: Applied to entire domain regions
- **Line source terms**: For 2D and 3D problems.

## Available benchmarks

To gain more insight into this process, you can investigate [heat conduction benchmarks]({{< ref "/docs/benchmarks/heatconduction" >}}).
