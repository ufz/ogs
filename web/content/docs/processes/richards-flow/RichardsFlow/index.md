+++
author = "Dmitri Naumov"
date = "2025-11-01"
title = "Richards Flow Process"
weight = 1
+++

## Introduction

The Richards Flow Process solves the Richards equation for variably saturated flow in porous media.
This is a nonlinear partial differential equation that models water flow in unsaturated soils, combining Darcy's law with mass conservation.
The process is applicable to phenomena such as:

- Infiltration and drainage in unsaturated soils
- Groundwater recharge and capillary rise
- Irrigation and drainage systems
- Vadose zone hydrology
- Contaminant transport in unsaturated media

### Key Features

The Richards Flow process is nonlinear due to the dependence of hydraulic conductivity and water content on pressure (or suction).
It accounts for gravity effects and can handle both saturated and unsaturated flow conditions.
Secondary variables include saturation and Darcy velocity, which can be output for visualization and analysis.
The process can include mass lumping for improved numerical stability in certain scenarios.

## Theoretical background

The Richards Flow Process solves the Richards equation in strong form:

$$
\left( s_p S^l + \phi S^l \frac{1}{\rho^l} \frac{\partial \rho^l}{\partial p} + \phi \frac{\partial S^l}{\partial p} \right) \frac{\partial p}{\partial t} = \nabla \cdot \left[ \frac{\mathbf{k} k_{rel}}{\mu} (\nabla p - \rho^l \mathbf{g}) \right] + Q,
$$

where:

- $p$ is the pore pressure (primary variable, negative in unsaturated zone) [M·L⁻¹·T⁻²]
- $\phi$ is the porosity [-]
- $S^l$ is the liquid saturation [-]
- $s_p$ is the specific storage coefficient [L·T²·M⁻¹]
- $\mathbf{k}$ is the intrinsic permeability tensor [L²]
- $k_{rel}$ is the relative permeability [-]
- $\mu$ is the liquid viscosity [M·L⁻¹·T⁻¹]
- $\rho^l$ is the density of liquid [M·L⁻³]
- $\mathbf{g}$ is the gravitational acceleration vector [L·T⁻²]
- $Q$ is the source/sink term [T⁻¹]

The liquid velocity $\mathbf{v}^l$ [L·T⁻¹] is described by Darcy's law as:

$$
\mathbf{v}^l = -\frac{\mathbf{k} k_{rel}}{\mu} (\nabla p - \rho^l \mathbf{g}).
$$

Note, that a part of the Laplace term is neglected, see [Note on the Laplace term with non-constant density](#note-on-the-laplace-term-with-non-constant-density).

### Finite Element Discretization

The Richards equation is discretized using the finite element method, resulting in a system of nonlinear equations.
The discrete system takes the form:

$$
\mathbf{M}_e \frac{d\mathbf{p}_e}{dt} + \mathbf{K}_e \mathbf{p}_e = \mathbf{b}_e
$$

where the element matrices are defined as follows.

#### Element Mass Matrix

The element mass matrix $\mathbf{M}_e$  [L⁴·T²·M⁻¹] accounts for storage effects:

$$
\mathbf{M}_e = \int_{\Omega^e} \mathbf{N}^T \left( s_p \cdot S^l + \phi S^l \frac{1}{\rho^l} \frac{\partial \rho^l}{\partial p} - \phi \frac{\partial S^l}{\partial p_c} \right) \mathbf{N} d\Omega,
$$

where $p_c = -p$ is the capillary pressure.

#### Element Stiffness Matrix

The element stiffness matrix $\mathbf{K}_e$ [L⁴·M⁻¹] represents the flow terms:

$$
\mathbf{K}_e = \int_{\Omega^e} (\nabla \mathbf{N})^T \frac{\mathbf{k} k_{rel}}{\mu} (\nabla \mathbf{N}) d\Omega.
$$

#### Element Load Vector

The load vector $\mathbf{b}_e$ [L³·T⁻¹] includes gravity terms:

$$
\mathbf{b}_e = \int_{\Omega^e} (\nabla \mathbf{N})^T \frac{\mathbf{k} k_{rel}}{\mu} \rho^l \mathbf{g} d\Omega
$$

## Definition in the project file

Richards Flow process has to be declared in the project file in the processes block.
For example in following way:

```xml
<process>
    <name>RichardsFlow</name>
    <type>RICHARDS_FLOW</type>
    <integration_order>2</integration_order>
    <specific_body_force>0 -9.81</specific_body_force>
    <mass_lumping>true</mass_lumping>
    <process_variables>
        <process_variable>pressure</process_variable>
    </process_variables>
    <secondary_variables>
        <secondary_variable name="saturation"/>
        <secondary_variable name="darcy_velocity"/>
    </secondary_variables>
</process>
```

For more detailed description of tags used in this snippet, please see [Processes]({{< ref "/docs/userguide/blocks/processes" >}}).

### Process variables

The Richards Flow process requires single scalar process variable `pressure` in the configuration above.
For more details, see [Process variables]({{< ref "process_variables" >}}).

## Media

Richards Flow requires properties for the aqueous liquid phase and for the medium.

### Medium properties

Required medium properties for each medium. See [medium properties]({{< ref "media#properties" >}}) for more details on defining them.

| Property name | Units | Notes |
|---|---|---|
| `permeability` | L² | - |
| `porosity` | - | - |
| `reference_temperature` | Θ | - |
| `relative_permeability` | - | RelPermBrooksCorey, RelativePermeabilityVanGenuchten *etc.*|
| `saturation` | - | SaturationBrooksCorey, SaturationVanGenuchten, *etc.* |
| `storage` | L T²/M | - |

### Liquid properties

Required liquid properties for each medium. See [phase properties]({{< ref "media#phases" >}}) for more details on defining them.

| Property name | Units | Notes |
|---|---|---|
| `density` | M·L⁻³ | $\frac{\partial \rho^l}{\partial p}$ is part of storage term. |
| `viscosity` | M·L⁻¹·T⁻¹ | - |

## Features

### Specific body force

The specific body force vector can be specified to account for gravity effects.

### Mass lumping

The diagonal lumping of the mass matrix can be enabled by adding the following tag to the `<process> </process>` block:

```xml
<mass_lumping>true</mass_lumping>
```

## Note on the Laplace term with non-constant density

If $\rho^l$ is not constant, the Laplace term for the density scaled volume balance equation is

$$
\frac{1}{\rho^l}\nabla\left(\rho^l\frac{{\mathbf k}k_{rel}}{\mu}(\nabla p-\rho^l\mathbf g)\right)
$$

and its corresponding weak form is

$$
\int\frac{1}{\rho^l}\nabla\left(\rho^l\frac{{\mathbf k}k_{rel}}{\mu}(\nabla p-\rho^l\mathbf g)\right)\psi\mathrm{d}\Omega,
$$

where $\psi$ the test function.
Denoting $-\left(\frac{{\mathbf k}k_{rel}}{\mu}(\nabla p-\rho^l\mathbf g)\right)$ as $\mathbf v$, that weak term can be expanded as

$$
\int\nabla(\mathbf{v}\psi)-\nabla(\frac{\psi}{\rho^l})\mathbf{v}\mathrm{d}\Omega =
\int(\mathbf{v}\cdot\mathbf{n}\psi)\mathrm{d}\Gamma-\int\nabla\psi\cdot\mathbf{v}\mathrm{d}\Omega
-\int\psi\nabla(\frac{1}{\rho^l})\mathbf{v}\mathrm{d}\Omega.
$$

We see that the third term above is an extra one to be computed if the volume balance equation and non constant density are considered.
The quantity of the extra term is usually tiny $\nabla(\frac{1}{\rho^l})=-\frac{1}{(\rho^l)^2}\frac{\partial \rho^l}{\partial p}\nabla p$ and it is neglected in the implementation.

## Available benchmarks

To gain more insight into Richards Flow process, you can investigate [Richards Flow benchmarks]({{< ref "/docs/benchmarks/richards-flow" >}}).
