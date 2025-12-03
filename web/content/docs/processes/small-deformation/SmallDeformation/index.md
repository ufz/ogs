+++
author = "Dmitri Naumov"
date = "2025-11-06"
title = "Small Deformation Process"
weight = 1
+++

## Introduction

This process simulates the mechanical behavior of solids under small deformations and is suitable for stress and deformation analyses.

It supports a variety of constitutive models to represent different deformation mechanisms, including elastic, elasto-plastic, viscoelastic, and creep behavior.

### Key Features

The major features of this process are:

- Integration with MFront for defining custom constitutive models.
- B-bar method for volumetric locking mitigation.
- Principal stress computation and output.
- Material force calculations and output.

The process handles both 2D and 3D problems with various element types and can account for initial stress conditions, body forces, and complex boundary conditions.

## Theoretical background

The small deformation process solves the equilibrium equations in strong form:

$$
\nabla \cdot \boldsymbol{\sigma} + \mathbf{b} = \mathbf{0}
$$

where:

- $\boldsymbol{\sigma}$ is the Cauchy stress tensor [M·L⁻¹·T⁻²]
- $\mathbf{b}$ is the body force vector [M·L⁻²·T⁻²].

The constitutive relationship is given by:

$$
\boldsymbol{\sigma} = \mathbf{C} : \boldsymbol{\varepsilon}
$$

where:

- $\mathbf{C}$ is the fourth-order stiffness tensor [M·L⁻¹·T⁻²].
- $\boldsymbol{\varepsilon}$ is the small strain tensor [-]

The strain-displacement relationship for small deformations is:

$$
\boldsymbol{\varepsilon} = \frac{1}{2} \left( \nabla \mathbf{u} + (\nabla \mathbf{u})^T \right)
$$

where:

- $\mathbf{u}$ is the displacement vector [L].

### Finite Element Discretization

The equilibrium equation is discretized using the finite element method, resulting in:

$$
\mathbf{K} \mathbf{u} = \mathbf{f}
$$

where the element matrices are defined as follows.

#### Element Stiffness Matrix

The element stiffness matrix $\mathbf{K}_e$ [M·L·T⁻²] represents the elastic energy:

$$
\mathbf{K}_e = \int_{\Omega^e} \mathbf{B}^T \mathbf{C} \mathbf{B}  d\Omega,
$$

where $\mathbf{B}$ is the strain-displacement matrix.

#### Element Load Vector

The load vector $\mathbf{f}_e$ [M·L·T⁻²] includes body forces and surface tractions:

$$
\mathbf{f}_e = \int_{\Omega^e} \mathbf{N}^T \mathbf{b}  d\Omega + \int_{\Gamma^e} \mathbf{N}^T \mathbf{\tau}  d\Gamma,
$$

where:

- $\mathbf{N}$ is the shape function matrix
- $\mathbf{\tau}=\boldsymbol{\sigma}\cdot\mathbf{n}$ is the surface traction vector [M·L⁻¹·T⁻²] with $\mathbf{n}$ the outer normal to the surface.

## Definition in the project file

The small deformation process has to be declared in the project file in the processes block.
For example in following way:

```xml
<process>
    <name>SmallDeformation</name>
    <type>SMALL_DEFORMATION</type>
    <integration_order>2</integration_order>
    <specific_body_force>0 -9.81</specific_body_force>
    <use_b_bar>true</use_b_bar>
    <constitutive_relation>
        <type>LinearElasticIsotropic</type>
        <youngs_modulus>E</youngs_modulus>
        <poissons_ratio>nu</poissons_ratio>
    </constitutive_relation>
    <process_variables>
        <process_variable>displacement</process_variable>
    </process_variables>
    <secondary_variables>
        <secondary_variable name="sigma"/>
        <secondary_variable name="epsilon"/>
    </secondary_variables>
</process>
```

For more detailed description of tags used in this snippet, please see [Processes]({{< ref "/docs/userguide/blocks/processes" >}}).

### Process variables

The small deformation process requires displacement process variable.
For 2D problems, the displacement variable should have 2 components ($u_x$, $u_y$).
For 3D problems, the displacement variable should have 3 components ($u_x$, $u_y$, $u_z$).
For more details, see [Process variables]({{< ref "process_variables" >}}).

## Media

The small deformation process requires properties for the solid phase for each medium.

### Solid properties

Required solid property

| Property name | Units | SI | Notes |
|---|---|---|---|
| `density` | M·L⁻³ | kg·m⁻³ | Mass density of the solid material |

See [solid properties]({{< ref "media#solids" >}}) for more details on defining them.

### Material Models

Material constitutive relations and their parameters are specified in the process section, not in the media section.
See the example in [Definition in the project file](#definition-in-the-project-file) above.

## Features

### Specific body force

The specific body force vector can be specified to account for gravity effects:

```xml
<specific_body_force>0 -9.81</specific_body_force>
```

### B-bar method

The B-bar method can be enabled to mitigate volumetric locking in nearly incompressible materials:

```xml
<use_b_bar>true</use_b_bar>
```

### Initial stress conditions

Initial stress conditions can be specified:

```xml
<initial_stress>parameter_name</initial_stress>
```

### Reference temperature

Reference temperature can be specified for temperature dependent material models to simulate thermo-mechanical coupling:

```xml
<reference_temperature>parameter_name</reference_temperature>
```

## Material Models

The small deformation process supports various constitutive models:

- **LinearElastic**: Linear elastic material behavior
- **MohrCoulomb**: Mohr-Coulomb failure criterion
- **DruckerPrager**: Drucker-Prager plasticity model
- **CamClay**: Modified Cam Clay model
- **Burgers**: Viscoelastic Burgers model
- **Ehlers**: Elastoplastic model with damage
- **CreepBGR**
- **Lubby2**
- **MFront**: Generic material model interface

Each model has specific parameters and capabilities.
Refer to the material model documentation for detailed information.

## Available benchmarks

To gain more insight into this process, you can investigate [small deformation benchmarks]({{< ref "/docs/benchmarks/small-deformations" >}}).
