+++
author = "Feliks Kiszkurno, Wenqing Wang"
date = "2023-01-10"
title = "Thermo-Richards-Mechanics Process"
weight = 2
+++

This page describes Thermo-Richards-Mechanics Process (TRM)

<div class="note">

### Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

</div>

## Introduction

## Theoretical background

Global assembler for the monolithic scheme of the non-isothermal Richards flow coupled with mechanics.

### Governing equations without vapor diffusion

The energy balance equation is given by:

$$
(\rho c_p)^{eff}\dot T -
\nabla (\mathbf{k}_T^{eff} \nabla T)+\rho^l c_p^l \nabla T \cdot
\mathbf{v}^l
= Q_T
$$

with $T$ the temperature, $(\rho c_p)^{eff}$ the effective volumetric heat capacity, $\mathbf{k}_T^{eff}$ the effective thermal conductivity, $\rho^l$ the density of liquid, $c_p^l$ the specific heat capacity of liquid, $\mathbf{v}^l$ the liquid velocity, and $Q_T$ the point heat source.

The effective volumetric heat can be considered as a composite of the contributions of solid phase and the liquid phase as
$$
(\rho c_p)^{eff} = (1-\phi) \rho^s c_p^s + S^l \phi \rho^l c_p^l
$$
with $\phi$ the porosity, $S^l$  the liquid saturation, $\rho^s$ the solid density, and $c_p^s$ the specific heat capacity of solid. Similarly, the effective thermal conductivity is given by:
$$
\mathbf{k}_T^{eff} = (1-\phi) \mathbf{k}_T^s + S^l \phi k_T^l \mathbf I
$$

where $\mathbf{k}_T^s$ is the thermal conductivity tensor of solid, $k_T^l$ is the thermal conductivity of liquid, and $\mathbf I$ is the identity tensor.

The mass balance equation is given by:
$$
\begin{eqnarray*}{
 \left(S^l\beta - \phi\frac{\partial S}{\partial p_c}\right) \rho^l\dot p - S \left( \frac{\partial \rho^l}{\partial T} + \rho^l(\alpha_B -S) \alpha_T^s \right)\dot T\\
+\nabla (\rho^l \mathbf{v}^l) + S \alpha_B \rho^l \nabla \cdot \dot{\mathbf u}= Q_H
 }
 \end{eqnarray*}
 $$

 where $p$ is the pore pressure, $p_c$ is the capillary pressure, which is $-p$ under the single phase assumption, $\beta$ is a composite coefficient by the liquid compressibility and solid compressibility, $\alpha_B$ is the Biot's constant, $\alpha_T^s$ is the linear thermal expansivity of solid, $Q_H$ is the point source or sink term,  $\mathbf u$ is the displacement, and $H(S-1)$ is the Heaviside function.
 The liquid velocity $\mathbf{v}^l$ is described by the Darcy's law as
 $$
 \mathbf{v}^l=-\frac{{\mathbf k} k_{ref}}{\mu} (\nabla p - \rho^l \mathbf g)
 $$

 with ${\mathbf k}$ the intrinsic permeability, $k_{ref}$ the relative permeability, $\mathbf g$ the gravitational force.

 The momentum balance equation takes the form of
 $$
 \nabla (\mathbf{\sigma}-b(S)\alpha_B p^l \mathbf I) +\mathbf f=0
 $$

 with $\mathbf{\sigma}$ the effective stress tensor, $b(S)$ the Bishop model, $\mathbf f$ the body force, and $\mathbf I$ the identity.
 The primary unknowns of the momentum balance equation are the displacement $\mathbf u$, which is associated with the stress by the generalized Hook's law as
 $$
 {\dot {\mathbf {\sigma}}} = C {\dot {\mathbf \epsilon}}^e
 = C ( {\dot {\mathbf \epsilon}} - {\dot {\mathbf \epsilon}}^T
 -{\dot {\mathbf \epsilon}}^p - {\dot {\mathbf \epsilon}}^{sw}-\cdots )
 $$

 with $C$ the forth order elastic tensor,
 ${\dot {\mathbf \epsilon}}$ the total strain rate,
 ${\dot {\mathbf \epsilon}}^e$ the elastic strain rate,
 ${\dot {\mathbf \epsilon}}^T$ the thermal strain rate,
 ${\dot {\mathbf \epsilon}}^p$ the plastic strain rate,
 ${\dot {\mathbf \epsilon}}^{sw}$ the swelling strain rate.

 The strain tensor is given by displacement vector as
 $$
 \mathbf \epsilon =
 \frac{1}{2} \left((\nabla \mathbf u)^{\text T}+\nabla \mathbf u\right)
 $$

 where the superscript ${\text T}$ means transpose.

## Implementation

### Supported phases

- Aqueous liquid
- Solid

## Input variables and parameters

List of medium properties required by TRM process.

### Medium phase properties

Those properties are defined on the phase level for each medium. See [phase properties]({{< ref "media#phases" >}}) for more details on defining them.

|Property name | Mandatory | Constant | Function | Linear | Parameter | Other |
|---|---|---|---|---|---|---|
| Bulk modulus | Yes | Yes | No | No | No | - |
| Density | Yes | Yes | Yes | No | No | - |
| Latent heat | Yes | No | No | No | No | LatentWaterVapourLatentHeat |
| Specific heat capacity | Yes | Yes | No | Yes | No | - |
| Storage | Yes | Yes | No | No | No | - |
| Swelling stress rate | Yes | No | No | No | No | SaturationDependentSwelling |
| Thermal conductivity | Yes | Yes | No | No | No | - |
| Thermal diffusion enhancement factor | Yes | Yes | No | No | No | - |
| Vapour density | Yes | No | No | No | No | WaterVapourDensity |
| Vapour diffusion | Yes | No | No | No | No | VapourDiffusionFEBEX |
| Thermal expansivity | No | Yes | No | No | Yes | - |
| Thermo-osmosis coefficient | No | Yes | No | No | No | - |

### Medium properties

Those properties are defined on medium level. See [medium properties]({{< ref "media#properties" >}}) for more details on defining them.

|Property name | Mandatory | Constant | Function | Curve | Parameter | Other |
|---|---|---|---|---|---|---|
| Biot coefficient | Yes | Yes | No | No | No | - |
| Bishop effective stress | Yes | No | No | No | No | BishopPowerLaw, BishopSturationCutoff |
| Permeability | Yes | Yes | Yes | No | Yes | PermeabilityOrthotropicPowerLaw |
| Porosity | Yes | Yes | No | No | No | PorosityFromMassBalance |
| Relative permeability | Yes | Yes | No | Yes | No | RelPermBrooksCorey, RelativePermeabilityVanGenuchten |
| Saturation | Yes | Yes | No | Yes | No | SaturationLiakopoulos, SaturationVanGenuchten |
| Storage | Yes | Yes | No | No | No | - |
| Thermal conductivity | Yes | Yes | Yes | Yes | Yes | [EffectiveThermalConductivityPorosityMixing](#thermal-porosity-mixing) |
| Transport porosity | Yes | No | No | No | No | TransportPorosityFromMassBalance |

## Input parameters in the project file

TRM process has to be declared in the project file in the processes block. For example in following way:

```xml
<processes>
    <process>
        <type>THERMO_RICHARDS_MECHANICS</type>
    </process>
</processes>
```

### Process variables

Following process variables are available in TRM process:

- `temperature`
- `pressure`
- `displacement`

For more details, see [Process variables]({{< ref "process_variables" >}}).

### Example of full section defining TRM process

For more detailed description of tags used in this snippet, please see [Processes]({{< ref "/docs/userguide/blocks/processes#process-variables" >}}).

```xml
  <processes>
    <process>
      <name>BodyForceTest</name>
      <type>THERMO_RICHARDS_MECHANICS</type>
      <mass_lumping>false</mass_lumping>
      <integration_order>3</integration_order>
      <constitutive_relation id="0">
            ...
      </constitutive_relation>
      <process_variables>
        <temperature>temperature</temperature>
        <pressure>pressure</pressure>
        <displacement>displacement</displacement>
      </process_variables>
      <secondary_variables>
            ...
      </secondary_variables>
      <specific_body_force>0 -9.81</specific_body_force>
      <apply_body_force_for_deformation>false</apply_body_force_for_deformation>
      <initial_stress>Initial_stress</initial_stress>
    </process>
  </processes>
```

For more information on tags `<apply_body_force_to_deformation>` and `<mass_lumping>` see section [Features](#features) at this page.

## Features

### Specific body force

### Mass lumping

The diagonal lumping of the mass matrix of the Richards equation can be using by adding following tag to the `<process> </process>` block:

```xml
<mass_lumping>true</mass_lumping>
```

### Applying body force on deformation

### Thermal porosity mixing

TRM can automatically obtain thermal conductivity for the medium based on thermal conductivities of phases and porosity.

See [Thermal conductivity: effective porosity mixing]({{% ref "effective-porosity-mixing" %}}) for more information.

#### Examples

## Available benchmarks

To gain more insight into TRM process, you can investigate [TRM benchmarks]({{< ref "thermo-richards-mechanics" >}}).

## References
