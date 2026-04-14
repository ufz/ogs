+++
author = "Thomas Fisher, Dmitri Naumov, Fabien Magri, Marc Walther, Tianyuan Zheng, Olaf Kolditz"
date = "2026-04-04"
title = "Hydro-Thermal Process (HT)"
weight = 3
+++

## Introduction

The Hydro-Thermal (HT) process models coupled groundwater flow and heat transport in porous media.
Both sub-processes are described by parabolic partial differential equations that are coupled through temperature-dependent fluid properties (density, viscosity) and the advective heat transport term.

### Key features

- Monolithic and staggered coupling schemes.
- Fluid compressibility and optional solid thermal expansion coupling.
- Hydrodynamic thermal dispersion (longitudinal and transversal dispersivities).
- Numerical stabilisation for advection-dominated problems.
- Fracture flow support via aperture size parameter.
- Optional solid thermal expansion with Biot constant.
- Surface flux calculation.

### Physical variables

- **Primary variables**: pressure $p$ and temperature $T$.
- **Secondary variable**: Darcy velocity $\mathbf{q}$.

## Theoretical background

### Balance equation framework

Both the flow and heat transport processes are derived from integral conservation laws. Starting with the conservation of mass (or energy), we write:

$$
\frac{\partial}{\partial t} \int_{\Omega} S(u) \, d\Omega = -\int_{\Gamma} \mathbf{J} \cdot \mathbf{n} \, d\sigma + \int_{\Omega} Q \, d\Omega
$$

where $S(u)$ is the volume density of the conserved quantity, $\mathbf{J}$ is the flux, $\mathbf{n}$ is the outward normal, and $Q$ represents sources/sinks.

Applying the divergence theorem (Gauss) to convert the boundary integral to a volume integral:

$$
\int_{\Omega} \left[ \frac{\partial S(u)}{\partial t} + \nabla \cdot \mathbf{J} - Q \right] d\Omega = 0
$$

Since this holds for any domain $\Omega$, the strong form follows:

$$
\frac{\partial S(u)}{\partial t} + \nabla \cdot \mathbf{J} - Q = 0
$$

Substituting constitutive laws for flux $\mathbf{J}$ (such as Darcy's law for groundwater and Fourier's law for heat) yields the specific governing equations.

### Groundwater flow equation

The Darcy velocity is given by:

$$
\mathbf{q} = -\frac{\boldsymbol{\kappa}}{\mu(T)} \left( \nabla p + \varrho_f(T) g \mathbf{e}_z \right)
$$

The pressure equation including fluid compressibility reads:

$$
\left( \frac{\phi}{\varrho_f} \frac{\partial \varrho_f}{\partial p} + S \right) \frac{\partial p}{\partial t} - \nabla \cdot \left[ \frac{\boldsymbol{\kappa}}{\mu(T)} \nabla \Psi \right] - Q = 0
$$

where $\Psi = p + \varrho_f(T) g z$ is the piezometric head, $S$ is the storage coefficient, $\phi$ is the porosity, and $Q$ is a source/sink term.

**Derivation of the storage term:** Applying the balance equation framework with $S(p) = \phi \varrho_f(p,T)$ (fluid mass per unit volume) yields:

$$
\frac{\partial \phi \varrho_f}{\partial t} - \nabla \cdot \left[ \frac{\boldsymbol{\kappa}}{\mu} \nabla \Psi \right] - Q = 0
$$

Assuming the solid matrix is incompressible ($\partial \phi/\partial t = 0$), the storage term expands as:

$$
\frac{\partial \phi \varrho_f}{\partial t} = \phi \left( \frac{\partial \varrho_f}{\partial p} \frac{\partial p}{\partial t} + \frac{\partial \varrho_f}{\partial T} \frac{\partial T}{\partial t} \right)
$$

Under the Boussinesq approximation, the temperature dependence of density in storage is neglected, and fluid compressibility $\partial \varrho_f/\partial p$ is assumed constant, allowing the definition of an effective storage coefficient $S = \phi \frac{\partial \varrho_f}{\partial p}$. Temperature dependence is retained only in the buoyancy term through $\varrho_f(T)$ in the piezometric head.

### Heat transport equation

$$
c_p \frac{\partial T}{\partial t} - \nabla \cdot (\boldsymbol{\lambda} \nabla T) + \varrho_f c_f \langle \mathbf{q}, \nabla T \rangle = 0
$$

where:

| Symbol | Definition |
|---|---|
| $c_p = \varrho_f \phi c_f + \varrho_s (1 - \phi) c_s$ | Volumetric heat capacity of the mixture |
| $\boldsymbol{\lambda} = \boldsymbol{\lambda}^{\mathrm{cond}} + \boldsymbol{\lambda}^{\mathrm{disp}}$ | Hydrodynamic thermo-dispersion tensor |
| $\boldsymbol{\lambda}^{\mathrm{cond}}$ | Effective thermal conductivity (from the medium's `thermal_conductivity` MPL property, e.g. `EffectiveThermalConductivityPorosityMixing`) |
| $\boldsymbol{\lambda}^{\mathrm{disp}} = \varrho_f c_f \left[ \alpha_T \lVert\mathbf{q}\rVert \mathbf{I} + (\alpha_L - \alpha_T) \frac{\mathbf{q} \mathbf{q}^T}{\lVert\mathbf{q}\rVert} \right]$ | Thermal dispersivity |
| $\alpha_L$, $\alpha_T$ | Longitudinal and transversal thermo-dispersivities |

### Coupling

The fluid density $\varrho_f(T, p)$ and viscosity $\mu(T)$ couple the hydraulic and thermal equations.
Fluid compressibility $\partial\varrho_f/\partial p$ enters the storage term, and temperature-dependent density drives buoyancy.

**Note:** The classical Boussinesq (Oberbeck--Boussinesq) approximation, where density variations appear only in the buoyancy term, can be recovered by choosing a density model with no pressure dependence (i.e. $\partial\varrho_f/\partial p = 0$, e.g. a constant or temperature-only dependent density) and omitting the `solid_thermal_expansion` configuration.

### Finite element discretization

Both equations are discretized into the standard form $\mathbf{M} \dot{\mathbf{u}} + \mathbf{K} \mathbf{u} = \mathbf{f}$.

#### Pressure equation

$$
\mathbf{M}^p_{ij} = \int_{\Omega} N_i \left( \frac{\phi}{\varrho_f} \frac{\partial \varrho_f}{\partial p} + S \right) N_j \, d\Omega, \qquad
\mathbf{K}^p_{ij} = \int_{\Omega} \nabla N_i^T \frac{\boldsymbol{\kappa}}{\mu} \nabla N_j \, d\Omega
$$

$$
\mathbf{f}^p_i = -\int_{\Omega} \nabla N_i^T \frac{\boldsymbol{\kappa} \varrho_f g}{\mu} \mathbf{e}_z \, d\Omega + \int_{\Omega} N_i \, Q \, d\Omega + \int_{\Gamma_N} N_i \, g_N \, d\sigma
$$

#### Weak formulation

The weak form is obtained by multiplying the strong form by a test function $v \in H_0^1(\Omega)$ and integrating over the domain:

$$
\int_{\Omega} v \left[ \left( \frac{\phi}{\varrho_f} \frac{\partial \varrho_f}{\partial p} + S \right) \frac{\partial p}{\partial t} - \nabla \cdot \left[ \frac{\boldsymbol{\kappa}}{\mu} \nabla \Psi \right] - Q \right] d\Omega + \int_{\Gamma_N} v g_N \, d\sigma = 0
$$

Integration by parts applied to the diffusion term:

$$
\int_{\Omega} v \nabla \cdot \left[ \frac{\boldsymbol{\kappa}}{\mu} \nabla \Psi \right] d\Omega = -\int_{\Omega} \nabla v^T \frac{\boldsymbol{\kappa}}{\mu} \nabla \Psi \, d\Omega + \int_{\Gamma} v \frac{\boldsymbol{\kappa}}{\mu} \nabla \Psi \cdot \mathbf{n} \, d\sigma
$$

Applying the Neumann boundary condition and setting $v=0$ on Dirichlet boundaries yields the weak form:

$$
\int_{\Omega} v \left( \frac{\phi}{\varrho_f} \frac{\partial \varrho_f}{\partial p} + S \right) \frac{\partial p}{\partial t} d\Omega + \int_{\Omega} \nabla v^T \frac{\boldsymbol{\kappa}}{\mu} \nabla \Psi \, d\Omega - \int_{\Omega} v Q \, d\Omega - \int_{\Gamma_N} v g_N \, d\sigma = 0
$$

Substituting the Galerkin discretization $p \approx \sum_j N_j a_j$ and choosing $v = N_i$ recovers the matrix form $\mathbf{M}^p \dot{\mathbf{a}} + \mathbf{K}^p \mathbf{a} = \mathbf{f}^p$ presented above.

When `solid_thermal_expansion` is configured, an additional thermal expansion source term is added to the load vector:

$$
\mathbf{f}^p_{\mathrm{therm},i} = \int_{\Omega} \left[ 3(\alpha_B - \phi)\alpha_s - \frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T} \right] \dot{T} \, N_i \, d\Omega
$$

where $\alpha_B$ is the Biot constant, $\alpha_s$ is the solid thermal expansion coefficient, and $\dot{T}$ is approximated by backward finite differences.

#### Temperature equation

$$
\mathbf{M}^T_{ij} = \int_{\Omega} N_i \, c_p \, N_j \, d\Omega
$$

$$
\mathbf{K}^T_{ij} = \int_{\Omega} \nabla N_i^T \boldsymbol{\lambda} \nabla N_j \, d\Omega + \int_{\Omega} N_i \, \varrho_f c_f \, \mathbf{q}^T \nabla N_j \, d\Omega
$$

## Definition in the project file

The HT process is declared in the `<processes>` block of the project file.

### Monolithic scheme (default)

```xml
<process>
    <name>HydroThermal</name>
    <type>HT</type>
    <integration_order>2</integration_order>
    <process_variables>
        <temperature>T</temperature>
        <pressure>p</pressure>
    </process_variables>
    <specific_body_force>0 -9.81</specific_body_force>
    <secondary_variables>
        <secondary_variable name="darcy_velocity"/>
    </secondary_variables>
</process>
```

### Staggered scheme

To use the staggered coupling scheme, add:

```xml
<coupling_scheme>staggered</coupling_scheme>
```

In the staggered scheme, the heat transport equation (process ID 0) and the hydraulic equation (process ID 1) are solved sequentially.

### Process variables

The HT process requires two process variables: `temperature` and `pressure`.
Each should have 1 component.
For more details, see [Process variables]({{< ref "process_variables" >}}).

## Media properties

The HT process requires properties for the porous medium, the liquid phase, and the solid phase.

### Medium properties

| Property name | Units | SI | Notes |
|---|---|---|---|
| `permeability` | [L$^2$] | [m$^2$] | Intrinsic permeability tensor |
| `porosity` | [-] | [-] | Porous medium porosity |
| `thermal_conductivity` | [M$\cdot$L/(T$^3\cdot\Theta$)] | [W/(m$\cdot$K)] | Effective thermal conductivity of the medium |
| `thermal_longitudinal_dispersivity` | [L] | [m] | Longitudinal thermo-dispersivity $\alpha_L$ |
| `thermal_transversal_dispersivity` | [L] | [m] | Transversal thermo-dispersivity $\alpha_T$ |

### Liquid phase properties

| Property name | Units | SI | Notes |
|---|---|---|---|
| `density` | [M/L$^3$] | [kg/m$^3$] | Fluid mass density $\large^{\star}$ |
| `viscosity` | [M/(L$\cdot$T)] | [Pa$\cdot$s] | Dynamic fluid viscosity $\large^{\star}$ |
| `specific_heat_capacity` | [L$^2$/(T$^2\cdot\Theta$)] | [J/(kg$\cdot$K)] | Specific heat capacity of the fluid |
| `thermal_conductivity` | [M$\cdot$L/(T$^3\cdot\Theta$)] | [W/(m$\cdot$K)] | Thermal conductivity of the fluid |

<small>$\large^{\star}$ Functional dependencies (e.g. on temperature) can be specified. See the [OGS User Guide]({{< ref "/docs/userguide/blocks/media#properties" >}}).</small>

### Solid phase properties

| Property name | Units | SI | Notes |
|---|---|---|---|
| `density` | [M/L$^3$] | [kg/m$^3$] | Solid mass density |
| `specific_heat_capacity` | [L$^2$/(T$^2\cdot\Theta$)] | [J/(kg$\cdot$K)] | Specific heat capacity of the solid |
| `thermal_conductivity` | [M$\cdot$L/(T$^3\cdot\Theta$)] | [W/(m$\cdot$K)] | Thermal conductivity of the solid |
| `storage` | [L$\cdot$T$^2$/M] | [1/Pa] | Storage coefficient |

## Features

### Specific body force

The gravity vector is specified as:

```xml
<specific_body_force>0 -9.81</specific_body_force>
```

### Solid thermal expansion

Optional coupling of thermal expansion effects on the pore pressure:

```xml
<solid_thermal_expansion>
    <thermal_expansion>alpha_s</thermal_expansion>
    <biot_constant>biot</biot_constant>
</solid_thermal_expansion>
```

where `alpha_s` and `biot` are parameter names defined in the `<parameters>` block.

### Aperture size

For lower-dimensional fracture elements, an aperture size parameter can be specified:

```xml
<aperture_size>
    <parameter>fracture_aperture</parameter>
</aperture_size>
```

### Surface flux calculation

Surface flux output can be configured:

```xml
<calculatesurfaceflux>
    ...
</calculatesurfaceflux>
```

### Numerical stabilisation

The HT process supports numerical stabilisation methods for advection-dominated transport.

## Benchmarks

See the OGS benchmark gallery for [Hydro-Thermal examples]({{< ref "/docs/benchmarks/hydro-thermal" >}}).
