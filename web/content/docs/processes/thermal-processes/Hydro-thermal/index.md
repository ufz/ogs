+++
author = "Thomas Fischer, Dmitri Naumov, Fabien Magri, Marc Walther, Tianyuan Zheng, Olaf Kolditz, Wenqing Wang"
date = "2026-04-04"
title = "Hydro-Thermal Process (HT)"
weight = 3
+++

## Introduction

The Hydro-Thermal (HT) process models coupled groundwater flow and heat transport in porous media.
Both subprocesses are governed by coupled parabolic partial differential equations, linked via temperature dependent fluid properties and advective heat transport.

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

Both the flow and heat transport processes are derived from integral conservation laws.

### Mass balance equation

The Darcy velocity is given by:

$$
\mathbf{q} = -\frac{\boldsymbol{\kappa}}{\mu} \left( \nabla p - \varrho_f  \mathbf g \right),
$$
where $\boldsymbol{\kappa}$ is the permeability, $p$ is the pore pressure, $\varrho_f$ is the liquid density, and $\mathbf{g}$ is the gravitational force.  

The mass balance equation reads:
$$
 \frac{\partial }{\partial t} (\phi \varrho_f) + \nabla \cdot ({\varrho_f} \mathbf{q}) = Q_H,
$$
where $Q_H$ is the sink or source term.

Assuming the solid matrix is compressible (i.e., $\frac{\partial \phi }{\partial t}\neq0$), the mass balance equation can be expanded as

$$
\begin{align*}
\left(\phi  \frac{\partial \varrho_f}{\partial p}  + S_s {\varrho_f}\right) \frac{\partial p}{\partial t}
&-\overbrace{\left(3(\alpha_B-\phi)\alpha_T^s-\frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T}\right)\varrho_f\frac{\partial T}{\partial t}}^{\text{Thermal expansion}} \\
&+ \nabla \cdot ({\varrho_f} \mathbf{q}) = Q_H
\end{align*},
$$

where $S_s$ is the solid compressibility, $\alpha_B$ is the Biot coefficient, $\alpha_T^s$ is the linear solid thermal expansivity, and $T$ is the temperature.

The value of $S_s$ can be computed as $(1 - \phi) / K_s$, where $K_s$ is the intrinsic bulk modulus of the solid phase.

In certain scenarios, such as far-field simulations, the fluid density is often assumed constant. Consequently, the volume balance equation can be derived by scaling the mass balance equation by the constant fluid density:
$$
\begin{align*}
\left(\phi \frac{1}{\varrho_f} \frac{\partial }{\partial p}  + S_s {\varrho_f}\right) \frac{\partial p}{\partial t}
&-\overbrace{\left(3(\alpha_B-\phi)\alpha_T^s-\frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T}\right)\frac{\partial T}{\partial t}}^{\text{Thermal expansion}} \\
&+ \nabla \cdot ( \mathbf{q}) = Q_H/{\varrho_f}
\end{align*},
$$

**Note:** In OGS, the thermal expansion term is optional.

### Heat transport equation

$$
c_p \frac{\partial T}{\partial t} - \nabla \cdot (\boldsymbol{\lambda} \nabla T) + \varrho_f c_f \langle \mathbf{q}, \nabla T \rangle = Q_T,
$$

where:

| Symbol | Definition |
|---|---|
| $c_p = \varrho_f \phi c_f + \varrho_s (1 - \phi) c_s$ | Volumetric heat capacity of the mixture |
| $\boldsymbol{\lambda} = \boldsymbol{\lambda}^{\mathrm{cond}} + \boldsymbol{\lambda}^{\mathrm{disp}}$ | Hydrodynamic thermo-dispersion tensor |
| $\boldsymbol{\lambda}^{\mathrm{cond}}$ | Effective thermal conductivity (from the medium's `thermal_conductivity` MPL property, e.g. `EffectiveThermalConductivityPorosityMixing`) |
| $\boldsymbol{\lambda}^{\mathrm{disp}} = \varrho_f c_f \left[ \alpha_T \lVert\mathbf{q}\rVert \mathbf{I} + (\alpha_L - \alpha_T) \frac{\mathbf{q} \mathbf{q}^T}{\lVert\mathbf{q}\rVert} \right]$ | Thermal dispersivity |
| $\alpha_L$, $\alpha_T$ | Longitudinal and transversal thermo-dispersivities |
| $Q_T$ | Source or sink term |

### Weak formulation

The weak form is obtained by multiplying the strong form by a test function $v \in H_0^1(\Omega)$ and integrating over the domain.

For the pressure field, the weak form is:
$$
\begin{align*}
\int_{\Omega}  \left(\phi  \frac{\partial \varrho_f}{\partial p}  + S_s {\varrho_f}\right) \frac{\partial p}{\partial t} v \mathrm{d}\Omega
&-\int_{\Omega}\left(3(\alpha_B-\phi)\alpha_T^s-\frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T}\right)\varrho_f\frac{\partial T}{\partial t}v \mathrm{d}\Omega \\
&- \int_{\Omega} {\varrho_f} \nabla v \cdot  \mathbf{q} \mathrm{d}\Omega+ \int_{\Gamma} {\varrho_f} \mathbf{q}\cdot\mathbf{n}v \partial\Omega = \int_{\Omega}Q_H v \mathrm{d}\Omega
\end{align*},
$$
for the mass balance, and
$$
\begin{align*}
\int_{\Omega}  \left(\phi \frac{1}{\varrho_f} \frac{\partial }{\partial p}  + S_s {\varrho_f}\right) \frac{\partial p}{\partial t} v \mathrm{d}\Omega
&-\int_{\Omega}\left(3(\alpha_B-\phi)\alpha_T^s-\frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T}\right)\frac{\partial T}{\partial t}v \mathrm{d}\Omega \\
&- \int_{\Omega}  \nabla v \cdot  \mathbf{q} \mathrm{d}\Omega+ \int_{\Gamma}  \mathbf{q}\cdot\mathbf{n}v \partial\Omega = \int_{\Omega}Q_H v/\varrho_f \mathrm{d}\Omega
\end{align*},
$$
for the volume balance, where $\mathbf{n}$ is the outward unit normal to the domain boundary $\Gamma$.

For the temperature field, the weak form is:
$$
\begin{align*}
\int_{\Omega}  c_p \frac{\partial T}{\partial t}v \mathrm{d}\Omega + &
\int_{\Omega}  \boldsymbol{\lambda} \nabla T \cdot \nabla  v\mathrm{d}\Omega
+\int_{\Gamma} \boldsymbol{\lambda} \nabla T \cdot \mathbf{n} v\partial\Omega \\
+& \int_{\Omega}  \varrho_f c_f \langle \mathbf{q}, \nabla T \rangle v\mathrm{d}\Omega = \int_{\Omega}  Q_T v\mathrm{d}\Omega.
\end{align*}
$$

### Finite element discretization

Both equations are discretized into the standard form $\mathbf{M} \dot{\mathbf{u}} + \mathbf{K} \mathbf{u} = \mathbf{f}$ by substituting the Galerkin discretization $p \approx \sum_j N_j a_j$ and choosing $v = N_i$.

#### Pressure equation

For the mass balance equation, the resulting discretized matrices and vectors are:
$$
\mathbf{M}^p_{ij} = \int_{\Omega}  \varrho_f\left( \frac{\phi}{\varrho_f} \frac{\partial \varrho_f}{\partial p} + S_s \right) N_i N_j \, \mathrm{d}\Omega, \qquad
\mathbf{K}^p_{ij} = \int_{\Omega} \varrho_f\nabla N_i^T \frac{\boldsymbol{\kappa}}{\mu} \nabla N_j \, \mathrm{d}\Omega
$$

$$
\mathbf{f}^p_i = \int_{\Omega} \varrho_f^2 \nabla N_i^T \frac{\boldsymbol{\kappa} }{\mu} \mathbf{g} \, \mathrm{d}\Omega + \int_{\Omega}  Q_H \, N_i \, \mathrm{d}\Omega + \int_{\Gamma} \varrho_f \mathbf{q}\cdot\mathbf{n}N_i \partial\Omega
$$

In OGS, the thermal expansion term is optional. When configured, it is added to the right hand vector for the staggered scheme as

$$
\mathbf{f}^p_{\mathrm{therm},i} = \int_{\Omega} \varrho_f \left[ 3(\alpha_B - \phi)\alpha_T^s - \frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T} \right] \dot{T} \, N_i \, d\Omega,
$$

and it is added to the mass matrix for the monolithic scheme as
$$
\mathbf{M}^{pT}_{ij} = -\int_{\Omega}\varrho_f \left[ 3(\alpha_B - \phi)\alpha_T^s - \frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T} \right] N_i  N_j \, \mathrm{d}\Omega.
$$

The corresponding terms for the volume balance equation are obtained by scaling all the above integrals by the fluid density.

#### Temperature equation

$$
\mathbf{M}^T_{ij} = \int_{\Omega} N_i \, c_p \, N_j \, d\Omega
$$

$$
\mathbf{K}^T_{ij} =
\begin{cases}
\int_{\Omega} \nabla N_i^T \boldsymbol{\lambda} \nabla N_j \, d\Omega + \int_{\Omega} N_i \, \varrho_f c_f \, \mathbf{q}^T \nabla N_j \, d\Omega, \quad\text{for the monolithic scheme}\\
\int_{\Omega} \nabla N_i^T \boldsymbol{\lambda} \nabla N_j \, d\Omega, \quad\text{for the staggered scheme}
\end{cases}
$$

The right hand side vector is given by
$$
\mathbf{f}^T_{i} =
\begin{cases}
\int_{\Omega} Q_T \, N_i \, d\Omega, \quad\text{for the monolithic scheme}\\
\int_{\Omega} Q_T \, N_i \, d\Omega -  \int_{\Omega}  \varrho_f c_f \langle \mathbf{q}, \nabla T \rangle \, N_i \, d\Omega, \quad\text{for the staggered scheme}
\end{cases}
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

### Balance equation type: mass or volume

Similar to other hydraulics-related processes, HT supports an optional `<equation_balance_type>` tag with values `mass` or `volume` to select the balance equation type for the hydraulics process. For example:

```xml
<process>
    <name>HydroThermal</name>
    <type>HT</type>
    <integration_order>2</integration_order>
    <equation_balance_type>mass</equation_balance_type>
    ...
</process>
```

By default, the value is `volume`. If the project file uses a nonlinear fluid density model without specifying this tag, OGS will issue a fatal error prompting the user to add the tag. In such cases, the input values for Neumann boundary conditions and source/sink terms must be adjusted accordingly, for example:

- Neumann condition: from volume rate per area (SI unit: $\text{m}\cdot\text{s}^{-1}$) to mass rate per area (SI unit: $\text{kg}\cdot\text{m}^{-2}\text{s}^{-1}$).
- source/sink: from volume rate (SI unit: $\text{m}^{3}\text{s}^{-1}$) to mass rate (SI unit: $\text{kg}\cdot\text{s}^{-1}$).

### Thermal expansion in the mass/volume balance equation

The computation of the thermal expansion term is enabled only if the linear solid thermal expansivity (property name: `thermal_expansivity`) and the Biot's coefficient (property name: `biot_coefficient`) are defined in the input project file. If linear solid thermal expansivity is anisotropic, the average of its components can be used.  

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
| `biot_coefficient` | - | - |Only for the thermal expansion computation |

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
| `thermal_expansivity` | [1/T] | [1/K] | Only  the thermal expansion computation |

## Features

### Specific body force

The gravity vector is specified as:

```xml
<specific_body_force>0 -9.81</specific_body_force>
```

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
