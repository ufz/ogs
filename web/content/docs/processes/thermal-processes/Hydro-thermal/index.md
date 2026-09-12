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
- Fluid compressibility and optional solid thermal expansion coupling via the Biot coefficient.
- Hydrodynamic thermal dispersion (longitudinal and transversal dispersivities).
- Numerical stabilisation for advection-dominated problems.
- Fracture flow support via aperture size parameter.
- Surface flux calculation.

### Physical variables

- **Primary variables**: pressure $p$ and temperature $T$.
- **Secondary variable**: Darcy velocity $\mathbf{q}$.

## Theoretical background

Both the flow and heat transport processes are derived from integral conservation laws.
At the moment there is no coupling by source or sink terms, i.e., the coupling is
 implemented only through density changes due to temperature changes in the buoyancy
 term of the groundwater flow. The coupling scheme is referred to as the Boussinesq approximation.

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
&+ \nabla \cdot ({\varrho_f} \mathbf{q}) = Q_H,
\end{align*}
$$

where $S_s$ is the specific storage of the solid phase, $\alpha_B$ is the Biot coefficient, $\alpha_T^s$ is the linear solid thermal expansivity, and $T$ is the temperature.

The value of $S_s$ can be computed as $(\alpha_B - \phi)(1-\alpha_B) / K$, where $K$ is the drained bulk modulus. Equivalently, $S_s = (\alpha_B - \phi) / K_s$ with $K_s = K / (1 - \alpha_B)$ the intrinsic bulk modulus of the solid phase.

In certain scenarios, such as far-field simulations, the fluid density is often assumed constant. Consequently, the volume balance equation can be derived by dividing the mass balance equation by the constant fluid density:
$$
\begin{align*}
\left(\frac{\phi}{\varrho_f} \frac{\partial \varrho_f}{\partial p}  + S_s \right) \frac{\partial p}{\partial t}
&-\overbrace{\left(3(\alpha_B-\phi)\alpha_T^s-\frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T}\right)\frac{\partial T}{\partial t}}^{\text{Thermal expansion}} \\
&+ \nabla \cdot ( \mathbf{q}) = Q_H/{\varrho_f},
\end{align*}
$$

**Note:** In OGS, the fluid part of the thermal expansion term is always
present; its solid part $3(\alpha_B-\phi)\alpha_T^s$ is optional.

### Heat transport equation

$$
c_p \frac{\partial T}{\partial t} - \nabla \cdot (\boldsymbol{\lambda} \nabla T) + \varrho_f c_f \langle \mathbf{q}, \nabla T \rangle = Q_T,
$$

where:

| Symbol | Definition |
| --- | --- |
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
&- \int_{\Omega} {\varrho_f} \nabla v \cdot  \mathbf{q} \mathrm{d}\Omega+ \int_{\Gamma} {\varrho_f} \mathbf{q}\cdot\mathbf{n}v \mathrm{d}\Gamma = \int_{\Omega}Q_H v \mathrm{d}\Omega,
\end{align*}
$$
for the mass balance, and
$$
\begin{align*}
\int_{\Omega}  \left(\frac{\phi}{\varrho_f} \frac{\partial \varrho_f}{\partial p}  + S_s \right) \frac{\partial p}{\partial t} v \mathrm{d}\Omega
&-\int_{\Omega}\left(3(\alpha_B-\phi)\alpha_T^s-\frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T}\right)\frac{\partial T}{\partial t}v \mathrm{d}\Omega \\
&- \int_{\Omega}  \nabla v \cdot  \mathbf{q} \mathrm{d}\Omega+ \int_{\Gamma}  \mathbf{q}\cdot\mathbf{n}v \mathrm{d}\Gamma = \int_{\Omega}Q_H v/\varrho_f \mathrm{d}\Omega,
\end{align*}
$$
for the volume balance, where $\mathbf{n}$ is the outward unit normal to the domain boundary $\Gamma$.

**Note:** The term $\frac{\phi}{\varrho_f} \frac{\partial \varrho_f}{\partial p}$ represents the fluid compressibility contribution to the storage coefficient.

For the temperature field, the weak form is:
$$
\begin{align*}
\int_{\Omega}  c_p \frac{\partial T}{\partial t}v \mathrm{d}\Omega + &
\int_{\Omega}  \boldsymbol{\lambda} \nabla T \cdot \nabla  v\mathrm{d}\Omega
+\int_{\Gamma} \boldsymbol{\lambda} \nabla T \cdot \mathbf{n} v \mathrm{d}\Gamma \\
+& \int_{\Omega}  \varrho_f c_f \langle \mathbf{q}, \nabla T \rangle v\mathrm{d}\Omega = \int_{\Omega}  Q_T v\mathrm{d}\Omega.
\end{align*}
$$

### Finite element discretization

Both equations are discretized into the standard form $\mathbf{M} \dot{\mathbf{u}} + \mathbf{K} \mathbf{u} = \mathbf{f}$ by substituting the Galerkin discretization $p \approx \sum_j N_j a_j$ and choosing $v = N_i$.

#### Pressure equation

For the mass balance equation, the resulting discretized matrices and vectors are:
$$
\mathbf{M}^p_{ij} = \int_{\Omega}  \varrho_f\left( \frac{\phi}{\varrho_f} \frac{\partial \varrho_f}{\partial p} + S_s \right) N_i N_j \, \mathrm{d}\Omega, \qquad
\mathbf{K}^p_{ij} = \int_{\Omega} \varrho_f\nabla N_i^T \frac{\boldsymbol{\kappa}}{\mu} \nabla N_j \, \mathrm{d}\Omega,
$$

$$
\mathbf{f}^p_i = \int_{\Omega} \varrho_f^2 \nabla N_i^T \frac{\boldsymbol{\kappa} }{\mu} \mathbf{g} \, \mathrm{d}\Omega + \int_{\Omega}  Q_H \, N_i \, \mathrm{d}\Omega + \int_{\Gamma} \varrho_f \mathbf{q}\cdot\mathbf{n}N_i \mathrm{d}\Gamma.
$$

In OGS, the thermal expansion term is always assembled; only its solid part $3(\alpha_B - \phi)\alpha_T^s$ is optional and drops out when the solid thermal expansivity is undefined. It is added to the right hand vector for the staggered scheme as

$$
\mathbf{f}^p_{\mathrm{therm},i} = \int_{\Omega} \varrho_f \left[ 3(\alpha_B - \phi)\alpha_T^s - \frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T} \right] \dot{T} \, N_i \, d\Omega,
$$

and it is added to the mass matrix for the monolithic scheme as
$$
\mathbf{M}^{pT}_{ij} = -\int_{\Omega}\varrho_f \left[ 3(\alpha_B - \phi)\alpha_T^s - \frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T} \right] N_i  N_j \, \mathrm{d}\Omega.
$$

The corresponding terms for the volume balance equation are obtained by dividing all the above integrals by the fluid density.

#### Temperature equation

$$
\mathbf{M}^T_{ij} = \int_{\Omega} N_i \, c_p \, N_j \, d\Omega,
$$

$$
\mathbf{K}^T_{ij} = \int_{\Omega} \nabla N_i^T \boldsymbol{\lambda} \nabla N_j \, d\Omega + \int_{\Omega} N_i \, \varrho_f c_f \, \mathbf{q}^T \nabla N_j \, d\Omega.
$$

The advection term is part of $\mathbf{K}^T$ in both schemes; they differ only in
where the Darcy velocity $\mathbf{q}$ comes from, namely from the pressure
solution of the same equation system in the monolithic scheme and from the
previous solution of the pressure equation in the staggered scheme.

The right hand side vector is given by
$$
\mathbf{f}^T_{i} = \int_{\Omega} Q_T \, N_i \, d\Omega.
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

By default, the value is `volume`. The volume balance requires the liquid phase `density` property to be of type `Constant`; with any other property type, for instance `Linear`, `Function` or `WaterDensityIAPWSIF97Region1`, OGS issues a fatal error prompting the user to add the tag. In such cases, the input values for Neumann boundary conditions and source/sink terms must be adjusted accordingly, for example:

- Neumann condition: from volume rate per area (SI unit: $\text{m}\cdot\text{s}^{-1}$) to mass rate per area (SI unit: $\text{kg}\cdot\text{m}^{-2}\text{s}^{-1}$).
- source/sink: from volume rate (SI unit: $\text{m}^{3}\text{s}^{-1}$) to mass rate (SI unit: $\text{kg}\cdot\text{s}^{-1}$).

The choice of the balance equation type is orthogonal to the material property combinations listed in *Material property restrictions* below: the mass balance scales $\mathbf{M}^p$, $\mathbf{K}^p$, $\mathbf{f}^p$ and the thermal expansion term by $\varrho_f$ uniformly, and it does not change which combination applies.

### Thermal expansion in the mass/volume balance equation

The fluid part of the thermal expansion term, $-\frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T}$, is always computed. The solid part, $3(\alpha_B - \phi)\alpha_T^s$, is added only if the linear solid thermal expansivity (property name: `thermal_expansivity`) and the Biot's coefficient (property name: `biot_coefficient`) are defined in the input project file. If linear solid thermal expansivity is anisotropic, the average of its components can be used.  

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
| --- | --- | --- | --- |
| `permeability` | [L$^2$] | [m$^2$] | Intrinsic permeability tensor |
| `porosity` | [-] | [-] | Porous medium porosity |
| `thermal_conductivity` | [M$\cdot$L/(T$^3\cdot\Theta$)] | [W/(m$\cdot$K)] | Effective thermal conductivity of the medium |
| `thermal_longitudinal_dispersivity` | [L] | [m] | Longitudinal thermo-dispersivity $\alpha_L$ |
| `thermal_transversal_dispersivity` | [L] | [m] | Transversal thermo-dispersivity $\alpha_T$ |
| `biot_coefficient` | - | - | Only for the thermal expansion computation |

### Material property restrictions for Biot's coefficient, porosity, storage, and thermal expansivity combinations

The Biot coefficient and the solid thermal expansivity are optional material parameters used to account for the solid phase contribution to thermal expansion. Once these parameters are set in the project file, they are associated with the specific storage of the solid phase $S_s$, which is input via `storage`, and the effective thermal expansivity $\beta_{\mathrm{eff}}=\left(3(\alpha_B-\phi)\alpha_T^s-\frac{\phi}{\varrho_f}\frac{\partial \varrho_f}{\partial T}\right)$. This can be treated as a parameter combination when both the Biot coefficient and the solid thermal expansivity are defined. The rules are:

1. $\alpha_B$ is read only when $\alpha_T^s$ is defined. Therefore:
   - $\alpha_T^s$ defined, $\alpha_B$ undefined → fatal error (undefined property).
   - $\alpha_B$ defined, $\alpha_T^s$ undefined → $\alpha_B$ is silently ignored (never read).
2. The specific storage $S_s$ given by the solid-phase `storage` property is always read.
3. The combination $\alpha_B = 1$ and $S_s \neq 0$ leads to a fatal error because
   $S_s=(\alpha_B-\phi)(1-\alpha_B)/K$. The guard is only reached when
   $\alpha_T^s$ is defined, since that is when $\alpha_B$ is read at all. The
   evaluated values are compared once per element, at its integration points,
   during the initialisation of the local assemblers, at $t=0$ and
   $\Delta t=0$. That is the literal time zero, not the initial time of the
   time loop, which may be a different value. This covers every property that depends on
   neither the primary variables nor the time, whatever its type. A property
   depending on the time is compared at $t=0$ only, which the simulation need
   not pass through at all, so a violation occurring at any later time is not
   caught here. A property that does depend on the primary variables evaluates
   to NaN at that point and is not rejected: a Biot coefficient of NaN makes the
   comparison $\alpha_B = 1$ false, and a specific storage of NaN is excluded
   explicitly. Both cases are decided by the same comparison, which is
   repeated at every integration point during the assembly, where the
   properties have values. With $\alpha_T^s$ undefined the combination is
   accepted silently, because $\alpha_B$ is then never read.

The ranges of these four parameters are listed in the following table:

| Quantity | Physical range | Code guard? |
| --- | --- | --- |
| Porosity $\phi$ | $(0,1)$ | none |
| Biot $\alpha_B$ | $[\phi,\,1]$ | only $\alpha_B=1 \Rightarrow S_s=0$, and only when $\alpha_T^s$ is defined |
| Specific storage of the solid phase $S_s$ | $\geq 0$ | none |
| Solid thermal expansivity $\alpha_T^s$ | $\geq 0$ | none |

This means that if these parameters are left unchecked and fall outside the stated ranges, only the combination $\alpha_B=1,\,S_s\neq0$ is detected, and only when $\alpha_T^s$ is defined.

The relation $S_s=(\alpha_B-\phi)(1-\alpha_B)/K$ also constrains the reverse
direction: $S_s=0$ with $\alpha_B<1$ implies $\alpha_B=\phi$ or $K=\infty$.
That combination is not guarded, so row 4 below runs without a warning even
though only $\alpha_B=1$ makes $S_s=0$ consistent with that relation. Row 1
carries no such inconsistency, since $\alpha_B$ is undefined there and never
read; its risk is the vanishing storage term alone.

The following table summarises the behaviour of the different combinations of the four parameters (Legend: **✓ defined**, **✗ undefined**):

| # | $\alpha_T^s$ | $\alpha_B$ | $\alpha_B$ value | $S_s$ | $\beta_\mathrm{eff}$ | Storage $\mathbf{M}^p$ | Outcome |
| --- | :---: | :---: | --- | --- | --- | --- | --- |
| 1 | ✗ | ✗ | — | $0$ | fluid only: $-\phi\,\frac{\partial \varrho_f}{\partial T}/\varrho_f$ | fluid compressibility only | Runs. **Risk:** if $\varrho_f$ depends on $T$ only (no $p$) and $S_s=0$, the storage term vanishes, which is physically inconsistent and may lead to numerical instability. |
| 2 | ✗ | ✗ | — | $>0$ | fluid only | fluid compressibility $+\,S_s$ | Valid. Solid thermal expansion ignored. |
| 3 | ✗ | ✓ | any | any | fluid only | fluid compressibility $+\,S_s$ | Runs. $\alpha_B$ is **silently ignored** (never read since $\alpha_T^s$ is unset). |
| 4 | ✓ | ✓ | $[\phi,1)$ | $0$ | fluid $+\ 3(\alpha_B-\phi)\alpha_T^s$ | fluid compressibility only | Runs, unguarded. Storage is carried entirely by fluid compressibility, but $S_s=0$ with $\alpha_B<1$ is only consistent for $\alpha_B=\phi$ or $K=\infty$. Carries the same **Risk** as row 1. |
| 5 | ✓ | ✓ | $[\phi,1)$ | $>0$ | fluid $+\ 3(\alpha_B-\phi)\alpha_T^s$ | fluid compressibility $+\,S_s$ | Valid. General poroelastic-consistent case. |
| 6 | ✓ | ✓ | $=1$ | $0$ | fluid $+\ 3(1-\phi)\alpha_T^s$ | fluid compressibility only | Valid. Incompressible solid grains ($\alpha_B=1 \Rightarrow S_s=0$). Carries the same **Risk** as row 1: with a temperature-only $\varrho_f$ the storage term vanishes identically. |

### Liquid phase properties

| Property name | Units | SI | Notes |
| --- | --- | --- | --- |
| `density` | [M/L$^3$] | [kg/m$^3$] | Fluid mass density $\large^{\star}$ |
| `viscosity` | [M/(L$\cdot$T)] | [Pa$\cdot$s] | Dynamic fluid viscosity $\large^{\star}$ |
| `specific_heat_capacity` | [L$^2$/(T$^2\cdot\Theta$)] | [J/(kg$\cdot$K)] | Specific heat capacity of the fluid |
| `thermal_conductivity` | [M$\cdot$L/(T$^3\cdot\Theta$)] | [W/(m$\cdot$K)] | Thermal conductivity of the fluid |

<small>$\large^{\star}$ Functional dependencies (e.g. on temperature) can be specified. See the [OGS User Guide]({{< ref "/docs/userguide/blocks/media#properties" >}}).</small>

### Solid phase properties

| Property name | Units | SI | Notes |
| --- | --- | --- | --- |
| `density` | [M/L$^3$] | [kg/m$^3$] | Solid mass density |
| `specific_heat_capacity` | [L$^2$/(T$^2\cdot\Theta$)] | [J/(kg$\cdot$K)] | Specific heat capacity of the solid |
| `thermal_conductivity` | [M$\cdot$L/(T$^3\cdot\Theta$)] | [W/(m$\cdot$K)] | Thermal conductivity of the solid |
| `storage` | [L$\cdot$T$^2$/M] | [1/Pa] | Storage coefficient |
| `thermal_expansivity` | [1/T] | [1/K] | Only for the thermal expansion computation |

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
