+++
date = "2026-08-02"
title = "Liquid Flow (LIQUID_FLOW)"
author = "Philipp Selzer, Tom Fischer, Olaf Kolditz"
weight = 1
+++

## Introduction

The Liquid-Flow-Process models saturated single-phase flow, which can be of variable density, in porous and fractured media.
It may be used for modeling water flow but can be equally used to model flow of other liquids and also gas as long as flow follows the Darcy regime.
However, the natural application of the Liquid-Flow-Process is in modeling groundwater flow.
As such, the documentation above specifically targets hydrogeologists.
In the context of groundwater modeling, it is important to note that the pressure-based variant of the groundwater flow equation is implemented in OGS, which may not be the most intuitive for hydrogeologists.
To this end, the documentation above aims, among other things, to give conversions between the pressure-based and the hydraulic-head based formulations of the groundwater flow equation as well as it outlines how to use the Liquid-Flow-Process to natively model the hydraulic head-based formulation of the groundwater flow equation commonly used in hydrogeology.

### Physical key features

- Fluid flow: liquid (incompressible or compressible) flow; gas (compressible) flow is also supported when it follows the Darcy regime.
- Primary variable: Pressure (in most applications liquid-phase pressure, i.e. pressure of the water phase), in alternative parameterization the pressure can be converted to hydraulic head.
- Secondary variable: Darcy velocity (aka "specific discharge").

## Theoretical background

The `LiquidFlow` process can solve the fluid-flow equation in two variants.

### Mass balance formulation

The more general formulation is the mass-balance version of the equation:

$$
\left( \rho_w \frac{\partial n_e}{\partial p_w} + n_e \frac{\partial \rho_w }{\partial p_w} \right) \frac{\partial p_w}{\partial t}  -\nabla \cdot \left( \rho_w \frac{\mathbf{K}_0}{\mu_w} \left( \nabla p_w - \rho_w \mathbf{g} \right)\right) = \rho_w W_0. \qquad \text{(mass)}
$$

This variant accounts for compressible (density-dependent) flow. It is selected by including the following tag in the project file:

```xml
<equation_balance_type>mass</equation_balance_type>
```

### Volume balance formulation

Dividing the mass balance equation by $\rho_w$ and assuming density is constant in space (so that $\nabla \cdot (\rho_w \, \cdot) \approx \rho_w \, \nabla \cdot (\cdot)$) gives the volume-balance form. Substituting the OGS storage coefficient $\beta_s = \partial n_e / \partial p_w$:

$$
\left(\beta_s + \frac{n_e}{\rho_w} \frac{\partial \rho_w }{\partial p_w} \right) \frac{\partial p_w}{\partial t}  -\nabla \cdot \left(\frac{\mathbf{K}_0}{\mu_w} \left( \nabla p_w - \rho_w \mathbf{g} \right)\right) = W_0.
$$

In the code, the volume balance mode **enforces a constant (pressure-independent) fluid density**, i.e.\ $\partial \rho_w / \partial p_w \equiv 0$. The equation therefore simplifies to:

$$
\beta_s \frac{\partial p_w}{\partial t}  -\nabla \cdot \left(\frac{\mathbf{K}_0}{\mu_w} \left( \nabla p_w - \rho_w \mathbf{g} \right)\right) = W_0. \qquad \text{(volume, default)} \quad (1)
$$

This is the **default** formulation. To explicitly select it, set the following tag in the project file:

```xml
<equation_balance_type>volume</equation_balance_type>
```

with the Darcy velocity as secondary variable:

$$
\begin{equation}
\mathbf{q} = - \frac{\mathbf{K_0}}{\mu_w} \left( \nabla p_w - \rho_w \mathbf{g} \right)
\end{equation}
$$

where the subscript $w$ denotes "water" and in the following table all variables used for the former and further equations are listed:

| Symbols | Units | SI | Definition |
|---|---|---|---|
| $S_0$ | [1/L] | [1/m] | specific storage coefficient in hydrogeology |
| $\rho_w$ , ($\rho^f$) | [M/L³] | [kg/m³] | fluid density |
| $g$ | [L/T²] | [m/s²] | gravitational acceleration coefficient |
| $\beta_s = \partial n_e / \partial p_w$ | [LT²/M] | [1/Pa] |storage coefficient of  OGS |
| $p_w$ , ($p^f$) | [M/LT²] | [Pa] | liquid phase (i.e. water) pressure |
| $t$ | [T] | [s] | time |
| $\nabla$ | [1/L] | [1/m] | divergence, gradient operator |
| $\mathbf{K_0}$ | [L²] | [m²] | intrinsic permeability |
| $\mathbf{g}$ | [L/T²] | [m/s²] | gravitational acceleration vector  (0,0,g) |
| $\mu_w$ | [M/LT] | [Pa s] | dynamic fluid viscosity |
| $n_e$ | [-] | [-] | effective porosity |
| $h$ | [L] | [m] | hydraulic head |
| $W_0$ | [1/T] | [1/s] | source/sink term |

Please note that in OGS the gravitational acceleration vector is defined as $\mathbf{b} = - \mathbf{g}$. Thus, if you want to work with the pressure-based variant of groundwater flow, the input for the gravitational acceleration in 3D is:

```xml
<specific_body_force> 0 0 -9.81</specific_body_force>
```

### Hydraulic head parameterisation

It is important to note that OGS also has a "storage" coefficient, which is defined as such: $\beta_s = \partial n_e / \partial p_w$. One can show that for equation (1) under the assumption $\partial \rho_w / \partial p_w = 0$ the following relation holds: $\beta_s = S_0/(\rho_w g)$, relating the specific storage of hydrogeology to the storage implemented in OGS.

One can use equation (1) to emulate the hydrogeological variant of the groundwater equation by applying the following measures:

- Ignore the partial derivative $\partial \rho_w / \partial p_w$, then OGS will set it to zero
- Give $(0,0,0)^T$ as gravitational vector $g$ using the tag "specific\_body\_force"
- Set $\rho_w = 1$
- Set $\mu_w=1$

Then:

- The liquid-phase pressure, $p_w$, will be equal to hydraulic head $h$
- The intrinsic permeability $\mathbf{K}_0$ will be the hydraulic conductivity $\mathbf{K}$
- The storage of OGS, $\beta_s$, will be equal to the specific storage $S_0$

Thus, you effectively model the groundwater flow equation in the following form:

$$
S_0 \frac{\partial h}{\partial t} - \nabla \cdot (\mathbf{K} \nabla h) = W_0. \qquad (2)
$$

More details how to use the liquid flow equation see the detailed theory documentation in: https://gitlab.opengeosys.org/ogs/documentation/liquidflow/-/jobs/artifacts/main/raw/main.pdf?job=build

### Finite Element Discretization

The Liquid-Flow equation is discretized using the finite element method.
The discrete system takes the form:

$$
\begin{equation}
\mathbf{M}_e \frac{d\mathbf{p}_e}{dt} + \mathbf{K}_e \mathbf{p}_e = \mathbf{b}_e
\end{equation}
$$

where the element matrices depend on the chosen equation balance type as follows.

#### Element mass matrix

|  | Volume balance (default) | Mass balance |
|---|---|---|
| $\mathbf{M}_e$ | $\displaystyle\int_{\Omega^e} \mathbf{N}^T \beta_s \mathbf{N} \, d\Omega$ | $\displaystyle\int_{\Omega^e} \mathbf{N}^T \!\left(\rho_w \beta_s + n_e \frac{\partial \rho_w}{\partial p_w}\right)\! \mathbf{N} \, d\Omega$ |

#### Element conductance matrix

|  | Volume balance (default) | Mass balance |
|---|---|---|
| $\mathbf{K}_e$ | $\displaystyle\int_{\Omega^e} (\nabla \mathbf{N})^T \frac{\mathbf{K_0}}{\mu} (\nabla \mathbf{N}) \, d\Omega$ | $\displaystyle\int_{\Omega^e} (\nabla \mathbf{N})^T \frac{\rho_w \mathbf{K_0}}{\mu} (\nabla \mathbf{N}) \, d\Omega$ |

#### Element load vector

|  | Volume balance (default) | Mass balance |
|---|---|---|
| $\mathbf{b}_e$ | $\displaystyle\int_{\Omega^e} (\nabla \mathbf{N})^T \frac{\mathbf{K_0}}{\mu} \rho_w \mathbf{g} \, d\Omega$ | $\displaystyle\int_{\Omega^e} (\nabla \mathbf{N})^T \frac{\rho_w^2 \mathbf{K_0}}{\mu} \mathbf{g} \, d\Omega$ |

## OGS project `prj` file

Liquid Flow process has to be declared in the project file `prj` in the processes block. For example in following way:

```xml
<process>
    <name>LiquidFlow</name>
    <type>LIQUID_FLOW</type>
    <integration_order>2</integration_order>
    <process_variables>
        <process_variable>pressure</process_variable>
    </process_variables>
    <secondary_variables>
        <secondary_variable internal_name="darcy_velocity" output_name="v"/>
    </secondary_variables>
    <specific_body_force>0.0 0.0</specific_body_force>
</process>
```

### Media properties

For Liquid-Flow fluid and porous medium properties need to be specified.

#### Fluid properties

| Property name | Units | SI | Notes |
|---|---|---|---|
| `density` | [M/L³] | [kg·m⁻³] | Fluid mass density: functional dependencies can be specified |
| `viscosity` | [M/LT] | [Pa·s] | Dynamic fluid viscosity: functional dependencies can be specified |

#### Medium properties

| Property name | Units | SI | Notes |
|---|---|---|---|
| `porosity` | [-] | [-] | Porous medium porosity: functional dependencies can be specified |
| `permeability` | [L²] | [m²] | Intrinsic permeability: functional dependencies can be specified |
| `reference_temperature` | [T] | [K] | Reference temperature for fluid properties |
| `storage` | [LT²/M] | [m·s²/kg] | Specific storage coefficient $\beta$ |

## Benchmarks

Liquid-Flow examples from the OGS benchmark gallery you find here: https://www.opengeosys.org/6.5.7/docs/benchmarks/liquid-flow/

## Workflows

Shortly we will provide complete workflows (Jupyter Notebooks) to set up Liquid-Flow examples from the scratch including geometric description, meshing, project generation, simulation and displaying results. See also OGSTools: https://ogstools.opengeosys.org/stable/auto_examples/index.html
