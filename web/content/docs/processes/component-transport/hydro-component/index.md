+++
date = "2022-03-21T12:00:13+01:00"
title = "ComponentTransport"
author = "Haibing Shao, Renchao Lu"
weight = 4

[menu]
  [menu.userguide]
    parent = "process-dependent-configuration"
+++

## Introduction

`ComponentTransport` process is widely used to predict the distribution of chemical components in the subsurface, which is controlled by the groundwater flow (advection), the hydrodynamic dispersion and diffusion, the sorption on the solid phase, as well as the decay of component.

## Mathematical framework

The governing equation implemented in OGS-6 is the so-called advective and diffusion equation (`ADE`), with the consideration of sorption and decay process. ADE is widely used to describe the concentration of chemical components in groundwater aquifer and porous media. The equations can be solved analytically in (simplified) 1D cases. For more complex geometry, especially with heterogeneous material properties, numerical solution is often preferred. The `ComponentTransport` process has the following assumptions.

* The model domain is a porous media and it is fully saturated with water.
* The fluid velocity in the porous media is assumed to be slow, and the flow process is regulated by Darcy's law.
* The sorption process redistributes the chemical component between the aqueous and solid phase.
* When a decay coefficient is given, the decay of the chemical component follows first-order reaction kinetics equation.

For the flow process, the continuity equation for flowing fluid in a saturated porous medium is as follows
$$
\frac{\partial \left(\phi \rho\right)}{\partial t} + \nabla \cdot \left(\textbf{q} \rho\right) + Q_{p} = 0,
$$
where $\phi$ [-] is the porosity, $\rho$ [kg/m$^3$] is the fluid density, $t$ [s] is the time, $\textbf{q}$ [m/s] is the Darcy flux with laminar flow assumptions, and $Q_{p}$ [kg/m$^3$/s] is the source-sink term. According to Darcy's Law, the flux $\textbf{q}$ is related to the pressure drop and gravitational body force through
$$
\textbf{q} = - \frac{\textbf{k}}{\mu} \left(\nabla p - \rho \textbf{g}\right),
$$
where $\textbf{k}$ [m$^2$] is the intrinsic permeability, $\mu$ [Pa$\cdot$s] is the fluid dynamic viscosity, and $\textbf{g}$ [m/s$^2$] is the gravity vector.

For each chemical component $\alpha = 1, .., N_p$, its corresponding advective and diffusion equation (`ADE`) reads,
$$
\begin{equation}
\frac{\partial \left(\phi R c_{\alpha}\right)}{\partial t} + \nabla \cdot \left( \textbf{q} c_{\alpha} - \textbf{D} \nabla c_{\alpha} \right) + Q_{c_{\alpha}} + \phi \lambda R c_{\alpha} = 0,
\end{equation}
$$
with the concentration $c_{\alpha}$ of the chemical component as the primary variable. $D$ [m$^2$/s] denotes the hydrodynamic dispersion tensor with the following relation

$$
D = (\phi D_{p} + \beta_T  \lVert \textbf{q} \rVert) \textbf{I} + （ \beta_L - \beta_T ） \frac{\textbf{q} \textbf{q}^{T}}{\lVert \textbf{q} \rVert}
$$

implemented, where $D_p$ [m$^2$/s] is the pore diffusion coefficient, $\beta_L$ and $\beta_T$ [m] are the longitudinal and transversal dispersion coefficients. $R$ [-] is the retardation factor defined as
$$
R = 1 + \rho_{b} K_{D} / \phi
$$
with the bulk density of the porous media $\rho_{b}$ [kg/m$^3$] and the distribution coefficient $K_{D}$ [m$^3$/kg], and $\lambda$ [1/s] is the first-order decay constant,
$$
\lambda = ln 2 / t_{1/2}
$$
where $t_{1/2}$ [s] is the half life of the decaying component.

It is worth noting that non-isothermal component transport process can be simulated by including the process variable `Temperature`. Currently, the non-isothermal component transport process is only available in staggered scheme.

The corresponding heat transport equation is given as follows.
$$
\begin{equation}
\left(\phi \rho_\text{f} c_\text{f} + (1-\phi) \rho_\text{s} c_\text{s} \right) \frac{\partial T}{\partial t} + \rho_\text{f} c_\text{f} \textbf{q} \cdot \nabla  T - \nabla \cdot \left (\Lambda \cdot \nabla T \right) = 0
\end{equation}
$$
where $c_\text{f}$ and $c_\text{s}$ [J/kg/K] refer to the specific heat capacity of fluid and solid, $\rho_\text{f}$ and $\rho_\text{s}$ [kg/m$^3$] denote the density of fluid and solid, $\phi$ [-] and $\Lambda$ [W/m/K] are the porosity and tensor of thermal dispersion.

## Input parameters

The following table shows an overview of all input parameters available in the ComponentTransport process.

| Parameter                  | Symbol      | Unit       | Doxygen and Example              |
| -------------------------- | ----------- | ---------- | ---------------------- |
| Porosity                   | $\phi$      |[-]       |[Link](https://doxygen.opengeosys.org/de/d8f/ogs_file_param__material__porous_medium__porous_medium__porosity.html),[Example](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/ComponentTransport/ConTracer/ConTracer_1d.prj) |
| Fluid density              | $\rho$      |[kg/m$^{3}$] |[Link](https://doxygen.opengeosys.org/d1/d47/ogs_file_param__material__fluid__density.html) |
| Intrinsic permeability     |$\textbf{k}$|[m$^{2}$]   | [Link](https://doxygen.opengeosys.org/d5/d06/ogs_file_param__material__porous_medium__permeability.html),[Example](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/ComponentTransport/ConTracer/ConTracer_1d.prj)  |
| Dynamic viscosity          | $\mu$  |[Pa$\cdot$s]|[Link](https://doxygen.opengeosys.org/da/d5d/ogs_file_param__material__fluid__viscosity.html)|
| Gravity vector (specific body force) | $\textbf{g}$|[m/s$^{2}$]  |                        | [Link](https://doxygen.opengeosys.org/db/d19/ogs_file_param__prj__processes__process__componenttransport__specific_body_force)
| Retardation factor         | $R$         |[-]         | [Example](https://doxygen.opengeosys.org/d0/d40/ogs_ctest_prj__parabolic__componenttransport__advectiondiffusionsorptiondecay__1d_advectiondiffusionsorptiondecay__prj) |
| First-order decay constant | $\lambda$   |[1/s]       | [Example](https://doxygen.opengeosys.org/d0/d40/ogs_ctest_prj__parabolic__componenttransport__advectiondiffusionsorptiondecay__1d_advectiondiffusionsorptiondecay__prj) |

## Input file definition

### Process definition

In the `ComponentTransport` process, the configuration is as follows.

* `<name>`: name of the chemical component.
* `<type>`: must be `ComponentTransport`.
* `<integration_order>`: This is the order of the integration method for element-wise integration. In common cases set to `2`.
* `<process_variables>`: The primary variables of the `ComponentTransport` process are either `<concentration>` or  `<pressure>`. For the variable concentration, the name of the chemical component is given. Like in the following example, there are 3 chemical components, i.e. Si, Al and Cl. The `<pressure>` process' variable is also named 'pressure', see `<process_variables>` section outside of process' definition.

```xml
<processes>
    <process>
        <name>hc</name>
        <type>ComponentTransport</type>
        <integration_order>2</integration_order>
        <process_variables>
            <concentration>Si</concentration>
            <concentration>Al</concentration>
            <concentration>Cl</concentration>
            <pressure>pressure</pressure>
        </process_variables>
        <specific_body_force>0 0</specific_body_force>
        <secondary_variables>
            <secondary_variable internal_name="darcy_velocity" output_name="darcy_velocity"/>
        </secondary_variables>
    </process>
</processes>
```

### Component definition

Under the keyword `<component>`, the properties of the transported chemical component are defined. The parameters here are the pore diffusion coefficient, the retardation factor, and the decay rate. Below is an example of the Si component with the corresponding transport parameters.

```xml
<component>
    <name>Si</name>
    <properties>
        <property>
            <name>pore_diffusion</name>
            <type>Constant</type>
            <value>1</value>
        </property>
        <property>
            <name>retardation_factor</name>
            <type>Constant</type>
            <value>0</value>
        </property>
        <property>
            <name>decay_rate</name>
            <type>Parameter</type>
            <parameter_name>decay</parameter_name>
        </property>
    </properties>
</component>
```

## Available benchmarks

* [Heterogeneous Saturated Mass Transport]({{< relref "hc_ogs6-vs-ogs5" >}})
* [Saturated Mass Transport]({{< relref "saturated-mass-transport" >}})
* [Saturated Variable-Density Flow and Mass Transport (Elder)](/docs/benchmarks/hydro-component/elder_jupyter)
* [Saturated Variable-Density Flow and Mass Transport (Goswami)]({{< relref "goswami" >}})
* [Theis solution for well pumping]({{< relref "hc_theis" >}})
* [Variable Dependent Boundary Condition]({{< relref "vdbc" >}})
* [Conservative tracer transport with time varying source (1D/2D)]({{< relref "contracer" >}})
* [(Advection-)diffusion-sorption-decay problem](https://www.opengeosys.org/docs/benchmarks/hydro-component/diffusionsorptiondecay/)
* [Two-layer diffusion problem](https://www.opengeosys.org/docs/benchmarks/hydro-component/multilayerdiffusion/)
