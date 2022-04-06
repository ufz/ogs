+++
date = "2022-03-21T12:00:13+01:00"
title = "Component Transport Process"
author = "Haibing Shao, Renchao Lu"
weight = 44

[menu]
  [menu.userguide]
    parent = "process-dependent-configuration"
+++

## Description of the ComponentTransport process

ComponentTransport process is widely used to predict the distribution of chemical components in the subsurface, which is controlled by the groundwater flow (advection), the hydrodynamic dispersion and diffusion, the sorption on the solid phase, as well as the decay of component. 

## Mathematical framework

The governing equation implemented in OGS-6 is the so-called advective and diffusion equation (`ADE`), with the consideration of sorption and decay process. ADE is widely used to describe the concentration of chemical components in groundwater aquifer and porous media. The equations can be solved analytically in (simplified) 1D cases. For more complex geometry, especially with heterogeneous material properties, numerical solution is often preferred. The ComponentTransport process has the following assumptions. 

* The model domain is a porous media and it is fully saturated with water.
* The fluid velocity in the porous media is assumed to be slow, and the flow process is regulated by Darcy's law.
* The sorption process redistributes the chemical component between the aqueous and solid phase.
* When a decay coefficient is given, the decay of the chemical component follows first-order reaction kinetics equation. 

For the flow process, the continuity equation for flowing fluid in a saturated porous medium is as follows
$$
S_{\textrm{s}} \frac{\partial p}{\partial t} + \nabla \cdot \textbf{q} + Q_{p} = 0,
$$
with the hydraulic pressure $p$ [Pa] of the fluid as the primary variable. In the above flow equation, $S_{\textrm{s}}$ [1/Pa] is the specific storage, $t$ [s] is the time, $\textbf{q}$ [m/s] is the Darcy flux with laminar flow assumptions, and $Q_{p}$ [1/s] is the source-sink term. According to Darcy's Law, the flux $\textbf{q}$ is related to the pressure drop and body forces through
$$
\textbf{q} = -\frac{\textbf{k}}{\mu} \left(\nabla p - \rho^{\textrm{l}} \textbf{g}\right),
$$
where $\textbf{k}$ [m$^2$] is the intrinsic permeability, $\mu$ [Pa$\cdot$s] is the fluid dynamic viscosity, and $\textbf{g}$ [m/s$^2$] is the gravity vector.

For each chemical component $\alpha = 1, .., N_p$, its corresponding advective and diffusion equation (`ADE`) reads,
$$
\begin{equation}
\frac{\partial \left(\phi R c_{\alpha}\right)}{\partial t} + \nabla \cdot \left( \textbf{q} c_{\alpha} - \phi \textbf{D} \nabla c_{\alpha} \right) + \phi \alpha R c_{\alpha} = 0,
\end{equation}
$$
with the concentration $c_{\alpha}$ of the chemical component as the primary variable. $D$ [m$^2$/s] is the diffusion/dispersion coefficient for the component, $R$ [-] is the retardation factor defined as
$$
R = 1 + \rho k_{d} / \phi
$$
with the dry density of the porous media $\rho$ [kg/m$^3$] and the distribution coefficient $k_d$ [m$^3$/kg], and $\alpha$ [1/s] is the first-order decay constant, i.e. 
$$
\alpha = ln 2 / t_{1/2}
$$
where $t_{1/2}$ [s] is the half life of the decaying component. 

## Input parameters

### Process definition

In the configuration of `ComponentTransport` process, it is generally configured as follows.

* `<name>`: name of the chemical component is given here.
* `<type>`: must be ComponentTransport.
* `<integration_order>`: It is the order of the integration method for element-wise integration, normally set to 2.
* `<process_variables>`: The primary variables of the ComponentTransport process are either `<concentration>` or  `<pressure>`. For the variable concentration, the name of the chemical component is given. Like in the following example, there are 3 chemical components, i.e. Si, Al and Cl. The `<pressure>` process' variable is usually also named 'pressure', see `<process_variables>` section outside of process' definition.

```bash
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

Under the keyword `<component>`, the properties of the transported chemical component must be given. Usually the parameters to be given are the pore diffusion coefficient, the retardation factor, and the decay rate. Below is an example of the Si component with its corresponding transport parameters.

```bash
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

* [Heterogeneous Saturated Mass Transport](https://www.opengeosys.org/docs/benchmarks/hydro-component/hc_ogs6-vs-ogs5/)
* [Saturated Mass Transport](https://www.opengeosys.org/docs/benchmarks/hydro-component/hydro-component/)
* [Saturated Variable-Density Flow and Mass Transport (Elder)](https://www.opengeosys.org/docs/benchmarks/hydro-component/elder/)
* [Saturated Variable-Density Flow and Mass Transport (Goswami)](https://www.opengeosys.org/docs/benchmarks/hydro-component/goswami/)
* [Theis solution for well pumping](https://www.opengeosys.org/docs/benchmarks/hydro-component/theis/hc_theis/)
* [Variable Dependent Boundary Condition](https://www.opengeosys.org/docs/benchmarks/hydro-component/vdbc/)
* [Conservative tracer transport with time varying source (1D/2D)](https://www.opengeosys.org/docs/benchmarks/hydro-component/contracer/contracer/)
* [(Advection-)diffusion-sorption-decay problem](https://www.opengeosys.org/docs/benchmarks/notebooks/diffusionsorptiondecay/)
* [Two-layer diffusion problem](https://www.opengeosys.org/docs/benchmarks/notebooks/multilayerdiffusion/)