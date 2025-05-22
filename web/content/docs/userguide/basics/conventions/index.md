+++
date = "2021-06-17T11:00:13+01:00"
title = "OGS Conventions"
author = "Joerg Buchwald"
weight = 71
+++

## Units

In general OpenGeoSys is not using any specified intrinsic units, i.e., OGS assumes that a self-consistent set of units (
SI-units for example) is used for all quantities.
However, there are some exceptions to that rule concerning empirical laws. For instance there are some expressions for the
vapor diffusion and latent heat implemented, which assume that temperature is given in Kelvin. Therefore, we recommend just
using SI units everywhere.

For more information on how to come up with a self-consistent unit scheme see [this PDF](units_ogs.pdf).
<!-- TODO: Consider updating and extending the pdf-file on unit conventions. -->

### Dimensional analysis of source and boundary terms based on their primary variables

Most processes in OpenGeoSys are derived and implemented with a 3D spatial setting in mind.
I.e., the conceptual model of the processes exists in a three-dimensional world and boundaries are two-dimensional.
The canonical (a) Neumann boundary conditions and (b) volumetric source/sink terms (STs) for those processes
are to (a) assign a flux (aka as flux density, i.e., a total flux $Q$ of some quantity (like volume of water per time) per a 2D
pre-defined area) to parts of the boundary of the domain or to
(b) add a <em>quantity rate density</em> (i.e., a total flux of some quantity $Q$ per 3D volume) to the simulation domain,
respectively.
Here, $Q$ in our formulation is an already integrated flux. Meaning it is typically a different quantity than the primary
unknown of the `process` under consideration. It may be e.g. a secondary variable being a function of the primary unknown. This
is for instance the case, if we consider prescribed fluxes of water in a porous medium on a Neumann boundary. Here, $Q$ would
be a prescribed Darcy velocity integrated over some area being a function of the pressure gradient.
For energy balances in OpenGeoSys the primary unknown is usually the temperature $T$,
while corresponding fluxes $Q$ are usually expressed in terms of energy.

**Note**: The connection between the primary unknown and these canonical BCs and STs is made via the implemented partial
differential equation (PDE).
If a PDE is implemented slightly differently (e.g. PDE (B) is PDE (A) multiplied by a factor of density),
the corresponding quantities in the BCs and STs generally also change (by a factor of density in this example).
<!-- TODO: Please explain why some PDEs are defined when differently. Are these changes made automatically, like a multipliction with density, or is it up to user to do that? -->

Lower-dimensional realizations of the mentioned processes are essentially 2D slices through a 3D world.
Hence, the parameterization does not change.
Despite solving a problem whose solution varies only in one or two dimensions, you still have to apply Neumann BCs and
volumetric STs
as integrated quantity per time per 2D area and quantity per time per 3D volume, respectively.
Sometimes, these lower dimensional realizations are referred to as <em>cross-sectional models</em>.

In contrast, for _real_ two-dimensional and one-dimensional problems the fluxes of Neumann boundary conditions and volumetric
source terms can be seen as normalized by their one- and two-dimensional boundaries, e.g., quantity per time per length<sup>d</sup>, where d stands for dimensionality, in the case of volumetric STs.

**Attention**: The specific units given below apply to 3D and mentioned cross-sectional models.
That applies to all of OpenGeoSys's processes; currently there are no processes in OpenGeoSys that are manifest 1D or 2D.

**Units:**

| Quantity name | Dimension symbol  | SI base units |
| ------------- | ----------------- | ------------- |
| time          | $\mathrm{T}$      | $\mathrm{s}$  |
| length        | $\mathrm{L}$      | $\mathrm{m}$  |
| mass          | $\mathrm{M}$      | $\mathrm{kg}$ |
| temperature   | $\mathrm{\Theta}$ | $\mathrm{K}$  |

**Temperature:**

* [Dirichlet BC](https://doxygen.opengeosys.org/d5/d71/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__dirichlet): Boundaries held at a fixed temperature (Units: $\mathrm{\Theta}$, SI: $\mathrm{K}$).
* [Neumann BC](https://doxygen.opengeosys.org/d1/d2e/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__neumann): Used to impose a heat flux through a boundary domain (Units: $\mathrm{M \cdot T^{-3}}$, SI: $\mathrm{W}\cdot \mathrm{m}^{-2}$).
* [Volumetric ST](https://doxygen.opengeosys.org/d0/d89/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__volumetric): Representing a volumetric heat source in the model domain (Units: $\mathrm{M \cdot L^{-1} \cdot T^{-3}}$, SI: $\mathrm{W} \cdot \mathrm{m}^{-3}$).
* [Nodal ST](https://doxygen.opengeosys.org/d0/d2c/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__nodal): Representing a heat source in the model domain (Units: $\mathrm{M \cdot L^{d-1} \cdot T^{-3}}$, SI: $\mathrm{W} \cdot \mathrm{m}^{d-3}$).
* [Line ST](https://doxygen.opengeosys.org/d9/d4a/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__line): Representing a heat source on a line shaped subdomain (Units: $\mathrm{M} \cdot \mathrm{L}^{d-2} \cdot \mathrm{T}^{-3}$, SI: $\mathrm{W} \cdot \mathrm{m}^{d-4}$).

**Displacement:**

* [Dirichlet BC](https://doxygen.opengeosys.org/d5/d71/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__dirichlet): Boundaries held at a fixed displacement (Units: $\mathrm{L}$, SI: $\mathrm{m}$).
* [Neumann BC](https://doxygen.opengeosys.org/d1/d2e/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__neumann): Used to apply an external traction at a boundary (Units: $\mathrm{M \cdot L^{-1} \cdot T^{-2}}$, SI: $\mathrm{N} \cdot \mathrm{m}^{-2}$).
* [Volumetric ST](https://doxygen.opengeosys.org/d0/d89/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__volumetric): A body force density acting in the model domain. One example is gravity, but gravity can be modeled by other means in OpenGeoSys (Units: $\mathrm{M \cdot L^{-2} \cdot T^{-2}}$, SI: $\mathrm{N} \cdot \mathrm{m}^{-3}$).

**Pressure:**

* [Dirichlet BC](https://doxygen.opengeosys.org/d5/d71/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__dirichlet): Boundaries held at a fixed pressure (Units: $\mathrm{M \cdot L^{-1} \cdot T^{-2}}$, SI: $\mathrm{Pa}$).
* [Neumann BC](https://doxygen.opengeosys.org/d1/d2e/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__neumann): Used to impose a mass flux through a domain boundary (Units: $\mathrm{M \cdot L^{-2} \cdot T^{-1}}$, SI: $\mathrm{kg} \cdot \mathrm{m}^{-2} \cdot \mathrm{s}^{-1}$).
* [Volumetric ST](https://doxygen.opengeosys.org/d0/d89/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__volumetric): Representing a volumetric mass source in the model domain (Units: $\mathrm{ M \cdot L^{-3} \cdot T^{-1}}$, SI: $\mathrm{kg} \cdot \mathrm{m}^{-3} \cdot \mathrm{s}^{-1}$).
* [Nodal ST](https://doxygen.opengeosys.org/d0/d2c/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__nodal): Representing a mass production rate in the model domain (Units: $\mathrm{M} \cdot \mathrm{L}^{d-3} \cdot \mathrm{T}^{-1}$, SI: $\mathrm{kg} \cdot \mathrm{m}^{d-3} \cdot \mathrm{s}^{-1}$).
* [Line ST](https://doxygen.opengeosys.org/d9/d4a/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__line): Representing a mass source on a line shaped subdomain (Units: $\mathrm{M} \cdot \mathrm{L}^{d-4} \cdot \mathrm{T}^{-1}$, SI:  $\mathrm{kg} \cdot \mathrm{m}^{d-4} \cdot \mathrm{s}^{-1}$).

**Pressure (Liquid Flow Process):**

The liquid flow process uses either a volume-based or a mass-based PDE, which
 can be selected by setting `equation_balance_type` to `volume` or `mass`
 correspondingly in the project file. For the mass-based PDE, the
 **pressure**-section above applies to the liquid flow process as well. For the
 volume-based formulation, you have to keep in mind that the density must be
 constant, and have to divide all units of natural boundary conditions and
 source terms by density. This leads to the following list:

* [Dirichlet BC](https://doxygen.opengeosys.org/d5/d71/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__dirichlet): Boundaries held at a fixed pressure (Units: $\mathrm{M \cdot L^{-1} \cdot T^{-2}}$, SI: $\mathrm{Pa}$).
* [Neumann BC](https://doxygen.opengeosys.org/d1/d2e/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__neumann): Used to impose a volume flux through a domain boundary (Units: $\mathrm{L} \cdot \mathrm{T}^{-1}$, SI: $\mathrm{m} \cdot \mathrm{s}^{-1} = \mathrm{m}^3 \cdot \mathrm{m}^{-2} \cdot \mathrm{s}^{-1}$).
* [Volumetric ST](https://doxygen.opengeosys.org/d0/d89/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__volumetric): Representing a volumetric volume source in the model domain (Units: $\mathrm{T}^{-1}$, SI: $\mathrm{s}^{-1} = \mathrm{m}^{3} \cdot \mathrm{m}^{-3} \cdot \mathrm{s}^{-1}$).
* [Nodal ST](https://doxygen.opengeosys.org/d0/d2c/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__nodal): Representing a volume production rate in the model domain (Units: $\mathrm{L}^{d} \cdot \mathrm{T}^{-1}$, SI: $\mathrm{m}^{d} \cdot \mathrm{s}^{-1}$).
* [Line ST](https://doxygen.opengeosys.org/d9/d4a/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__line): Representing a volume source on a one-dimensional subdomain (Units: $\mathrm{L}^{d-1} \cdot \mathrm{T}^{-1}$, SI:  $\mathrm{m}^{d-1} \cdot \mathrm{s}^{-1}$).

## Stress sign

OpenGeoSys uses the positive sign convention for tensile stresses, whereas compressional stresses carry a negative sign.

## Stress-strain assumption for 2D simulations

OpenGeoSys assumes plane strain conditions as default for 2D simulations.

## Order of process variables/global components

OpenGeoSys uses internally a vector of primary vector components ("global component vector").
Although, which set of components is used, varies from process to process, the order of primary variable components for most
THM-based processes is: $T$, $p$, $u_{x}$, $u_{y}$, $u_{z}$.
However, there are some exceptions and as there are also different primary variables for each process, you will find a list
containing the process variables as they appear in the global component vector.
This order is used, e.g., to display the per component convergence of the non-linear solver in the output.

<!-- TODO: The considered processes beneath need a more detailed explanation on which physical processes are exactly considered how. -->

* [ComponentTransport](https://doxygen.opengeosys.org/de/d0d/namespaceProcessLib_1_1ComponentTransport.html#processvariablescomponenttransport)
* [HeatConduction](https://doxygen.opengeosys.org/df/d77/namespaceProcessLib_1_1HeatConduction.html#processvariablesheatconduction)
* [HeatTransportBHE](https://doxygen.opengeosys.org/d9/d4d/namespaceProcessLib_1_1HeatTransportBHE.html#processvariablesbhe)
* [HydroThermal](https://doxygen.opengeosys.org/dd/d60/namespaceProcessLib_1_1HT.html#processvariablesht)
* [HydroMechanics](https://doxygen.opengeosys.org/d3/da3/namespaceProcessLib_1_1HydroMechanics.html#processvariableshm)
* [LIE HM](https://doxygen.opengeosys.org/de/d36/namespaceProcessLib_1_1LIE_1_1HydroMechanics.html#processvariablesliehm)
* [LIE SmallDeformation](https://doxygen.opengeosys.org/de/df3/namespaceProcessLib_1_1LIE_1_1SmallDeformation.html#processvariablesliesd)
* [LiquidFlow](https://doxygen.opengeosys.org/dc/dd7/namespaceProcessLib_1_1LiquidFlow.html#processvariableslf)
* [PhaseField](https://doxygen.opengeosys.org/d2/dc2/namespaceProcessLib_1_1PhaseField.html#processvariablespf)
* [RichardsComponentTransport](https://doxygen.opengeosys.org/d6/d3b/namespaceProcessLib_1_1RichardsComponentTransport.html#processvariablesrct)
* [RichardsMechanics](https://doxygen.opengeosys.org/d6/d4a/namespaceProcessLib_1_1RichardsMechanics.html#processvariablesrm)
* [SmallDeformation](https://doxygen.opengeosys.org/da/d84/namespaceProcessLib_1_1SmallDeformation.html#processvariablessd)
* [SteadyStateDiffusion](https://doxygen.opengeosys.org/d8/d59/namespaceProcessLib_1_1SteadyStateDiffusion.html#processvariablesssd)
* [TH2M](https://doxygen.opengeosys.org/d6/de7/namespaceProcessLib_1_1TH2M.html#processvariablesth2m)
* [ThermalTwoPhaseFlowWithPP](https://doxygen.opengeosys.org/dd/d48/namespaceProcessLib_1_1ThermalTwoPhaseFlowWithPP.html#processvariables)
* [ThermoMechanics](https://doxygen.opengeosys.org/d5/dd4/namespaceProcessLib_1_1ThermoMechanics.html#processvariablestm)
* [ThermoMechanicalPhaseField](https://doxygen.opengeosys.org/de/d90/namespaceProcessLib_1_1ThermoMechanicalPhaseField.html#processvariablestmpf)
* [ThermoHydroMechanics](https://doxygen.opengeosys.org/db/d5f/namespaceProcessLib_1_1ThermoHydroMechanics.html#processvariablesthm)
* [ThermoRichardsFlow](https://doxygen.opengeosys.org/d8/dd2/namespaceProcessLib_1_1ThermoRichardsFlow.html#processvariablestrf)
* [ThermoRichardsMechanics](https://doxygen.opengeosys.org/d9/de9/namespaceProcessLib_1_1ThermoRichardsMechanics.html#processvariablestrm)
* [TwoPhaseFlow with PP](https://doxygen.opengeosys.org/d0/d3f/namespaceProcessLib_1_1TwoPhaseFlowWithPP.html#processvariablestpfwpp)
* [TwoPhaseFlow with Prho](https://doxygen.opengeosys.org/d4/de4/namespaceProcessLib_1_1TwoPhaseFlowWithPrho.html#processvariablestpfwprho)

## <a name="symmetric-tensors"></a>  Symmetric tensors and Kelvin mapping

To map the elasticity/stiffness tensor OpenGeoSys internally uses a Kelvin mapping with an adapted component ordering for
computational reasons \[[1](https://arxiv.org/abs/1605.09606)\].
For 2D, the Kelvin-Vector of the stress tensor looks like $\sigma=(\sigma_{xx},\sigma_{yy},\sigma_{zz},\sqrt{2}\sigma_{xy})$
whereas the 3D version reads as $\sigma=(\sigma_{xx},\sigma_{yy},\sigma_{zz},\sqrt{2}\sigma_{xy}, \sqrt{2}\sigma_{yz},\sqrt{2}\sigma_{xz})$.

For Kelvin mapping also see the [conversion function documentation](https://doxygen.opengeosys.org/d6/dce/namespacemathlib_1_1kelvinvector#ad78b122c10e91732e95181b6c9a92486).

The input and output of symmetric tensors consists of the full (symmetric) tensor elements (without the factor $\sqrt{2}$),
retaining the same order.
I.e., Input and output components are $\sigma=(\sigma_{xx},\sigma_{yy},\sigma_{zz},\sigma_{xy})$
and $\sigma=(\sigma_{xx},\sigma_{yy},\sigma_{zz},\sigma_{xy},\sigma_{yz},\sigma_{xz})$, respectively.

## Staggered Scheme

A staggered scheme solves coupled problems by alternating on separate physical domains (e.g. thermal and mechanical) in
contrast to monolithic schemes which solve all domains simultaneously (e.g., thermomechanical).
Thus staggered schemes add another level of iterations, however this may pay off since the subproblems are smaller and they
enable a finer tuning of the specific solvers.

### Fixed-stress Split for Hydro-mechanical Processes

For hydro-mechanical processes the fixed-stress split has been implemented, since it turned out advantageous [[1]](#1).
For the sake of brevity, we do not describe the scheme itself, but intend to provide guidance for its stabilization parameter.
On this parameter depends how many coupling iterations are needed and thus how long it takes to obtain a solution.
The optimal value of this parameter is not a-priori known, only that it lies within a certain interval is known (see [[2]](#2)):

$$
\frac{1}{2}\frac{\alpha^2}{K_\mathrm{1D}} \le \beta_\mathrm{FS} \le \frac{\alpha^2}{K_\mathrm{ph}},
$$

where $\alpha$ denotes Biot's coefficient, and $K_\mathrm{1D}$ and $K_\mathrm{ph}$ are the bulk moduli described next.
By $K_\mathrm{ph}$ we mean the bulk modulus adjusted to the spatial dimension, so in three dimensions it coincides with the
drained bulk modulus, whereas in lower dimensions it becomes a constrained bulk modulus (2D plane strain, 1D uniaxial strain).
Assuming isotropic, linear elasticity we have

$$\begin{eqnarray}
K_\mathrm{3D} &=& \lambda + \frac{2}{3}\mu, \\
K_\mathrm{2D} &=& \lambda + \frac{2}{2}\mu, \\
K_\mathrm{1D} &=& \lambda + \frac{2}{1}\mu. \\
\end{eqnarray}$$

OGS sets the stabilization parameter, which corresponds to a coupling compressibility,
 via an optional tag `fixed_stress_stabilization_parameter` inside the tag `coupling_scheme`.

$$
\beta_\mathrm{FS} = p_\mathrm{FS} \frac{\alpha^2}{K_\mathrm{3D}},
$$

by default to $p_\mathrm{FS}=\frac{1}{2}$.
For isotropic, linear elasticity we provide the interval [[2]](#2) and the recommended value [[3]](#3) in dependence on
Poisson's ratio $\nu$ (note $\frac{\lambda}{\mu}=\frac{2\nu}{1-2\nu}$).

| | 2D  | 3D |
| ------ | ------ | ------ |
| $p_\mathrm{FS}^\mathrm{min}$          | $\frac{1}{6}\frac{1+\nu}{1-\nu}$      | same as 2D  |
| $p_\mathrm{FS}^\mathrm{MW}$        | $\frac{1+\nu}{3}$      | $\frac{1}{2}$  |
| $p_\mathrm{FS}^\mathrm{max}$          | $\frac{2(1+\nu)}{3}$      | $1$ |

For more information about the algorithms of the fixed stress splitting, please
 visit the page about a HM benchmark:
   [Injection and Production in 1D Linear Poroelastic
  Medium with the Staggered Scheme]({{< relref "InjectionProduction#staggered-scheme-fixed-stress-splitting" >}}).

## References

<a id="1">[1]</a>
{{< bib "kimtchjua2009" >}}

<a id="2">[2]</a>
{{< bib "stonor2019" >}}

<a id="3">[3]</a>
{{< bib "mikwhe2013" >}}
