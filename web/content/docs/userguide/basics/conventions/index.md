+++
date = "2021-06-17T11:00:13+01:00"
title = "OGS Conventions"
author = "Joerg Buchwald"
weight = 6


[menu]
  [menu.userguide]
    parent = "basics"
+++

## Units

In general OpenGeoSys is not using any intrinsic units, i.e. OGS assumes that a self-consistent set of units (SI for example) is used for all quantities.
However, there are some exceptions to the rule as for instance some empirical laws like some expressions for the vapor diffusion and latent heat that assume the temperature to be given in Kelvin.
Therefore, we recommend using SI base units.


### Dimensional analysis of source and boundary terms based on their primary variables

__Units:__

| Quantity name | Dimension symbol  | SI base units |
| ------------- | ----------------- | ------------- |
| time          | $\mathrm{T}$      | $\mathrm{s}$  |
| length        | $\mathrm{L}$      | $\mathrm{m}$  |
| mass          | $\mathrm{M}$      | $\mathrm{kg}$ |
| temperature   | $\mathrm{\Theta}$ | $\mathrm{K}$  |

__Temperature:__
* [Dirichlet BC](https://doxygen.opengeosys.org/d5/d71/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__dirichlet): Boundaries held at a fixed temperature (Units: $\mathrm{\Theta}$, SI: $\mathrm{K}$).
* [Neumann BC](https://doxygen.opengeosys.org/d1/d2e/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__neumann): Used to impose a heat flux through a boundary domain (Units: $\mathrm{M \cdot L^{3-d} \cdot T^{-3}}$, SI: $\mathrm{W}\cdot \mathrm{m}^{d-1}$).
* [Nodal ST](https://doxygen.opengeosys.org/d0/d2c/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__nodal): Representing a heat source in the model domain (Units: $\mathrm{M \cdot L^{2} \cdot T^{-3}}$, SI: $\mathrm{W}$).
* [Volumetric ST](https://doxygen.opengeosys.org/d0/d89/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__volumetric): Representing a volumetric heat source in the model domain (Units: $\mathrm{M \cdot L^{2-d} \cdot T^{-3}}$, SI: $\mathrm{W} \cdot \mathrm{m}^{-d}$).
* [Line ST](https://doxygen.opengeosys.org/d9/d4a/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__line): Representing a heat source on a line shaped subdomain (Units: $\mathrm{M \cdot L^{1} \cdot T^{-3}}$, SI: $\mathrm{W} \cdot \mathrm{m}^{-1}$).

__Displacement__
* [Dirichlet BC](https://doxygen.opengeosys.org/d5/d71/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__dirichlet): Boundaries held at a fixed displacement (Units: $\mathrm{L}$, SI: $\mathrm{m}$).
* [Neumann BC](https://doxygen.opengeosys.org/d1/d2e/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__neumann): Used to apply an external traction at a boundary (Units: $\mathrm{M \cdot L^{2-d} \cdot T^{-2}}$, SI: $\mathrm{N} \cdot \mathrm{m}^{d-1}$).

__Pressure__
* [Dirichlet BC](https://doxygen.opengeosys.org/d5/d71/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__dirichlet): Boundaries held at a fixed pressure (Units: $\mathrm{M \cdot L^{-1} \cdot T^{-2}}$, SI: $\mathrm{Pa}$).
* [Neumann BC](https://doxygen.opengeosys.org/d1/d2e/ogs_file_param__prj__process_variables__process_variable__boundary_conditions__boundary_condition__neumann): Used to impose a mass flux through a boundary domain (Units: $\mathrm{M \cdot L^{1-d} \cdot T^{-1}}$, SI: $\mathrm{kg} \cdot \mathrm{m}^{1-d} \cdot \mathrm{s}^{-1}$).
* [Nodal ST](https://doxygen.opengeosys.org/d0/d2c/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__nodal): Representing a mass production rate in the model domain (Units: $\mathrm{M \cdot T^{-1}}$, SI: $\mathrm{kg} \cdot \mathrm{s}^{-1}$).
* [Volumetric ST](https://doxygen.opengeosys.org/d0/d89/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__volumetric): Representing a volumetric mass source in the model domain (Units: $\mathrm{ M \cdot L^{-d} \cdot T^{-1}}$, SI: $\mathrm{kg} \cdot \mathrm{m}^{-d} \cdot \mathrm{s}^{-1}$).
* [Line ST](https://doxygen.opengeosys.org/d9/d4a/ogs_file_param__prj__process_variables__process_variable__source_terms__source_term__line): Representing a mass source on a line shaped subdomain (Units: $\mathrm{ M \cdot L^{-1} \cdot T^{-1}}$, SI:  $\mathrm{kg} \cdot \mathrm{m}^{-1} \cdot \mathrm{s}^{-1}$).


## Stress sign

OpenGeoSys uses the positive sign convention for tensile stresses, whereas compressional stresses carry a negative sign.

## Stress-strain assumption for 2D simulations

OpenGeoSys assumes plane strain conditions as default for 2D simulations.

## Order of process variables/global components

OpenGeoSys uses internally a vector of primary vector components ("global component vector").
Although, which set of components is used, varies from process to process, the order of primary variable components for most THM-based processes is: $T$, $p$, $u_{x}$, $u_{y}$, $u_{z}$.
However, there are some exceptions and as there are also different primary variables for each process, you will find a list containing the process variables as they appear in the global component vector.
This order is used, e.g., to display the per component convergence of the non-linear solver in the output.

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
* [SmallDeformation Nonlocal](https://doxygen.opengeosys.org/d9/d9a/namespaceProcessLib_1_1SmallDeformationNonlocal.html#processvariablessdnl)
* [SteadyStateDiffusion](https://doxygen.opengeosys.org/d8/d59/namespaceProcessLib_1_1SteadyStateDiffusion.html#processvariablesssd)
* [StokesFlow](https://doxygen.opengeosys.org/d7/d18/namespaceProcessLib_1_1StokesFlow.html#processvariablessf)
* [TES](https://doxygen.opengeosys.org/d5/d3a/namespaceProcessLib_1_1TES.html#processvariablestes)
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

To map the elasticity/stiffness tensor OpenGeoSys internally uses a Kelvin mapping with an adapted component ordering for computational reasons \[[1](https://arxiv.org/abs/1605.09606)\].
For 2D, the Kelvin-Vector of the stress tensor looks like $\sigma=(\sigma_{xx},\sigma_{yy},\sigma_{zz},\sqrt{2}\sigma_{xy})$ whereas the 3D version reads as $\sigma=(\sigma_{xx},\sigma_{yy},\sigma_{zz},\sqrt{2}\sigma_{xy}, \sqrt{2}\sigma_{yz},\sqrt{2}\sigma_{xz})$.

For Kelvin mapping also see the [conversion function documentation](https://doxygen.opengeosys.org/d6/dce/namespacemathlib_1_1kelvinvector#ad78b122c10e91732e95181b6c9a92486).

The input and output of symmetric tensors consists of the full (symmetric) tensor elements (without the factor $\sqrt{2}$), retaining the same order.
I.e., Input and output components are be $\sigma=(\sigma_{xx},\sigma_{yy},\sigma_{zz},\sigma_{xy})$
and $\sigma=(\sigma_{xx},\sigma_{yy},\sigma_{zz},\sigma_{xy},\sigma_{yz},\sigma_{xz})$, respectively.

## Staggered Scheme

A staggered scheme solves coupled problems by alternating on separate physical domains (e.g. thermal and mechanical) in contrast to monolithic schemes which solve all domains simultaneously (e.g. thermomechanical).
Thus staggered schemes add another level of iterations, however this may pay off since the subproblems are smaller and they enable a finer tuning of the specific solvers.

### Fixed-stress Split for Hydro-mechanical Processes

For hydro-mechanical processes the fixed-stress split has been implemented, since it turned out advantageous [[1]](#1).
For sake of brevity, we do not describe the scheme itself, but intend to provide guidance for its stabilization parameter.
On this parameter depends how many coupling iterations are needed and thus how long it takes to obtain a solution.
The optimal value of this parameter is not a-priori known, only that it lies within a certain interval [[2]](#2)
\begin{equation}
\frac{1}{2}\frac{\alpha^2}{K_\mathrm{1D}} \le \beta_\mathrm{FS} \le \frac{\alpha^2}{K_\mathrm{ph}},
\end{equation}
where $\alpha$ denotes Biot's coefficient, and $K_\mathrm{1D}$ and $K_\mathrm{ph}$ are the bulk moduli described next.
By $K_\mathrm{ph}$ we mean the bulk modulus adjusted to the spatial dimension, so in three dimensions it coincides with the drained bulk modulus, whereas in lower dimensions it becomes a constrained bulk modulus (2D plane strain, 1D uniaxial strain).
Assuming isotropic, linear elasticity we have
\begin{eqnarray}
K_\mathrm{3D} &=& \lambda + \frac{2}{3}\mu, \\
K_\mathrm{2D} &=& \lambda + \frac{2}{2}\mu, \\
K_\mathrm{1D} &=& \lambda + \frac{2}{1}\mu. \\
\end{eqnarray}
OGS sets the stabilization parameter, which corresponds to a coupling compressibility, via the `coupling_scheme_parameter`
\begin{equation}
\beta_\mathrm{FS} = p_\mathrm{FS} \frac{\alpha^2}{K_\mathrm{3D}},
\end{equation}
by default to $p_\mathrm{FS}=\frac{1}{2}$.
For isotropic, linear elasticity we provide the interval [[2]](#2) and the recommended value [[3]](#3) in dependence on Poisson's ratio $\nu$ (note $\frac{\lambda}{\mu}=\frac{2\nu}{1-2\nu}$).

| | 2D  | 3D |
| ------ | ------ | ------ |
| $p_\mathrm{FS}^\mathrm{min}$          | $\frac{1}{6}\frac{1+\nu}{1-\nu}$      | same as 2D  |
| $p_\mathrm{FS}^\mathrm{MW}$        | $\frac{1+\nu}{3}$      | $\frac{1}{2}$  |
| $p_\mathrm{FS}^\mathrm{max}$          | $\frac{2(1+\nu)}{3}$      | $1$ |

## References
<a id="1">[1]</a>
{{< bib "kimtchjua2009" >}}

<a id="2">[2]</a>
{{< bib "stonor2019" >}}

<a id="3">[3]</a>
{{< bib "mikwhe2013" >}}
