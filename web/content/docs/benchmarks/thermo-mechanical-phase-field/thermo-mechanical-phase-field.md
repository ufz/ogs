+++
project = "ThermoMechanicalPhaseField/cube_1e0.prj"
author = "Xing-Yuan Miao"
date = "2018-03-11T09:10:33+01:00"
title = "Non-isothermal test"
weight = 158

[menu]
  [menu.benchmarks]
    parent = "thermo-mechanical-phase-field"

+++

{{< data-link >}}

## Problem description

We solve a homogeneous damage under non-isothermal conditions. The size of the cubic model is 1$\times$1$\times1$ mm. Detailed model description can refer the latest published benchmark book "Thermo-Hydro-Mechanical-Chemical Processes in Fractured Porous Media: Modelling and Benchmarking" (Chapter 9.7 -- A Phase-Field Model for Brittle Fracturing of Thermo-Elastic Solids).

A unconfined compression test was applied as a comparison.
The thermal expansion test was implemented by imposing a temperature increase to the domain while the top surface of the model was held in place. The temperature loading was chosen to achieve the same compressive load as that imposed in the unconfined compression test.

## Models, results, and evaluation

Results show phase-field evolution in the thermo-mechanical case can follow the mechanical case, and both solutions correspond to the analytical solution:

{{< img src="../uncon_com_bc.png" >}}
{{< img src="../therm_exp_bc.png" >}}
{{< img src="../t_pf.png" >}}

The analytical solution is: $$d = \dfrac{G\textrm{c}}{G\textrm{c}+4\epsilon \psi_\textrm{e}^+}$$
where due to negative (elastic) volume strains only the deviatoric energy drives the phase field. 

$$
\begin{equation}
\psi_\textrm{e}^+ = \mu \mathbf{\epsilon}^\textrm{D} : \mathbf{\epsilon}^\textrm{D} = \frac{2\mu}{3} \left(\frac{u(1+\nu)}{L} \right)^2
\end{equation}
$$

for mechanical case, and

$$
\begin{equation}
\psi_\textrm{e}^+ = \mu \mathbf{\epsilon}^\textrm{D} : \mathbf{\epsilon}^\textrm{D} = \frac{2\mu}{3} [\alpha \Delta T(1+\nu)]^2
\end{equation}
$$

for thermo-mechanical case.
Where the Poisson's ratio is evolving with the degradation of the shear modulus 

$$
\begin{equation}
\nu(d) = \frac{3K-2Gd^2}{2(3K+Gd^2)}
\end{equation}
$$
