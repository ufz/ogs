+++
project = ["ThermoMechanics/cube_1e3.prj"]
author = "Xing-Yuan Miao"
date = "2017-05-08T15:10:33+01:00"
title = "Thermoelastic Stress"
image = "stress.png"
+++

{{< data-link >}}

## Problem description

We solve a thermo-mechanical homogeneous model in cube domain. The dimensions of
this cube model are 1 m in all directions. The boundary conditions and
temperature loadings, as well as the material can refer Chapter 14 in Kolditz et
al. \cite Kolditz2012 for detailed problem description.

## Results and evaluation

Result showing temperature and stresses development with time in the centre node
of the model:

{{< figure src="temperature.png" >}}
{{< figure src="stress.png" >}}

The analytical solution of stresses after heating is:
$$\begin{equation}
\sigma_{xx} = \sigma_{yy} = \sigma_{zz} = - \frac{\alpha \Delta T E}{1 - 2 \nu}
= - 3.260869\, \textrm{MPa}
\end{equation}$$

The relative error between the numerical simulation and the analytical solution
is 9.2<span class="math inline">⋅10<sup>-13</sup></span>.

## References

{{< bib "kolditz2012" >}}
