+++
project = "https://github.com/ufz/ogs-data/blob/master/ThermoMechanics/cube_1e3.prj"
author = "Xing-Yuan Miao"
date = "2017-05-08T15:10:33+01:00"
title = "Thermoelastic Stress"
weight = 156

[menu]
  [menu.benchmarks]
    parent = "thermo-mechanics"

+++

{{< project-link >}}

## Problem description

We solve a thermo-mechanical homogeneous model in cube domain. The dimensions of this cube model are 1\,m in all directions. The boundary conditions and temperature loadings, as well as the material can refer Chapter 14 in Kolditz et al. for detailed problem description.

## Results and evaluation

Result showing temperature and stresses development with time in the centre node of the model:

{{< img src="../temperature.png" >}}
{{< img src="../stress.png" >}}

The analytical solution of stresses after heating is:
$$
\begin{equation}
\sigma_{xx} = \sigma_{yy} = \sigma_{zz} = - \frac{\alpha \Delta T E}{1 - 2 \nu} = - 3.260869\, \mathrm{MPa}
\end{equation}
$$

The relative error between the numerical simulation and the analytical solution is 9.2 $\cdot$ 10\uexp{-13}.

Kolditz, Olaf, Uwe-Jens Görke, Hua Shao, and Wenqing Wang, eds. Thermo-hydro-mechanical-chemical processes in porous media: benchmarks and examples. Vol. 86. Springer Science & Business Media, 2012.



