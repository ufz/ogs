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

We solve a thermo-mechanical homogeneous model in cube domain. The dimensions of this cube model are 1\,m in all directions. The boundary conditions and temperature loadings, as well as the material can refer [this PDF](../Thermo-Mechanics.pdf) for detailed problem description.

## Results and evaluation

Result showing temperature and stresses development with time in the centre node of the model:

{{< img src="../temperature".pdf >}}
{{< img src="../stress".pdf >}}

The analytical solution of stresses after heating is:
$$
\begin{equation}
\sigma_{xx} = \sigma_{yy} = sigma_{zz} = - \frac{\alpha \Delta T E}{1 - 2 \nu}\end{equation}
$$

