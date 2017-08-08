+++
project = "https://github.com/ufz/ogs-data/blob/master/PhaseField/beam3d.prj"
author = "Xing-Yuan Miao"
date = "2017-05-19T09:10:33+01:00"
title = "Crack beam under tension"
weight = 158

[menu]
  [menu.benchmarks]
    parent = "phase-field"

+++

{{< project-link >}}

## Problem description

We solve a homogeneous beam model under a given displacement loading. The length of the beam is 2\,mm. Detailed model description can refer [this PDF](../Miao_Biot2017.pdf).
## Results and evaluation

Results show crack phase-field and displacement field distributions through the length of the beam:

{{< img src="../beam.png" >}}
{{< img src="../beam_d.png" >}}
{{< img src="../beam_u.png" >}}

For highlighting the deviation between the analytical and numerical solution, we provide the absolute error of the analytical solution and numerical simulation as follows:
{{< img src="../error_u.png" >}}

The analytical solution is:
$$
\begin{equation}
d (x) = 1 - {\mathrm{e}}^{\frac{- |x|}{2 \varepsilon}}
\end{equation}
$$
$$
\begin{equation}
u (x) = \dfrac{\sigma}{E} \int_0^x \dfrac{1}{d (x)^2 + k} \mathrm{d}x
\end{equation}
$$
with
$$
\begin{equation}
\sigma = \dfrac{E u_0}{I (\varepsilong, k)}
\end{equation}
$$
\begin{equation}
I (\varepsilong, k) =  \int_0^1  \dfrac{1}{d (x)^2 + k} \mathrm{d}x
\end{equation}
$$