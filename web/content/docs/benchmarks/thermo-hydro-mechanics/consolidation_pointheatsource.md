+++
project = "ThermoHydroMechanics/Linear/Point_injection/pointheatsource_quadratic-mesh.prj"

author = "Jörg Buchwald"
date = "2019-08-05T12:39:58+01:00"
title = "Consolidation around a point heat source"
weight = 70

[menu]
  [menu.benchmarks]
    parent = "thermo-hydro-mechanics"

+++

{{< data-link >}}

## Problem description

The problem describes a heat source embedded in a fluid-saturated porous medium.
The spherical symmetry is modeled using a 10 m x 10 m disc with a point heat source ($Q=150\\;\mathrm{W}$) placed at one corner ($r=0$) and a curved boundary at $r=10\\;\mathrm{m}$. Applying rotational axial symmetry at one of the linear boundaries, the model region transforms into a half-space configuration of the spherical symmetrical problem.
The initial temperature and the pore pressure are 273.15 K and 0 Pa, respectively.
The axis-normal displacements along the symmetry (inner) boundaries were set to zero, whereas the pore pressure, as well as the temperature, are set to their initial values along the outer (curved) boundary.
The heat coming from the point source is propagated through the medium, causing the fluid and the solid to expand at different rates.
The resulting pore pressure (gradient) is triggering a thermally driven consolidation process caused by the fluid flow away from the heat source until equilibrium is reached.
The corresponding derivation of the analytical solution can be found in the works cited below.
The main project input file is `square_1e2.prj`. Geometry and mesh are stored in `square_1x1.gml` and `quarter_002_2nd.vtu`.

## Equations

The problem equations can be found in the original work of Booker and Savvidou (1985) or Chaudhry et al. (2019).
The analytical solution of the coupled THM consolidation problem can be expressed in terms of some derived parameters:

\begin{equation}
    \kappa = \dfrac{K}{m}
\end{equation}
\begin{equation}\label{eq:consolidation}
    c = \dfrac{k_\text{s}}{\eta}\left(\lambda + 2G\right)
\end{equation}
\begin{equation}
    X = a_\text{u}\left(\lambda+2G\right)-b^{\prime}
\end{equation}
\begin{equation}
    Y = \dfrac{1}{\lambda+2G}\left(\dfrac{X}{\left(1-\dfrac{c}{\kappa}\right)a_\text{u}}+\dfrac{b^{\prime}}{a_\text{u}}\right)
\end{equation}
\begin{equation}
    Z = \dfrac{1}{\lambda+2G}\left(\dfrac{X}{\left(1-\dfrac{c}{\kappa}\right)a_\text{u}}\right)
\end{equation}
\begin{equation}
    r =\sqrt{x_{1}^{2}+x_{2}^{2}+x_{3}^{2}}
\end{equation}
\begin{equation}
    f^{A}=\text{erfc}\left(\dfrac{r}{2\sqrt{At}}\right),\quad A=\kappa,c
\end{equation}
\begin{equation}
    g^{A}=\dfrac{At}{r^{2}}+\left(\frac{1}{2}-\dfrac{At}{r^{2}}\right)f^{A}-\sqrt{\dfrac{At}{\pi r^{2}}} \exp\left(-\dfrac{r^{2}}{4At}\right)
\end{equation}
\begin{equation}
    g^{\ast} = Yg^{\kappa}-Zg^{c}
\end{equation}

and

\begin{equation}
    g^{A}\_{,i} = \dfrac{2x_{i}At}{r^{4}}\left(f^{A}-1+\dfrac{r}{\sqrt{\pi At}}\exp\left(-\dfrac{r^{2}}{4At}\right)\right),\quad i=1,2,3
\end{equation}
\begin{equation}
    g^{\ast}\_{,i} = Yg^{\kappa}\_{,i}-Zg^{c}_{,i}
\end{equation}

For the temperature, porepressure and displacements, the correct solution can be found in the original work:
\begin{equation}
    \Delta T = \dfrac{Q}{4\pi Kr}f^{\kappa}
\end{equation}
\begin{equation}
    p = \dfrac{X\,Q}{\left(1-\dfrac{c}{\kappa}\right)4\pi Kr}\left(f^{\kappa}-f^{c}\right)
\end{equation}
\begin{equation}
    u_{i} = \dfrac{Q a_\text{u}x_{i}}{4\pi Kr}\;g^{\ast}
\end{equation}

For the stress components the corrected expressions can be found in the work of Chaudhry et al. (2019):

\begin{equation}
    \sigma^{\prime}\_{ij\,|\,j=i} = \dfrac{Q a_\text{u}}{4\pi Kr}\left( 2G\left[g^{\ast}\left(1-\dfrac{x^{2}_{i}}{r^{2}}\right)+x_{i}g^{\ast}_{,i}\right]+\lambda \left[x_{i}g^{\ast}_{,i}+2g^{\ast}\right]\right)-b^{\prime}\Delta T
\end{equation}
\begin{equation}
    \sigma^{\prime}\_{ij\,|\,j \neq i} = \dfrac{Q a_\text{u}}{4\pi Kr}\left( G\left[x_{i}g^{\ast}_{,j}+x_{j}g^{\ast}_{,i}-2g^{\ast}\dfrac{x_{i}x_{j}}{r^{2}}\right]\right)
\end{equation}

## Results and evaluation

The analytical expressions (12-16) together with the numerical model can now be evaluated at different points as a function of time or for a given time as a function of their spatial coordinates.
The results below were taken from the benchmark published in Chaudhry et al. (2019) and might slightly differ from the benchmark in the OGS6 repo.

{{< img src="../images/resp_vs_t_square.png" >}}

In the pictures above, the analytical and numerical results for temperature ($T$), pressure ($p$), displacement ($u_i$) and stress ($\sigma_{ij}$) are plotted as function of time ($t$) at point $P=(1.3,0.682,0.0)$ and along the radial coordinate ($r$ ) at time $t=5\cdot 10^5$ (below).

{{< img src="../images/resp_vs_x_square.png" >}}

(Figures were taken from Chaudhry et al. (2019).)

The absolute errors between OGS6 and the analytical solution for temperature, pressure, and displacement are depicted below. For all three response variables, one observes that the error reaches its maximum around the same time when also the slope of the response variable is maximal.

{{< img src="../images/errorpT_vs_t.png" >}}

{{< img src="../images/errordispl_vs_t.png" >}}

## References

[1] Booker, J. R.; Savvidou, C. (1985), Consolidation around a point heat source. International Journal for Numerical and Analytical Methods in Geomechanics, 1985, 9. Jg., Nr. 2, S. 173-184.

[2] Chaudhry, A. A.; Buchwald, J.; Kolditz, O. and Nagel, T. (2019), Consolidation around a point heatsource (correction & verification). International Journal for Numerical and Analytical Methods in Geomechanics, 2019, <https://doi.org/10.1002/nag.2998>.
