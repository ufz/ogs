+++
author = "Chaofan Chen, Haibing Shao"
date = "2020-02-24T13:44:00+01:00"
title = "Wellbore Heat Transport - EUBHE"
weight = 123
project = "Parabolic/T/BHE_1P/BHE_1P.prj"

[menu]
  [menu.benchmarks]
    parent = "heat-transport-bhe"

+++

{{< data-link >}}

## Problem description

Ramey (Ramey et al. (1962)) proposed the analytical solution concerning the wellbore heat transmission, which can be used to quantify the fluid temperature change in the wellbore. In order to verify the single pipe flow model in the OGS, the numerical results was compared with the [Ramey's analytical solution](Analytical_wellbore_heat_transport.zip). The detailed calculation of the Ramey's analytical solution is given below.

## Model Setup

In this benchmark, the length of the wellbore is 30 m as shown in Figure 1 and the cold water is injected into the inlet point of the wellbore with temperature of 20 $^{\circ}$C. The initial temperature of the fluid and grout in the wellbore is 20 $^{\circ}$C, and temperature of the surrounding rock is 55 $^{\circ}$C. The wellbore and pipe diameter are 0.28 m and 0.25826 m, respectively. And the flow rate is 0.0002 $m^3/s$.

{{< img src="pipe_flow_3d_model.png" width="80">}}

Figure 1: Single pipe flow model

## Ramey's analytical solution

In Ramey's analytical solution (Ramey et al. (1962)), the outlet temperature of the pipe inside the wellbore can be calculated by,

\begin{equation}
    T_o(t) = T_{s} + (T_i(t) - T_{s})\exp(-\Delta z/X)
\end{equation}

\noindent where, $q$ is the flow rate of the fluid in the wellbore and coefficient $X$ is determined by,

\begin{equation}
    X = \frac{q\rho_fc_{p,f}(\lambda_{s}+r_pUf(t))}{2\pi r_pU \lambda_{s}}
\end{equation}

with dimensionless time $t_D = \frac{\lambda_{s}t}{(\rho_{s}c_{p,s}r_b)}$, the time function $f(t)$ can be calculated by,

\begin{align}
    &f(t) = [0.4063+0.5\ln(t_D)][1+\frac{0.6}{t_D}], & t_D > 1.5
    \\
    &f(t) = 1.1281\sqrt{t_D}(1-0.3\sqrt{t_D}), & t_D \leqslant 1.5
\end{align}

\noindent and the overall heat transfer coefficient $U$ is written as follows,

\begin{equation}
    U = [\frac{r_{pi}+t_{pi}}{r_{pi}h}+(r_{pi}+t_{pi})(\frac{\ln{\frac{r_{pi}+t_{pi}}{r_{pi}}}}{\lambda_{pi}}+\frac{\ln{\frac{r_b}{r_{pi}+t_{pi}}}}{\lambda_{grout}})]^{-1}
\end{equation}

\begin{equation}
    h = \frac{\lambda_f Nu}{2r_{pi}}
\end{equation}

The Nusselt number can be determined according to the Gnielinski's equation (Gnielinski et al. (1975)),

\begin{align}
    & Nu = 4.364, & Re < 2300 \\
    & Nu = \frac{\frac{f}{8}(Re - 1000)Pr}{1+12.7\sqrt{\frac{f}{8}}(Pr^{\frac{2}{3}}-1)}, &  2300\leqslant Re < 5 \times 10^6
\end{align}

Pr is the Prandtl number, and the friction factor $f$, is evaluated by Churchill correlation (Churchill et al. (1977)),

\begin{equation}
    f = \frac{1}{(\frac{1}{[((\frac{8}{Re})^{10}+(\frac{Re}{36500})^{20})]^{1/2}}+[2.21(\ln{\frac{Re}{7}})]^{10})^{1/5}}
\end{equation}

The Prandtl and Reynolds number can be calculated as follows,

\begin{align}
    & Pr = \frac{\mu_f c_{p,f}}{\lambda_f}
    & Re = \frac{\rho_f v d_{pi}}{\mu_f}
\end{align}
\noindent where, $\mu_f, \rho_f$ and $\lambda_f$ is the fluid viscosity, density and thermal conductivity.

## Results and discussion

The outlet temperature change over time was compared against analytical solution and presented in Figure 2. After 30 days, the fluid temperature distribution in the wellbore is shown in Figure 3. The maximum relative error between the numerical model and Ramey's analytical solution is less than 0.15 \%.

In numerical model, the outlet temperature at beginning stage is affected by the initial temperature in the pipe inside the wellbore. The initial fluid temperature set in the benchmark means there is water with 20 $^{\circ}$C filled in the wellbore already before injecting water into the wellbore. But in the analytical solution, no initial temperature is set and the temperature keeps equilibrium state at every moment. The impact of initial temperature condition in numerical model is decreasing with increasement of the operational time as shown in Figure 2.

{{< img src="T_out_comparison.png" width="120">}}

Figure 2: Comparison with analytical solution results

{{< img src="absolute_error_fluid_T_30d.png" width="200">}}

Figure 3: Distributed temperature of fluid and absolute error.

## References

[1] Ramey Jr, H. J. (1962). Wellbore heat transmission. Journal of petroleum Technology, 14(04), 427-435.

[2] Gnielinski, V. (1975). New equations for heat and mass transfer in the turbulent flow in pipes and channels. NASA STI/recon technical report A, 75, 8-16.

[3] Churchill, S. W. (1977). Comprehensive correlating equations for heat, mass and momentum transfer in fully developed flow in smooth tubes. Industrial & Engineering Chemistry Fundamentals, 16(1), 109-116.
