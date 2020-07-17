+++
date = "2017-02-15T11:46:49+01:00"
title = "Robin boundary condition"
project = "Elliptic/line_1_SteadyStateDiffusion/line_1e1_robin_left_picard.prj"
author = "Thomas Fischer, Dmitri Naumov"
weight = 104

[menu]
  [menu.benchmarks]
    parent = "elliptic"

+++

{{< data-link >}}

## Equations

We start with simple linear homogeneous elliptic problem:
$$
\begin{equation*}
k\; \Delta h = 0 \quad \text{in }\Omega
\end{equation*}$$
w.r.t boundary conditions
$$
\eqalign{
h = g_D &\quad \text{on }\Gamma_D,\cr
k{\partial h \over \partial n} = g_N &\quad \text{on }\Gamma_N,\cr
{\partial h \over \partial n} = \alpha (h_0 - h(x))  &\quad \text{on }\Gamma_R,
}$$
where $h$ could be hydraulic head, pressure, or temperature and $k$ is
the diffusion tensor (hydraulic conductivity, permeability divided by dynamic
viscosity, or heat conductivity). The subscripts $D,$ $N,$ and $R$ denote the
Dirichlet-, Neumann-, and Robin-type boundary conditions, $n$ is the normal
vector pointing outside of $\Omega$, and $\Gamma = \Gamma_D \cup \Gamma_N \cup
\Gamma_R$ and $\Gamma_D \cap \Gamma_N \cap \Gamma_R = \emptyset$.

## First benchmark: Problem specification

We solve the Laplace equation on a line domain $[0, 1]$ with $k = 1$
w.r.t. the specific boundary conditions:
$$
\eqalign{
{\partial h \over \partial n} = \alpha (h_0 - h(x)) &\quad \text{for } x=0,\cr
h(x) = g_D &\quad \text{for } x=1,
}$$
see
[`line_1e1_robin_left_picard.prj`](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/Tests/Data/Elliptic/line_1_SteadyStateDiffusion/line_1e1_robin_left_picard.prj).

### Analytical solution

One particular solution is
$$
\begin{equation*}
h(x) = A x + B.
\end{equation*}
$$

The normal direction is facing out of the bulk domain. The Robin-type boundary
condition in this example is set on the left side of the line domain.
Consequently, in this case the directional derivative is the negative derivative
$$
\begin{equation*}
\left.\frac{\partial h}{\partial n}\right\rvert_{x=0} = -h'(x)|_{x=0}.
\end{equation*}
$$
From the evaluation of the the Robin-type boundary condition it follows
$$
\begin{equation*}
\left.\frac{\partial h}{\partial n}\right\rvert_{x=0} = -A = \alpha (h_0 - h(0)) = \alpha (h_0 - B).
\end{equation*}
$$
Using the expression for $A$ in the Dirichlet-type boundary condition
$$
\begin{equation*}
h(x)|_{x=1} = A + B = -\alpha (h_0 - B) + B = -\alpha h_0 + (1+\alpha) B = g_D
\end{equation*}
$$
yields for $\alpha \not= -1$:
$$
\begin{align*}
B &= \frac{g_D + \alpha h_0}{1 + \alpha} \quad \textrm{and}\\
A &= -\alpha \left( h_0 - \frac{g_D + \alpha h_0}{1 + \alpha} \right)
= -\alpha \left( \frac{h_0 + \alpha h_0 - g_D - \alpha h_0}{1 + \alpha} \right)
= -\alpha \left( \frac{h_0 - g_D}{1 + \alpha} \right).
\end{align*}
$$
The particular solution is
$$
\begin{equation*}
h(x) = \frac{\alpha (g_D - h_0)}{1 + \alpha} x + \frac{g_D + \alpha h_0}{1 + \alpha}.
\end{equation*}
$$
Using the values from the project file $\alpha = -2,$ $h_0 = 1.5$, $g_D = 2$
results in
$$
\begin{equation*}
h(x) = \frac{-2 (2 - 1.5)}{1 + (-2)} x + \frac{2+(-2) \times 1.5}{1 + (-2)}
    = x + 1.
\end{equation*}
$$

### Results and evaluation

The left figure shows the pressure along the line, in the right figure the
difference between the analytical solution and the numerical calculated solution
is plotted.

{{< img src="../line_1e1_robin_left.png" >}}

## Second benchmark: Problem specification and analytical solution

We solve the Laplace equation on a line domain $[0, 1]$ with $k = 1$
w.r.t. the specific boundary conditions:
$$
\eqalign{
h(x) = g_D &\quad \text{for } x=0,\cr
{\partial h \over \partial n} = \alpha (h_0 - h(x)) &\quad \text{for } x=1,
}$$
see
[`line_1e1_robin_right_picard.prj`](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/Tests/Data/Elliptic/line_1_SteadyStateDiffusion/line_1e1_robin_right_picard.prj).

One particular solution is
$$
\begin{equation*}
h(x) = A x + B.
\end{equation*}
$$
Due to the Dirichlet boundary condition it follows:
$$
\begin{equation*}
h(0) = g_D = B \quad \Rightarrow \quad h(x) = A x + g_D.
\end{equation*}
$$
From the Robin-type boundary condition we get
$$
\begin{equation*}
h'(x)|_{x=1} = A = \alpha \left(h_0 - h(x)|_{x=1} \right)
    = \alpha \left.\left(h_0 - (Ax+g_D)\right)\right\rvert_{x=1}
    = \alpha (h_0 - g_D) - \alpha A.
\end{equation*}
$$
$$
\begin{equation*}
\Rightarrow A = \frac{\alpha (h_0 - g_D)}{1+\alpha}
\end{equation*}
$$
$$
\begin{equation*}
h(x) = \frac{\alpha (h_0 - g_D)}{1+\alpha} x + g_D.
\end{equation*}
$$
The values from the project file are: $\alpha = -2,$ $h_0 = 1.5$, $g_D = 1$ yielding
$$
\begin{equation*}
h(x) = x + 1.
\end{equation*}
$$
