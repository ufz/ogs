+++
date = "2017-02-15T11:46:49+01:00"
title = "Groundwater Flow (Robin)"
project = "Elliptic/line_1_GroundWaterFlow/line_1e1_robin_left_picard.prj"
author = "Dmitri Naumov"
weight = 103

[menu]

  [menu.benchmarks]
    #weight = 103
    parent = "elliptic"

+++

{{< project-link >}}

## Equations

We start with simple linear homogeneous elliptic problem:
$$
\begin{equation}
k\; \Delta h = 0 \quad \text{in }\Omega
\end{equation}$$
w.r.t boundary conditions
$$
\eqalign{
h(x) = g_D(x) &\quad \text{on }\Gamma_D,\cr
k{\partial h(x) \over \partial n} = g_N(x) &\quad \text{on }\Gamma_N,
}$$
where $h$ could be hydraulic head, the subscripts $D$ and $N$ denote the Dirichlet- and Neumann-type boundary conditions, $n$ is the normal vector pointing outside of $\Omega$, and $\Gamma = \Gamma_D \cup \Gamma_N$ and $\Gamma_D \cap \Gamma_N = \emptyset$.

## Problem specification and analytical solution

We solve the Laplace equation on a line domain $[0\times 1]^2$ with $k = 1$ w.r.t. the specific boundary conditions:

$$
\eqalign{
h(x,y) = 1 &\quad \text{on } (x=0,y) \subset \Gamma_D,\cr
h(x,y) = 1 &\quad \text{on } (x,y=0) \subset \Gamma_D,\cr
k {\partial h(x,y) \over \partial n} = 1 &\quad \text{on } (x=1,y) \subset \Gamma_N,\cr
k {\partial h(x,y) \over \partial n} = 0 &\quad \text{on } (x,y=1) \subset \Gamma_N.
}$$

The solution of this problem is
$$
\begin{equation}
h(x,y) = 1 + \sum_{k=1}^\infty A_k \sin\bigg(C_k y\bigg) \sinh\bigg(C_k x\bigg),
\end{equation}
$$
