+++
author = "Thomas Fisher, Dmitri Naumov, Fabien Magri, Marc Walther, Tianyuan Zheng, Olaf Kolditz"
date = "2025-10-13"
title = "Hydro-Thermal Process"
weight = 3
+++

The content of this webpage was taken from [this document](/docs/benchmarks/hydro-thermal/constant-viscosity/HT-Process.pdf).

## $1 H T$ process

The hydro-thermal ( $H T$ ) process in porous media consists of the coupling of groundwater and thermal flow processes.
Both processes are described using partial differential equations (PDE) of parabolic type.

Sec. 1.1 is a general motivation to parabolic PDEs which is very similar to 1 .
Further reading: (4) [5] [2]

### 1.1 Balance Equations

Let $\Omega$ be a domain, $\Gamma$ the boundary of the domain and let $u$ be an intrinsic quantity (for instance mass or heat) and the volume density is described by a function $S(u)$.
The amount of the quantity in the domain can vary within time by two reasons.
Firstly, new quantity can accumulate by flow over $\Gamma$ or secondly it can be generated due to the presence of sources or sinks within $\Omega$.
Consequently, the balance reads

$$
\begin{equation*}
\frac{\partial}{\partial t} \int_{\Omega} S(u(x, t)) \mathrm{d} x=-\int_{\Gamma}\langle J(x, t) \mid n(x)\rangle \mathrm{d} \sigma+\int_{\Omega} Q(x, t) \mathrm{d} x \tag{1.1}
\end{equation*}
$$

where $J(x, t)$ is the flow over the boundary, $n$ is normal vector pointing outside of $\Omega, \mathrm{d} \sigma$ is an infinitesimal small surface element and $Q(x, t)$ describes sources and sinks within $\Omega$.
Further mathematical manipulations leads to

$$
\begin{equation*}
\int_{\Omega} \frac{\partial S(u(x, t))}{\partial t} \mathrm{~d} x+\int_{\Gamma}\langle J(x, t) \mid n(x)\rangle \mathrm{d} \sigma-\int_{\Omega} Q(x, t) \mathrm{d} x=0 \tag{1.2}
\end{equation*}
$$

Applying the theorem of Gauss yields to

$$
\begin{equation*}
\int_{\Omega} \frac{\partial S(u(x, t))}{\partial t} \mathrm{~d} x+\int_{\Omega} \operatorname{div} J(x, t) \mathrm{d} x-\int_{\Omega} Q(x, t) \mathrm{d} x=0 \tag{1.3}
\end{equation*}
$$

Finally,

$$
\begin{equation*}
\int_{\Omega}\left[\frac{\partial S(u(x, t))}{\partial t}+\operatorname{div} J(x, t)-Q(x, t)\right] \mathrm{d} x=0 \tag{1.4}
\end{equation*}
$$

Since the domain is arbitrary it holds:

$$
\begin{equation*}
\frac{\partial S(u(x, t))}{\partial t}+\operatorname{div} J(x, t)-Q(x, t)=0 \tag{1.5}
\end{equation*}
$$

Depending on the constitutive law that describes the flow $J$, we obtain the balance equation of the considered process.
Important practical laws are

$$
\begin{equation*}
J^{(1)}=-\boldsymbol{K} \boldsymbol{g} \mathbf{r a d} u=-\boldsymbol{K} \nabla u \tag{1.6}
\end{equation*}
$$

which describes diffusive flow and

$$
\begin{equation*}
J^{(2)}=c u \quad \text { (where } c \text { is a velocity vector) } \tag{1.7}
\end{equation*}
$$

which describes advective flow or a combination of (1.6) and (1.7).
For instance, substituting (1.6) in (1.5) leads to the following parabolic partial differential equation:

$$
\begin{equation*}
\frac{\partial S(u(x, t))}{\partial t}-\nabla \cdot[\boldsymbol{K}(x, t) \nabla u(x, t)]-Q(x, t)=0 \tag{1.8}
\end{equation*}
$$

while the description of the flow by a combination of (1.6) and (1.7) yields to

$$
\begin{equation*}
\frac{\partial S(u(x, t))}{\partial t}-\nabla \cdot[\boldsymbol{K}(x, t) \nabla u(x, t)-c u(x, t)]-Q(x, t)=0 . \tag{1.9}
\end{equation*}
$$

### 1.2 Groundwater Flow

### 1.2.1 Constitutive Law: Darcy Flow

In the modeling of groundwater flow we assume the validity of the Darcy law:

$$
\begin{equation*}
q=-\frac{\boldsymbol{\kappa}}{\mu(T)}\left(\boldsymbol{\operatorname { g r a d }} p+\varrho_{f}(T) \cdot g \cdot e_{z}\right)=-\frac{\boldsymbol{\kappa}}{\mu(T)}(\boldsymbol{\operatorname { g r a d }}(\underbrace{p+\varrho_{f}(T) \cdot g \cdot z}_{\Psi}))=-\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi \tag{1.10}
\end{equation*}
$$

where

- $T$ is the temperature in $[\Theta]$
- $q(x, t, T)$ is the Darcy velocity in $\left[\frac{L}{T}\right]$
- $p(x, t)$ is the pressure $\left[\frac{M L}{T^{2}} \frac{1}{L^{2}}\right]$,
- $\kappa$ is the anisotropic intrinsic permeability tensor of the porous medium (that can depend on the saturation $S\left[L^{2}\right]$ which will lead to Richards flow)
- $\mu(T, p)$ is the temperature and pressure dependent dynamic viscosity $\left[\frac{M L}{T^{2}} \cdot T\right]$,
- $\varrho_{f}(x, t, T, p)$ is the temperature and pressure dependent mass density of the fluid $\left[\frac{M}{L^{3}}\right]$,
- $g$ the gravitation constant $\left[\frac{L}{T^{2}}\right]$

### 1.2.2 Balance Equation

In the groundwater flow the function $S$ in the balance equations (1.8) or (1.9) is replaced by $\phi \rho(x, t, p)$ :

$$
\begin{equation*}
\frac{\partial \phi \rho(p, T)}{\partial t}-\operatorname{div} \frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi-Q(x, t)=0 \tag{1.11}
\end{equation*}
$$

where

- $\phi$ is the porosity of the solid
- $Q(x, t)$ describes the inner sources or sinks, in coupled processes sources and sinks $Q$ can also result from changes of the other primary variable, for instance through the changing of temperature sources and sinks can arise

For the implementation it is assumed that the medium is incompressible, i.e., the porosity does not change and thus $\frac{\partial \phi}{\partial t}=0$. Therefore, the first term of 1.12 is

$$
\frac{\partial \phi \rho(p, T)}{\partial t}=\underbrace{\frac{\partial \phi}{\partial t}}_{=0} \varrho(p, T)+\phi \frac{\partial \varrho(p, T)}{\partial t}=\phi\left(\frac{\partial \varrho}{\partial p} \frac{\partial p}{\partial t}+\frac{\partial \varrho}{\partial T} \frac{\partial T}{\partial t}\right)
$$

As a part of the Boussinesq approximation it is assumed that the last term of the above equation $\frac{\partial \varrho}{\partial T} \frac{\partial T}{\partial t}$ vanishes.
Furthermore, it is assumed that the density depends linearly on the pressure, i.e. $\frac{\partial \varrho}{\partial p}$ is constant.
Under this assumptions it is possible to summarize the (constant) porosity and the constant derivation $\frac{\partial \varrho}{\partial p}$ into a new constant $S$ and (1.11) changes to

$$
\begin{equation*}
S \frac{\partial p}{\partial t}-\operatorname{div}\left[\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi\right]-Q(x, t)=0 \tag{1.12}
\end{equation*}
$$

### 1.2.3 Boundary Conditions

$$
\begin{array}{rll}
p-g_{D, p}=0 & \text { on } \quad \Gamma_{D} & \text { (Dirichlet type boundary conditions) } \\
\left\langle\left.\frac{\boldsymbol{\kappa}}{\mu(T)} \operatorname{grad} \Psi \right\rvert\, n\right\rangle+g_{N}=0 \quad \text { on } \quad \Gamma_{N} & \text { (Neumann type boundary conditions) } \tag{1.14}
\end{array}
$$

### 1.2.4 Weak Formulation

Multiplying (1.12) with -1 and summing up with (1.14) leads to

$$
\begin{equation*}
-S \frac{\partial p}{\partial t}+\operatorname{div}\left[\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi\right]+Q(x, t)+\left\langle\left.\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi \right\rvert\, n\right\rangle+g_{N}=0 . \tag{1.15}
\end{equation*}
$$

Since (1.15) holds true for arbitrary points of the domain, the equation stays valid if it is multiplied by test functions $v, \bar{v} \in H_{0}^{1}(\Omega)$ and the integration over the domain $\Omega$ and the Neumann boundary $\Gamma_{N, p}$, respectively:

$$
\begin{equation*}
\int_{\Omega} v\left(-S \frac{\partial p}{\partial t}+\operatorname{div}\left[\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi\right]+Q(x, t)\right) \mathrm{d} x+\int_{\Gamma_{N}} \bar{v}\left(\left\langle\left.\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi \right\rvert\, n\right\rangle+g_{N}\right) \mathrm{d} \sigma=0 \tag{1.16}
\end{equation*}
$$

or equivalently
(1.17)
$-\int_{\Omega} v S \frac{\partial p}{\partial t} \mathrm{~d} x+\int_{\Omega} v\left(\operatorname{div}\left[\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi\right]\right) \mathrm{d} x+\int_{\Omega} v Q(x, t) \mathrm{d} x+\int_{\Gamma_{N}} \bar{v}\left\langle\left.\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi \right\rvert\, n\right\rangle \mathrm{d} \sigma+\int_{\Gamma_{N}} \bar{v} g_{N} \mathrm{~d} \sigma=0$
Integration by parts of the second term of 1.17 results in

$$
\begin{equation*}
\int_{\Omega} v\left(\operatorname{div}\left[\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi\right]\right) \mathrm{d} x=-\int_{\Omega}\left\langle\boldsymbol{\operatorname { g r a d }} v \left\lvert\, \frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi\right.\right\rangle \mathrm{d} x+\int_{\Omega} \operatorname{div}\left(v\left[\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi\right]\right) \mathrm{d} x \tag{1.18}
\end{equation*}
$$

Using Green's formula for the last term of the above expression

$$
\begin{array}{rl}
\int_{\Omega} \operatorname{div}\left(v\left[\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi\right]\right) \mathrm{d} & x=\oint_{\Gamma}\left\langle\left. v \frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi \right\rvert\, n\right\rangle \mathrm{d} \sigma  \tag{1.19}\\
& =\int_{\Gamma_{D}}\left\langle\left. v \frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi \right\rvert\, n\right\rangle \mathrm{d} \sigma+\int_{\Gamma_{N}}\left\langle\left. v \frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi \right\rvert\, n\right\rangle \mathrm{d} \sigma
\end{array}
$$

and the integral on the Dirichlet boundary $\Gamma_{D}$ vanishes because $v=0$ holds.
Finally, the expression (1.18) takes the form

$$
\begin{equation*}
\int_{\Omega} v\left[\operatorname{div}\left(\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi\right)\right] \mathrm{d} x=-\int_{\Omega}\langle\boldsymbol{\operatorname { g r a d }} v \mid \boldsymbol{K} \boldsymbol{\operatorname { g r a d }} \Psi\rangle \mathrm{d} x+\int_{\Gamma_{N}}\left\langle\left. v \frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi \right\rvert\, n\right\rangle \mathrm{d} \sigma \tag{1.20}
\end{equation*}
$$

Putting (1.20) in (1.17) yields to

$$
\begin{aligned}
0= & -\int_{\Omega} v S \frac{\partial p}{\partial t} \mathrm{~d} x-\int_{\Omega}\left\langle\boldsymbol{\operatorname { g r a d }} v \left\lvert\, \frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi\right.\right\rangle \mathrm{d} x+\int_{\Gamma_{N}}\left\langle\left. v \frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi \right\rvert\, n\right\rangle \mathrm{d} \sigma \\
& +\int_{\Omega} v Q(x, t) \mathrm{d} x+\int_{\Gamma_{N}} \bar{v}\left\langle\left.\frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi \right\rvert\, n\right\rangle \mathrm{d} \sigma+\int_{\Gamma_{N}} \bar{v} g_{N} \mathrm{~d} \sigma
\end{aligned}
$$

Since the test functions are arbitrary, by setting $v=-\bar{v}$ the second and fourth term cancel each other.
Multiplying by -1 results in

$$
\begin{equation*}
0=\int_{\Omega} v S \frac{\partial p}{\partial t} \mathrm{~d} x+\int_{\Omega}\left\langle\boldsymbol{\operatorname { g r a d }} v \left\lvert\, \frac{\boldsymbol{\kappa}}{\mu(T)} \boldsymbol{\operatorname { g r a d }} \Psi\right.\right\rangle \mathrm{d} x-\int_{\Omega} v Q(x, t) \mathrm{d} x-\int_{\Gamma_{N}} v g_{N} \mathrm{~d} \sigma \tag{1.21}
\end{equation*}
$$

### 1.2.5 Finite Element Discretization

$$
\begin{equation*}
p \approx \sum N_{j} a_{j}=N a, \quad \Psi \approx \sum N_{j} a_{j}+\varrho_{f} g z \tag{1.22}
\end{equation*}
$$

where $N_{i}(x, y, z)$ are the shape functions and $a_{i}$ are coefficients. Galerkin principle:

$$
\begin{equation*}
v=N_{i} \tag{1.23}
\end{equation*}
$$

Substituting (1.22) and (1.23) in (1.21) leads to

$$
\begin{equation*}
\left[\int_{\Omega} N_{i} S N_{j} \mathrm{~d} x\right] \frac{\partial a_{j}}{\partial t}+\left[\int_{\Omega} \nabla^{T} N_{i} \frac{\boldsymbol{\kappa}}{\mu(T)} \nabla N \mathrm{~d} x\right] a+\int_{\Omega} \nabla^{T} N_{i} \frac{\boldsymbol{\kappa} \varrho_{f} g}{\mu(T)} e_{z} \mathrm{~d} x-\int_{\Omega} N_{i} Q(x, t) \mathrm{d} x-\int_{\Gamma_{N}} N_{i} g_{N} \mathrm{~d} \sigma=0, \tag{1.24}
\end{equation*}
$$

$i, j=1, \ldots, n$, which is a set of linear equations of the form

$$
\begin{equation*}
\boldsymbol{C} \dot{a}+\boldsymbol{K} a+f=0 \tag{1.25}
\end{equation*}
$$

with

$$
\begin{align*}
\boldsymbol{C}_{i j} & =\int_{\Omega} N_{i} S N_{j} \mathrm{~d} x  \tag{1.26}\\
\boldsymbol{K}_{i j} & =\int_{\Omega} \nabla^{T} N_{i} \frac{\boldsymbol{\kappa}}{\mu(T)} \nabla N_{j} \mathrm{~d} x  \tag{1.27}\\
f_{i} & =-\int_{\Omega} N_{i} Q(x, t) \mathrm{d} x-\int_{\Gamma_{N}} N_{i} g_{N} \mathrm{~d} \sigma+\int_{\Omega} \nabla^{T} N_{i} \frac{\boldsymbol{\kappa} \varrho_{f} g}{\mu(T)} e_{z} \mathrm{~d} x \tag{1.28}
\end{align*}
$$

### 1.3 Heat Conduction Equation
<!-- vale off -->
### 1.3.1 Fick's Law
<!-- vale on -->
$$
\begin{equation*}
J=-\lambda \boldsymbol{g r a d} T \tag{1.29}
\end{equation*}
$$

where $T$ is the temperature, $\lambda$ is the hydrodynamic thermo-dispersion tensor and $J$ the heat flow

### 1.3.2 Balance Equation

The heat conduction equation is also known as the advection-diffusion or convection-diffusion equation.

$$
\begin{equation*}
\underbrace{\left[\varrho_{f}(x) \phi c_{f}+\varrho_{s}(x)(1-\phi) c_{s}\right]}_{:=c_{p}} \cdot \frac{\partial T}{\partial t}-\operatorname{div}(\lambda \boldsymbol{\operatorname { g r a d }} T)+\varrho_{f} \cdot c_{f}(x) \cdot\langle q \mid \boldsymbol{\operatorname { g r a d }} T\rangle=0 \tag{1.30}
\end{equation*}
$$

where

- $q$ is the Darcy velocity defined by
- $\varrho_{s}$ is the solid density $\left[\frac{M}{L^{3}}\right]$
- $\lambda$ hydrodynamic thermo-dispersion tensor,

$$
\lambda=\lambda^{\text {cond }}+\lambda^{\text {disp }}
$$

- $\lambda^{\text {cond }}(\phi)=\phi \lambda_{f}+(1-\phi) \lambda_{s}$ is the thermal conductivity
- $\lambda^{\operatorname{disp}}\left(\alpha_{T}, \alpha_{L}\right)=\varrho_{f} c_{f}\left[\alpha_{T}\|q\|_{2} \boldsymbol{I}+\left(\alpha_{L}-\alpha_{T}\right) \frac{q q^{T}}{\|q\|_{2}}\right]$ is the thermal dispersivity, where $\alpha_{T}$ is the transverse thermo-dispersivity and $\alpha_{L}$ is the longitudinal thermo-dispersivity of the fluid
- $\phi$ is the porosity

### 1.3.3 Boundary Conditions

$$
\begin{array}{rlll}
T & =g_{D}^{T} & \text { on } & \Gamma_{D}
\end{array} \quad \text { (Dirichlet type boundary conditions) }
$$

### 1.3.4 Weak Formulation

The integration of the reformulated Neumann type boundary condition, i.e., $\langle\lambda \boldsymbol{\operatorname { g r a d }} T \mid n\rangle+g_{N}^{T}=0$, into (1.30), multiplying with arbitrary test functions $v, \bar{v} \in H_{0}^{1}(\Omega)$ and integration over $\Omega$ results in

$$
\begin{align*}
& -\int_{\Omega} v \cdot \operatorname{div}(\lambda \boldsymbol{\operatorname { g r a d }} T) \mathrm{d} \Omega+\int_{\Omega} v \cdot \varrho_{f} \cdot c_{f}(x)\langle q \mid \boldsymbol{\operatorname { g r a d }} T\rangle \mathrm{d} \Omega  \tag{1.33}\\
& +\int_{\Omega} v \cdot c_{p}(x) \cdot \frac{\partial T}{\partial t} \mathrm{~d} \Omega+\int_{\Gamma_{N}} \bar{v} \cdot\left[\langle\lambda \boldsymbol{\operatorname { g r a d }} T \mid n\rangle+g_{N}^{T}\right] \mathrm{d} \sigma=0
\end{align*}
$$

Integration by parts of the first term in the above equation yields:

$$
\begin{equation*}
\int_{\Omega} v \operatorname{div}[\lambda \boldsymbol{\operatorname { g r a d }} T] \mathrm{d} \Omega=-\int_{\Omega}\langle\boldsymbol{\operatorname { g r a d }} v \mid[\lambda \boldsymbol{\operatorname { g r a d }} T]\rangle \mathrm{d} \Omega+\int_{\Omega} \operatorname{div}[v \lambda \boldsymbol{\operatorname { g r a d }} T] \mathrm{d} \Omega \tag{1.34}
\end{equation*}
$$

Using Green's formulae for the last term of the above expression

$$
\begin{equation*}
\int_{\Omega} \operatorname{div}[v \lambda \boldsymbol{\operatorname { g r a d }} T] \mathrm{d} \Omega=\oint_{\Gamma}\langle v \lambda \boldsymbol{\operatorname { g r a d }} T \mid n\rangle \mathrm{d} \sigma=\int_{\Gamma_{D}}\langle v \lambda \boldsymbol{\operatorname { g r }} \boldsymbol{a} \boldsymbol{d} T \mid n\rangle \mathrm{d} \sigma+\int_{\Gamma_{N}}\langle v \lambda \boldsymbol{\operatorname { g r }} \boldsymbol{a} \boldsymbol{d} T \mid n\rangle \mathrm{d} \sigma \tag{1.35}
\end{equation*}
$$

and since $v$ vanishes on $\Gamma_{D}$ the integral over $\Gamma_{D}$ also vanishes, this leads to

$$
\begin{equation*}
\int_{\Omega} v[\operatorname{div}(\lambda \boldsymbol{\operatorname { g r a d }} T)] \mathrm{d} \Omega=-\int_{\Omega}\langle\boldsymbol{\operatorname { g r a d }} v \mid \lambda \boldsymbol{\operatorname { g r a d }} T\rangle \mathrm{d} \Omega+\int_{\Gamma_{N}}\langle v \lambda \boldsymbol{\operatorname { g r a d }} T \mid n\rangle \mathrm{d} \sigma \tag{1.36}
\end{equation*}
$$

Thus (1.33) reads:

$$
\begin{align*}
0=- & \int_{\Omega}\langle\boldsymbol{\operatorname { g r a d }} v \mid \lambda \boldsymbol{\operatorname { g r a d }} T\rangle \mathrm{d} \Omega+\int_{\Gamma_{N}}\langle v \lambda \boldsymbol{\operatorname { g r a d }} T \mid n\rangle \mathrm{d} \sigma+\int_{\Omega} v \cdot \varrho_{f} \cdot c_{f}(x)\langle q \mid \boldsymbol{\operatorname { g r a d }} T\rangle \mathrm{d} \Omega \\
& +\int_{\Omega} v \cdot c_{p}(x) \cdot \frac{\partial T}{\partial t} \mathrm{~d} \Omega+\int_{\Gamma_{N}} \bar{v} \cdot\langle\lambda \boldsymbol{\operatorname { g r a d }} T \mid n\rangle \mathrm{d} \sigma+\int_{\Gamma_{N}} \bar{v} \cdot g_{N}^{T} \mathrm{~d} \sigma \tag{1.37}
\end{align*}
$$

Setting $v=-\bar{v}$ and multiplying by -1 :
(1.38)
$\int_{\Omega}\langle\boldsymbol{\operatorname { g r a d }} v \mid \lambda \boldsymbol{\operatorname { g r a d }} T\rangle \mathrm{d} \Omega-\int_{\Omega} v \cdot \varrho_{f} \cdot c_{f}(x)\langle q \mid \boldsymbol{\operatorname { g r a d }} T\rangle \mathrm{d} \Omega-\int_{\Omega} v \cdot c_{p}(x) \cdot \frac{\partial T}{\partial t} \mathrm{~d} \Omega+\int_{\Gamma_{N}} v \cdot g_{N}^{T} \mathrm{~d} \sigma=0$

### 1.3.5 Finite Element Discretization

Analogously to the approximation (1.22) the temperature is approximated by:

$$
\begin{equation*}
T \approx \sum N_{i} b_{i}=N b \tag{1.39}
\end{equation*}
$$

using the same shape functions $N_{i}$ and time dependent coefficients $b_{i}$.
Using the shape functions again as test functions (Galerkin principle (1.23) the discretization of (1.38)) takes the following form

$$
\begin{align*}
& 0=\left[\int_{\Omega} \nabla^{T} N_{i} \lambda \nabla N \mathrm{~d} \Omega-\int_{\Omega} N_{i} \cdot \varrho_{f} \cdot c_{f}(x) \cdot q^{T} \cdot \nabla N \mathrm{~d} \Omega\right] b  \tag{1.40}\\
& +\left[\int_{\Omega} N_{i} \cdot c_{p}(x) N \mathrm{~d} \Omega\right] \frac{\mathrm{d} b}{\mathrm{~d} t}+\int_{\Gamma_{N}} N_{i} g_{N}^{T} \mathrm{~d} \sigma \quad(i=1, \ldots, n)
\end{align*}
$$

This is again a set of equations of the form

$$
\begin{equation*}
\boldsymbol{C} \dot{b}+\boldsymbol{K} b+f=0 \tag{1.41}
\end{equation*}
$$

with

$$
\begin{align*}
\boldsymbol{K}_{i j} & =\int_{\Omega} \nabla^{T} N_{i} \lambda \nabla N_{j} \mathrm{~d} \Omega-\int_{\Omega} N_{i} \varrho_{f} c_{f}\left\langle q \mid \nabla N_{j}\right\rangle \mathrm{d} \Omega  \tag{1.42}\\
f_{i} & =-\int_{\Omega} N_{i} Q(x, t) \mathrm{d} \Omega-\int_{\Gamma_{N}} N_{i} g_{N}^{T} \mathrm{~d} \sigma  \tag{1.43}\\
\boldsymbol{C}_{i j} & =\int_{\Omega} N_{i} c_{p} N_{j} \mathrm{~d} \Omega \tag{1.44}
\end{align*}
$$

### 1.4 Coupling the Processes

The heat conduction process is coupled with the confined groundwater flow process by the advective term in 1.30.

The fluid density $\varrho_{f}$ as well as the viscosity $\mu$ used in (1.10) (and respectively (1.12)) are coupled to the heat conduction process by their temperature dependencies.

For the fluid density the following linear dependency

$$
\begin{equation*}
\varrho_{f}(T)=\varrho_{\mathrm{ref}}\left(1-\beta(x)\left(T-T_{\mathrm{ref}}\right)\right) \tag{1.45}
\end{equation*}
$$

is implemented, where the fluid thermal expansion coefficient $\beta(x)\left[K^{-1}\right]$ depends on the medium and $T_{\text {ref }}$ is the reference temperature.

The temperature dependency of the fluid viscosity is implemented by the function:

$$
\begin{equation*}
\mu(T)=\mu_{0} \mathrm{e}^{-\frac{T-T_{c}}{T_{\nu}}} \tag{1.46}
\end{equation*}
$$

There is not implemented any coupling by source and sink terms in OGS-6, i.e., the density changes due to temperature changes affects only the buoyancy term $\varrho_{f}(T) \cdot g \cdot z$ in (1.10).
The currently implemented coupling schema is referred to as the Boussinesq approximation.

These simplified EOS are those used in the large-scale benchmark allowing to apply linear stability analysis of HT problem [6].

The IAPWS EOS for fluid density and viscosity are also implemented into OGS.

## Bibliography
<!-- vale off -->
[1] Lutz Angermann and Peter Knabner. Numerical Methods for Elliptic and Parabolic Partial Differential Equations, volume 44 of Texts in Applied Mathematics. Springer New York, 2003.

[2] H.-J.G. Diersch and O. Kolditz. Coupled groundwater flow and transport: 2. thermohaline and 3d convection systems. Advances in Water Resources, 21(5):401-425, 1998.

[3] H.-J.G. Diersch and O. Kolditz. Variable-density flow and transport in porous media: approaches and challenges. Advances in Water Resources, 25(8-12):899-944, 2002.

[4] P. Knabner, Ch. Tapp, and K. Thiele. Adaptivity in the finite volume discretization of variable density flows in porous media. Physics and Chemistry of the Earth, Part B: Hydrology, Oceans and Atmosphere, 26(4):319-324, 2001.

[5] Olaf Kolditz, Rainer Ratke, Hans-Jörg G. Diersch, and Werner Zielke. Coupled groundwater flow and transport: 1. verification of variable density flow and transport models. Advances in Water Resources, 21(1):27-46, 1998.

[6] V.I. Malkovsky and F. Magri. Thermal convection of temperature-dependent viscous fuids within three-dimensional faulted geothermal systems: Estimation from linear and numerical analyses. Water Resources Research, (52):2855-2867, 2016.
