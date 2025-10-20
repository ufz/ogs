+++
date = "2021-06-17T11:00:13+01:00"
title = "OGS Units"
author = "Thomas Nagel"
weight = 72
+++

<!-- Content taken from docs/userguide/basic/conventions/units_ogs.pdf -->

## 1 Preliminary remarks

OGS does not have any module that automatically takes care of unit consistency.
This is the user's job.
In many cases, units are obvious, but in some cases one might have to burn some sugar in one's brain. What helps, is looking at the implemented equations and making sure that indeed all units are consistent.
The PDE alone may not be sufficient for this exercise; instead, the weak form provides the most insight into the implemented situation.
We'll illustrate this here using a simple example.

This document touches but a minor part of the conventions used in OGS.
More important definitions are listed in our documentation.

Typical questions by users regarding units are:

- In what units should I give the Neumann boundary conditions?
- In what units should I give the source terms?
- What signs do the boundary conditions and source terms have?
- Do the units change when I switch from 3D, to 2D, or to axisymmetric?
- ...

## 2 Example: fluid phase mass balance

### 2.1 The PDE

Most of our processes involve a hydraulic process of some sort.
We take a simple saturated flow case as an example here, departing from a general source-free mass balance.

$$
\begin{gather*}
\frac{\mathrm{d}_{\alpha} \varrho_{\alpha}}{\mathrm{d} t}+\varrho_{\alpha} \operatorname{div} \mathbf{v}_{\alpha}=0  \tag{1}\\
\frac{\partial \varrho_{\alpha}}{\partial t}+\operatorname{div}\left(\varrho_{\alpha} \mathbf{v}_{\alpha}\right)=0 \tag{2}
\end{gather*}
$$

From the mass balance of an intrinsically compressible solid we can obtain:

$$
\begin{equation*}
\frac{\mathrm{d}_{\mathrm{S}} \phi}{\mathrm{~d} t}=\left(\alpha_{\mathrm{B}}-\phi\right)\left[\operatorname{div} \mathbf{v}_{\mathrm{S}}+\beta_{p, \mathrm{SR}} \frac{\mathrm{~d}_{\mathrm{S}} p}{\mathrm{~d} t}\right] \tag{3}
\end{equation*}
$$

Allowing the fluid-phase to be compressible as well, we arrive at a typical formulation:

$$
\begin{equation*}
\varrho_{\mathrm{LR}}\left[\phi \beta_{p, \mathrm{LR}}+\left(\alpha_{\mathrm{B}}-\phi\right) \beta_{p, \mathrm{SR}}\right] \frac{\mathrm{d}_{\mathrm{S}} p_{\mathrm{LR}}}{\mathrm{~d} t}+\operatorname{div}\left(\varrho_{\mathrm{LR}} \widetilde{\mathbf{w}}_{\mathrm{LS}}\right)+\alpha_{\mathrm{B}} \varrho_{\mathrm{LR}} \operatorname{div} \mathbf{v}_{\mathrm{S}}=0 \tag{4}
\end{equation*}
$$

Typically, Darcy's laws will be used:

$$
\begin{equation*}
\widetilde{\mathbf{w}}_{\mathrm{LS}}=-\frac{\mathbf{k}}{\mu_{\mathrm{LR}}}\left(\operatorname{grad} p_{\mathrm{LR}}-\varrho_{\mathrm{LR}} \mathbf{b}\right) \tag{5}
\end{equation*}
$$

### 2.2 The weak form

The finite element method rests on weighted residuals.
More specifically, OGS uses a Bubnov-Galerkin scheme where test functions are introduced to transform the strong form of the PDE into its weak form.
Using standard function spaces we arrive at the following form after a few manipulations

$$
\begin{align*}
& \int_{\Omega} v^{p} \varrho_{\mathrm{LR}}\left[\phi \beta_{p, \mathrm{LR}}+\left(\alpha_{\mathrm{B}}-\phi\right) \beta_{p, \mathrm{SR}}\right] \frac{\mathrm{d}_{\mathrm{S}} p_{\mathrm{LR}}}{\mathrm{~d} t} \mathrm{~d} \Omega-\int_{\Omega} \varrho_{\mathrm{LR}} \widetilde{\mathbf{w}}_{\mathrm{LS}} \cdot \operatorname{grad} v^{p} \mathrm{~d} \Omega+  \tag{6}\\
& +\int_{\Omega} \alpha_{\mathrm{B}} \varrho_{\mathrm{LR}} \operatorname{div} \mathbf{v}_{\mathrm{S}} \mathrm{~d} \Omega=\int_{\partial \Omega_{\mathrm{N}}} v^{p} \overline{\dot{m}}_{\mathrm{LS}} \mathrm{~d} \Gamma \tag{7}
\end{align*}
$$

which also give rise to the natural (Neumann) boundary conditions.

### 2.3 Boundary conditions

Now, the boundary conditions can be clearly defined

$$
\begin{align*}
p & =\bar{p} \quad \text { on } \partial \Omega_{\mathrm{D}}^{p}  \tag{8}\\
-\varrho_{\mathrm{LR}} \widetilde{\mathbf{w}}_{\mathrm{LS}} \cdot \mathbf{n} & =\overline{\dot{m}}_{\mathrm{LS}} \quad \text { on } \partial \Omega_{\mathrm{N}}^{p} \tag{9}
\end{align*}
$$

where $\partial \Omega_{\mathrm{D}}$ and $\partial \Omega_{\mathrm{N}}$ are complementary, $\partial \Omega_{\mathrm{D}} \cup \partial \Omega_{\mathrm{N}}=\Gamma$ and $\partial \Omega_{\mathrm{D}} \cap \partial \Omega_{\mathrm{N}}=0$.

### 2.4 Discretization

Space (finite element) and time (backward Euler) discretization finally yield

$$
\begin{align*}
& \int_{\Omega} \mathbf{N}^{\mathrm{T}} \varrho_{\mathrm{LR}}\left[\phi \beta_{p, \mathrm{LR}}+\left(\alpha_{\mathrm{B}}-\phi\right) \beta_{p, \mathrm{SR}}\right] \mathbf{N}^{\mathrm{T}} \mathrm{~d} \Omega \frac{\hat{\mathbf{p}}_{\mathrm{LR}}^{t+1}-\hat{\mathbf{p}}_{\mathrm{LR}}^{t}}{\Delta t}+\int_{\Omega} \nabla \mathbf{N}^{\mathrm{T}} \varrho_{\mathrm{LR}} \frac{\mathbf{k}}{\mu_{\mathrm{LR}}} \nabla \mathbf{N} \mathrm{~d} \Omega \hat{\mathbf{p}}_{\mathrm{LR}}- \\
& -\int_{\Omega} \nabla \mathbf{N}^{\mathrm{T}} \varrho_{\mathrm{LR}}^{2} \frac{\mathbf{k}}{\mu_{\mathrm{LR}}} \mathbf{b} \mathrm{~d} \Omega+\int_{\Omega} \mathbf{N}^{\mathrm{T}} \alpha_{\mathrm{B}} \varrho_{\mathrm{LR}} \mathbf{I}^{\mathrm{T}} \mathbf{B}_{u} \mathrm{~d} \Omega \frac{\hat{\mathbf{u}}^{t+1}-\hat{\mathbf{u}}^{t}}{\Delta t}=\int_{\partial \Omega_{\mathrm{N}}} \mathbf{N}^{\mathrm{T}} \overline{\dot{m}}_{\mathrm{LS}} \mathrm{~d} \Gamma \tag{10}
\end{align*}
$$

Source terms $r$ are added in a general manner to the RHS:

$$
\begin{align*}
& \int_{\Omega} \mathbf{N}^{\mathrm{T}} \varrho_{\mathrm{LR}}\left[\phi \beta_{p, \mathrm{LR}}+\left(\alpha_{\mathrm{B}}-\phi\right) \beta_{p, \mathrm{SR}}\right] \mathbf{N}^{\mathrm{T}} \mathrm{~d} \Omega \frac{\hat{\mathbf{p}}_{\mathrm{LR}}^{t+1}-\hat{\mathbf{p}}_{\mathrm{LR}}^{t}}{\Delta t}+\int_{\Omega} \nabla \mathbf{N}^{\mathrm{T}} \varrho_{\mathrm{LR}} \frac{\mathbf{k}}{\mu_{\mathrm{LR}}} \nabla \mathbf{N} \mathrm{~d} \Omega \hat{\mathbf{p}}_{\mathrm{LR}} \\
& -\int_{\Omega} \nabla \mathbf{N}^{\mathrm{T}} \varrho_{\mathrm{LR}}^{2} \frac{\mathbf{k}}{\mu_{\mathrm{LR}}} \mathbf{b} \mathrm{~d} \Omega+\int_{\Omega} \mathbf{N}^{\mathrm{T}} \alpha_{\mathrm{B}} \varrho_{\mathrm{LR}} \mathbf{I}^{\mathrm{T}} \mathbf{B}_{u} \mathrm{~d} \Omega \frac{\hat{\mathbf{u}}^{t+1}-\hat{\mathbf{u}}^{t}}{\Delta t}=\int_{\partial \Omega_{\mathrm{N}}} \mathbf{N}^{\mathrm{T}} \overline{\dot{m}}_{\mathrm{LS}} \mathrm{~d} \Gamma+\int_{\Omega} \mathbf{N}^{\mathrm{T}} r \mathrm{~d} \Omega=\hat{\mathbf{f}}_{p} \tag{11}
\end{align*}
$$

## 3 Observations

- The integration domain $\Omega$ (mesh) of dimension $N$ is associated with a boundary domain $\partial \Omega$ of codimension 1 (i.e. area elements $\Gamma$ are of dimension $N-1$ ).

- The Neumann bc is here positive for an inflow of mass. Likewise the source term definition on the RHS fixes the sign. Neumann boundary conditions for vectorial quantities (e.g. tractions) defined with respect to the Cartesian reference system are positive (negative) when pointing in the positive (negative) coordinate direction defined by the mesh [^0].

- Each integral of the weak form has to yield the same dimension. This fact helps to fix the units and also represents the units of the nodal fluxes.
- The (generalized) nodal forces $\hat{\mathbf{f}}_{p}$ represent the residuals before the addition of source and Neumann terms. This helps to evaluate mass flow rates, heat flow rates, reaction forces, etc. in response to external input.
- The user typically chooses a unit system convenient for the problem, for example: m N s K . Based on this choice, there needs to be consistency between
- The units of the physical properties and variables (density, compressibility, permeability, pressure, ...)
- The spatial units of the mesh
- The temporal units given in the time stepping scheme

## 4 Example

Basic choices are made already when specifying the mesh and the time stepping scheme, fixing length and time units, respectively.
Primary and secondary variables follow next.
Boundary and source terms are of particular interest and may seem non-trivial when entities with different dimensions are linked, e.g. 2D sources in 3D domains.
The axisymmetric case integrates over $r \mathrm{~d} \Omega$, thus solves a 3D problem on a 2D mesh.
Table 1 gives some examples in different unit systems.
Some rest on base units ( kg m s ), others make use of derived units for practical convenience (mm N d).

The source terms $r$ referenced in Tab. 1 are defined in OGS as volumetrically distributed source terms over subdomains of varying dimension.
OGS also offers the possibility to use nodal source terms, i.e. directly specified nodal forces.
This corresponds to an already performed integration on the discrete system, i.e.

$$
\begin{equation*}
\mathbf{r}_{p}=\int_{\Omega} \mathbf{N}^{\mathrm{T}} r \mathrm{~d} \Omega \tag{12}
\end{equation*}
$$

where directly the value for $\mathbf{r}_{p}$ will be specified at certain nodes.
Remark on the term (generalized) nodal forces: The term generalized force stems from the fact that in solid mechanical problems nodal forces in the proper sense are assembled to obtain mechanical equilibrium.
Once we move to thermo-dynamical systems with generalized degrees of freedom, we still speak of the associated thermodynamic forces by analogy.
Another useful structural similarity can be seen when recognizing that a force is a momentum rate: $1 \mathrm{~N}=1 \mathrm{~kg} \mathrm{~m} \mathrm{~s}^{-1} \mathrm{~s}^{-1}$.
While the mechanical system solves a momentum balance by means of equilibrium of nodal forces (momentum rates), transport processes solve mass balances by balancing nodal mass rates in $\mathrm{kg} \mathrm{s}^{-1}$, heat transport processes solve energy balances by means of nodal energy rates in $\mathrm{J} \mathrm{s}^{-1}$, and so on.
Thus, one can speak of generalized forces or, when normalized to an area (i.e. per square meter), of generalized fluxes.

**Table 1.** The following is based on Eq. (11).

| Object                | Symbol                     | Unit in 3D (m Ns) | Unit in 3D (mm Nᵈ) | Unit in 2D (m Ns) | Unit in 3D (kg ms) | Unit in 3D (kg mm s) |
|-----------------------|----------------------------|-------------------|---------------------|-------------------|--------------------|----------------------|
| **Domain mesh**       | $\Omega$                   | $\mathrm{m}^3$    | $\mathrm{mm}^3$     | $\mathrm{m}^2$    | $\mathrm{m}^3$     | $\mathrm{mm}^3$      |
| **Boundary mesh**     | $\Gamma$                   | $\mathrm{m}^2$    | $\mathrm{mm}^2$     | m                 | $\mathrm{m}^2$     | $\mathrm{mm}^2$      |
| **Time discretization** | $\Delta t$                | s                 | d                   | s                 | s                  | s                    |
| **Generalized nodal forces** | $\hat{\mathbf{f}}_p$ | $\mathrm{kg}\,\mathrm{s}^{-1}$ | $\mathrm{kg}\,\mathrm{d}^{-1}$ | $\mathrm{kg}\,\mathrm{s}^{-1}\,\mathrm{m}^{-1}$ | $\mathrm{kg}\,\mathrm{s}^{-1}$ | $\mathrm{kg}\,\mathrm{s}^{-1}$ |
| **Pore pressure**    | $p_{\text{LR}}$            | Pa                | MPa                 | Pa                | Pa                 | kPa                  |
| **Displacement**      | $u$                        | m                 | mm                  | m                 | m                  | mm                   |
| **Density**           | $\varrho_{\text{LR}}$      | $\mathrm{kg}\,\mathrm{m}^{-3}$ | $\mathrm{kg}\,\mathrm{mm}^{-3}$ | $\mathrm{kg}\,\mathrm{m}^{-3}$ | $\mathrm{kg}\,\mathrm{m}^{-3}$ | $\mathrm{kg}\,\mathrm{mm}^{-3}$ |
| **Permeability**      | $k$                        | $\mathrm{m}^2$    | $\mathrm{mm}^2$     | $\mathrm{m}^2$    | $\mathrm{m}^2$     | $\mathrm{mm}^2$      |
| **Viscosity**         | $\mu_{\text{LR}}$          | Pa s              | MPa d               | Pa s              | Pa s               | kPa s                |
| **Specific body force** | $b$                     | $\mathrm{N}\,\mathrm{kg}^{-1}= \mathrm{m}\,\mathrm{s}^{-2}$ | $\mathrm{N}\,\mathrm{kg}^{-1}= \mathrm{m}\,\mathrm{s}^{-2}$ | $\mathrm{N}\,\mathrm{kg}^{-1}= \mathrm{m}\,\mathrm{s}^{-2}$ | $\mathrm{N}\,\mathrm{kg}^{-1}= \mathrm{m}\,\mathrm{s}^{-2}$ | $\mathrm{mN}\,\mathrm{kg}^{-1}= \mathrm{mm}\,\mathrm{s}^{-2}$ |
| **Compressibility**   | $\beta_{p,\alpha R}$       | $\mathrm{Pa}^{-1}$ | $\mathrm{MPa}^{-1}$ | $\mathrm{Pa}^{-1}$ | $\mathrm{Pa}^{-1}$ | $\mathrm{kPa}^{-1}$ |
| **Boundary flux**     | $\dot{m}_{\text{LS}}$      | $\mathrm{kg}\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}$ | $\mathrm{kg}\,\mathrm{mm}^{-2}\,\mathrm{d}^{-1}$ | $\mathrm{kg}\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}$ | $\mathrm{kg}\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}$ | $\mathrm{kg}\,\mathrm{mm}^{-2}\,\mathrm{s}^{-1}$ |
| **Source term on 3D** | $r$                        | $\mathrm{kg}\,\mathrm{m}^{-3}\,\mathrm{s}^{-1}$ | $\mathrm{kg}\,\mathrm{mm}^{-3}\,\mathrm{d}^{-1}$ | – | $\mathrm{kg}\,\mathrm{m}^{-3}\,\mathrm{s}^{-1}$ | $\mathrm{kg}\,\mathrm{mm}^{-3}\,\mathrm{s}^{-1}$ |
| **Source term on 2D** | $r$                        | $\mathrm{kg}\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}$ | $\mathrm{kg}\,\mathrm{mm}^{-2}\,\mathrm{d}^{-1}$ | $\mathrm{kg}\,\mathrm{m}^{-3}\,\mathrm{s}^{-1}$[^1] | $\mathrm{kg}\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}$ | $\mathrm{kg}\,\mathrm{mm}^{-2}\,\mathrm{s}^{-1}$ |
| **Source term on 1D** | $r$                        | $\mathrm{kg}\,\mathrm{m}^{-1}\,\mathrm{s}^{-1}$ | $\mathrm{kg}\,\mathrm{mm}^{-1}\,\mathrm{d}^{-1}$ | $\mathrm{kg}\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}$[^2] | $\mathrm{kg}\,\mathrm{m}^{-1}\,\mathrm{s}^{-1}$ | $\mathrm{kg}\,\mathrm{mm}^{-1}\,\mathrm{s}^{-1}$ |

[^0]: The source term here is assigned as a function defined on the entire domain $\Omega$.
Practically, the source term is often assigned to a subdomain $\Omega_{\mathrm{s}} \subseteq \Omega$.
It can also be assigned to lower-dimensional entities (surfaces, lines) with an according change in physical dimension.

[^1]: While in a 3D simulation, a 2D source term lives on an entity of codimension 1, it is simply a subdomain with the same dimension in 2D simulations.
It is thus subjected to the same modelling assumption allowing a 2D representation and does not change physical dimension.
This can be read as $\mathrm{kg m}^{-2}$ $\mathrm{s}^{-1}$ $\mathrm{m}^{-1}$.
In other words, the volumetric source term is represented on the 2D domain per running meter.
In a 2D simulation, the 2D source term represents a 3D source condensed in the out-of-plane direction.

[^2]: In a 2D simulation, a line-source actually represents an area source, but condensed in the out-of-plane direction.
It could thus also be read as $\mathrm{kg}$ $\mathrm{m}^{-1}$ $\mathrm{s^{-1}}$ $\mathrm{m}^{-1}$.
