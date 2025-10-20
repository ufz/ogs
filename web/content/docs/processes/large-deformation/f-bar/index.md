+++
author = "Wenqing Wang, Thomas Nagel"
date = "2025-10-13"
title = "Large deformation - F bar method"
weight = 5
+++


According to the references $[1,2]$, a modified deformation gradient, $\overline{\mathbf{F}}$, is introduced to compute stresses in order to alleviate the spurious locking exhibited by the standard bi-linear and tri-linear elements near the incompressible limit.
<!-- vale off -->
The deformation gradient $\mathbf{F}$ can be expressed as a composition of dilatational change and deviatoric change as
<!-- vale on -->
$$
\mathbf{F}=\mathbf{F}_{d} \mathbf{F}_{v}
$$

With $\left(\mathbf{F}_{0}\right)_{v}=\operatorname{det}\left[\mathbf{F}_{0}\right]^{1 / n} \mathbf{I}$ and $\mathbf{F}_{d}=\operatorname{det}[\mathbf{F}]^{-1 / n} \mathbf{F}_{d}, \overline{\mathbf{F}}$ is defined as the composition of the deviatoric component of $\mathbf{F}$ with the volumetric component of $\mathbf{F}_{0}$ as

$$
\overline{\mathbf{F}}:=\mathbf{F}_{d}\left(\mathbf{F}_{0}\right)_{v}=\left(\frac{\operatorname{det}\left[\mathbf{F}_{0}\right]}{\operatorname{det}[\mathbf{F}]}\right)^{\frac{1}{n}} \mathbf{F}
$$

with $\mathbf{F}_{0}$ the value of $\mathbf{F}$ at the element center, $n$ the space dimension.
Alternatively, $\mathbf{F}_{0}$ can be computed as a average value as

$$
\mathbf{F}_{0}=\frac{\int_{\Omega_{e}} \mathbf{F} \mathrm{~d} \Omega}{\int_{\Omega_{e}} \mathrm{~d} \Omega}
$$

Hereafter, we denote $\left(\frac{\operatorname{det}\left[\mathbf{F}_{0}\right]}{\operatorname{det}[\mathbf{F}]}\right)^{\frac{1}{n}}$ as $\alpha$.

## 1 Equilibrium equations

We assume that $\left\{\Omega^{t}: \mathbf{x} \in \mathbb{R}^{n}\right\}$ is the current deformed configuration, $\left\{\Omega: \mathbf{X} \in \mathbb{R}^{n}\right\}$ is the reference configuration, and $\phi(\mathbf{X})$ is the coordinate mapping $\left\{\phi(\mathbf{X}): \Omega \rightarrow \Omega^{t}\right\}$ such as $\mathbf{x}=\phi(\mathbf{X})$. The displacement and its gradient can be written as

$$
\begin{align*}
& \mathbf{u}=\mathbf{x}-\mathbf{X}=\phi(\mathbf{X})-\mathbf{X}  \tag{1}\\
& \mathbf{F}=\nabla \phi(\mathbf{X}) \tag{2}
\end{align*}
$$

Let $\mathcal{S}$ is the space of admissible deformations defined by

$$
\begin{equation*}
\mathcal{S}=\phi: \Omega \rightarrow \mathbb{R}^{n}\left|\operatorname{det}\left(\nabla_{X} \phi\right)>0, \phi\right|_{\partial \Omega_{\phi}}=\phi_{b} \tag{3}
\end{equation*}
$$

and $\mathcal{V}_{\phi}$ is the tangent space to $\mathcal{S}$ at $\phi$ as

$$
\begin{equation*}
\mathcal{V}_{\phi}=d \mathbf{w} \circ \phi: \Omega \rightarrow \mathbb{R}^{n}\left|\operatorname{det}\left(\nabla_{X} \phi\right)>0, d \mathbf{w} \circ \phi\right|_{\partial \Omega_{\phi}}=0 \tag{4}
\end{equation*}
$$

In the total Lagrangian formulation, the equilibrium equations are derived from the principle of virtual work in the reference configuration $\Omega$. This leads to find $\phi \in \mathcal{S}$ such that $\forall d \mathbf{w} \in \mathcal{V}_{\phi}$ it satisfies

$$
\begin{equation*}
\int_{\Omega_{e}} \mathbf{S}: d \overline{\mathbf{E}}(d \mathbf{w}) \mathrm{d} \Omega=\int_{\Omega_{e}} \mathbf{f} \cdot d \mathbf{w} \mathrm{~d} \Omega+\int_{\left.\partial \Omega\right|_{\tau}} d \mathbf{w} \cdot \tau \mathrm{~d} \Gamma \tag{5}
\end{equation*}
$$

where $\mathbf{S}$ is the second Piola-Kirchhoff stress tensor, $\overline{\mathbf{E}}$ is the Green-Lagrange strain calculated with $\overline{\mathbf{F}}$ and $\mathbf{f}$ is the external load vector, The Green-Lagrange strain is defined by

$$
\begin{equation*}
\mathbf{E}=\frac{1}{2}\left(\mathbf{F}^{\mathrm{T}} \mathbf{F}-\mathbf{I}\right), \tag{6}
\end{equation*}
$$

Using the F bar method, the modified Green-Lagrange strain $\overline{\mathbf{E}}$ is

$$
\begin{align*}
\overline{\mathbf{E}} & =\frac{1}{2}\left(\overline{\mathbf{F}}^{\mathrm{T}} \overline{\mathbf{F}}-\mathbf{I}\right),  \tag{7}\\
& =\frac{1}{2}\left(\alpha^{2} \mathbf{F}^{\mathrm{T}} \mathbf{F}-\alpha^{2} \mathbf{I}+\alpha^{2} \mathbf{I}-\mathbf{I}\right),  \tag{8}\\
& =\alpha^{2} \mathbf{E}+\frac{1}{2}\left(\alpha^{2}-1\right) \mathbf{I} . \tag{9}
\end{align*}
$$

Consequently, the stress is computed as

$$
\begin{equation*}
\mathbf{S}:=\mathbf{S}(\overline{\mathbf{E}}) . \tag{10}
\end{equation*}
$$

## 2 Linearization

### 2.1 Basic derivatives

To linearize the equilibrium equations, we first derive some fundamental derivatives.

### 2.1.1 Directional derivatives

The directional derivative of a multivariable differentiable (scalar) function along a given vector $\mathbf{v}$ at a given point $\mathbf{x}$ intuitively represents the instantaneous rate of change of the function, moving through $\mathbf{x}$ with a velocity specified by $\mathbf{v}$ : The directional derivative of a scalar function $f(\mathbf{x}), \mathbf{x} \in \mathbb{R}^{n}$ along a vector $\mathbf{v} \in \mathbb{R}^{n}$ is the function $\nabla_{\mathbf{v}} f(\mathbf{x})$ defined by the limit:

$$
\begin{equation*}
\nabla_{\mathbf{v}} f(\mathbf{x})=\lim _{h \rightarrow 0} \frac{f(\mathbf{x}+h \mathbf{v})-f(\mathbf{x})}{h}=\left.\frac{\partial}{\partial h} f(\mathbf{x}+h \mathbf{v})\right|_{\lim _{h \rightarrow 0}} \tag{11}
\end{equation*}
$$

If $f(\mathbf{x})$ is differentiable at $\mathbf{x}$, the following equation holds after applying the first order Taylor approximation to $f(\mathbf{x}+h \mathbf{v})$ in the above definition

$$
\begin{equation*}
\nabla_{\mathbf{v}} f(\mathbf{x})=\nabla f(\mathbf{x}) \mathbf{v}=\left.\frac{\partial}{\partial h} f(\mathbf{x}+h \mathbf{v})\right|_{\lim _{h \rightarrow 0}} \tag{12}
\end{equation*}
$$

We will use the definition of directional derivatives to simplify the linearization.

### 2.1.2 Virtual deformation gradient $d \mathbf{F}$

$$
\begin{equation*}
d \mathbf{F}=\nabla_{\phi} \mathbf{F} d \mathbf{w}=\left.\frac{\partial}{\partial h} \mathbf{F}(\phi+h d \mathbf{w})\right|_{\lim _{h \rightarrow 0}}=\nabla d \mathbf{w} \tag{13}
\end{equation*}
$$

### 2.1.3 Virtual strain $d \mathbf{E}$

$$
\begin{equation*}
d \mathbf{E}=\nabla_{\phi} \mathbf{E} d \mathbf{w}=\frac{1}{2}\left((\nabla d \mathbf{w})^{\mathrm{T}} \mathbf{F}+\mathbf{F}^{\mathrm{T}} \nabla d \mathbf{w}\right) \tag{14}
\end{equation*}
$$

### 2.1.4 Virtual strain $d \overline{\mathbf{E}}$

Therefore, the variation of the modified Green-Lagrange strain gives

$$
\begin{equation*}
d \overline{\mathbf{E}}=\alpha^{2} d \mathbf{E}+\alpha(2 \mathbf{E}+\mathbf{I}) d \alpha \tag{15}
\end{equation*}
$$

Note that the derivative of the determinant of a matrix with respect to the matrix itself is used to obtain the above derivatives, which is

$$
\begin{equation*}
\frac{\partial}{\partial \mathbf{A}}(\operatorname{det}(\mathbf{A}))=\operatorname{det}(\mathbf{A}) \mathbf{A}^{-\mathrm{T}} \tag{16}
\end{equation*}
$$

This gives

$$
\begin{align*}
d \alpha & =\frac{\alpha}{\mathrm{n}}\left(\mathbf{F}_{0}^{-\mathrm{T}}: d \mathbf{F}_{0}-\mathbf{F}^{-\mathrm{T}}: d \mathbf{F}\right)  \tag{17}\\
& =\frac{\alpha}{\mathrm{n}}\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right) \tag{18}
\end{align*}
$$

Consequently,

$$
\begin{align*}
d \overline{\mathbf{E}} & =\alpha^{2} d \mathbf{E}+\alpha(2 \mathbf{E}+\mathbf{I})\left(\frac{\alpha}{\mathrm{n}}\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right)\right)  \tag{20}\\
& =\alpha^{2}\left(d \mathbf{E}+\frac{1}{\mathrm{n}}(2 \mathbf{E}+\mathbf{I})\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right)\right) \tag{21}
\end{align*}
$$

### 2.1.5 Jacobian

Assume that the body force $\mathbf{f}$ and the traction $\tau$ are independent of the displacement, the Jacobian for the Newton-Raphson method can be obtained by deriving the variation of the virtual strain energy as

$$
\begin{equation*}
d_{u} \int_{\Omega_{e}} \mathbf{S}(\overline{\mathbf{E}}): d \overline{\mathbf{E}} \mathrm{~d} \Omega=\int_{\Omega_{e}} \nabla_{u}(\mathbf{S}(\overline{\mathbf{E}}): d \overline{\mathbf{E}}) \delta \mathbf{u} \mathrm{d} \Omega \tag{22}
\end{equation*}
$$

with $\delta \mathbf{u} \in \mathcal{V}_{\phi}$
According to the directional derivative rule,

$$
\begin{align*}
\nabla_{u}(\mathbf{S}(\overline{\mathbf{E}}): d \overline{\mathbf{E}}) \delta \mathbf{u}= & \left.\frac{\partial}{\partial h}(\mathbf{S}(\overline{\mathbf{E}}(\phi+h \delta \mathbf{u})): d \overline{\mathbf{E}}(\phi+h \delta \mathbf{u}))\right|_{\lim _{h \rightarrow 0}}  \tag{23}\\
= & \frac{\partial \mathbf{S}(\overline{\mathbf{E}}(\phi+h \delta \mathbf{u})}{\partial \overline{\mathbf{E}}(\phi+h \delta \mathbf{u})}: \frac{\partial \overline{\mathbf{E}}(\phi+h \delta \mathbf{u})}{\partial h}:\left.d \overline{\mathbf{E}}(\phi+h \delta \mathbf{u})\right|_{\lim _{h \rightarrow 0}}  \tag{24}\\
& +\mathbf{S}(\overline{\mathbf{E}}(\phi+h \delta \mathbf{u})):\left.\frac{\partial}{\partial h}(d \overline{\mathbf{E}}(\phi+h \delta \mathbf{u}))\right|_{\lim _{h \rightarrow 0}} \tag{25}
\end{align*}
$$

where $\partial \mathbf{S}\left(\overline{\mathbf{E}}(\phi+h \delta \mathbf{u}) /\left.\partial \overline{\mathbf{E}}(\phi+h \delta \mathbf{u})\right|_{\lim _{h \rightarrow 0}}=\partial \mathbf{S}(\overline{\mathbf{E}}(\phi) / \partial \overline{\mathbf{E}}(\phi)\right.$ is the material tangential, which is a forth order tensor, hereafter we denote it as $\mathbf{C}$.

Since

$$
\begin{equation*}
\left.\frac{\partial \overline{\mathbf{E}}(\phi+h \delta \mathbf{u})}{\partial h}\right|_{\lim _{h \rightarrow 0}}=\nabla_{\phi} \overline{\mathbf{E}}(\phi) \delta \mathbf{u}=\delta \overline{\mathbf{E}}(\phi), \tag{26}
\end{equation*}
$$

we have

$$
\begin{equation*}
\delta \overline{\mathbf{E}}=\alpha^{2}\left(\delta \mathbf{E}+\frac{1}{\mathrm{n}}(2 \mathbf{E}+\mathbf{I})\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla \delta \mathbf{u}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}\right)\right) \tag{27}
\end{equation*}
$$

by the same way deriving $d \overline{\mathbf{E}}$.
Expanding $\partial\left(d \overline{\mathbf{E}}(\phi+h \delta \mathbf{u}) /\left.\partial h\right|_{\lim _{h \rightarrow 0}}\right.$ leads to

$$
\begin{align*}
\left.\frac{\partial}{\partial h}(d \overline{\mathbf{E}}(\phi+h \delta \mathbf{u}))\right|_{\lim _{h \rightarrow 0}=} & \left.\frac{\partial}{\partial h}\left(\alpha^{2}\left(d \mathbf{E}+\frac{1}{\mathrm{n}}(2 \mathbf{E}+\mathbf{I})\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right)\right)\right)\right|_{\lim _{h \rightarrow 0}}  \tag{28}\\
= & \left.2 \alpha\left(d \mathbf{E}+\frac{1}{\mathrm{n}}(2 \mathbf{E}+\mathbf{I})\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right)\right) \frac{\partial \alpha}{\partial h}\right|_{\lim _{h \rightarrow 0}}  \tag{29}\\
& +\left.\alpha^{2} \frac{\partial d \mathbf{E}}{\partial h}\right|_{\lim _{h \rightarrow 0}}+\left.\frac{2 \alpha^{2}}{\mathrm{n}} \frac{\partial \mathbf{E}}{\partial h}\right|_{\lim _{h \rightarrow 0}}\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right)  \tag{30}\\
& +\frac{\alpha^{2}}{\mathrm{n}}(2 \mathbf{E}+\mathbf{I})\left(\left.\frac{\partial}{\partial h}\left(\mathbf{F}_{0}^{-\mathrm{T}}\right)\right|_{\lim _{h \rightarrow 0}}: \nabla d \mathbf{w}_{0}-\left.\frac{\partial}{\partial h}\left(\mathbf{F}^{-\mathrm{T}}\right)\right|_{\lim _{h \rightarrow 0}}: \nabla d \mathbf{w}\right) \tag{31}
\end{align*}
$$

where

$$
\begin{gather*}
\left.\frac{\partial \alpha}{\partial h}\right|_{\lim _{h \rightarrow 0}}=\left.\frac{\partial \alpha(\phi+h \delta \mathbf{u})}{\partial h}\right|_{\lim _{h \rightarrow 0}}=\nabla_{\mathbf{u}} \alpha \delta \mathbf{u}=\delta \alpha=\frac{\alpha}{\mathrm{n}}\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla \delta \mathbf{u}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}\right),  \tag{32}\\
\left.\frac{\partial d \mathbf{E}}{\partial h}\right|_{\lim _{h \rightarrow 0}}=\left.\frac{1}{2} \frac{\partial}{\partial h}\left((\nabla d \mathbf{w})^{\mathrm{T}} \mathbf{F}+\mathbf{F}^{\mathrm{T}} \nabla d \mathbf{w}\right)\right|_{\lim _{h \rightarrow 0}}  \tag{33}\\
=\frac{1}{2}\left((\nabla d \mathbf{w})^{\mathrm{T}} \nabla \delta \mathbf{u}+(\nabla \delta \mathbf{u})^{t r} \nabla d \mathbf{w}\right),  \tag{34}\\
\left.\frac{\partial \mathbf{E}}{\partial h}\right|_{\lim _{h \rightarrow 0}}=\left.\frac{\partial \mathbf{E}(\phi+h \delta \mathbf{u})}{\partial h}\right|_{\lim _{h \rightarrow 0}}=\frac{1}{2}\left((\nabla \delta \mathbf{u})^{\mathrm{T}} \mathbf{F}+\mathbf{F}^{\mathrm{T}} \nabla \delta \mathbf{u}\right),  \tag{35}\\
\left.\frac{\partial \mathbf{F}_{0}^{-\mathrm{T}}}{\partial h}\right|_{\lim _{h \rightarrow 0}}=\left.\frac{\partial F_{0}^{-\mathrm{T}}}{\partial \mathbf{F}_{0}} \frac{\partial \mathbf{F}_{0}(\phi+h \delta \mathbf{u})}{\partial h}\right|_{\lim _{h \rightarrow 0}}=-\mathbf{F}_{0}^{-\mathrm{T}} \otimes \mathbf{F}_{0}^{-\mathrm{T}}: \nabla \delta \mathbf{u}_{0},  \tag{36}\\
\left.\frac{\partial \mathbf{F}^{-\mathrm{T}}}{\partial h}\right|_{\lim _{h \rightarrow 0}}=\left.\frac{\partial F^{-\mathrm{T}}}{\partial \mathbf{F}} \frac{\partial \mathbf{F}(\phi+h \delta \mathbf{u})}{\partial h}\right|_{\lim _{h \rightarrow 0}}=-\mathbf{F}^{-\mathrm{T}} \otimes \mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u} . \tag{37}
\end{gather*}
$$

We have

$$
\begin{align*}
\left.\frac{\partial}{\partial h}(d \overline{\mathbf{E}}(\phi+h \delta \mathbf{u}))\right|_{\lim _{h \rightarrow 0}}= & \delta(d \overline{\mathbf{E}})  \tag{38}\\
= & \frac{2 \alpha^{2}}{\mathrm{n}}\left(d \mathbf{E}+\frac{1}{\mathrm{n}}(2 \mathbf{E}+\mathbf{I})\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right)\right)  \tag{39}\\
& \quad\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla \delta \mathbf{u}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}\right)  \tag{40}\\
& +\left.\alpha^{2} \frac{\partial d \mathbf{E}}{\partial h}\right|_{\lim _{h \rightarrow 0}}  \tag{41}\\
& +\frac{\alpha^{2}}{\mathrm{n}}\left((\nabla \delta \mathbf{u})^{\mathrm{T}} \mathbf{F}+\mathbf{F}^{\mathrm{T}} \nabla \delta \mathbf{u}\right)\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right)  \tag{42}\\
& \quad-\frac{\alpha^{2}}{\mathrm{n}}(2 \mathbf{E}+\mathbf{I})\left(\mathbf{F}_{0}^{-\mathrm{T}} \otimes \mathbf{F}_{0}^{-\mathrm{T}}: \nabla \delta \mathbf{u}_{0}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}} \otimes \mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}: \nabla d \mathbf{w}\right) \tag{43}
\end{align*}
$$

Finally, we obtain the expression the variation of the virtual strain energy as

$$
\begin{align*}
\int_{\Omega_{e}} \nabla_{u}(\mathbf{S}(\overline{\mathbf{E}}): d \overline{\mathbf{E}}) \delta \mathbf{u} \mathrm{d} \Omega= & \int_{\Omega_{e}} \mathbf{C}: \delta \overline{\mathbf{E}}: d \overline{\mathbf{E}} \mathrm{~d} \Omega  \tag{44}\\
& +\int_{\Omega_{e}} \mathbf{S}(\overline{\mathbf{E}}): \delta(d \overline{\mathbf{E}}(\phi+h \delta \mathbf{u})) \mathrm{d} \Omega \tag{45}
\end{align*}
$$

## 3 Finite element

After the discretization with the Galerkin approach, we have

$$
\begin{equation*}
\delta \mathbf{u}=\mathbf{N} \delta \hat{\mathbf{u}}, \quad d \mathbf{u}=\mathbf{N} \hat{\mathbf{u}} \tag{46}
\end{equation*}
$$

with $\mathbf{N}$ the shape functions, $\hat{\delta \mathbf{u}}$ and $\hat{d \mathbf{u}}$ the arbitrary virtual nodal displacements. This gives

$$
\begin{equation*}
d \mathbf{E}=\frac{\partial \mathbf{E}}{\partial \hat{\mathbf{u}}}: \hat{d \mathbf{u}}, \quad \delta \mathbf{E}=\frac{\partial \mathbf{E}}{\partial \hat{\mathbf{u}}}: \hat{\delta \mathbf{u}} \tag{47}
\end{equation*}
$$

Assuming that stress tensor and strain tensor are symmetry, and considering the matrix form of $\partial \mathbf{E} / \partial \hat{\mathbf{u}}$ gives

$$
\begin{equation*}
\frac{\partial \mathbf{E}}{\partial \hat{\mathbf{u}}}:=\mathrm{B} . \tag{48}
\end{equation*}
$$

where $\mathrm{B}=\mathrm{B}_{l}+\mathrm{B}_{n l}$ with $\mathrm{B}_{n}$ and $\mathrm{B}_{n l}$ the linear and non-linear parts of the B matrix, respectively.

### 3.1 Modified Jacobian

Additional B matrix from $d \overline{\mathbf{E}}$ and $\delta \overline{\mathbf{E}}$
As for $d \overline{\mathbf{E}}$, let's denote $\alpha^{2}(2 \mathbf{E}+\mathbf{I})\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d w_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d w\right) / \mathrm{n}$ as $d \bar{\varepsilon}$, which gives

$$
\begin{equation*}
d \overline{\mathbf{E}}=\alpha^{2} d \mathbf{E}+d \bar{\varepsilon} \tag{49}
\end{equation*}
$$

Expanding $\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right)$ gives:

$$
\begin{equation*}
\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right)=\left(\mathbf{F}_{0}^{-T}\right)_{i j} d \hat{w}_{i}^{L} N_{, j}^{L}\left(\xi_{0}\right)-\left(\mathbf{F}^{-T}\right)_{i j} d \hat{w}_{i}^{L} N_{, j}^{L}=\left(\mathbf{q}_{0}-\mathbf{q}\right) d \hat{\mathbf{w}} \tag{50}
\end{equation*}
$$

where $\mathbf{q}$ is a $1 \times \mathrm{n} \cdot$ NE matrix or a transposed vector as

$$
\begin{gather*}
\mathbf{q}=\left(\mathbf{q}_{\mathrm{col} 1}, \mathbf{q}_{\mathrm{col} 2}\right), \text { or } \mathbf{q}=\left(\mathbf{q}_{\mathrm{col} 1}, \mathbf{q}_{\mathrm{col} 2}, \mathbf{q}_{\mathrm{col} 3}\right)  \tag{51}\\
\mathbf{q}_{\mathrm{col} \mathrm{i}}=\left(\left(\mathbf{F}^{-1}\right)_{j i} N_{, j}^{1}, \quad\left(\mathbf{F}^{-1}\right)_{j i} N_{, j}^{2}, \quad \cdots,\left(\mathbf{F}^{-1}\right)_{j i} N_{, j}^{\mathrm{NE}}\right), \quad i, j=1, \cdots, \operatorname{dim} . \tag{52}
\end{gather*}
$$

with NE the number of nodes, and $\mathbf{q}_{0}$ the value of $\mathbf{q}$ at the element center $\xi_{0}$.
This can be written into a matrix form as

$$
\begin{equation*}
[d \bar{\varepsilon}]=\hat{\mathrm{B}} d \hat{\mathbf{w}} \tag{53}
\end{equation*}
$$

with $[d \bar{\varepsilon}]$ the vector form of $d \bar{\varepsilon}$, and

$$
\begin{equation*}
\hat{\mathrm{B}}=\frac{\alpha^{2}}{\mathrm{n}}[2 \mathbf{E}+\mathbf{I}]\left(\mathbf{q}_{0}-\mathbf{q}\right) \tag{54}
\end{equation*}
$$

where $[2 \mathbf{E}+\mathbf{I}]=\left(2 E_{11}+1,2 E_{22}+1,2 E_{33}+1,2 E_{12}\right)^{\mathrm{T}}$ for plane strain problems, and $[2 \mathbf{E}+\mathbf{I}]=\left(2 E_{11}+1,2 E_{22}+1,2 E_{33}+1,2 E_{12}, 2 E_{23}, 2 E_{13}\right)^{\mathrm{T}}$ for 3D problems, respectively.

Note: the shear strain $E_{12}, E_{23}, E_{13}$ are assumed being scaled with $\sqrt{2}$ for the computation with the Kelvin vector.

The same for $\delta \overline{\mathbf{E}}$, we have $\delta \bar{\varepsilon}=\hat{\mathrm{B}} \delta \hat{\mathbf{u}}$. Since $d \overline{\mathbf{E}}=\alpha^{2} d \mathbf{E}+d \bar{\varepsilon}$,

$$
\begin{equation*}
[d \overline{\mathbf{E}}]=\left(\alpha^{2} \mathrm{~B}+\hat{\mathrm{B}}\right) d \hat{\mathbf{w}} \tag{55}
\end{equation*}
$$

where $[d \overline{\mathbf{E}}]$ means the vector form of $d \overline{\mathbf{E}}$.
We denote $\alpha^{2} \mathrm{~B}+\hat{\mathrm{B}}$ as $\overline{\mathrm{B}}$, which simplifies the expression of the Jacobian from $\int_{\Omega_{e}} \mathbf{C}$ : $\delta \overline{\mathbf{E}}: d \overline{\mathbf{E}} \mathrm{~d} \Omega$ as

$$
\begin{equation*}
\int_{\Omega_{e}} \mathbf{C}:(\delta \overline{\mathbf{E}}): d \overline{\mathbf{E}} \mathrm{~d} \Omega \xrightarrow{\text { matrix-vector form }} \int_{\Omega_{e}} \overline{\mathrm{~B}}^{\mathrm{T}}[\mathbf{C}] \overline{\mathrm{B}} \mathrm{~d} \Omega \tag{56}
\end{equation*}
$$

where $[\mathbf{C}]$ is matrix from of $\mathbf{C}$.

### 3.2 Additional contributions to Jacobian from $\int_{\Omega_{e}} \delta(d \overline{\mathbf{E}}) \mathrm{d} \Omega$

### 3.2.1 Term $\int_{\Omega_{e}} \alpha^{2} \mathbf{S}: \delta(d \mathbf{E}) \mathrm{d} \Omega=\int_{\Omega_{e}} \mathbf{S}: \alpha^{2} \partial d \mathbf{E} /\left.\partial h \mathrm{~d} \Omega\right|_{\lim _{h \rightarrow 0}}$

From this term, we obtain the standard G matrix related Jacobian contribution as

$$
\begin{equation*}
\int_{\Omega_{e}} \alpha^{2} \mathbf{G}^{\mathrm{T}}[[\mathbf{S}]] \mathbf{G} \mathrm{d} \Omega \tag{57}
\end{equation*}
$$

with $[[\mathbf{S}]]$ for a matrix with stress matrix as diagonal blocks.
3.2.2 Term with $\frac{2 \alpha^{2}}{\mathbf{n}}\left(d \mathbf{E}+\frac{1}{\mathbf{n}}(2 \mathbf{E}+\mathbf{I})\left(\mathbf{F}_{0}^{-\mathbf{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathbf{T}}: \nabla d \mathbf{w}\right)\right)$

$$
\left(\mathbf{F}_{0}^{-\mathbf{T}}: \nabla \delta \mathbf{u}_{0}-\mathbf{F}^{-\mathbf{T}}: \nabla \delta \mathbf{u}\right)
$$

The corresponding term in the linearized weak form is

$$
\begin{array}{r}
\int_{\Omega_{e}} \mathbf{S}: \frac{2 \alpha^{2}}{\mathrm{n}}\left(d \mathbf{E}+\frac{1}{\mathrm{n}}(2 \mathbf{E}+\mathbf{I})\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right)\right) \\
\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla \delta \mathbf{u}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}\right) \mathrm{d} \Omega \tag{59}
\end{array}
$$

which can be written:

$$
\begin{equation*}
\frac{2}{\mathrm{n}} \int_{\Omega_{e}} \mathbf{S}: d \overline{\mathbf{E}}\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla \delta \mathbf{u}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}\right) \mathrm{d} \Omega \tag{60}
\end{equation*}
$$

Note that

$$
\begin{equation*}
\mathbf{S}: d \overline{\mathbf{E}}=(\overline{\mathrm{B}} d \hat{\mathbf{w}})^{\mathrm{T}}[\mathbf{S}] \tag{61}
\end{equation*}
$$

with $[\mathbf{S}]$ the stress in vector type, e.g. the stress in the Kelvin vector.
While

$$
\begin{equation*}
\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla \delta \mathbf{u}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}\right)=\left(\mathbf{q}_{0}-\mathbf{q}\right) \hat{\delta \mathbf{u}} \tag{62}
\end{equation*}
$$

Therefore the additional Jacobian obtained from this term is

$$
\begin{equation*}
\frac{2}{\mathrm{n}} \int_{\Omega_{e}}(\overline{\mathrm{~B}})^{\mathrm{T}}[\mathbf{S}]\left(\mathbf{q}_{0}-\mathbf{q}\right) \mathrm{d} \Omega \tag{63}
\end{equation*}
$$

### 3.2.3 Term with $\left.\frac{2 \alpha^{2}}{\mathbf{n}} \frac{\partial \mathbf{E}}{\partial h}\right|_{\lim _{h \rightarrow 0}}\left(\mathbf{F}_{0}^{-\mathbf{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathbf{T}}: \nabla d \mathbf{w}\right)$

Note that $\left.\frac{\partial \mathbf{E}}{\partial h}\right|_{\lim _{h \rightarrow 0}}=\delta \mathbf{E}$ in that term, the term corresponding integration term is

$$
\begin{equation*}
\int_{\Omega_{e}} \frac{2 \alpha^{2}}{\mathrm{n}} \mathbf{S}: \delta \mathbf{E}\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right) \mathrm{d} \Omega \tag{64}
\end{equation*}
$$

We see that

$$
\begin{equation*}
\mathbf{S}: \delta \mathbf{E}=[\mathbf{S}]^{t r} \mathrm{~B} \hat{\delta} \hat{\mathbf{u}} \tag{65}
\end{equation*}
$$

with $[\mathbf{S}]$ the stress in vector type.
The same for $\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla \delta \mathbf{u}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}\right)$, the discretized of $\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right)$ takes the form

$$
\begin{equation*}
\left(\mathbf{F}_{0}^{-\mathrm{T}}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}}: \nabla d \mathbf{w}\right):=\left(\mathbf{q}_{0}-\mathbf{q}\right) \hat{d \mathbf{w}}=\hat{d \mathbf{w}^{t r}}\left(\mathbf{q}_{0}^{t r}-\mathbf{q}^{t r}\right) \tag{66}
\end{equation*}
$$

Therefore the integration can be written as

$$
\begin{equation*}
\frac{2}{\mathrm{n}} \int_{\Omega_{e}} \hat{d \mathbf{w}} \hat{\alpha}^{t r}\left[\mathbf{q}_{0}^{t r}-\mathbf{q}^{t r}\right](\mathbf{S})^{t r} \mathrm{~B} \hat{\delta} \hat{\mathbf{u}} \mathrm{~d} \Omega \tag{67}
\end{equation*}
$$

This Jacobian contribution from this integration is

$$
\begin{equation*}
\frac{2}{\mathrm{n}} \int_{\Omega_{e}} \alpha^{2}\left[\mathbf{q}_{0}^{t r}-\mathbf{q}^{t r}\right](\mathbf{S})^{t r} \operatorname{Bd} \Omega \tag{68}
\end{equation*}
$$

### 3.2.4 Term with $-\frac{\alpha^{2}}{\mathbf{n}}(2 \mathbf{E}+\mathbf{I})\left(\mathbf{F}_{0}^{-\mathbf{T}} \otimes \mathbf{F}_{0}^{-\mathbf{T}}: \nabla \delta \mathbf{u}_{0}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathbf{T}} \otimes \mathbf{F}^{-\mathbf{T}}: \nabla \delta \mathbf{u}: \nabla d \mathbf{w}\right)$

The corresponding integration is

$$
\begin{equation*}
-\frac{1}{\mathrm{n}} \int_{\Omega_{e}} \alpha^{2} \mathbf{S}:(2 \mathbf{E}+\mathbf{I})\left(\mathbf{F}_{0}^{-\mathrm{T}} \otimes \mathbf{F}_{0}^{-\mathrm{T}}: \nabla \delta \mathbf{u}_{0}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}} \otimes \mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}: \nabla d \mathbf{w}\right) \mathrm{d} \Omega \tag{69}
\end{equation*}
$$

Note that $\mathbf{F}^{-\mathrm{T}} \otimes \mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}: \nabla d \mathbf{w}=\left(\mathbf{F}^{-\mathrm{T}}: d \mathbf{w}\right)\left(\mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}\right)$

From the above description, we know that $\mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}:=\mathbf{q} \hat{\delta \mathbf{u}}$ and $\mathbf{F}_{0}^{-\mathrm{T}}: \nabla \delta \mathbf{u}:=\mathbf{q}_{0} \hat{\delta \mathbf{u}}$ after discretization. Therefore, the integration can be written as

$$
\begin{align*}
& -\frac{1}{\mathrm{n}} \int_{\Omega_{e}} \alpha^{2} \mathbf{S}:(2 \mathbf{E}+\mathbf{I})\left(\mathbf{F}_{0}^{-\mathrm{T}} \otimes \mathbf{F}_{0}^{-\mathrm{T}}: \nabla \delta \mathbf{u}_{0}: \nabla d \mathbf{w}_{0}-\mathbf{F}^{-\mathrm{T}} \otimes \mathbf{F}^{-\mathrm{T}}: \nabla \delta \mathbf{u}: \nabla d \mathbf{w}\right) \mathrm{d} \Omega  \tag{70}\\
& \quad=-\frac{1}{\mathrm{n}} \int_{\Omega_{e}} \alpha^{2} \mathbf{S}:(2 \mathbf{E}+\mathbf{I})\left(\left(\mathbf{q}_{0} d \mathbf{w}\right)^{\mathrm{T}} \mathbf{q}_{0} \delta \mathbf{u}-(\mathbf{q} d \mathbf{w})^{\mathrm{T}} \mathbf{q} \delta \mathbf{u}\right) \mathrm{d} \Omega  \tag{71}\\
& \left.\quad=-\frac{1}{\mathrm{n}} \int_{\Omega_{e}} \alpha^{2} \mathbf{S}:(2 \mathbf{E}+\mathbf{I})(d \mathbf{w})^{\mathrm{T}} \mathbf{q}_{0}^{t r} \mathbf{q}_{0} \delta \mathbf{u}-d \mathbf{w}^{\mathrm{T}} \mathbf{q}^{t r} \mathbf{q} \delta \mathbf{u}\right) \mathrm{d} \Omega \tag{72}
\end{align*}
$$

Therefore the Jacobian contribution from this term is

$$
\begin{equation*}
-\frac{1}{n} \int_{\Omega_{e}} \alpha^{2} \mathbf{S}:(2 \mathbf{E}+\mathbf{I})\left(\mathbf{q}_{0}^{t r} \mathbf{q}_{0}-\mathbf{q}^{t r} \mathbf{q}\right) \mathrm{d} \Omega \tag{73}
\end{equation*}
$$

Note that $\mathbf{S}:(2 \mathbf{E}+\mathbf{I})$ is a scalar, and it can be computed by the dot product of stress vectors as $[\mathbf{S}] \cdot[(2 \mathbf{E}+\mathbf{I})]$.

### 3.3 Jacobian and residual

Finally, we obtain the Jacobian for the total Lagrange formulation with the F bar method:

$$
\begin{align*}
& \int_{\Omega_{e}} \overline{\mathrm{~B}}^{\mathrm{T}}[\mathbf{C}] \overline{\mathrm{B}} \mathrm{~d} \Omega+\int_{\Omega_{e}} \alpha^{2} \mathbf{G}^{\mathrm{T}}[[\mathbf{S}]] \mathbf{G} \mathrm{d} \Omega+\frac{2}{\mathrm{n}} \int_{\Omega_{e}}(\overline{\mathrm{~B}})^{\mathrm{T}}[\mathbf{S}]\left(\mathbf{q}_{0}-\mathbf{q}\right) \mathrm{d} \Omega  \tag{74}\\
& +\frac{2}{\mathrm{n}} \int_{\Omega_{e}} \alpha^{2}\left[\mathbf{q}_{0}^{t r}-\mathbf{q}^{t r}\right](\mathbf{S})^{t r} \mathrm{Bd} \Omega-\frac{1}{\mathrm{n}} \int_{\Omega_{e}} \alpha^{2} \mathbf{S}:(2 \mathbf{E}+\mathbf{I})\left(\mathbf{q}_{0}^{t r} \mathbf{q}_{0}-\mathbf{q}^{t r} \mathbf{q}\right) \mathrm{d} \Omega \tag{75}
\end{align*}
$$

With the equilibrium equation (5), the discretized residual is

$$
\begin{equation*}
\mathbf{R}=\int_{\Omega_{e}} \overline{\mathrm{~B}}[\mathbf{S}] \mathrm{d} \Omega-\int_{\Omega_{e}} \mathbf{f} \mathrm{~N} \mathrm{~d} \Omega-\int_{\left.\partial \Omega\right|_{\tau}} \tau \mathrm{N} \mathrm{~d} \Gamma=\mathbf{0} \tag{76}
\end{equation*}
$$

with N the shape function matrix.

## References
<!-- vale off -->
[1] EA de Souza Neto, D Perić, M Dutko, and DRJ1400785 Owen. Design of simple low order finite elements for large strain analysis of nearly incompressible solids. International Journal of Solids and Structures, 33(20-22):3277-3296, 1996.

[2] Thomas Elguedj, Yuri Bazilevs, Victor M Calo, and Thomas JR Hughes. $\overline{\mathrm{B}}$ and $\overline{\mathrm{F}}$ projection methods for nearly incompressible linear and non-linear elasticity and plasticity using higher-order nurbs elements. Computer methods in applied mechanics and engineering, 197(33-40):2732-2762, 2008.
