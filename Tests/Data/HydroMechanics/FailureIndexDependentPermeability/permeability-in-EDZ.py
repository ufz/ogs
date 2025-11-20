# %% [markdown]
# +++
# date = "2020-06-11"
# title = "Permeability in EDZ"
# weight = 151
# project = ["HydroMechanics/FailureIndexDependentPermeability/quad_with_half_hole.prj"]
# author = "Wenqing Wang, Dmitri Naumov, Thomas Nagel, Olaf Kolditz, Noor Hasan"
# web_subsection = "hydro-mechanics"
# +++


# %% [markdown]
# ## Theoretical Background
# This implementation report presents the theoretical basis of the failure
# index dependent permeabilitymodel for the excavation damaged zone (EDZ)
# by taking into account HM response in rock due to excavation [2] as well
# as the implementation in OGS-6.
# The type of that permeability model in OGS-6 is
# PermeabilityMohrCoulombFailureIndexModel.
#
# The failure index dependent permeability model according to [2] is
# defined as
#
# $$
# \mathbf{k} = \mathbf{k}_0 + H(f - 1) k_r e^{bf} \mathbf{I} \tag{1}
# $$
#
# where $k_0$ is the intrinsic permeability of the undamaged material,
# $H$ is the Heaviside step function, $f$ is the failure index, $k_r$
# is a reference permeability, $b$ is a fitting parameter.
# $k_r$ and $b$ can be calibrated by experimental data.
#
# The failure index $f$ is calculated from the Mohr-Coulomb failure
# criterion comparing an acting shearstress.
# With the conventional mechanics notations, which mean that tensile
# stress is positive, the Mohr-Coulomb failure criterion [1] takes
# the form
#
# $$
# \tau_f(\sigma) = c - \sigma \tan\phi \tag{2}
# $$
#
# with $\tau$ the shear strength, $c$ the cohesion, $\sigma$ the normal
# stress, and $\phi$ the internal friction angle.
# We further introduce the maximum shear stress
# $\tau_{\mathrm{m}}=\left(\sigma_3-\sigma_1\right) / 2$
# and the mean stress $\sigma_{\mathrm{m}}=\left(\sigma_1+\sigma_3\right) / 2$,
# where $\sigma_1$ and $\sigma_3$ are the minimum and maximum shear
# stress, respectively.
# A second criterion is implemented similar to a tension cut-off.
# Let $\sigma_{\mathrm{m}}^{\max } \in(0, c / \tan \phi)$ be a
# limit value related (but not equivalent) to tensile strength
# of the material.
# Then, the failure index is determined by
#
# $$
# f= \begin{cases}\frac{\left|\tau_{\mathrm{m}}\right|}{\cos (\phi) \tau_{\mathrm{f}}\left(\sigma_{\mathrm{m}}\right)} & \text { if } \sigma_{\mathrm{m}} \leq \sigma_{\mathrm{m}}^{\max } \\ \max \left\{\frac{\left|\tau_{\mathrm{m}}\right|}{\cos (\phi) \tau_{\mathrm{f}}\left(\sigma_{\mathrm{m}}\right)}, \frac{\sigma_{\mathrm{m}}}{\sigma_{\mathrm{m}}^{\max }}\right\} & \text { if } \sigma_{\mathrm{m}}>\sigma_{\mathrm{m}}^{\max }\end{cases} \tag{3}
# $$
#
# The computed permeability components are restricted with an upper
# bound, i.e. $\mathbf{k} := k_{ij} < k_{max}$
#
# The material properties for the test example are given in Table 1.
#
# **Table 1: Material properties**
# | Property | Value | Unit |
# | :--- | :--- | :--- |
# | Fluid |  |  |
# | Density | 1000 | $\mathrm{~kg} / \mathrm{m}^3$ |
# | Fluid viscosity | $10^{-3}$ | Pa s |
# | Solid |  |  |
# | Density | 2650 | $\mathrm{~kg} / \mathrm{m}^3$ |
# | Porous medium |  |  |
# | Porosity | 0.15 | - |
# | Intrinsic permeability | the EDZ model | $\mathrm{m}^2$ |
# | Elasticity |  |  |
# | Young's modulus | $6 \cdot 10^9$ | Pa |
# | Poisson's ratio | 0.3 | - |
# | Biot's coefficient | 0.6 | - |
#
# The parameters of the EDZ permeability with the Mohr Coulumb
# failure index are
#
# $$
# \begin{aligned}
# & \mathbf{k}_0 = \{10^{-20}\} \mathbf{m}^2, \\
# & k_r = 10^{-19} \mathbf{m}^2, b=3.0, c=1 \mathbf{MPa}, \phi=15^\circ, \\
# & k_{\max} = 10^{-6} \mathbf{m}^2, \sigma_{\mathrm{m}}^{\max} = 0.8 \frac{c}{\tan \phi} = 2.985640646055102e6 \mathbf{MPa}
# \end{aligned} \tag{4}
# $$
#
# This geometry of this example is a square of
# [0,50] $\times$ [-25,25] $m²$ with a half circle hole with a radius
# of 2.3 m and a center at (0,0).
#
# The initial pore pressure is 4.7 MPa.
# The initial displacement and stress components are all zero.
#
# **Table 2: Boundary**
# | Boundary | Mass balance equation | Momentum balance equation |
# | :--- | :--- | :--- |
# | Left | No flux | $u_x=0, \tau_y=0$ |
# | Right | No flux | $\tau_x=-15 \mathrm{MPa}, \tau_y=0$ |
# | Bottom | No flux | $u_y=0, \tau_x=0$ |
# | Top | No flux | $\tau_x=0, \tau_y=-12 \mathrm{MPa}$ |
# | Hole surface | $p=0.1 \mathrm{MPa}$ | $\tau_n=0$ |

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import ogstools as ot
from ogstools.meshlib import MeshSeries

# %%
# creating output path if it doesn't exist already.
out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

# %%
model = ot.Project(
    input_file="quad_with_half_hole.prj", output_file="quad_with_half_hole.prj"
)

# %%
model.write_input()

# %%
model.run_model(logfile=Path(out_dir) / "log.txt")

# %% [markdown]
# The computed permeability is shown in Fig. 1, in which one can
# see that the permeability near the hole is increased with a
# reasonable distribution pattern.
# This implies that the permeability model can describe the
# permeability change in EDZ.

# %%
ms = MeshSeries("quad_with_half_hole.pvd").mesh(1)

fig_perm, ax = plt.subplots(figsize=(8, 6))
ms.plot_contourf("permeability", fig=fig_perm, ax=ax, fontsize=15)
ax.set_title("Figure 1: Calculated permeability distribution in [m$^2$].")
plt.tight_layout()
plt.show()

# %% [markdown]
# The distributions of the horizontal stress, $\sigma_{xx}$ , and
# the pore pressure are illustrated in Fig. 2.

# %%
pressure = ot.variables.pressure.replace(data_name="pressure_interpolated")
fig_stress, ax = plt.subplots(nrows=1, ncols=2, figsize=(20, 8))

ms.plot_contourf(ot.variables.stress["xx"], fig=fig_stress, ax=ax[0], fontsize=15)
ms.plot_contourf(pressure, fig=fig_stress, ax=ax[1], fontsize=15)
plt.suptitle(
    "Figure 2: Calculated distributions of horizontal stress (left) and pore pressure (right), in [MPa] respectively",
    fontsize=20,
)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## References
# [1] J.F. Labuz and A. Zang. Mohr–Coulomb failure criterion.
# _Rock Mechanics and Rock Engineering_,45(6):975–979, 2012.
#
# [2] W.Q. Wang, H. Shao, Th. Nagel, and O. Kolditz. Analysis
# of coupled thermal-hydro-mechanical processes during small
# scale in-situ heater experiment in Callovo-Oxfordian clay
# rock introducing a failure-index permeability model. _International
# Journal of Rock Mechanics and Mining Sciences_,
# revised manuscript under review, 2020.
# %%
