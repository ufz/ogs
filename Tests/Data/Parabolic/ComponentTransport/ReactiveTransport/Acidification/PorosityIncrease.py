# ---
# jupyter:
#   jupytext:
#     notebook_metadata_filter: note
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.1
#   kernelspec:
#     display_name: Python (.venv)
#     language: python
#     name: venv
#   note: 'TODO: remove after 6.5.4 release: ci skip'
# ---

# %% [raw]
# +++
# title = "Matrix acidification in a calcite-containing interlayer"
# date = "2025-10-20"
# author = "Mostafa Mollaali"
# image = "figures/schematic.png"
# web_subsection = "reactive-transport"
# weight = 3
# +++

# %%
# %%

from __future__ import annotations

import os
from pathlib import Path
from shutil import copy2

import gmsh
import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pyvista as pv
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from ogstools.meshlib import MeshSeries

ot.plot.setup.show_region_bounds = False

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Matrix acidification in a calcite-containing interlayer
#
# ## Problem overview
#
# This study models how the injection of sulfuric acid alters the porosity and permeability of a layered porous medium. The domain is water-saturated and consists of three distinct layers:
#
# - Upstream quartz sand (inert)
# - Reactive sand-calcite interlayer
# - Downstream quartz sand (inert)
#
# <figure>
#   <img src="./figures/schematic.png" alt="Schematic" style="max-width: 1200px;">
#   <figcaption>Figure: Schematic</figcaption>
# </figure>
#
# The quartz sand is chemically inert, while the key process is the dissolution of calcite in the middle layer by acid, which increases porosity. This porosity change modifies permeability according to the Kozeny-Carman relationship ($
# k = k_0 \frac{(1 - \phi_0)^2}{\phi_0^3} \frac{\phi^3}{(1 - \phi)^2}
# $), thereby influencing subsequent fluid flow and solute transport. Aqueous chemical speciation is computed using the [LLNL thermodynamic database](https://search.r-project.org/CRAN/refmans/phreeqc/html/llnl.dat.html).
#
#
# Boundary conditions, initial conditions, and specific parameter values are provided in the Input Data section. The geochemical reactions are implemented using PHREEQC, tracking five primary components ($\mathrm{S(6)}$, $\mathrm{C}$, $\mathrm{Ca}$, $\mathrm{Na}$, $\mathrm{Cl}$) with calcite dissolution kinetics following a three-mechanism rate law (acid, neutral, and $\mathrm{CO_2}$-promoted pathways).
#
# ## Governing equations: transport and chemical reaction
#
# The model couples fluid flow with reactive solute transport through the advection-dispersion-reaction equation:
#
# $$
# \frac{\partial c_{T\alpha}}{\partial t}
# +\nabla\!\cdot\!\big(\mathbf{q}\,c_{T\alpha}-\phi\,\mathbf{D}\,\nabla c_{T\alpha}\big)
# +R_{\min,\alpha}=0
# $$
#
# $$
# \mathbf{q}=-\frac{k}{\mu}\big(\nabla p-\rho_\ell\,\mathbf{g}\big),\qquad
# \mathbf{D}=D_p\,\mathbf{I}+\alpha_T\,|\mathbf{v}|\,\mathbf{I}
# +(\alpha_L-\alpha_T)\,\frac{\mathbf{v}\mathbf{v}^\top}{|\mathbf{v}|},\qquad
# \mathbf{v}=\frac{\mathbf{q}}{\phi}.
# $$
#
# - $c_{T\alpha}$: total concentration of primary species $\alpha$.
# - $\mathbf{q}$: Darcy flux; body force term $-\rho_\ell\,\mathbf{g}$ included.
# - $\mathbf{D}$: hydrodynamic dispersion tensor with pore diffusion $D_p$ and dispersivities $\alpha_L,\alpha_T$ (using $\mathbf{v}=\mathbf{q}/\phi$).
# - $R_{\min,\alpha}$: heterogeneous reaction term (mineral dissolution/precipitation).
#
# Although Darcy's law is written here with the gravitational body-force term for generality, in this benchmark we neglect gravity by setting the specific body force to zero, so gravity does not contribute to the flow in the simulations.
#
# <!-- ### Chemical components
#
# The model tracks five primary aqueous components:
#
# | Transported Component | Symbol | Primary Role in the System |
# | :--- | :--- | :--- |
# | Sulfate | $\mathrm{S(6)}$ | Tracks the injected acid ($\mathrm{H_2SO_4}$) anion |
# | Total Inorganic Carbon | $\mathrm{C}$ | Represents all carbonate species ($\mathrm{CO_2(aq)}$, $\mathrm{H_2CO_3}$, $\mathrm{HCO_3^-}$, $\mathrm{CO_3^{2-}}$) |
# | Calcium | $\mathrm{Ca}$ | Product of calcite dissolution; controls mineral saturation |
# | Sodium & Chloride | $\mathrm{Na}$, $\mathrm{Cl}$ | Background ions for charge balance and ionic strength |
#
# ## System chemistry
#
# ### Aqueous speciation and reactions
#
# The core chemistry involves calcite dissolution and carbonate equilibrium:
#
# $$
# \begin{aligned}
# &\mathrm{CaCO_3(s) \rightleftharpoons Ca^{2+} + CO_3^{2-}} \quad &\mathrm{(Calcite Dissolution)} \\
# &\mathrm{CO_3^{2-} + H^+ \rightleftharpoons HCO_3^-} \quad &\mathrm{(Carbonate Protonation)} \\
# &\mathrm{HCO_3^- + H^+ \rightleftharpoons CO_2(aq) + H_2O} \quad &\mathrm{(Formation of Dissolved CO$_2$)} \\
# &\mathrm{H^+ + OH^- \rightleftharpoons H_2O} \quad &\mathrm{(Water Dissociation)}
# \end{aligned}
# $$
#
# The carbonate ion released from calcite reacts with protons in acidic conditions, driving dissolution forward.
#
# ### Thermodynamic driving force: saturation state
#
# The tendency for calcite dissolution is determined by the **Saturation Ratio** ($SR$):
#
# $$
# SR = \frac{a(\mathrm{Ca^{2+}}) \cdot a(\mathrm{CO_3^{2-}})}{K_{eq}}
# $$
#
# where $K_{eq} = 10^{-8.48}$ is the calcite equilibrium constant at 25°C, and **activity** ($a$) is the effective concentration accounting for ionic interactions.
#
# - $SR < 1$: Solution is undersaturated → dissolution occurs
# - $SR \approx 1$: System at equilibrium → no net reaction
# - $SR > 1$: Solution is supersaturated → precipitation may occur
#
# ### Kinetic rate of dissolution
#
# The calcite dissolution rate follows a multi-mechanism rate law:
#
# $$
# r_{\mathrm{cal}} = \left( k_a \cdot a(\mathrm{H^+}) + k_n + k_c \cdot \mathrm{SR_{CO_2(g)}} \right) \cdot \left( 1 - SR \right)
# $$
#
# **Mechanisms:**
# 1. **Acid mechanism** ($k_a \cdot a(\mathrm{H^+})$): Dominant at low pH
# 2. **Neutral mechanism** ($k_n$): pH-independent water-mediated dissolution
# 3. **CO$_2$ mechanism** ($k_c \cdot \mathrm{SR_{CO_2(g)}}$): Dependent on CO$_2$ saturation state
#
# The affinity term $(1 - SR)$ ensures the rate approaches zero at equilibrium.
#
# **Rate constants at 298.15 K:**
# $$
# \begin{aligned}
# &\log_{10} k_a = -0.30 \\
# &\log_{10} k_n = -5.81 \\
# &\log_{10} k_c = -3.48
# \end{aligned}
# $$
#
# Temperature dependence for each mechanism $i \in \{a, n, c\}$:
# $$
# k_i(T_K) = k_i(298.15) \exp\left[-\frac{E_{a,i}}{R} \left(\frac{1}{T_K} - \frac{1}{298.15}\right)\right]
# $$
# with activation energies $E_{a,a} = 14.4$ kJ/mol, $E_{a,n} = 23.5$ kJ/mol, $E_{a,c} = 35.4$ kJ/mol, and gas constant $R = 8.314$ J/mol·K.
#
# ## Coupling chemistry and physical change
#
# ### From surface rate to porosity change
#
# The volumetric reaction rate is calculated from the surface rate:
#
# $$
# R_{\mathrm{vol}} = A \cdot r_{\mathrm{cal}}
# $$
#
# where the total reactive surface area $A = s_s \cdot V_{\mathrm{bulk}} = s_s \cdot (V_m \cdot n_0)$ combines:
# - Specific surface area $s_s = 20$ m²/m³
# - Molar volume $V_m = 3.692 \times 10^{-5}$ m³/mol
# - Initial calcite amount $n_0 = \mathrm{PARM}(1)$
#
# The moles of calcite dissolved over time interval $\Delta t$:
# $$
# \Delta n_{\mathrm{cal}} = R_{\mathrm{vol}} \cdot \Delta t
# $$
#
# The resulting porosity change:
# $$
# \Delta \phi = \frac{\Delta n_{\mathrm{cal}} \cdot V_m}{V}
# $$
# where $V$ is the volume of the computational cell.
#
# **Mass conservation safeguards:**
# - Mineral exhaustion: reaction stops when calcite depleted ($M \leq 0$)
# - Dissolution only: precipitation prevented ($r \geq 0$)
#
# ### Permeability evolution
#
# Permeability in the reactive layer is updated using the Kozeny-Carman equation:
#
# $$
# k = k_0 \frac{(1 - \phi_0)^2}{\phi_0^3} \frac{\phi^3}{(1 - \phi)^2}
# $$
#
# where $k_0$ and $\phi_0$ are initial permeability and porosity, and $k$ and $\phi$ are updated values.
#
#
#
# All thermodynamic calculations use the LLNL database with temperature correction of equilibrium constants. -->
#
# ---
#

# %% [markdown]
# # PHREEQC geochemical modeling
#
# ## PHREEQC configuration
#
# ```xml
# <chemical_system chemical_solver="Phreeqc">
#     <mesh>domain</mesh>
#     <linear_solver>general_linear_solver</linear_solver>
#     <database>llnl.dat</database>
# ```
#
# - **Geochemical Solver**: PHREEQC computes aqueous speciation and activity coefficients
# - **Spatial Coupling**: Domain mesh synchronizes chemical and transport grids
# - **Thermodynamic Database**: LLNL database provides temperature-corrected equilibrium constants
#
# ## Aqueous system definition
#
# ```xml
# <solution>
#     <temperature>25</temperature>
#     <pressure>1</pressure>
#     <pe>4</pe>
#     <components>
#         <component>S(6)</component>
#         <component>C</component>
#         <component>Ca</component>
#         <component>Na</component>
#         <component>Cl</component>
#     </components>
# </solution>
# ```
#
# - **Reference Conditions**: 25°C, 1 atm for thermodynamic calculations
# - **Redox Environment**:
#  $$
#   \mathrm{pe} = -\log(a_{e^-}) = 4 \quad \mathrm{(mildly\;oxidizing\;conditions)}
# $$
# - **Basis Components**: Five master species define the chemical system:
#   - **S(6)**: Total sulfate species $\mathrm{SO}_4^{2-}, \mathrm{HSO}_4^-, \mathrm{NaSO}_4^-, \mathrm{CaSO}_4^0$
#   - **C**: Inorganic carbon system $\mathrm{CO}_2(aq), \mathrm{H}_2\mathrm{CO}_3, \mathrm{HCO}_3^-, \mathrm{CO}_3^{2-}$
#   - **Ca**: Calcium aqueous complexes $\mathrm{Ca}^{2+}, \mathrm{CaCO}_3^0, \mathrm{CaSO}_4^0, \mathrm{CaHCO}_3^+$
#   - **Na/Cl**: Conservative electrolytes for ionic strength and charge balance
#
# - **Initial and inflow pH:**
#   The initial pore water in the quartz layers has pH $5.6$, while the reactive interlayer starts at pH $6.12$. The injected $0.005\,\mathrm{M}$ $\mathrm{H_2SO_4}$ solution corresponds to an inflow pH of $2.13$.
#
# - **pH calculation in PHREEQC:**
#   In the reactive-transport coupling, pH is not prescribed explicitly in PHREEQC. Instead, PHREEQC computes the proton activity $a(\mathrm{H}^+)$ and thus
#   $$
#     \mathrm{pH} = -\log_{10} a(\mathrm{H}^+)
#   $$
#   from the specified total concentrations of the master components (S(6), C, Ca, Na, Cl) together with the imposed redox condition ($\mathrm{pe} = 4$), by enforcing electroneutrality and aqueous speciation equilibria. The initial and inflow pH values reported above correspond to these PHREEQC-equilibrated solutions and are used to define the Dirichlet data for the transported proton component $H$ in the ComponentTransport process.
#
# ## Kinetic reaction framework
#
# ```xml
# <kinetic_reactants>
#     <kinetic_reactant>
#         <name>Calcite</name>
#     </kinetic_reactant>
# </kinetic_reactants>
# ```
#
# Calcite dissolution modeled as time-dependent kinetic process rather than instantaneous equilibrium.
#
# ## Calcite dissolution kinetics implementation
#
# ### Fundamental constants
# ```xml
# <statement>R = 8.314472</statement>
# <statement>deltaT = 1/TK - 1/298.15</statement>
# ```
#
# $$
# R = 8.314472\,\mathrm{J\,mol^{-1}\,K^{-1}}, \qquad
# \Delta\!\left(\frac{1}{T_K}\right)
# = \frac{1}{T_K} - \frac{1}{298.15}\ \mathrm{K^{-1}} .
# $$
#
#
# Calcite dissolution requires three mechanisms to accurately model its pH-dependent behavior: acid mechanism for low pH conditions, neutral mechanism for circumneutral pH, and carbonate mechanism for CO$_2$-rich systems. This multi-mechanism approach ensures thermodynamic consistency across environmental conditions.
#
# #### Acid mechanism (H$^+$-Promoted)
# ```xml
# <statement>logK25 = -0.30</statement>
# <statement>Ea = 14400</statement>
# <statement>ny = 1.0</statement>
# <statement>mech_a = (10^logK25) * ACT("H+")^ny * EXP(- Ea / R * deltaT)</statement>
# ```
# $$
# r_a = k_{a,25} \cdot (a_{\mathrm{H^+}})^{n_y} \cdot \exp\left[-\frac{E_a}{R} \left(\frac{1}{T_K} - \frac{1}{298.15}\right)\right]
# $$
# where:
# - $k_{a,25} = 10^{-0.30}$ mol/m²/s (rate constant at 25°C)
# - $E_a = 14,400$ J/mol (activation energy)
# - $n_y = 1.0$ (reaction order with respect to H⁺ activity)
#
# #### Neutral mechanism (Water-Promoted)
# ```xml
# <statement>logK25 = -5.81</statement>
# <statement>Ea = 23500</statement>
# <statement>mech_b = (10^logK25) * EXP(- Ea / R * deltaT)</statement>
# ```
# $$
# r_n = k_{n,25} \cdot \exp\left[-\frac{E_a}{R} \left(\frac{1}{T_K} - \frac{1}{298.15}\right)\right]
# $$
# where:
# - $k_{n,25} = 10^{-5.81}$ mol/m²/s
# - $E_a = 23,500$ J/mol
#
# #### Carbonate mechanism (CO$_2$-Promoted)
# ```xml
# <statement>logK25 = -3.48</statement>
# <statement>Ea = 35400</statement>
# <statement>ny = 1.0</statement>
# <statement>mech_c = (10^logK25) * SR("CO2(g)")^ny * EXP(- Ea / R * deltaT)</statement>
# ```
# $$
# r_c = k_{c,25} \cdot (\mathrm{SR}_{\mathrm{CO_2(g)}})^{n_y} \cdot \exp\left[-\frac{E_a}{R} \left(\frac{1}{T_K} - \frac{1}{298.15}\right)\right]
# $$
# where:
# - $k_{c,25} = 10^{-3.48}$ mol/m²/s
# - $E_a = 35,400$ J/mol
# - $\mathrm{SR}_{\mathrm{CO_2(g)}}$ = saturation ratio relative to CO$_2$(g)
#
# ### Combined surface rate
# ```xml
# <statement>K_temp = mech_a + mech_b + mech_c</statement>
# ```
# $$
# r_{\mathrm{surface}} = r_a + r_n + r_c \quad \mathrm{[mol·m}^{-2}\mathrm{·s}^{-1}\mathrm{]}
# $$
#
#
# ## **Thermodynamic driving force**
#
# ### **Saturation Ratio**
# The **Saturation Ratio (SR)** quantifies how far the solution is from mineral equilibrium:
#
# ```xml
# <statement>SR("Calcite")</statement>
# ```
# $$
# \mathrm{SR} = \frac{\mathrm{Ionic\;Activity\;Product}}{\mathrm{Equilibrium\;Constant}} = \frac{a(\mathrm{Ca^{2+}}) \cdot a(\mathrm{CO_3^{2-}})}{K_{eq}}
# $$
# where $K_{eq} = 10^{-8.48}$ for calcite at 25°C.
#
#
# The affinity term $(1 - SR)$ provides the **thermodynamic driving force**:
#
# - **$\mathbf{SR < 1}$**: Solution is **undersaturated** → dissolution occurs
# - **$\mathbf{SR \approx 1}$**: System at **equilibrium** → zero net reaction
# - **$\mathbf{SR > 1}$**: Solution is **supersaturated** → precipitation may occur
#
# **Physical Meaning:** Without $(1 - SR)$, dissolution would continue indefinitely even at equilibrium. This term ensures the reaction rate **approaches zero** as the system reaches chemical equilibrium, making the kinetics **thermodynamically consistent**.
#
# **In the rate law:**
# $$
# r_{\mathrm{vol}} = A \cdot r_{\mathrm{surface}} \cdot \underbrace{\left[1 - (\mathrm{SR})^\theta\right]^\eta}_{\mathrm{Driving\;force}}
# $$
#
#
# This is the mathematical equivalent of sugar stopping dissolution when tea becomes saturated. It's the "fullness brake" that makes the model physically realistic.
#
#
# ### Affinity term parameters
# ```xml
# <statement>teta = 1</statement>
# <statement>eta = 1</statement>
# ```
# $$
# \theta = 1, \quad \eta = 1 \quad \mathrm{(affinity\;term\;exponents)}
# $$
#
# ### Surface area calculation
# ```xml
# <statement>molar_volume = 3.692e-05</statement>
# <statement>bulk_volume = molar_volume * PARM(1)</statement>
# <statement>specific_surface_area = 20</statement>
# <statement>surface_area = specific_surface_area * bulk_volume</statement>
# ```
# $$
# A = s_s \cdot V_{\mathrm{bulk}} = s_s \cdot (V_m \cdot n_0)
# $$
# where:
# - $s_s = 20$ m²/m³ (specific surface area)
# - $V_m = 3.692 \times 10^{-5}$ m³/mol (calcite molar volume)
# - $n_0 = \mathrm{PARM}(1)$ (initial calcite moles)
#
# This calculation determines how much mineral surface is exposed to water. Just as crushed ice melts faster than a solid cube, higher surface area dramatically increases dissolution rates by providing more contact points for chemical reactions.
#
# ### Final volumetric rate
# ```xml
# <statement>rate = surface_area * K_temp * (1 - SR("Calcite")^teta)^eta</statement>
# ```
# $$
# r_{\mathrm{vol}} = A \cdot r_{\mathrm{surface}} \cdot \left[1 - (\mathrm{SR}_{\mathrm{Calcite}})^\theta\right]^\eta \quad \mathrm{[mol·s}^{-1}\mathrm{]}
# $$
# where the saturation ratio is:
# $$
# \mathrm{SR}_{\mathrm{Calcite}} = \frac{a_{\mathrm{Ca^{2+}}} \cdot a_{\mathrm{CO_3^{2-}}}}{K_{\mathrm{eq}}}, \quad K_{\mathrm{eq}} = 10^{-8.48}
# $$
#
# ## Mass conservation implementation
#
# ```xml
# <statement>if (M &lt;= 0) then goto 25</statement>
# <statement>if rate &lt; 0 then rate = 0</statement>
# <statement>moles = rate * time</statement>
# <statement>save moles</statement>
# ```
#
# - **Mineral Exhaustion**: Terminates reaction when calcite depleted ($M \leq 0$)
# - **Dissolution Only**: Prevents precipitation ($r_{\mathrm{vol}} \geq 0$)
# - **Mass Change**: $\Delta n = r_{\mathrm{vol}} \cdot \Delta t$ (moles reacted)
# - **State Update**: Saves mass change for system update
#
# ## Numerical controls
#
# ```xml
# <knobs>
#     <max_iter>100</max_iter>
#     <relative_convergence_tolerance>1e-12</relative_convergence_tolerance>
#     <tolerance>1e-15</tolerance>
#     <step_size>100</step_size>
#     <scaling>0</scaling>
# </knobs>
# ```
#
# - **Convergence**: Relative tolerance = $10^{-12}$, absolute tolerance = $10^{-15}$
# - **Iterations**: Maximum 100 iterations
# - **Time Step**: $\Delta t = 100$ seconds
# - **Scaling**: No automatic convergence scaling
#
# ## Variable
#
# | Symbol   | Code Variable          | Meaning                      | Units                     | Value                                  |
# |----------|------------------------|------------------------------|---------------------------|----------------------------------------|
# | $R$      | `R`                    | Gas constant                 | $\mathrm{J\,mol^{-1}\,K^{-1}}$ | $8.314472$                          |
# | $T_K$      | `TK`                   | Temperature                  | $\mathrm{K}$              | –                                      |
# | $\Delta T$ | `deltaT`             | Inverse temperature difference | $\mathrm{K^{-1}}$       | –                                      |
# | $k_{a,25}$ | `10^logK25`          | Acid mechanism rate constant     | $\mathrm{mol\,m^{-2}\,s^{-1}}$ | $10^{-0.30}$                     |
# | $k_{n,25}$ | `10^logK25`          | Neutral mechanism rate constant  | $\mathrm{mol\,m^{-2}\,s^{-1}}$ | $10^{-5.81}$                     |
# | $k_{c,25}$ | `10^logK25`          | Carbonate mechanism rate constant | $\mathrm{mol\,m^{-2}\,s^{-1}}$ | $10^{-3.48}$                     |
# | $E_a$    | `Ea`                  | Activation energy            | $\mathrm{J\,mol^{-1}}$    | $14{,}400 / 23{,}500 / 35{,}400$       |
# | $n_y$    | `ny`                  | Reaction order               | –                         | $1.0$                                 |
# | $\theta$ | `teta`                | Saturation ratio exponent    | –                         | $1$                                   |
# | $\eta$   | `eta`                 | Affinity term exponent       | –                         | $1$                                   |
# | $s_s$    | `specific_surface_area` | Specific surface area       | $\mathrm{m^{2}\,m^{-3}}$  | $20$                                  |
# | $V_m$    | `molar_volume`        | Molar volume                 | $\mathrm{m^{3}\,mol^{-1}}$ | $3.692\times10^{-5}$                 |
# | $n_0$    | `PARM(1)`             | Initial calcite moles        | $\mathrm{mol}$            | User input                             |
#

# %% [markdown]
# # Input data
#
# ### Geometry & Layers
#
# | Item                       | Symbol | Unit |    Value |
# | :------------------------- | :----- | :--: | -------: |
# | Column height              | $L$    |   m  |   $0.07$ |
# | Inner diameter             | —      |   m  | $0.0099$ |
# | Upstream layer length      | —      |   m  |  $0.024$ |
# | Reactive interlayer length | —      |   m  |  $0.022$ |
# | Downstream layer length    | —      |   m  |  $0.024$ |
#
# ### Fluid & Transport
#
# | Item                            | Symbol       |         Unit          |                               Value |
# | :------------------------------ | :----------- | :-------------------: | ----------------------------------: |
# | Fluid density                   | $\rho_\ell$  | $\mathrm{kg\,m^{-3}}$ |                   $1.0\times10^{3}$ |
# | Dynamic viscosity               | $\mu$        |   $\mathrm{Pa\,s}$    |                  $1.0\times10^{-3}$ |
# | Pore diffusion (all components) | $D_p$        |   $\mathrm{m^{2}\,s^{-1}}$ |                  $1.0\times10^{-9}$ |
# | Longitudinal dispersivity       | $\alpha_L$   |       $\mathrm{m}$    |                  $1.0\times10^{-4}$ |
# | Transverse dispersivity         | $\alpha_T$   |       $\mathrm{m}$    | $1.0\times10^{-4}$ *(1D: not used)* |
# | Volumetric flow (reference)     | —            |   $\mathrm{m^{3}\,s^{-1}}$ |                $1.67\times10^{-10}$ |
#
#
# ### Media: Initial Porosity $,\phi_0,$ & Permeability $,k_0$
#
# | Layer               | $\phi_0$ [–] |        $k_0$ [m$^2$] |
# | :------------------ | -----------: | -------------------: |
# | Upstream quartz     |     $0.4277$ | $2.95\times10^{-12}$ |
# | Reactive interlayer |     $0.3255$ | $1.52\times10^{-12}$ |
# | Downstream quartz   |     $0.4685$ | $3.21\times10^{-12}$ |
#
#
# ### Solid Phase (Calcite)
#
# | Item                                        | Symbol |        Unit        |                Value |
# | :------------------------------------------ | :----- | :----------------: | -------------------: |
# | Molar volume                                | $V_m$  | $\mathrm{m^3\,mol^{-1}}$ | $3.693\times10^{-5}$ |
# | Specific surface area                       | $s_s$  | $\mathrm{m^2\,m^{-3}}$   |                 $20$ |
# | Calcite volume fraction (interlayer, init.) | —      |         —          |             $0.2185$ |
#
#
# ### Transported Components
#
# | Transported Component  | Symbol                           | Primary Role in the System                                                                |
# | :--------------------- | :------------------------------- | :---------------------------------------------------------------------------------------- |
# | Sulfate                | $\mathrm{S(6)}$                  | Tracks the injected acid anion $\mathrm{H_2SO_4}$ (conservative w.r.t. mineral reactions) |
# | Proton                 | $\mathrm{H^+}$                   | Drives acidity (low pH) and dissolution rates                                             |
# | Total Inorganic Carbon | $\mathrm{C}$                     | Aggregates $\mathrm{CO_2(aq)}$, $\mathrm{HCO_3^-}$, $\mathrm{CO_3^{2-}}$                  |
# | Calcium                | $\mathrm{Ca^{2+}}$               | Product of calcite dissolution; controls saturation                                       |
# | Sodium & Chloride      | $\mathrm{Na^+}$, $\mathrm{Cl^-}$ | Background ions for charge balance / ionic strength                                       |
#
# ### Initial Pore-Water Composition (by Zone)
#
# | Zone                    | $c_C$ [mol L$^{-1}$] | $c_{\mathrm{Ca}}$ [mol L$^{-1}$] | $c_{\mathrm{S(6)}}$ [mol L$^{-1}$] |         pH | $[\mathrm{H}^+]$ [mol L$^{-1}$] |
# | :---------------------- | -------------------: | -------------------------------: | ---------------------------------: | ---------: | ------------------------------: |
# | Quartz layers (up/down) |   $1.6\times10^{-5}$ |               $1.0\times10^{-8}$ |                 $1.0\times10^{-8}$ |      $5.6$ |           $2.5119\times10^{-6}$ |
# | Reactive interlayer     | $3.793\times10^{-2}$ |            $7.9567\times10^{-3}$ |                 $1.0\times10^{-8}$ | **$6.12$** |         **$7.59\times10^{-7}$** |
#
# ### Boundary / Inlet Chemistry
#
# | Item              | Symbol                         |      Unit       |                                              Value |
# | :---------------- | :----------------------------- | :-------------: | -------------------------------------------------: |
# | Infiltrant acid   | —                              |        —        | $0.005\,\mathrm{M}\ \mathrm{H_2SO_4}$ (pH $=2.128$) |
# | Inlet sulfate     | $c_{\mathrm{S(6)},\mathrm{in}}$ | $\mathrm{mol\,L^{-1}}$ |                                 $5.0\times10^{-3}$ |
# | Inlet proton      | $[\mathrm{H}^+]_{\mathrm{in}}$  | $\mathrm{mol\,L^{-1}}$ |                 $10^{-2.1281}=7.4456\times10^{-3}$ |
# | Inlet carbon      | $c_{C,\mathrm{in}}$            | $\mathrm{mol\,L^{-1}}$ |                                 $1.6\times10^{-5}$ |
# | Inlet calcium     | $c_{\mathrm{Ca},\mathrm{in}}$  | $\mathrm{mol\,L^{-1}}$ |                                 $1.0\times10^{-8}$ |
# | Inlet NaCl        | —                              | $\mathrm{mol\,L^{-1}}$ |                                 $3.0\times10^{-8}$ |
#
#
# **Notes**
#
# * 1D setups don't use $\alpha_T$. setting $\alpha_T=\alpha_L$ is harmless and keeps 2D ready.
# * For the interlayer pH = $6.12$, $[\mathrm{H}^+] = 10^{-6.12} = 7.59\times10^{-7}\ \mathrm{mol,L^{-1}}$.
#
# ---
#

# %% [markdown]
# ## Input Data

# %%
dimension_str = "1d"  # "1d" or "2d"

# %%
_dim_map = {"1d": 1, "2d": 2}
dim_norm = str(dimension_str).strip().lower()
if dim_norm not in _dim_map:
    msg = (
        f"Unsupported dimension '{dimension_str}'. Use one of: {list(_dim_map.keys())}"
    )
    raise ValueError(msg)
target_dim = _dim_map[dim_norm]
dim_suffix = dim_norm

MESH_NAME = "rect_band_default"
MESH_ROOT = out_dir / "meshes"
MESH_ROOT.mkdir(parents=True, exist_ok=True)
(MESH_ROOT / dim_suffix).mkdir(parents=True, exist_ok=True)


# geometry / mesh
GEOMETRY = {
    "L": 0.07,  # m
    "H": 0.0099,  # m
    "left_L": 0.024,  # m
    "band_L": 0.022,  # m
    "right_L": 0.024,  # m
    "h": 0.002,  # m
}
MESH = {
    "mesh_name": MESH_NAME,
    "msh_file": MESH_ROOT / dim_suffix / f"{MESH_NAME}.msh",
    "vtu_file": MESH_ROOT / dim_suffix / "domain.vtu",
}

# chemical fields
FIELDS = {
    "C_main": 1.6e-5,
    "Ca_main": 1.0e-8,
    "H_main": 2.5119e-6,  # [H⁺] = 10^(-pH) = 10^(-5.6) = 2.5119e-6 mol/L
    "C_mid": 3.7934e-2,
    "Ca_mid": 7.9567e-3,
    "H_mid": 7.59e-7,  # [H⁺] = 10^(-pH) = 10^(-6.12) ≈ 7.59e-7 mol/L
}

# material properties
KAPPA = {
    "kappa_upstream": 2.95e-12,
    "kappa_reactive_init": 1.52e-12,
    "kappa_downstream": 3.21e-12,
}
FLUID_PHYSICS = {
    "rho_liq": 1000.0,  # kg/m^3
    "mu_liq": 1.0e-3,  # Pa·s
    "alpha_L": 1e-4,  # m
    "alpha_T": 1e-4,  # m
    "pore_diffusion": 1.0e-9,  # m^2/s
    "retardation_factor": 1.0,
    "decay_rate": 0.0,
    "components": ["S(6)", "C", "Ca", "Na", "Cl", "H"],
}
SOLID = {
    "molar_volume": 3.693e-05,  # m^3/mol
    "volume_fraction_default": 0.0,
    "volume_fraction_by_mid": {"1": 0.2185},
}

# initials / boundaries (scalars)
PARAMS = {
    "porosity_upstream": 0.4277,
    "porosity_reactive": 0.3255,
    "porosity_downstream": 0.4685,
    "p_initial": 101325,
    "p_outlet": 101325,
    "q_inlet": 2.16948e-3,  # mass-flux = rho*q if rho=1000 and q≈2.16948e-6 m/s
    "S6_initial": 1.0e-8,  # Initial sulfate in pore water: 10⁻⁵ mmol/L = 1e-8 mol/L
    "S6_inlet": 5.0e-3,  # Inlet sulfate from 0.005M H₂SO₄ = 0.005 mol/L SO₄²⁻
    "H_inlet": 7.4456e-3,  # [H⁺] = 10^(-pH) = 10^(-2.1281) ≈ 0.0074456
    "C_inlet": 1.6e-5,
    "Ca_inlet": 1.0e-8,
    "NaCl_initial": 3.0e-8,
    "NaCl_inlet": 3.0e-8,
}

# project / io paths
PROJECT = {
    "basename": "porosityIncrease_hc",
    "pvd_path": Path(f"{out_dir}/hc_porosity_change_{dim_suffix}.pvd"),
    "chem_db_src": Path.cwd() / "llnl.dat",
    "chem_db_dst": out_dir / "llnl.dat",
    "workdir": Path.cwd(),
    "out_dir": out_dir,
}


#  time stepping values
time_settings = {
    "t_initial": 0,
    "t_end": 360000,
    "initial_dt": 120,
    "minimum_dt": 120,
    "maximum_dt": 120,
}

# # copy chemistry database
if not PROJECT["chem_db_src"].is_file():
    msg = f"Source file not found: {PROJECT['chem_db_src']}"
    raise FileNotFoundError(msg)
copy2(PROJECT["chem_db_src"], PROJECT["chem_db_dst"])
# %% [markdown]
# ## Mesh generation using GMSH


# %%
def _compute_nx(
    left_L: float, band_L: float, right_L: float, h: float
) -> tuple[int, int, int]:
    """Compute element counts in x for left, band, right regions."""
    nx_left = max(1, round(left_L / h))
    nx_band = max(2, round(band_L / h))
    nx_right = max(1, round(right_L / h))

    return nx_left, nx_band, nx_right


# --------------------------------------------------------------------------------------
# 2D structured rectangle
# --------------------------------------------------------------------------------------
def create_rect_structured(
    filepath: str | Path,
    L: float = 0.07,
    H: float = 0.0099,
    left_L: float = 0.024,
    band_L: float = 0.022,
    right_L: float = 0.024,
    h: float = 0.0005,
    generate_mesh: bool = True,
) -> Path:
    filepath = Path(filepath)
    model_name = filepath.stem

    if abs((left_L + band_L + right_L) - L) > 1e-12:
        msg = f"Invalid lengths: left_L + band_L + right_L = {left_L + band_L + right_L} != L = {L}"
        raise ValueError(msg)
    if not (0.0 < band_L < L):
        msg = f"Invalid band_L={band_L} for L={L}."
        raise ValueError(msg)
    if H <= 0:
        msg = f"Invalid height H={H}."
        raise ValueError(msg)

    # element counts
    nx_left, nx_band, nx_right = _compute_nx(left_L, band_L, right_L, h)
    ny = max(1, round(H / h))

    # geometry coordinates
    xL = -L / 2.0
    x1 = xL + left_L
    x2 = x1 + band_L
    xR = L / 2.0
    yB = -H / 2.0
    yT = H / 2.0

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Verbosity", 1)
        gmsh.option.setNumber("Mesh.RecombineAll", 1)  # global recombine

        gmsh.model.add(model_name)
        geo = gmsh.model.geo
        mgeo = gmsh.model.geo.mesh

        # -----------------------
        # Points (counter-clockwise)
        # -----------------------
        geo.addPoint(xL, yB, 0, h, tag=1)
        geo.addPoint(x1, yB, 0, h, tag=2)
        geo.addPoint(x2, yB, 0, h, tag=3)
        geo.addPoint(xR, yB, 0, h, tag=4)
        geo.addPoint(xR, yT, 0, h, tag=5)
        geo.addPoint(x2, yT, 0, h, tag=6)
        geo.addPoint(x1, yT, 0, h, tag=7)
        geo.addPoint(xL, yT, 0, h, tag=8)

        # -----------------------
        # Lines
        # -----------------------
        geo.addLine(1, 2, tag=1)  # bottom left
        geo.addLine(2, 3, tag=2)  # bottom band
        geo.addLine(3, 4, tag=3)  # bottom right
        geo.addLine(4, 5, tag=4)  # right vertical
        geo.addLine(5, 6, tag=5)  # top right
        geo.addLine(6, 7, tag=6)  # top band
        geo.addLine(7, 8, tag=7)  # top left
        geo.addLine(8, 1, tag=8)  # left vertical
        geo.addLine(2, 7, tag=9)  # internal left vertical
        geo.addLine(3, 6, tag=10)  # internal right vertical

        # -----------------------
        # Surfaces
        # -----------------------
        def plane(loop_id: int, edges: list[int], surf_id: int) -> None:
            geo.addCurveLoop(edges, loop_id)
            geo.addPlaneSurface([loop_id], surf_id)

        plane(1, [1, 9, 7, 8], 1)
        plane(2, [2, 10, 6, -9], 2)
        plane(3, [3, 4, 5, -10], 3)

        # -----------------------
        # Transfinite definitions
        # -----------------------
        mgeo.setTransfiniteCurve(1, nx_left + 1)
        mgeo.setTransfiniteCurve(7, nx_left + 1)

        mgeo.setTransfiniteCurve(2, nx_band + 1)
        mgeo.setTransfiniteCurve(6, nx_band + 1)

        mgeo.setTransfiniteCurve(3, nx_right + 1)
        mgeo.setTransfiniteCurve(5, nx_right + 1)

        mgeo.setTransfiniteCurve(8, ny + 1)
        mgeo.setTransfiniteCurve(4, ny + 1)
        mgeo.setTransfiniteCurve(9, ny + 1)
        mgeo.setTransfiniteCurve(10, ny + 1)

        for s in (1, 2, 3):
            mgeo.setTransfiniteSurface(s, "Alternate")
            mgeo.setRecombine(2, s)

        geo.synchronize()

        # -----------------------
        # Physical groups
        # -----------------------
        pg = gmsh.model.addPhysicalGroup
        pn = gmsh.model.setPhysicalName

        # boundaries (lines, dim = 1)
        right = pg(1, [4])
        pn(1, right, "Right")

        left = pg(1, [8])
        pn(1, left, "Left")

        # domains (surfaces, dim = 2)
        dom1 = pg(2, [1])
        pn(2, dom1, "domain1")

        mid = pg(2, [2])
        pn(2, mid, "Reactive_domain")

        dom2 = pg(2, [3])
        pn(2, dom2, "domain2")

        comp = pg(2, [1, 2, 3])
        pn(2, comp, "Computational_domain")
        if generate_mesh:
            gmsh.model.mesh.generate(2)

        out = filepath if filepath.suffix == ".msh" else filepath.with_suffix(".msh")
        gmsh.write(str(out))
        return out

    finally:
        gmsh.finalize()


# --------------------------------------------------------------------------------------
# 1D structured line
# --------------------------------------------------------------------------------------
def create_line_structured(
    filepath: str | Path,
    L: float = 0.07,
    left_L: float = 0.024,
    band_L: float = 0.022,
    right_L: float = 0.024,
    h: float = 0.0005,
    generate_mesh: bool = True,
) -> Path:
    filepath = Path(filepath)
    model_name = filepath.stem

    if abs((left_L + band_L + right_L) - L) > 1e-12:
        msg = f"Invalid lengths: left_L + band_L + right_L = {left_L + band_L + right_L} != L = {L}"
        raise ValueError(msg)
    if not (0.0 < band_L < L):
        msg = f"band_L must be in (0, L). Got band_L={band_L} for L={L}."
        raise ValueError(msg)

    nx_left, nx_band, nx_right = _compute_nx(left_L, band_L, right_L, h)

    xL = -L / 2.0
    x1 = xL + left_L
    x2 = x1 + band_L
    xR = L / 2.0

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Verbosity", 1)
        gmsh.model.add(model_name)
        geo = gmsh.model.geo
        mgeo = gmsh.model.geo.mesh

        # points
        geo.addPoint(xL, 0.0, 0.0, h, tag=1)
        geo.addPoint(x1, 0.0, 0.0, h, tag=2)
        geo.addPoint(x2, 0.0, 0.0, h, tag=3)
        geo.addPoint(xR, 0.0, 0.0, h, tag=4)

        # lines
        geo.addLine(1, 2, tag=1)
        geo.addLine(2, 3, tag=2)
        geo.addLine(3, 4, tag=3)

        # transfinite along the line
        mgeo.setTransfiniteCurve(1, nx_left + 1)
        mgeo.setTransfiniteCurve(2, nx_band + 1)
        mgeo.setTransfiniteCurve(3, nx_right + 1)

        geo.synchronize()

        pg = gmsh.model.addPhysicalGroup
        pn = gmsh.model.setPhysicalName

        # boundaries
        right = pg(0, [4])
        pn(0, right, "Right")

        left = pg(0, [1])
        pn(0, left, "Left")

        # domains (lines, dim = 1)
        dom1 = pg(1, [1])
        pn(1, dom1, "domain1")

        mid = pg(1, [2])
        pn(1, mid, "Reactive_domain")

        dom2 = pg(1, [3])
        pn(1, dom2, "domain2")

        comp = pg(1, [1, 2, 3])
        pn(1, comp, "Computational_domain")

        if generate_mesh:
            gmsh.model.mesh.generate(1)

        out = filepath if filepath.suffix == ".msh" else filepath.with_suffix(".msh")
        gmsh.write(str(out))
        return out

    finally:
        gmsh.finalize()


def generate_mesh(dim: str) -> Path:
    d = dim.strip().lower()
    mesh_dir = MESH_ROOT / d
    mesh_dir.mkdir(parents=True, exist_ok=True)
    fp = mesh_dir / f"{MESH_NAME}.msh"

    if d == "2d":
        return create_rect_structured(
            filepath=fp,
            L=GEOMETRY["L"],
            H=GEOMETRY["H"],
            left_L=GEOMETRY["left_L"],
            band_L=GEOMETRY["band_L"],
            right_L=GEOMETRY["right_L"],
            h=GEOMETRY["h"],
            generate_mesh=True,
        )
    if d == "1d":
        return create_line_structured(
            filepath=fp,
            L=GEOMETRY["L"],
            left_L=GEOMETRY["left_L"],
            band_L=GEOMETRY["band_L"],
            right_L=GEOMETRY["right_L"],
            h=GEOMETRY["h"],
            generate_mesh=True,
        )

    msg = "dim must be '1d' or '2d'"
    raise ValueError(msg)


# %%
msh_path = generate_mesh(dim_suffix)
# %% [markdown]
# ## Convert msh to vtu for running OpenGeoSys simulations


# %%
def plot_2d_meshes(meshes_2d):
    figures = []
    if not meshes_2d:
        return figures

    domain_mesh = meshes_2d.get("domain", next(iter(meshes_2d.values())))

    try:
        fig = plt.figure(figsize=(10, 3.0), dpi=150)
        ax = fig.add_subplot(111)

        domain_mesh.plot_contourf(
            "MaterialIDs", ax=ax, cbar=True, linewidth=0.2, cmap="viridis"
        )

        ax.set_aspect("equal")
        ax.set_xlabel("x / m", fontsize=11)
        ax.set_ylabel("y / m", fontsize=11)
        ax.tick_params(labelsize=10)
        ax.set_title("2D Mesh - Material IDs", fontsize=13, pad=15)

        boundary_colors = {"Left": "red", "Right": "blue"}
        for name, mesh in meshes_2d.items():
            for boundary_key, color in boundary_colors.items():
                if boundary_key.lower() in name.lower():
                    try:
                        points = mesh.points
                        if len(points) > 1:
                            if boundary_key in ["Left", "Right"]:
                                sorted_idx = np.argsort(points[:, 1])
                            else:
                                sorted_idx = np.argsort(points[:, 0])
                            ax.plot(
                                points[sorted_idx, 0],
                                points[sorted_idx, 1],
                                color=color,
                                linewidth=5,
                                alpha=0.8,
                                label=name,
                            )
                    except Exception:
                        continue

        if any(
            any(
                boundary_key.lower() in name.lower() for boundary_key in boundary_colors
            )
            for name in meshes_2d
        ):
            ax.legend(fontsize=9, loc="upper right", frameon=True)

        fig.tight_layout()
        figures.append(fig)

    except Exception as e:
        print(f"Could not create 2D contour plot: {e}")

    return figures


def plot_1d_meshes(meshes_1d):
    figures = []
    if not meshes_1d:
        return figures

    mesh = meshes_1d.get("domain", next(iter(meshes_1d.values())))
    fig, ax = plt.subplots(figsize=(10, 2.5), dpi=150)

    points = mesh.points
    if len(points) == 0:
        return figures

    cell_data = None
    for key in ("MaterialIDs", "MaterialID", "MaterialId"):
        if key in mesh.cell_data and len(mesh.cell_data[key]) > 0:
            cell_data = np.asarray(mesh.cell_data[key]).flatten()
            break

    segments = []
    if hasattr(mesh, "cells"):
        cells = np.asarray(mesh.cells)
        if len(cells) > 0:
            if cells.ndim == 1:
                i = 0
                while i < len(cells):
                    n_points = cells[i]
                    cell_points = cells[i + 1 : i + 1 + n_points]
                    for j in range(len(cell_points) - 1):
                        x1 = points[cell_points[j]][0]
                        x2 = points[cell_points[j + 1]][0]
                        segments.append([[x1, 0], [x2, 0]])
                    i += 1 + n_points
            else:
                for cell in cells:
                    if len(cell) >= 2:
                        x1 = points[cell[0]][0]
                        x2 = points[cell[1]][0]
                        segments.append([[x1, 0], [x2, 0]])

    if not segments and len(points) > 1:
        for i in range(len(points) - 1):
            x1 = points[i][0]
            x2 = points[i + 1][0]
            segments.append([[x1, 0], [x2, 0]])
        if cell_data is None:
            cell_data = np.ones(len(segments), dtype=int)

    if not segments:
        return figures

    segments = np.array(segments)

    lc = LineCollection(
        segments,
        linewidth=8.0,
        capstyle="round",
        joinstyle="round",
        alpha=0.9,
        zorder=2,
        cmap="viridis",
    )

    if cell_data is not None and len(cell_data) == len(segments):
        lc.set_array(cell_data)
        cbar = plt.colorbar(lc, ax=ax, pad=0.02, shrink=0.8)
        cbar.set_label("MaterialIDs", fontsize=11, labelpad=8)
        cbar.ax.tick_params(labelsize=9)
        if np.issubdtype(cell_data.dtype, np.integer):
            unique_vals = np.unique(cell_data)
            cbar.set_ticks(unique_vals)
            cbar.set_ticklabels([str(int(v)) for v in unique_vals])

    ax.add_collection(lc)

    x_coords = points[:, 0]
    ax.plot(
        x_coords,
        np.zeros_like(x_coords),
        "o",
        markersize=4,
        color="white",
        markeredgecolor="black",
        markeredgewidth=1.2,
        alpha=0.95,
        zorder=3,
        label="Nodes",
    )

    boundary_colors = {"Left": "red", "Right": "blue"}
    boundary_markers = []

    for name, boundary_mesh in meshes_1d.items():
        for boundary_key, color in boundary_colors.items():
            if boundary_key.lower() in name.lower():
                try:
                    boundary_points = boundary_mesh.points
                    if len(boundary_points) > 0:
                        if boundary_key == "Left":
                            x_pos = boundary_points[0][0]
                            marker = ">"
                        else:
                            x_pos = (
                                boundary_points[0][0]
                                if len(boundary_points) == 1
                                else boundary_points[-1][0]
                            )
                            marker = "<"

                        ax.plot(
                            x_pos,
                            0,
                            marker=marker,
                            markersize=10,
                            color=color,
                            markeredgecolor="black",
                            markeredgewidth=1.0,
                            zorder=4,
                        )
                        boundary_markers.append(True)
                except Exception:
                    continue

    if len(x_coords) > 0:
        x_min, x_max = x_coords.min(), x_coords.max()
        if x_max == x_min:
            x_min -= 0.5
            x_max += 0.5
        else:
            padding = (x_max - x_min) * 0.05
            x_min -= padding
            x_max += padding
        ax.set_xlim(x_min, x_max)

    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])

    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color("#404040")
    ax.spines["bottom"].set_linewidth(0.8)

    ax.tick_params(axis="x", which="both", bottom=True, labelsize=10)
    ax.set_xlabel("x / m", fontsize=11, labelpad=8)
    ax.set_ylabel("")

    ax.grid(True, axis="x", alpha=0.2, linestyle="-", linewidth=0.5)
    ax.set_title("1D Mesh - Material IDs", fontsize=13, pad=15, weight="normal")

    legend_handles = []
    legend_handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            color="white",
            markeredgecolor="black",
            markersize=6,
            linestyle="None",
            label="Nodes",
        )
    )

    if boundary_markers:
        legend_handles.append(
            Line2D(
                [0],
                [0],
                marker=">",
                color="red",
                markeredgecolor="black",
                markersize=8,
                linestyle="None",
                label="Left",
            )
        )
        legend_handles.append(
            Line2D(
                [0],
                [0],
                marker="<",
                color="blue",
                markeredgecolor="black",
                markersize=8,
                linestyle="None",
                label="Right",
            )
        )

    if legend_handles:
        ax.legend(
            handles=legend_handles,
            loc="upper right",
            frameon=True,
            framealpha=0.9,
            edgecolor="#CCCCCC",
            fontsize=9,
            handletextpad=0.5,
        )

    ax.set_facecolor("#FAFAFA")
    fig.patch.set_facecolor("white")
    fig.set_constrained_layout(True)

    figures.append(fig)
    return figures


def generate_vtu_meshes_from_Gmsh(msh_file: Path, target_dim: int):
    if not msh_file.is_file():
        raise FileNotFoundError(msh_file)
    out_dir = msh_file.parent
    figures = []
    if target_dim == 2:
        meshes = ot.meshes_from_gmsh(msh_file, dim=[2], log=False)
        for name, mesh in meshes.items():
            pv.save_meshio(out_dir / f"{name}.vtu", mesh)
        if meshes:
            figures.extend(plot_2d_meshes(meshes))
    elif target_dim == 1:
        meshes = ot.meshes_from_gmsh(msh_file, dim=[1], log=False)
        for name, mesh in meshes.items():
            pv.save_meshio(out_dir / f"{name}.vtu", mesh)
        if meshes:
            figures.extend(plot_1d_meshes(meshes))
    if figures:
        plt.show()
    return figures


def run_mesh(dim: str):
    d = dim.strip().lower()
    td = {"1d": 1, "2d": 2}[d]
    base_msh = MESH["msh_file"]
    base_vtu = MESH["vtu_file"]
    msh = base_msh.parent.parent / d / base_msh.name
    vtu = base_vtu.parent.parent / d / base_vtu.name
    if not msh.is_file():
        raise FileNotFoundError(msh)
    generate_vtu_meshes_from_Gmsh(msh, target_dim=td)
    return vtu


def assert_structured_mesh_size(mesh_path: Path, target_dim: int) -> None:
    mesh = pv.read(str(mesh_path))

    nx_left, nx_band, nx_right = _compute_nx(
        GEOMETRY["left_L"],
        GEOMETRY["band_L"],
        GEOMETRY["right_L"],
        GEOMETRY["h"],
    )

    if target_dim == 1:
        expected_cells = nx_left + nx_band + nx_right
        expected_points = expected_cells + 1
    elif target_dim == 2:
        ny = max(1, round(GEOMETRY["H"] / GEOMETRY["h"]))
        nx_total = nx_left + nx_band + nx_right
        expected_cells = nx_total * ny
        expected_points = (nx_total + 1) * (ny + 1)
    else:
        msg = f"Unsupported target_dim={target_dim!r} (expected 1 or 2)."
        raise ValueError(msg)

    if mesh.n_cells != expected_cells or mesh.n_points != expected_points:
        msg = (
            "Unexpected mesh size:\n"
            f"  cells : {mesh.n_cells} (expected {expected_cells})\n"
            f"  points: {mesh.n_points} (expected {expected_points})"
        )
        raise RuntimeError(msg)

    print(
        f"Mesh check passed: {mesh.n_cells} cells, {mesh.n_points} points "
        f"(dim = {target_dim}D)."
    )


# %%
vtu_path = run_mesh(dim_suffix)
assert_structured_mesh_size(vtu_path, target_dim)
# %% [markdown]
# ## Pre-processing


# %%
def pre_processing_chemo_by_params(
    mesh_path: str | Path,
    L: float,
    H: float,
    left_L: float,
    band_L: float,
    right_L: float,
    C_main=1.6e-5,
    Ca_main=1.0e-8,
    H_main=2.461e-6,
    C_mid=3.7934e-2,
    Ca_mid=7.9567e-3,
    H_mid=7.485e-7,
    tol: float = 1e-12,
    save_reactive_submesh: bool = True,
):
    mesh_path = Path(mesh_path)
    if not mesh_path.is_file():
        raise FileNotFoundError(mesh_path)

    if abs((left_L + band_L + right_L) - L) > 1e-12 or not (0 < band_L < L):
        msg = f"Invalid geometry parameters: left_L={left_L}, band_L={band_L}, right_L={right_L}, L={L}"
        raise ValueError(msg)

    xL = -L / 2.0
    x1 = xL + left_L
    x2 = x1 + band_L

    m = pv.read(str(mesh_path))
    if m.n_points == 0:
        msg = f"No points in mesh: {mesh_path}"
        raise RuntimeError(msg)

    x = m.points[:, 0]
    is_mid = (x >= x1 - tol) & (x <= x2 + tol)

    n = m.n_points
    C = np.full(n, float(C_main))
    Ca = np.full(n, float(Ca_main))
    H = np.full(n, float(H_main))

    C[is_mid] = float(C_mid)
    Ca[is_mid] = float(Ca_mid)
    H[is_mid] = float(H_mid)

    m.point_data["C"] = C
    m.point_data["Ca"] = Ca
    m.point_data["H"] = H
    m.save(str(mesh_path))

    if save_reactive_submesh:
        reactive_mesh_path = mesh_path.with_name("physical_group_Reactive_domain.vtu")
        if reactive_mesh_path.exists():
            r = pv.read(str(reactive_mesh_path))
            if r.n_points > 0:
                nr = r.n_points
                r.point_data["C"] = np.full(nr, float(C_mid))
                r.point_data["Ca"] = np.full(nr, float(Ca_mid))
                r.point_data["H"] = np.full(nr, float(H_mid))
                r.save(str(reactive_mesh_path))
                return str(mesh_path), str(reactive_mesh_path)
        else:
            print("Reactive submesh not found; nothing to update.")

    return str(mesh_path)


# %%
pre_processing_chemo_by_params(
    mesh_path=vtu_path,
    L=GEOMETRY["L"],
    H=GEOMETRY["H"],
    left_L=GEOMETRY["left_L"],
    band_L=GEOMETRY["band_L"],
    right_L=GEOMETRY["right_L"],
    C_main=FIELDS["C_main"],
    Ca_main=FIELDS["Ca_main"],
    H_main=FIELDS["H_main"],
    C_mid=FIELDS["C_mid"],
    Ca_mid=FIELDS["Ca_mid"],
    H_mid=FIELDS["H_mid"],
)


# %%
def identify_subdomains(dim: str):
    d = dim.strip().lower()
    path = MESH_ROOT / d
    bulk = path / "domain.vtu"
    if not bulk.is_file():
        return

    names = [
        "domain.vtu",
        "physical_group_Computational_domain.vtu",
        "physical_group_domain1.vtu",
        "physical_group_domain2.vtu",
        "physical_group_Left.vtu",
        "physical_group_Reactive_domain.vtu",
        "physical_group_Right.vtu",
    ]

    cwd = Path.cwd()
    os.chdir(path)
    try:
        ot.cli().identifySubdomains("-m", "domain.vtu", "--", *names)
    finally:
        os.chdir(cwd)


identify_subdomains(dim_suffix)
# %% [markdown]
# ### prepare project file

# %%
prj_in = Path(f"{PROJECT['basename']}.prj")
prj_out = out_dir / f"{PROJECT['basename']}_updated_{dim_suffix}.prj"
prj = ot.Project(input_file=prj_in, output_file=prj_out)


# %%
def _set_values(name: str, vec4):
    """Write a 4-vector into <values> (or <value>) for parameter `name`"""
    txt = " ".join(f"{v:.16g}" for v in vec4)
    n = prj.replace_text(txt, xpath=f".//parameters/parameter[name='{name}']/values")
    if n == 0:
        n = prj.replace_text(txt, xpath=f".//parameters/parameter[name='{name}']/value")
    if n == 0:
        msg = f"Parameter '{name}' not found (neither <values> nor <value>)."
        raise KeyError(msg)


def _set_value(name: str, val):
    """Write a scalar into <value> (or <values>) for parameter `name`"""
    sval = f"{val:.16g}" if isinstance(val, int | float) else str(val)
    n = prj.replace_text(sval, xpath=f".//parameters/parameter[name='{name}']/value")
    if n == 0:
        n = prj.replace_text(
            sval, xpath=f".//parameters/parameter[name='{name}']/values"
        )
    if n == 0:
        msg = f"Parameter '{name}' not found (neither <value> nor <values>)."
        raise KeyError(msg)


def _set_all(xpath: str, val):
    """Replace all matches at `xpath` with `val`"""
    sval = f"{val:.16g}" if isinstance(val, int | float) else str(val)
    n = prj.replace_text(sval, xpath=xpath)
    if n == 0:
        msg = f"No nodes matched xpath: {xpath}"
        raise KeyError(msg)


def add_kappa_parameters(prj, dim: int, kappa: dict):
    if dim == 1:
        for name, K in kappa.items():
            prj.parameters.add_parameter(
                name=name, type="Constant", value=f"{float(K):.16g}"
            )
    elif dim == 2:
        for name, K in kappa.items():
            prj.parameters.add_parameter(
                name=name,
                type="Constant",
                values=f"{float(K):.16g} 0 0 {float(K):.16g}",
            )
    elif dim == 3:
        for name, K in kappa.items():
            prj.parameters.add_parameter(
                name=name,
                type="Constant",
                values=" ".join(f"{v:.16g}" for v in [K, 0, 0, 0, K, 0, 0, 0, K]),
            )
    else:
        msg = f"Invalid dim: {dim}. Must be 1, 2, or 3."
        raise ValueError(msg)


# %%
add_kappa_parameters(prj, dim=target_dim, kappa=KAPPA)

for k, v in PARAMS.items():
    _set_value(k, v)

for param, value in time_settings.items():
    _set_all(f".//time_stepping/{param}", value)

# dispersion (all media)
_set_all(
    ".//media/medium/properties/property[name='longitudinal_dispersivity']/value",
    FLUID_PHYSICS["alpha_L"],
)
_set_all(
    ".//media/medium/properties/property[name='transversal_dispersivity']/value",
    FLUID_PHYSICS["alpha_T"],
)

# fluid properties (AqueousLiquid in all media)
_set_all(
    ".//phases/phase[type='AqueousLiquid']/properties/property[name='density']/value",
    FLUID_PHYSICS["rho_liq"],
)
_set_all(
    ".//phases/phase[type='AqueousLiquid']/properties/property[name='viscosity']/value",
    FLUID_PHYSICS["mu_liq"],
)

# component transport (each listed component in all media)
for comp in FLUID_PHYSICS["components"]:
    base = f".//phases/phase[type='AqueousLiquid']/components/component[name='{comp}']/properties/property"
    _set_all(base + "[name='pore_diffusion']/value", FLUID_PHYSICS["pore_diffusion"])
    _set_all(
        base + "[name='retardation_factor']/value", FLUID_PHYSICS["retardation_factor"]
    )
    _set_all(base + "[name='decay_rate']/value", FLUID_PHYSICS["decay_rate"])

# solid Calcite (all media)
_set_all(
    ".//phases/phase[type='Solid']"
    "/components/component[name='Calcite']"
    "/properties/property[name='molar_volume']/value",
    SOLID["molar_volume"],
)

_set_all(
    ".//phases/phase[type='Solid']"
    "/components/component[name='Calcite']"
    "/properties/property[name='volume_fraction']/value",
    SOLID["volume_fraction_default"],
)

for mid, vf in SOLID["volume_fraction_by_mid"].items():
    n = prj.replace_text(
        f"{vf:.16g}",
        xpath=(
            f".//media/medium[@id='{mid}']"
            f"/phases/phase[type='Solid']"
            f"/components/component[name='Calcite']"
            f"/properties/property[name='volume_fraction']/value"
        ),
    )
    if n == 0:
        msg = f"Calcite volume_fraction node not found for medium @id='{mid}. Check your PRJ structure and ids."
        raise KeyError(msg)

force_vec = {1: "0", 2: "0 0", 3: "0 0 0"}[target_dim]


prj.replace_text(force_vec, xpath=".//processes/process/specific_body_force")
prj.replace_text(f"hc_porosity_change_{dim_suffix}", xpath=".//output/prefix")

prj.write_input()
print("Wrote:", prj_out)


# %% [markdown]
# ## Run simulation

# %%
prj.run_model(
    logfile=Path(out_dir, f"out_{dim_suffix}.txt"),
    args=f"-o {out_dir} -m {MESH_ROOT / dim_suffix}",
)
# %% [markdown]
# ## Post-processing

# %%
VARS = ["Calcite_avg", "porosity_avg"]
PVD_PATH = PROJECT["pvd_path"]


def _var_tex(var: str) -> str:
    v = var.lower()
    if "porosity" in v:
        return r"\varphi_\mathrm{avg}"
    if "calcite" in v:
        return r"\mathrm{Ca}_\mathrm{avg}"
    return var


def get_centerline_cell_indices(mesh, y_tol: float = 1e-8):
    """Return cell indices and x-coordinates along the x-centerline."""
    xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
    ext = np.array([xmax - xmin, ymax - ymin, zmax - zmin])
    nz_dims = int(np.sum(ext > 1e-10))

    if nz_dims not in (1, 2):
        msg = "Only 1D or 2D meshes are supported for centerline extraction."
        raise ValueError(msg)

    centers = mesh.cell_centers().points
    xs = centers[:, 0]

    if nz_dims == 1:
        order = np.argsort(xs)
        return order, xs[order], nz_dims, ymin

    ys = centers[:, 1]
    center_y = 0.5 * (ymin + ymax)

    y_unique = np.unique(np.round(ys / y_tol) * y_tol)
    row_y = y_unique[np.argmin(np.abs(y_unique - center_y))]

    row_mask = np.abs(ys - row_y) < y_tol
    if not np.any(row_mask):
        msg = "No centerline row found in y."
        raise RuntimeError(msg)

    row_ids = np.nonzero(row_mask)[0]
    row_xs = xs[row_mask]

    order = np.argsort(row_xs)
    centerline_ids = row_ids[order]
    x_coords = row_xs[order]

    return centerline_ids, x_coords, nz_dims, row_y


def plot_centerline_profiles(
    pvd_path: Path,
    vars_to_plot,
    out_dir: Path | None = None,
):
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )

    ms = ot.MeshSeries(pvd_path)
    n = len(ms.timevalues)
    if n == 0:
        msg = f"No timesteps in PVD file: {pvd_path}"
        raise RuntimeError(msg)

    idxs = sorted({round(f * (n - 1)) for f in (0.0, 0.25, 0.5, 0.75, 1.0)})
    times = [ms.timevalues[i] for i in idxs]

    ref_mesh = ms[0]
    xmin, xmax, *_ = ref_mesh.bounds
    centerline_ids, x_coords, *_ = get_centerline_cell_indices(ref_mesh)
    nx = len(x_coords)

    npz_data: dict[str, np.ndarray] = {
        "axis_label": np.array("x"),
        "coord": x_coords,
        "times": np.array(times, dtype=float),
    }

    cmap = plt.get_cmap("tab10")
    colors = [cmap(i % 10) for i in range(len(idxs))]
    markevery = max(1, nx // 25)

    for var in vars_to_plot:
        fig, ax = plt.subplots(figsize=(4.8, 3.2), dpi=200)
        plotted_any = False
        var_arr = np.full((len(idxs), nx), np.nan, dtype=float)

        for k, (i, t, c) in enumerate(zip(idxs, times, colors, strict=False)):
            mesh = ms[i]

            var_name = None
            for avail_var in mesh.array_names:
                if avail_var.lower() == var.lower():
                    var_name = avail_var
                    break
            if var_name is None:
                continue
            if var_name not in mesh.cell_data:
                continue

            all_cell_values = np.asarray(mesh.cell_data[var_name])
            if all_cell_values.shape[0] <= centerline_ids.max():
                msg = (
                    f"Cell data for variable '{var}' has insufficient length "
                    f"at t={t:g} s."
                )
                raise RuntimeError(msg)

            values = all_cell_values[centerline_ids]
            if values.size != nx:
                continue

            var_arr[k, :] = values
            ax.plot(
                x_coords,
                values,
                "-",
                marker="o",
                markevery=markevery,
                markersize=2.5,
                lw=1.2,
                alpha=0.95,
                color=c,
                markerfacecolor="none",
                label=f"t={t:g} s",
            )
            plotted_any = True

        if not plotted_any:
            plt.close(fig)
            print(f"[skip] '{var}' not found as cell data in any timestep.")
            continue

        npz_data[var] = var_arr

        var_tex = _var_tex(var)

        ax.set_xlim(xmin, xmax)
        ax.set_xlabel(r"$x\,/\,\mathrm{m}$")
        ax.set_ylabel(rf"${var_tex}$")
        ax.tick_params(axis="both", labelsize=9)
        ax.grid(True, which="both", alpha=0.3, linestyle=":")
        ax.legend(loc="best", frameon=False, fontsize=7)
        ax.set_title(rf"${var_tex}$")

        fig.tight_layout()
        plt.show()

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        npz_path = out_dir / f"{pvd_path.stem}_centerline.npz"
        np.savez_compressed(npz_path, **npz_data)
        print(f"[saved] centerline profiles -> {npz_path}")


# %%
plot_centerline_profiles(
    pvd_path=PVD_PATH,
    vars_to_plot=VARS,
    out_dir=Path(PROJECT["out_dir"]),
)


# %% [markdown]
# ## Comparison between 1D and 2D simulations


# %%
def _swap_dim_suffix(stem: str) -> tuple[str, str, str]:
    s = stem.lower()
    if s.endswith("1d"):
        curr = "1d"
        other = "2d"
        other_stem = stem[:-2] + "2d"
    elif s.endswith("2d"):
        curr = "2d"
        other = "1d"
        other_stem = stem[:-2] + "1d"
    else:
        msg = (
            f"Cannot infer dimension from stem '{stem}'. "
            "Expected it to end with '1d' or '2d'."
        )
        raise ValueError(msg)
    return curr, other, other_stem


def _align_other_to_base(
    base_x: np.ndarray, other_x: np.ndarray, other_vals: np.ndarray
) -> np.ndarray:
    if np.allclose(base_x, other_x, rtol=1e-10, atol=1e-12):
        return other_vals

    aligned = np.empty((other_vals.shape[0], base_x.shape[0]))
    for it in range(other_vals.shape[0]):
        aligned[it, :] = np.interp(base_x, other_x, other_vals[it, :])
    return aligned


def _var_tex(var: str) -> str:
    v = var.lower()
    if "porosity" in v:
        return r"\varphi_\mathrm{avg}"
    if "calcite" in v:
        return r"\mathrm{Ca}_\mathrm{avg}"
    return var


def compare_centerline_npz(
    base_npz: Path,
    other_npz: Path,
    vars_to_plot,
    use_relative: bool = False,
    base_label: str = "BASE",
    other_label: str = "OTHER",
):
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )

    base_data = np.load(base_npz)
    other_data = np.load(other_npz)

    base_x = base_data["coord"]
    other_x = other_data["coord"]
    base_t = base_data["times"]
    other_t = other_data["times"]

    if base_t.shape != other_t.shape or not np.allclose(base_t, other_t):
        print("[warn] Times differ between base and other.")
    times = base_t

    print(f"[compare_centerline_npz] {base_label} file : {base_npz}")
    print(f"[compare_centerline_npz] {other_label} file: {other_npz}")
    print(
        "[compare_centerline_npz] x-base:  "
        f"{base_x[0]:.6g} .. {base_x[-1]:.6g} ({len(base_x)} pts)"
    )
    print(
        "[compare_centerline_npz] x-other: "
        f"{other_x[0]:.6g} .. {other_x[-1]:.6g} ({len(other_x)} pts)"
    )

    nt = len(times)
    cmap = plt.get_cmap("tab10")
    colors = [cmap(i % 10) for i in range(nt)]
    markevery = max(1, len(base_x) // 25)

    nvars = len(vars_to_plot)
    if nvars == 0:
        print("[compare_centerline_npz] No variables to plot.")
        return

    ncols = nvars
    nrows = 3 if use_relative else 2

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.7 * ncols, 2.3 * nrows * 1.2),
        dpi=200,
        sharex="col",
    )

    if nvars == 1:
        axes = np.atleast_2d(axes)

    for j, var in enumerate(vars_to_plot):
        if var not in base_data or var not in other_data:
            print(f"[skip] '{var}' not in both files.")
            continue

        base_vals = np.asarray(base_data[var])
        other_vals = np.asarray(other_data[var])

        other_vals_aligned = _align_other_to_base(base_x, other_x, other_vals)
        abs_err = np.abs(base_vals - other_vals_aligned)

        if use_relative:
            eps = 1e-14
            denom = np.maximum(np.abs(base_vals), eps)
            rel_err = abs_err / denom
        else:
            rel_err = None

        var_tex = _var_tex(var)

        ax_prof = axes[0, j]
        ax_abs = axes[1, j]
        ax_rel = axes[2, j] if use_relative else None

        for it, t in enumerate(times):
            c = colors[it]
            ax_prof.plot(
                base_x,
                base_vals[it, :],
                "-",
                marker="o",
                markevery=markevery,
                markersize=2.5,
                lw=1.2,
                alpha=0.95,
                color=c,
                markerfacecolor="none",
                label=f"{base_label}  t={t:g} s",
            )
            ax_prof.plot(
                base_x,
                other_vals_aligned[it, :],
                "--",
                marker="s",
                markevery=markevery,
                markersize=2.5,
                lw=1.2,
                alpha=0.95,
                color=c,
                markerfacecolor="none",
                label=f"{other_label} t={t:g} s",
            )

        ax_prof.set_ylabel(rf"${var_tex}$")
        ax_prof.grid(True, which="both", alpha=0.3, linestyle=":")
        if j == 0:
            ax_prof.legend(loc="best", frameon=False, fontsize=7)
        ax_prof.set_title(rf"${var_tex}$")

        # absolute error
        max_abs_err = np.max(abs_err, axis=0)
        for it, t in enumerate(times):
            c = colors[it]
            ax_abs.plot(
                base_x,
                abs_err[it, :],
                "-",
                marker="^",
                markevery=markevery,
                markersize=4,
                lw=1.0,
                alpha=0.9,
                color=c,
                markerfacecolor="none",
                label=f"t={t:g} s",
            )

        ax_abs.fill_between(
            base_x,
            0.0,
            max_abs_err,
            color="black",
            alpha=0.18,
        )

        ax_abs.set_ylabel(rf"$|\Delta {var_tex}|$")
        ax_abs.grid(True, which="both", alpha=0.3, linestyle=":")
        if j == 0:
            ax_abs.legend(loc="best", frameon=False, fontsize=7)

        # relative error
        if ax_rel is not None and rel_err is not None:
            max_rel_err = np.max(rel_err, axis=0)
            for it, t in enumerate(times):
                c = colors[it]
                ax_rel.plot(
                    base_x,
                    rel_err[it, :],
                    "-",
                    marker="v",
                    markevery=markevery,
                    lw=1.0,
                    alpha=0.9,
                    color=c,
                    markerfacecolor="none",
                    label=f"t={t:g} s",
                )

            ax_rel.fill_between(
                base_x,
                0.0,
                max_rel_err,
                color="black",
                alpha=0.18,
            )

            ax_rel.set_ylabel(rf"$|\Delta {var_tex}| / |{var_tex}|$")
            ax_rel.grid(True, which="both", alpha=0.3, linestyle=":")
            if j == 0:
                ax_rel.legend(loc="best", frameon=False, fontsize=7)

    for j in range(ncols):
        axes[-1, j].set_xlabel(r"$x\,/\,\mathrm{m}$")

    fig.tight_layout()
    plt.show()


def compare_with_other_dim(vars_to_plot, use_relative: bool = True):
    stem = PVD_PATH.stem
    curr_suffix, other_suffix, other_stem = _swap_dim_suffix(stem)

    out_dir = Path(PROJECT["out_dir"])
    exp_dir = Path("expected")

    current_npz = out_dir / f"{stem}_centerline.npz"
    other_npz = exp_dir / f"{other_stem}_centerline.npz"

    if not current_npz.is_file():
        msg = (
            f"Current dim NPZ not found: {current_npz}.\n"
            "Run plot_centerline_profiles first for the current dimension."
        )
        raise FileNotFoundError(msg)
    if not other_npz.is_file():
        msg = (
            f"Expected other-dim NPZ not found: {other_npz}.\n"
            "Check 'expected' folder and the naming convention."
        )
        raise FileNotFoundError(msg)

    base_label = other_suffix.upper()
    other_label = curr_suffix.upper()

    compare_centerline_npz(
        base_npz=other_npz,
        other_npz=current_npz,
        vars_to_plot=vars_to_plot,
        use_relative=use_relative,
        base_label=base_label,
        other_label=other_label,
    )


# %%
compare_with_other_dim(
    vars_to_plot=VARS,
    use_relative=False,
)

# %% [markdown]
# ## Verification of mesh and porosity against reference solution

# %%
ms = MeshSeries(PVD_PATH)
mesh_new = ms[-1]
mesh_ref = pv.read(
    Path("expected") / "hc_porosity_change_1d_ts_3000_t_360000.000000.vtu"
)

if mesh_new.n_points != mesh_ref.n_points or mesh_new.n_cells != mesh_ref.n_cells:
    msg = "Mesh points/cells counts differ."
    raise ValueError(msg)
if not np.array_equal(mesh_new.cells, mesh_ref.cells) or not np.array_equal(
    mesh_new.celltypes, mesh_ref.celltypes
):
    msg = "Cell connectivity/types differ."
    raise ValueError(msg)


def _assert_cell(name, atol=1e-6, rtol=0):
    np.testing.assert_allclose(
        np.asarray(mesh_new.cell_data[name]),
        np.asarray(mesh_ref.cell_data[name]),
        rtol=rtol,
        atol=atol,
    )
    print(f"OK CELL: {name}")


for n, atol, rtol in [("porosity_avg", 1e-6, 0)]:
    if n in mesh_new.cell_data and n in mesh_ref.cell_data:
        _assert_cell(n, atol, rtol)

print("All requested cell arrays matched within tolerances.")

# %% [markdown]
# # References
#
# 1. Lu, Renchao, Thomas Nagel, Jenna Poonoosamy, Dmitri Naumov, Thomas Fischer, Vanessa Montoya, Olaf Kolditz, and Haibing Shao. **A New Operator-Splitting Finite Element Scheme for Reactive Transport Modeling in Saturated Porous Media.** *Computers & Geosciences* 163 (2022): 105106. [https://doi.org/10.1016/j.cageo.2022.105106](https://doi.org/10.1016/j.cageo.2022.105106)
#
# 2. Mollaali, Mostafa, Keita Yoshioka, Renchao Lu, Vanessa Montoya, Victor Vilarrasa, and Olaf Kolditz. **Variational Phase-Field Fracture Approach in Reactive Porous Media.** *International Journal for Numerical Methods in Engineering* 126, no. 1 (2025): e7621. [https://doi.org/10.1002/nme.7621](https://doi.org/10.1002/nme.7621)

# %%
