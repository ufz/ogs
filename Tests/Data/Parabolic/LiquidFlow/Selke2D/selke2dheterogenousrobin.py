# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "Simulating the Selke 2D process with heterogeneous aquifer and Robin boundary conditions"
# date = "2026-07-29"
# author = "Erik Nixdorf, Niklas Ritter"
# web_subsection = "liquid-flow"
# weight = 175
# +++

# %% [markdown]
# # Simulating the Selke 2D process with heterogeneous aquifer and Robin boundary conditions
# The River Selke is a fourth-order stream with a length of 64 km draining a part of the forested Harz Mountains and the agriculturally dominated
# Northern Harz foreland in Central Germany. Due to different water management targets such as the design of infiltration wells along a former pit
# mine the area is of particular interest for hydrogeological research (cf. [Nixdorf et al., 2025](#references))
# In this notebook we show how to append the [previous notebook on the Selke 2D process with Neumann boundary conditions](../selke2ddirichletneumann) to model a heterogeneous aquifer and Robin boundary constraints

# %%
from pathlib import Path

import numpy as np
import ogstools as ot
import pyvista as pv

ot.plot.setup.show_region_bounds = False
# %% [markdown]
# ## 1 Load data and run Dirichlet/Neumann simulation
# We begin by loading the necessary meshes and the project file for the Selke 2D process with Dirichlet and Neumann boundary conditions which we want to modify.  \
# Note, that we use a different domain mesh than before which contains two instead of one MaterialID which allows us to differentiate between Quaternary and non-Quaternary Materials

# %%
project_robin = ot.Project("Selke_Basin.prj").copy()
# We need to add the changes we've made in the previous notebook. To do this we add an xml patch during execution of the models
execution_patch = ot.Execution(
    args="-p " + str(Path.cwd().resolve() / "neumann_patch.xml")
)
input_meshes = [
    "Selke_Basin_Domain_2ids.vtu",
    "Selke_Basin_PG_Concordia.vtu",
    "Selke_Basin_PG_Koenigsauer.vtu",
    "Selke_Basin_PG_Wilsleber.vtu",
    "Selke_Basin_PL_Bode.vtu",
    "Selke_Basin_PL_Wipper.vtu",
    "Selke_Basin_Pnt_Wells.vtu",
    "Selke_Basin_PL_Selke.vtu",
    "Selke_Basin_PL_Upstream.vtu",
    "Selke_Basin_PG_recharge.vtu",
]
meshes_robin = ot.Meshes.from_files(input_meshes)
output_names_dict = {
    "Selke_Basin_Domain_2ids": "Selke Basin",
    "Selke_Basin_PG_Concordia": "Lake Concordia",
    "Selke_Basin_PG_Koenigsauer": "Lake Königsauer",
    "Selke_Basin_PG_Wilsleber": "Lake Wilsleber",
    "Selke_Basin_PL_Bode": "River Bode",
    "Selke_Basin_PL_Wipper": "River Wipper",
    "Selke_Basin_Pnt_Wells": "Wells",
    "Selke_Basin_PL_Selke": "River Selke",
    "Selke_Basin_PL_Upstream": "Upstream",
    "Selke_Basin_PG_recharge": "Selke Basin Recharge",
}
meshes_robin.rename_output_names(output_names_dict)


# %% [markdown]
# ## 2 Modify the project file

# %% [markdown]
# ## 2.1 Changing to heterogeneous aquifer
# We begin by replacing the reference to the old homogeneous mesh with the new 2 MaterialID mesh \
# Material ID 0 represents the Quaternary medium, MaterialID 1 the non-Quaternary
# %%
project_robin.replace_mesh("Selke_Basin_Domain_homo.vtu", "Selke_Basin_Domain_2ids.vtu")
# Visualise MaterialIDs of new domain mesh
fig_material_ids = ot.plot.contourf(
    meshes_robin.domain, "MaterialIDs", fontsize=20, show_edges=False
)
fig_material_ids.axes[0].set_title(
    "MaterialIDs in Selke Basin",
    fontsize=24,
    loc="center",
)
# %% [markdown]
# After this, we add all required properties to the new medium of id 1
# %%
# When we add a property to a previously non-existent medium_id and phase type, OGSTools creates the new medium and phase automatically.
project_robin.media.add_property(
    medium_id="1",
    phase_type="AqueousLiquid",
    name="viscosity",
    type="Constant",
    value=0.0011373,
)
project_robin.media.add_property(
    medium_id="1",
    phase_type="AqueousLiquid",
    name="density",
    type="Constant",
    value=1000,
)
project_robin.media.add_property(
    medium_id="1", name="permeability", type="Constant", value=9.5e-11
)
project_robin.media.add_property(
    medium_id="1", name="reference_temperature", type="Constant", value=288.15
)
project_robin.media.add_property(
    medium_id="1", name="porosity", type="Constant", value=0.05
)
project_robin.media.add_property(
    medium_id="1", name="storage", type="Constant", value=0
)
# %% [markdown]
# We also set the permeability of the Quaternary medium to $5\cdot 10^{-9}$

# %%
project_robin.replace_medium_property_value(0, "permeability", 5e-9)

# %% [markdown]
# ## 2.2 Adding Robin boundary conditions
# The current Dirichlet boundary conditions we use to characterise the rivers overvalues groundwater. \
# We therefore want to change the project file to use Robin boundary conditions, which are given by
# $$\begin{equation} Q_{ex}=K\cdot W\cdot L\cdot \frac{H_{SW}-H_{GW}}{M}\end{equation}$$
# In OGS, the value $\frac{K\cdot W\cdot L}{M}$ is combined into the single parameter bed_conductance

# %%
project_robin.parameters.add_parameter(
    name="bed_conductance",
    type="Constant",
    value=5e-9,
)
# We first have to remove the obsolete boundary conditions
project_robin.remove_element(
    ".//boundary_conditions/boundary_condition[mesh='Selke_Basin_PL_Selke']"
)
project_robin.remove_element(
    ".//boundary_conditions/boundary_condition[mesh='Selke_Basin_PL_Wipper']"
)
project_robin.remove_element(
    ".//boundary_conditions/boundary_condition[mesh='Selke_Basin_PL_Bode']"
)
# Now we can readd the improved boundary conditions
project_robin.process_variables.add_bc(
    process_variable_name="pressure",
    type="Robin",
    mesh="Selke_Basin_PL_Bode",
    u_0="p_river_bode",
    alpha="bed_conductance",
)
project_robin.process_variables.add_bc(
    process_variable_name="pressure",
    type="Robin",
    mesh="Selke_Basin_PL_Selke",
    u_0="p_river_selke",
    alpha="bed_conductance",
)
project_robin.process_variables.add_bc(
    process_variable_name="pressure",
    type="Robin",
    mesh="Selke_Basin_PL_Wipper",
    u_0="p_river_wipper",
    alpha="bed_conductance",
)

# %% [markdown]
# # 3 Running the model and analysis
# ## 3.1 Running the heterogeneous Robin model
# Having finished the setup, we now build and run the model and show the hydraulic head at the final timestep

# %%
model_robin = ot.Model(
    project=project_robin, meshes=meshes_robin, execution=execution_patch
)
sim_robin = model_robin.run()
print(f"Simulation status: {sim_robin.status_str}")
assert sim_robin.status == ot.Simulation.Status.done
# %%
sim_robin.meshseries.scale(spatial="km", time="d")
sim_robin.meshseries.point_data["hydraulic head h / m"] = (
    sim_robin.meshseries.point_data["pressure"] / (1000 * 9.81)
)
for mesh in sim_robin.meshseries:
    mesh.point_data["depth to groundwater h / m"] = (
        mesh.points[:, 2] * 1000 - mesh.point_data["hydraulic head h / m"]
    )
fig_robin = ot.plot.contourf(
    sim_robin.meshseries[-1], "hydraulic head h / m", fontsize=20
)
fig_robin.axes[0].set_title(
    "Hydraulic Head at simulation end with \n  Robin boundary conditions and heterogeneous aquifer",
    fontsize=24,
    loc="center",
)
fig_robin = ot.plot.contourf(
    sim_robin.meshseries[-1], "depth to groundwater h / m", fontsize=20
)
fig_robin.axes[0].set_title(
    "Depth to groundwater at simulation end with \n  Robin boundary conditions and heterogeneous aquifer",
    fontsize=24,
    loc="center",
)


# %% [markdown]
# ## 3.2 Comparison with Dirichlet + Neumann Model
# We set up and run the previous simulation as a point of comparison. \
# Having done that, we can show the differences both in hydraulic head as well as depth to groundwater between the two simulations

# %%
project_neumann = ot.Project("Selke_Basin.prj").copy()
# Load Meshes
input_meshes_neumann = [
    "Selke_Basin_Domain_homo.vtu",
    "Selke_Basin_PG_Concordia.vtu",
    "Selke_Basin_PG_Koenigsauer.vtu",
    "Selke_Basin_PG_Wilsleber.vtu",
    "Selke_Basin_PL_Bode.vtu",
    "Selke_Basin_PL_Wipper.vtu",
    "Selke_Basin_Pnt_Wells.vtu",
    "Selke_Basin_PL_Selke.vtu",
    "Selke_Basin_PL_Upstream.vtu",
    "Selke_Basin_PG_recharge.vtu",
]
meshes_neumann = ot.Meshes.from_files(input_meshes_neumann)
# Set up and run Dirichlet simulation
model_neumann = ot.Model(
    project=project_neumann, meshes=meshes_neumann, execution=execution_patch
)
sim_neumann = model_neumann.run()
assert sim_neumann.status == ot.Simulation.Status.done
sim_neumann.meshseries.scale(spatial="km", time="d")
sim_neumann.meshseries.point_data["hydraulic head h / m"] = (
    sim_neumann.meshseries.point_data["pressure"] / (1000 * 9.81)
)
for mesh in sim_neumann.meshseries:
    mesh.point_data["depth to groundwater h / m"] = (
        mesh.points[:, 2] * 1000 - mesh.point_data["hydraulic head h / m"]
    )
series_diff = ot.MeshSeries.difference(
    ms_a=sim_robin.meshseries, ms_b=sim_neumann.meshseries
)
fig_difference_head = ot.plot.contourf(
    series_diff[-1],
    "hydraulic head h / m",
    fontsize=20,
)
fig_difference_head.axes[0].set_title(
    "Hydraulic Head difference at end of simulation between \n heterogeneous Robin and Neumann+Dirichlet boundary conditions",
    fontsize=24,
    loc="center",
)
fig_difference_groundwater = ot.plot.contourf(
    series_diff[-1],
    "depth to groundwater h / m",
    fontsize=20,
)
fig_difference_groundwater.axes[0].set_title(
    "Depth to groundwater difference at end of simulation between \n heterogeneous Robin and Neumann+Dirichlet boundary conditions",
    fontsize=24,
    loc="center",
)

# %%
# Check that simulations have produced expected outputs
# The same test as at the end of the previous notebook. Since we had to run the Neumann simulation again for comparison, we check its outputs again.
sim_neumann_expected = pv.read("expected_neumann.vtu")
np.testing.assert_allclose(
    sim_neumann.meshseries[-1]["pressure"], sim_neumann_expected["pressure"], rtol=4e-3
)
# Test of the new simulation
sim_robin_expected = pv.read("expected_robin.vtu")
np.testing.assert_allclose(
    sim_robin.meshseries[-1]["pressure"], sim_robin_expected["pressure"], rtol=4e-3
)


# %% [markdown]
# ## OGS links
# - Project files:
#   - Original Project file: [Selke_Basin.prj](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/Selke_Basin.prj)
#   - XML Patch for Neumann boundary conditions: [neumann_patch.xml](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/neumann_patch.xml)
# - Mesh files:
#   - [Selke_Basin_Domain_2ids.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/Selke_Basin_Domain_2ids.vtu)
#   - [Selke_Basin_PG_Concordia.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/Selke_Basin_PG_Concordia.vtu)
#   - [Selke_Basin_PG_Koenigsauer.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/Selke_Basin_PG_Koenigsauer.vtu)
#   - [Selke_Basin_PG_Wilsleber.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/Selke_Basin_PG_Wilsleber.vtu)
#   - [Selke_Basin_PL_Bode.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/Selke_Basin_PL_Bode.vtu)
#   - [Selke_Basin_PL_Wipper.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/Selke_Basin_PL_Wipper.vtu)
#   - [Selke_Basin_Pnt_Wells.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/Selke_Basin_Pnt_Wells.vtu)
#   - [Selke_Basin_PL_Selke.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/Selke_Basin_PL_Selke.vtu)
#   - [Selke_Basin_PL_Upstream.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/Selke_Basin_PL_Upstream.vtu)
#   - [Selke_Basin_PG_recharge.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/Selke_Basin_PG_recharge.vtu)
#   - [Selke_Basin_Domain_homo.vtu (only for the Dirichlet+Neumann comparison model)](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Parabolic/LiquidFlow/Selke2D/Selke_Basin_Domain_homo.vtu)


# %% [markdown]
# <a id='references'></a>
# ## References
# - Nixdorf, E., Sips, M., & Morstein, P. (2025). Development of a digital workflow for layer reduction strategies and visual-analytical outcome analysis to enhance geo-hydraulic modelling efficiency. Digital Water, 3(1), 1-24. [doi:10.1080/28375807.2025.2460825](https://doi.org/10.1080/28375807.2025.2460825)
