/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateEquilibriumReactants.h"

#include <optional>

#include "BaseLib/ConfigTree.h"
#include "EquilibriumReactant.h"
#include "MeshLib/Mesh.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::vector<EquilibriumReactant> createEquilibriumReactants(
    std::optional<BaseLib::ConfigTree> const& config, MeshLib::Mesh& mesh)
{
    if (!config)
    {
        return {};
    }

    std::vector<EquilibriumReactant> equilibrium_reactants;
    for (
        auto const& equilibrium_reactant_config :
        //! \ogs_file_param{prj__chemical_system__equilibrium_reactants__phase_component}
        config->getConfigSubtreeList("phase_component"))
    {
        auto name =
            //! \ogs_file_param{prj__chemical_system__equilibrium_reactants__phase_component__name}
            equilibrium_reactant_config.getConfigParameter<std::string>("name");

        double const saturation_index =
            //! \ogs_file_param{prj__chemical_system__equilibrium_reactants__phase_component__saturation_index}
            equilibrium_reactant_config.getConfigParameter<double>(
                "saturation_index");

        auto reaction_irreversibility =
            //! \ogs_file_param{prj__chemical_system__equilibrium_reactants__phase_component__reaction_irreversibility}
            equilibrium_reactant_config.getConfigParameter<std::string>(
                "reaction_irreversibility", "");

        if (!reaction_irreversibility.empty() &&
            (reaction_irreversibility != "dissolve_only" &&
             reaction_irreversibility != "precipitate_only"))
        {
            OGS_FATAL(
                "{:s}: reaction direction only allows `dissolve_only` or "
                "`precipitate_only`",
                name);
        }

        auto molality = MeshLib::getOrCreateMeshProperty<double>(
            mesh, name, MeshLib::MeshItemType::IntegrationPoint, 1);

        auto volume_fraction = MeshLib::getOrCreateMeshProperty<double>(
            mesh, "phi_" + name, MeshLib::MeshItemType::IntegrationPoint, 1);

        auto mesh_prop_molality = MeshLib::getOrCreateMeshProperty<double>(
            mesh, name + "_avg", MeshLib::MeshItemType::Cell, 1);
        mesh_prop_molality->resize(mesh.getNumberOfElements());

        equilibrium_reactants.emplace_back(std::move(name),
                                           molality,
                                           volume_fraction,
                                           mesh_prop_molality,
                                           saturation_index,
                                           std::move(reaction_irreversibility));
    }

    return equilibrium_reactants;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
