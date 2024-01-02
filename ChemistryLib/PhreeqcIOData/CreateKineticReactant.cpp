/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateKineticReactant.h"

#include <optional>

#include "BaseLib/ConfigTree.h"
#include "KineticReactant.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::vector<KineticReactant> createKineticReactants(
    std::optional<BaseLib::ConfigTree> const& config, MeshLib::Mesh& mesh)
{
    if (!config)
    {
        return {};
    }

    std::vector<KineticReactant> kinetic_reactants;
    for (
        auto const& reactant_config :
        //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant}
        config->getConfigSubtreeList("kinetic_reactant"))
    {
        //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__name}
        auto name = reactant_config.getConfigParameter<std::string>("name");

        auto chemical_formula =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__chemical_formula}
            reactant_config.getConfigParameter<std::string>("chemical_formula",
                                                            "");

        auto parameters =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__parameters}
            reactant_config.getConfigParameter<std::vector<double>>(
                "parameters", {});

        bool const fix_amount =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__fix_amount}
            reactant_config.getConfigParameter<bool>("fix_amount", false);

        auto molality = MeshLib::getOrCreateMeshProperty<double>(
            mesh, name, MeshLib::MeshItemType::IntegrationPoint, 1);

        auto molality_prev = MeshLib::getOrCreateMeshProperty<double>(
            mesh, name + "_prev", MeshLib::MeshItemType::IntegrationPoint, 1);

        auto volume_fraction = MeshLib::getOrCreateMeshProperty<double>(
            mesh, "phi_" + name, MeshLib::MeshItemType::IntegrationPoint, 1);

        auto volume_fraction_prev = MeshLib::getOrCreateMeshProperty<double>(
            mesh,
            "phi_" + name + "_prev",
            MeshLib::MeshItemType::IntegrationPoint,
            1);

        auto mesh_prop_molality = MeshLib::getOrCreateMeshProperty<double>(
            mesh, name + "_avg", MeshLib::MeshItemType::Cell, 1);
        mesh_prop_molality->resize(mesh.getNumberOfElements());

        kinetic_reactants.emplace_back(std::move(name),
                                       std::move(chemical_formula),
                                       molality,
                                       molality_prev,
                                       volume_fraction,
                                       volume_fraction_prev,
                                       mesh_prop_molality,
                                       std::move(parameters),
                                       fix_amount);
    }

    return kinetic_reactants;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
