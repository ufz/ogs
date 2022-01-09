/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateKineticReactant.h"

#include <optional>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "KineticReactant.h"
#include "MeshLib/Mesh.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
std::unique_ptr<Kinetics> createKineticReactants(
    std::optional<BaseLib::ConfigTree> const& config, MeshLib::Mesh const& mesh)
{
    if (!config)
    {
        return nullptr;
    }

    std::vector<KineticReactant> kinetic_reactants;
    for (
        auto const& reactant_config :
        //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant}
        config->getConfigSubtreeList("kinetic_reactant"))
    {
        //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__name}
        auto name = reactant_config.getConfigParameter<std::string>("name");

        double const initial_amount =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__initial_amount}
            reactant_config.getConfigParameter<double>("initial_amount");

        auto amount = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh),
            name,
            MeshLib::MeshItemType::Node,
            1);
        std::fill(std::begin(*amount), std::end(*amount), initial_amount);

        kinetic_reactants.emplace_back(name, initial_amount);
    }

    return std::make_unique<Kinetics>(kinetic_reactants);
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
