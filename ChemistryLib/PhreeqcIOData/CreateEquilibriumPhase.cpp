/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <boost/optional/optional.hpp>

#include "BaseLib/ConfigTree.h"
#include "CreateEquilibriumPhase.h"
#include "EquilibriumPhase.h"
#include "MeshLib/Mesh.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::vector<EquilibriumPhase> createEquilibriumPhases(
    boost::optional<BaseLib::ConfigTree> const& config,
    MeshLib::Mesh const& mesh)
{
    if (!config)
    {
        return {};
    }

    std::vector<EquilibriumPhase> equilibrium_phases;
    for (
        auto const& equilibrium_phase_config :
        //! \ogs_file_param{prj__chemical_system__equilibrium_phases__equilibrium_phase}
        config->getConfigSubtreeList("equilibrium_phase"))
    {
        auto name =
            //! \ogs_file_param{prj__chemical_system__equilibrium_phases__equilibrium_phase__name}
            equilibrium_phase_config.getConfigParameter<std::string>("name");

        double const initial_amount =
            //! \ogs_file_param{prj__chemical_system__equilibrium_phases__equilibrium_phase__initial_amount}
            equilibrium_phase_config.getConfigParameter<double>(
                "initial_amount");

        double const saturation_index =
            //! \ogs_file_param{prj__chemical_system__equilibrium_phases__equilibrium_phase__saturation_index}
            equilibrium_phase_config.getConfigParameter<double>(
                "saturation_index");

        auto amount = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh),
            name,
            MeshLib::MeshItemType::Node,
            1);
        std::fill(amount->begin(), amount->end(), initial_amount);

        equilibrium_phases.emplace_back(
            std::move(name), amount, saturation_index);
    }

    return equilibrium_phases;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
