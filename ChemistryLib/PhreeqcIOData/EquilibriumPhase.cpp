/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EquilibriumPhase.h"
#include "BaseLib/ConfigTreeUtil.h"

namespace ChemistryLib
{
std::vector<EquilibriumPhase> createEquilibriumPhases(
    boost::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
        return {};

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

        equilibrium_phases.emplace_back(
            std::move(name), initial_amount, saturation_index);
    }

    return equilibrium_phases;
}
}  // namespace ChemistryLib
