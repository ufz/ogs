/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>

#include "KineticReactant.h"
#include "BaseLib/ConfigTreeUtil.h"

namespace ChemistryLib
{
std::vector<KineticReactant> createKineticReactants(
    boost::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
        return {};

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

        auto parameters =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__parameters}
            reactant_config.getConfigParameterOptional<std::vector<double>>(
                "parameters");

        kinetic_reactants.emplace_back(
            std::move(name), initial_amount, std::move(parameters));
    }

    return kinetic_reactants;
}

std::ofstream& operator<<(std::ofstream& out,
                          std::vector<KineticReactant> const& kinetic_reactants)
{
    for (auto const& kinetic_reactant : kinetic_reactants)
    {
        out << kinetic_reactant.name << "\n";

        out << "-m  " << kinetic_reactant.amount << "\n";

        if (!kinetic_reactant.parameters.empty())
        {
            out << "-parms";
            for (auto const& parameter : kinetic_reactant.parameters)
                out << " " << parameter;
            out << "\n";
        }
    }

    return out;
}
}  // namespace ChemistryLib
