/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "ReactionRate.h"

namespace ChemistryLib
{
std::vector<ReactionRate> createReactionRates(
    boost::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
    {
        return {};
    }

    std::vector<ReactionRate> reaction_rates;
    for (auto const& rate_config :
         //! \ogs_file_param{prj__chemical_system__rates__rate}
         config->getConfigSubtreeList("rate"))
    {
        auto kinetic_reactant =
            //! \ogs_file_param{prj__chemical_system__rates__rate__kinetic_reactant}
            rate_config.getConfigParameter<std::string>("kinetic_reactant");

        std::vector<std::string> expression_statements;
        auto const expression_config =
            //! \ogs_file_param{prj__chemical_system__rates__rate__expression}
            rate_config.getConfigSubtree("expression");
        for (
            auto const& expression_statement :
            //! \ogs_file_param{prj__chemical_system__rates__rate__expression__statement}
            expression_config.getConfigParameterList<std::string>("statement"))
        {
            expression_statements.push_back(expression_statement);
        }

        reaction_rates.emplace_back(std::move(kinetic_reactant),
                                    std::move(expression_statements));
    }

    return reaction_rates;
}
}  // namespace ChemistryLib