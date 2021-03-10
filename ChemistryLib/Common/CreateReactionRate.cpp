/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <optional>

#include "BaseLib/ConfigTree.h"
#include "ChemistryLib/PhreeqcIOData/ReactionRate.h"
#include "ChemistryLib/PhreeqcKernelData/ReactionRate.h"

namespace ChemistryLib
{
template <typename ReactionRate>
std::vector<ReactionRate> createReactionRates(
    std::optional<BaseLib::ConfigTree> const& config)
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

        auto const expression_config =
            //! \ogs_file_param{prj__chemical_system__rates__rate__expression}
            rate_config.getConfigSubtree("expression");
        auto const& statements =
            //! \ogs_file_param{prj__chemical_system__rates__rate__expression__statement}
            expression_config.getConfigParameterList<std::string>("statement");

        std::vector<std::string> expression_statements;
        expression_statements.reserve(statements.size());
        std::copy(begin(statements),
                  end(statements),
                  back_inserter(expression_statements));

        reaction_rates.emplace_back(std::move(kinetic_reactant),
                                    std::move(expression_statements));
    }

    return reaction_rates;
}

template std::vector<PhreeqcIOData::ReactionRate>
createReactionRates<PhreeqcIOData::ReactionRate>(
    std::optional<BaseLib::ConfigTree> const& config);

template std::vector<PhreeqcKernelData::ReactionRate>
createReactionRates<PhreeqcKernelData::ReactionRate>(
    std::optional<BaseLib::ConfigTree> const& config);
}  // namespace ChemistryLib
