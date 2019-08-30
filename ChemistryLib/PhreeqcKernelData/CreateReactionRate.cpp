/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <stdlib.h>
#include <boost/optional/optional.hpp>

#include "BaseLib/ConfigTree.h"
#include "CreateReactionRate.h"

#include "ThirdParty/iphreeqc/src/src/phreeqcpp/global_structures.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
std::tuple<rate*, int> createReactionRates(
    boost::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
    {
        return {};
    }

    int count_rates = 0;
    struct rate* rates = (struct rate*)malloc(sizeof(struct rate));
    for (auto const& rate_config :
         //! \ogs_file_param{prj__chemical_system__rates__rate}
         config->getConfigSubtreeList("rate"))
    {
        if (count_rates > 0)
            rates = (struct rate*)realloc(
                rates, (std::size_t)(count_rates + 1) * sizeof(struct rate));

        std::string* kinetic_reactant = new std::string(
            //! \ogs_file_param{prj__chemical_system__rates__rate__kinetic_reactant}
            rate_config.getConfigParameter<std::string>("kinetic_reactant"));
        rates[count_rates].name = kinetic_reactant->c_str();

        auto const expression_config =
            //! \ogs_file_param{prj__chemical_system__rates__rate__expression}
            rate_config.getConfigSubtree("expression");

        std::string statements;
        int line_number = 1;
        for (
            auto const& statement :
            //! \ogs_file_param{prj__chemical_system__rates__rate__expression__statement}
            expression_config.getConfigParameterList<std::string>("statement"))
        {
            statements += std::to_string(line_number) + " " + statement + "; ";
            ++line_number;
        }

        rates[count_rates].commands = new char[statements.size() + 1];
        std::copy(
            statements.begin(), statements.end(), rates[count_rates].commands);
        rates[count_rates].commands[statements.size()] = '\0';

        rates[count_rates].new_def = 1;
        ++count_rates;
    }

    return std::make_tuple(rates, count_rates);
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
