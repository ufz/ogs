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
#include "Exchange.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::vector<Exchange> createExchange(
    std::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
    {
        return {};
    }

    std::vector<Exchange> exchange;
    for (auto const& site_config :
         //! \ogs_file_param{prj__chemical_system__exchange_assemblage}
         config->getConfigSubtreeList("exchange_assemblage"))
    {
        //! \ogs_file_param{prj__chemical_system__exchange_assemblage__name}
        auto name = site_config.getConfigParameter<std::string>("name");

        auto const molality =
            //! \ogs_file_param{prj__chemical_system__exchange_assemblage__molality}
            site_config.getConfigParameter<double>("molality");

        exchange.emplace_back(
            std::move(name), molality);
    }

    return exchange;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
