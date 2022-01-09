/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "Function.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Function> createFunction(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Function");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create Function property {:s}.", property_name);

    std::vector<std::string> value_expressions;
    //! \ogs_file_param{properties__property__Function__value}
    auto const& value_config = config.getConfigSubtree("value");

    //! \ogs_file_param{properties__property__Function__value__expression}
    for (auto const& p : value_config.getConfigSubtreeList("expression"))
    {
        value_expressions.emplace_back(p.getValue<std::string>());
    }

    // For each derivative a name of the variable and the list of expressions.
    std::vector<std::pair<std::string, std::vector<std::string>>>
        dvalue_expressions;
    //! \ogs_file_param{properties__property__Function__dvalue}
    for (auto const& dvalue_config : config.getConfigSubtreeList("dvalue"))
    {
        auto variable_name =
            //! \ogs_file_param{properties__property__Function__dvalue__variable_name}
            dvalue_config.getConfigParameter<std::string>("variable_name");

        std::vector<std::string> expressions;
        //! \ogs_file_param{properties__property__Function__dvalue__expression}
        for (auto const& p : dvalue_config.getConfigSubtreeList("expression"))
        {
            expressions.emplace_back(p.getValue<std::string>());
        }
        dvalue_expressions.emplace_back(std::move(variable_name),
                                        std::move(expressions));
    }

    return std::make_unique<MaterialPropertyLib::Function>(
        std::move(property_name), value_expressions, dvalue_expressions);
}
}  // namespace MaterialPropertyLib
