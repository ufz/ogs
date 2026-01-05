// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <unordered_set>

#include "BaseLib/ConfigTree.h"
#include "Linear.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Linear> createLinear(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Linear");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create Linear property {:s}.", property_name);
    auto const reference_value =
        //! \ogs_file_param{properties__property__Linear__reference_value}
        config.getConfigParameter<double>("reference_value");

    std::vector<MaterialPropertyLib::IndependentVariable> ivs;
    for (auto const& independent_variable_config :
         //! \ogs_file_param{properties__property__Linear__independent_variable}
         config.getConfigSubtreeList("independent_variable"))
    {
        auto const& variable_name =
            //! \ogs_file_param{properties__property__Linear__independent_variable__variable_name}
            independent_variable_config.getConfigParameter<std::string>(
                "variable_name");
        auto const reference_condition =
            //! \ogs_file_param{properties__property__Linear__independent_variable__reference_condition}
            independent_variable_config.getConfigParameter<double>(
                "reference_condition");
        auto const slope =
            //! \ogs_file_param{properties__property__Linear__independent_variable__slope}
            independent_variable_config.getConfigParameter<double>("slope");

        auto const min =
            //! \ogs_file_param{properties__property__Linear__independent_variable__min}
            independent_variable_config.getConfigParameterOptional<double>(
                "min");

        auto const max =
            //! \ogs_file_param{properties__property__Linear__independent_variable__max}
            independent_variable_config.getConfigParameterOptional<double>(
                "max");

        static const std::unordered_set<std::string> filter_not_variables = {
            "t", "x", "y", "z"};
        MaterialPropertyLib::StringOrVariable ivt;
        if (!filter_not_variables.contains(variable_name))
        {
            ivt = MaterialPropertyLib::convertStringToVariable(variable_name);
        }
        else
        {
            ivt = variable_name;
        }

        MaterialPropertyLib::IndependentVariable iv{ivt, reference_condition,
                                                    slope, min, max};

        ivs.push_back(std::move(iv));
    }

    return std::make_unique<MaterialPropertyLib::Linear>(
        std::move(property_name), reference_value, ivs);
}
}  // namespace MaterialPropertyLib
