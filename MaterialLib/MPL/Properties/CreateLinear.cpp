/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 10, 2019
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "Linear.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Linear> createLinear(
    BaseLib::ConfigTree const& config)
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
    for (
        auto const& independent_variable_config :
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

        MaterialPropertyLib::Variable ivt =
            MaterialPropertyLib::convertStringToVariable(variable_name);

        MaterialPropertyLib::IndependentVariable iv{ivt, reference_condition,
                                                    slope};

        ivs.push_back(std::move(iv));
    }

    return std::make_unique<MaterialPropertyLib::Linear>(
        std::move(property_name), reference_value, ivs);
}
}  // namespace MaterialPropertyLib
