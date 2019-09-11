/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 10, 2019
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "ExponentialProperty.h"

namespace MaterialPropertyLib
{
std::unique_ptr<ExponentialProperty> createExponentialProperty(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "Exponential");
    DBUG("Create Exponential property");
    auto const reference_value =
        //! \ogs_file_param{properties__property__ExponentialProperty__reference_value}
        config.getConfigParameter<double>("reference_value");

    auto const& exponent_data_config =
        //! \ogs_file_param{properties__property__ExponentialProperty__exponent}
        config.getConfigSubtree("exponent");

    auto const& variable_name =
        //! \ogs_file_param{properties__property__ExponentialProperty__exponent__variable_name}
        exponent_data_config.getConfigParameter<std::string>("variable_name");
    auto const reference_condition =
        //! \ogs_file_param{properties__property__ExponentialProperty__exponent__reference_condition}
        exponent_data_config.getConfigParameter<double>("reference_condition");
    auto const factor =
        //! \ogs_file_param{properties__property__ExponentialProperty__exponent__factor}
        exponent_data_config.getConfigParameter<double>("factor");

    MaterialPropertyLib::Variable exp_data_type =
        MaterialPropertyLib::convertStringToVariable(variable_name);

    MaterialPropertyLib::ExponentData const exp_data{
        exp_data_type, reference_condition, factor};

    return std::make_unique<MaterialPropertyLib::ExponentialProperty>(
        reference_value, exp_data);
}
}  // namespace MaterialPropertyLib