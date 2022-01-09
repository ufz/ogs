/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 10, 2019
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "Exponential.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Exponential> createExponential(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Exponential");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create Exponential property {:s}.", property_name);
    auto const reference_value =
        //! \ogs_file_param{properties__property__Exponential__reference_value}
        config.getConfigParameter<double>("reference_value");

    auto const& exponent_data_config =
        //! \ogs_file_param{properties__property__Exponential__exponent}
        config.getConfigSubtree("exponent");

    auto const& variable_name =
        //! \ogs_file_param{properties__property__Exponential__exponent__variable_name}
        exponent_data_config.getConfigParameter<std::string>("variable_name");
    auto const reference_condition =
        //! \ogs_file_param{properties__property__Exponential__exponent__reference_condition}
        exponent_data_config.getConfigParameter<double>("reference_condition");
    auto const factor =
        //! \ogs_file_param{properties__property__Exponential__exponent__factor}
        exponent_data_config.getConfigParameter<double>("factor");

    auto const offset =
        //! \ogs_file_param{properties__property__Exponential__offset}
        config.getConfigParameter<double>("offset");

    MaterialPropertyLib::Variable exp_data_type =
        MaterialPropertyLib::convertStringToVariable(variable_name);

    MaterialPropertyLib::ExponentData const exp_data{
        exp_data_type, reference_condition, factor};

    return std::make_unique<MaterialPropertyLib::Exponential>(
        std::move(property_name), offset, reference_value, exp_data);
}
}  // namespace MaterialPropertyLib
