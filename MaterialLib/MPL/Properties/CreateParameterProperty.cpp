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
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"

#include "ParameterProperty.h"

namespace MaterialPropertyLib
{
std::unique_ptr<ParameterProperty> createParameterProperty(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    config.checkConfigParameter("type", "Parameter");
    DBUG("Create Parameter property");

    std::string const& parameter_name =
        //! \ogs_file_param{properties__property__Parameter__parameter_name}
        config.getConfigParameter<std::string>("parameter_name");
    auto const& parameter = ParameterLib::findParameter<double>(
        parameter_name, parameters, 0, nullptr);
    return std::make_unique<MaterialPropertyLib::ParameterProperty>(parameter);
}
}  // namespace MaterialPropertyLib