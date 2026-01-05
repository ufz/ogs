// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateParameter.h"

#include "BaseLib/ConfigTree.h"
#include "Parameter.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Parameter> createParameterProperty(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Parameter");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create Parameter property {:s}.", property_name);

    std::string const& parameter_name =
        //! \ogs_file_param{properties__property__Parameter__parameter_name}
        config.getConfigParameter<std::string>("parameter_name");
    auto const& parameter = ParameterLib::findParameter<double>(
        parameter_name, parameters, 0, nullptr);
    return std::make_unique<MaterialPropertyLib::Parameter>(
        std::move(property_name), parameter);
}
}  // namespace MaterialPropertyLib
