/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ParameterLib/Utils.h"

#include "DupuitPermeability.h"

namespace BaseLib
{
class ConfigTree;
}

namespace ParameterLib
{
struct ParameterBase;
}

namespace MaterialPropertyLib
{
std::unique_ptr<DupuitPermeability> createDupuitPermeability(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Dupuit");
    DBUG("Create Dupuit permeability.");

    std::string const& parameter_name =
        //! \ogs_file_param{properties__property__DupuitPermeability__parameter_name}
        config.getConfigParameter<std::string>("parameter_name");
    auto const& parameter = ParameterLib::findParameter<double>(
        parameter_name, parameters, 0, nullptr);
    return std::make_unique<MaterialPropertyLib::DupuitPermeability>(parameter);
}
}  // namespace MaterialPropertyLib
