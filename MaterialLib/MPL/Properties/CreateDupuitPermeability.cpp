/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateDupuitPermeability.h"

#include "BaseLib/ConfigTree.h"
#include "DupuitPermeability.h"
#include "ParameterLib/Utils.h"

namespace MaterialPropertyLib
{
std::unique_ptr<DupuitPermeability> createDupuitPermeability(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Dupuit");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create DupuitPermeability property {:s}.", property_name);

    std::string const& parameter_name =
        //! \ogs_file_param{properties__property__DupuitPermeability__parameter_name}
        config.getConfigParameter<std::string>("parameter_name");
    auto const& parameter = ParameterLib::findParameter<double>(
        parameter_name, parameters, 0, nullptr);
    return std::make_unique<MaterialPropertyLib::DupuitPermeability>(
        std::move(property_name), parameter);
}
}  // namespace MaterialPropertyLib
