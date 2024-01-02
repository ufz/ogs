/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateCubicLawPermeability.h"

#include "BaseLib/ConfigTree.h"
#include "CubicLawPermeability.h"
#include "ParameterLib/Utils.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createCubicLawPermeability(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "CubicLawPermeability");
    DBUG("Create CubicLaw Permeability model");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    auto const& b = ParameterLib::findParameter<double>(
        //! \ogs_file_param{properties__property__CubicLawPermeability__fracture_aperture}
        config.getConfigParameter<std::string>("fracture_aperture"), parameters,
        0, nullptr);

    return std::make_unique<CubicLawPermeability>(std::move(property_name), b);
}
}  // namespace MaterialPropertyLib
