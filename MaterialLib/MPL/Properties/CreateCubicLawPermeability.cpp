// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateCubicLawPermeability.h"

#include <limits>
#include <string>

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

    auto const fracture_aperture_name =
        //! \ogs_file_param{properties__property__CubicLawPermeability__fracture_aperture}
        config.getConfigParameter<std::string>("fracture_aperture", "");

    ParameterLib::Parameter<double>* b =
        fracture_aperture_name == ""
            ? nullptr
            : &ParameterLib::findParameter<double>(
                  fracture_aperture_name, parameters, 0, nullptr);

    return std::make_unique<CubicLawPermeability>(std::move(property_name), b);
}
}  // namespace MaterialPropertyLib
