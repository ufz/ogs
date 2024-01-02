/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 17, 2021, 3:47 PM
 */

#include "CreateSoilThermalConductivitySomerton.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createSoilThermalConductivitySomerton(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SoilThermalConductivitySomerton");

    OGS_FATAL(
        "The MPL property SoilThermalConductivitySomerton is "
        "deprecated. Please use SaturationWeightedThermalConductivity "
        "instead.");
}
}  // namespace MaterialPropertyLib
