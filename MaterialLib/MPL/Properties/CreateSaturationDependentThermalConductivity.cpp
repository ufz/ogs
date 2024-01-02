/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
namespace MaterialPropertyLib
{
std::unique_ptr<Property> createSaturationDependentThermalConductivity(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type",
                                "SaturationDependentThermalConductivity");

    OGS_FATAL(
        "The MPL property SaturationDependentThermalConductivity is "
        "deprecated. Please use SaturationWeightedThermalConductivity "
        "instead.");
}
}  // namespace MaterialPropertyLib
