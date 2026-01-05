// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
