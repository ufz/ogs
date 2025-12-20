// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
