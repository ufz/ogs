/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on November 1, 2016, 10:06 AM
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace PorousMedium
{
class CapillaryPressureSaturation;
/// Create a capillary pressure model
/// \param config  ConfigTree object has a tag of `<capillary_pressure>`
std::unique_ptr<CapillaryPressureSaturation> createCapillaryPressureModel(
    BaseLib::ConfigTree const& config);
}
}  // namespace MaterialLib
