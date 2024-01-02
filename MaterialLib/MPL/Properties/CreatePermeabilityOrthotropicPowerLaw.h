/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}
namespace ParameterLib
{
struct CoordinateSystem;
}
namespace MaterialPropertyLib
{
class Property;
}

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createPermeabilityOrthotropicPowerLaw(
    BaseLib::ConfigTree const& config,
    ParameterLib::CoordinateSystem const* const local_coordinate_system);
}  // namespace MaterialPropertyLib
