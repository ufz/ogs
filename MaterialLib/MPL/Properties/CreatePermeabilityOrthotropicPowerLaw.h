// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
