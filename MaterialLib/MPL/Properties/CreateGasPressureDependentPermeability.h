// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ParameterLib
{
struct CoordinateSystem;
struct ParameterBase;
}  // namespace ParameterLib

namespace MaterialPropertyLib
{
class Property;
}

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createGasPressureDependentPermeability(
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system);
}  // namespace MaterialPropertyLib
