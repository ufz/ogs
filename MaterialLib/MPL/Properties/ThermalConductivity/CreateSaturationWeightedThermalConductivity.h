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
struct ParameterBase;
}  // namespace ParameterLib

namespace MaterialPropertyLib
{
class Property;

std::unique_ptr<Property> createSaturationWeightedThermalConductivity(
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters);
}  // namespace MaterialPropertyLib
