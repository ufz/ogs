// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}
namespace MaterialPropertyLib
{
class Medium;
}
namespace ParameterLib
{
struct ParameterBase;
struct CoordinateSystem;
}  // namespace ParameterLib
namespace MathLib
{
class PiecewiseLinearInterpolation;
}

namespace MaterialPropertyLib
{
/// This function parses the "phases" and "properties" subtrees of the config
/// tree and calls create methods for the phase vector and the properties array.
/// Medium properties are optional. If not defined, default properties are
/// assigned.
std::unique_ptr<Medium> createMedium(
    int const material_id,
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);
}  // namespace MaterialPropertyLib
