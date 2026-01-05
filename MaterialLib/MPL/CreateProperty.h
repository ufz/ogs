// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <memory>
#include <optional>
#include <vector>

#include "PropertyType.h"

namespace BaseLib
{
class ConfigTree;
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
class Property;

/// This data type is based on a std::array. It can hold pointers to objects of
/// class Property or its inheritors. The size of this array is determined by
/// the number of entries of the PropertyType enumerator.
using PropertyArray =
    std::array<std::unique_ptr<Property>, PropertyType::number_of_properties>;

/// The method reads the 'properties' tag in the prj-file and creates component
/// properties accordingly.
///
/// First, a new property iy created based on the specified property type.
/// Then, the property name is evaluated and the property is copied into the
/// properties array.
std::unique_ptr<PropertyArray> createProperties(
    int const geometry_dimension,
    std::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);

}  // namespace MaterialPropertyLib
