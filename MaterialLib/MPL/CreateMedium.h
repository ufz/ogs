/**
 * \file
 * \author Norbert Grunwald
 * \date   07.09.2017
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
}
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
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);
}  // namespace MaterialPropertyLib
