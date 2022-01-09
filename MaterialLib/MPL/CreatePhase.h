/**
 * \file
 * \author Norbert Grunwald
 * \date   07.09.2017
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}
namespace ParameterLib
{
struct ParameterBase;
struct CoordinateSystem;
}
namespace MaterialPropertyLib
{
class Phase;
}
namespace MathLib
{
class PiecewiseLinearInterpolation;
}

namespace MaterialPropertyLib
{
/// A method that parses the phase details and stores them in the private
/// phases_ member.
///
/// This method creates the phases of the medium. Unlike a medium, a phase may
/// have a name. However, this is silly at the moment since this name still has
/// no effect (except of some benefits in regard of readability).
/// Phase components are required (a phase consists of at least one component).
/// Phase properties are optional. If not given, default properties are
/// assigned. These default properties average the component properties,
/// weighted by mole fraction.
std::vector<std::unique_ptr<Phase>> createPhases(
    int const geometry_dimension,
    std::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);
}  // namespace MaterialPropertyLib
