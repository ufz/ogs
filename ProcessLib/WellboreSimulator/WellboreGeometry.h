/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}
namespace ParameterLib
{
template <typename T>
struct Parameter;
struct ParameterBase;
}  // namespace ParameterLib

namespace ProcessLib
{
namespace WellboreSimulator
{
struct WellboreGeometry
{
    double length;
    ParameterLib::Parameter<double> const& diameter;
    ParameterLib::Parameter<double> const& casing_thickness;
    ParameterLib::Parameter<double> const& pipe_thickness;
    ParameterLib::Parameter<double> const& roughness;
};

WellboreGeometry createWellboreGeometry(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);
}  // namespace WellboreSimulator
}  // namespace ProcessLib
