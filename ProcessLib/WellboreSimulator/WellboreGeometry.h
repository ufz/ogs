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
