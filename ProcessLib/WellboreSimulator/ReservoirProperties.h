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
struct ReservoirProperties
{
    ParameterLib::Parameter<double> const& temperature;
    ParameterLib::Parameter<double> const& thermal_conductivity;
    ParameterLib::Parameter<double> const& density;
    ParameterLib::Parameter<double> const& specific_heat_capacity;
    ParameterLib::Parameter<double> const& pressure;
};

ReservoirProperties createReservoirProperties(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);
}  // namespace WellboreSimulator
}  // namespace ProcessLib
