/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
