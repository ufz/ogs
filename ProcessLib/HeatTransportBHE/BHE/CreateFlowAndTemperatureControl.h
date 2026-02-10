// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <memory>
#include <string>

#include "FlowAndTemperatureControl.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{
class PiecewiseLinearInterpolation;
}

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct RefrigerantProperties;

FlowAndTemperatureControl createFlowAndTemperatureControl(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    RefrigerantProperties const& refrigerant);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
