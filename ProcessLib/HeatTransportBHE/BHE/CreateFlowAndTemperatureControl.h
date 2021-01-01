/**
 * \file
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
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    RefrigerantProperties const& refrigerant);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
