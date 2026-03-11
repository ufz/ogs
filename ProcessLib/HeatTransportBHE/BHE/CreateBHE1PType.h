// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "ParameterLib/Parameter.h"

namespace BaseLib
{
class ConfigTree;
}
namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
template <typename T_BHE>
T_BHE createBHE1PType(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
