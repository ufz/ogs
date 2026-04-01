// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace BaseLib
{
class ConfigTree;
}
namespace MeshLib
{
class Node;
}
namespace ParameterLib
{
struct ParameterBase;
}
namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
template <typename T_BHE>
T_BHE createBHECoaxial(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::vector<MeshLib::Node*> const& bhe_nodes);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
