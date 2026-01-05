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
}  // namespace BaseLib

namespace MeshLib
{
class Mesh;
class Node;
}  // namespace MeshLib

namespace ParameterLib
{
struct ParameterBase;
}

namespace ProcessLib
{
struct DeactivatedSubdomain;
}
namespace ProcessLib
{
std::vector<DeactivatedSubdomain> createDeactivatedSubdomains(
    BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);

}  // namespace ProcessLib
