/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   DeactivatedSubdomain.h
 *
 * Created on November 29, 2018, 10:50 AM
 */
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
std::vector<std::unique_ptr<DeactivatedSubdomain const>>
createDeactivatedSubdomains(
    BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);

}  // namespace ProcessLib
