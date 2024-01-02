/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <memory>
#include <string>

#include "BHECommonCoaxial.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

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
T_BHE createBHECoaxial(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
