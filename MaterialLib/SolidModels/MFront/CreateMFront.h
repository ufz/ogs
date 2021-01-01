/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
namespace MFront
{
template <int DisplacementDim>
std::unique_ptr<MechanicsBase<DisplacementDim>> createMFront(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

extern template std::unique_ptr<MechanicsBase<2>> createMFront<2>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
extern template std::unique_ptr<MechanicsBase<3>> createMFront<3>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
