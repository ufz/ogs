// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

#include "BaseLib/ConfigTree-fwd.h"
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
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

extern template std::unique_ptr<MechanicsBase<2>> createMFront<2>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
extern template std::unique_ptr<MechanicsBase<3>> createMFront<3>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
