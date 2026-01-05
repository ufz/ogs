// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "LinearElasticTransverseIsotropic.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::unique_ptr<LinearElasticTransverseIsotropic<DisplacementDim>>
createLinearElasticTransverseIsotropic(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config, const bool skip_type_checking);
}  // namespace Solids
}  // namespace MaterialLib
