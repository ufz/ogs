// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "FractureModelBase.h"

namespace MaterialLib
{
namespace Fracture
{
template <int DisplacementDim>
std::unique_ptr<FractureModelBase<DisplacementDim>> createCoulomb(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

}  // namespace Fracture
}  // namespace MaterialLib
