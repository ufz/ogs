// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
template <int GlobalDim>
Eigen::Matrix<double, GlobalDim, GlobalDim> formEffectiveThermalConductivity(
    MaterialPropertyLib::PropertyDataType const& solid_thermal_conductivity,
    const double fluid_thermal_conductivity, const double porosity);
}  // namespace MaterialPropertyLib
