/*
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 31, 2019, 12:10 PM
 */

#pragma once

#include <Eigen/Dense>

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
template <int GlobalDim>
Eigen::Matrix<double, GlobalDim, GlobalDim> formEffectiveThermalConductivity(
    MaterialPropertyLib::PropertyDataType const& solid_thermal_conductivity,
    const double fluid_thermal_conductivity, const double porosity);
}  // namespace MaterialPropertyLib
