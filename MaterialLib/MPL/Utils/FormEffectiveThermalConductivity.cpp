/*
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 31, 2019, 12:10 PM
 */

#include "FormEffectiveThermalConductivity.h"

#include "FormEigenTensor.h"

namespace MaterialPropertyLib
{
template <int GlobalDim>
Eigen::Matrix<double, GlobalDim, GlobalDim> formEffectiveThermalConductivity(
    MaterialPropertyLib::PropertyDataType const& solid_thermal_conductivity,
    const double fluid_thermal_conductivity, const double porosity)
{
    return (1.0 - porosity) *
               formEigenTensor<GlobalDim>(solid_thermal_conductivity) +
           porosity * fluid_thermal_conductivity *
               Eigen::Matrix<double, GlobalDim, GlobalDim>::Identity();
}

template Eigen::Matrix<double, 1, 1> formEffectiveThermalConductivity<1>(
    MaterialPropertyLib::PropertyDataType const& solid_thermal_conductivity,
    const double fluid_thermal_conductivity, const double porosity);
template Eigen::Matrix<double, 2, 2> formEffectiveThermalConductivity<2>(
    MaterialPropertyLib::PropertyDataType const& solid_thermal_conductivity,
    const double fluid_thermal_conductivity, const double porosity);
template Eigen::Matrix<double, 3, 3> formEffectiveThermalConductivity<3>(
    MaterialPropertyLib::PropertyDataType const& solid_thermal_conductivity,
    const double fluid_thermal_conductivity, const double porosity);

}  // namespace MaterialPropertyLib
