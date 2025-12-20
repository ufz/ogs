// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
