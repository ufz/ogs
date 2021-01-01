/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

#include <memory>
#include <utility>

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}  // namespace MaterialLib
namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace ThermoMechanicalPhaseField
{
template <int DisplacementDim>
struct ThermoMechanicalPhaseFieldProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    std::map<int, std::unique_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    ParameterLib::Parameter<double> const& residual_stiffness;
    ParameterLib::Parameter<double> const& crack_resistance;
    ParameterLib::Parameter<double> const& crack_length_scale;
    ParameterLib::Parameter<double> const& kinetic_coefficient;
    ParameterLib::Parameter<double> const& solid_density;
    ParameterLib::Parameter<double> const& linear_thermal_expansion_coefficient;
    ParameterLib::Parameter<double> const& specific_heat_capacity;
    ParameterLib::Parameter<double> const& thermal_conductivity;
    ParameterLib::Parameter<double> const& residual_thermal_conductivity;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double const reference_temperature =
        std::numeric_limits<double>::quiet_NaN();
};

}  // namespace ThermoMechanicalPhaseField
}  // namespace ProcessLib
