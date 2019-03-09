/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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
    ThermoMechanicalPhaseFieldProcessData(
        MeshLib::PropertyVector<int> const* const material_ids_,
        std::map<int,
                 std::unique_ptr<
                     MaterialLib::Solids::MechanicsBase<DisplacementDim>>>&&
            solid_materials_,
        ParameterLib::Parameter<double> const& residual_stiffness_,
        ParameterLib::Parameter<double> const& crack_resistance_,
        ParameterLib::Parameter<double> const& crack_length_scale_,
        ParameterLib::Parameter<double> const& kinetic_coefficient_,
        ParameterLib::Parameter<double> const& solid_density_,
        ParameterLib::Parameter<double> const&
            linear_thermal_expansion_coefficient_,
        ParameterLib::Parameter<double> const& specific_heat_capacity_,
        ParameterLib::Parameter<double> const& thermal_conductivity_,
        ParameterLib::Parameter<double> const& residual_thermal_conductivity_,
        double const reference_temperature_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_)
        : material_ids(material_ids_),
          solid_materials{std::move(solid_materials_)},
          residual_stiffness(residual_stiffness_),
          crack_resistance(crack_resistance_),
          crack_length_scale(crack_length_scale_),
          kinetic_coefficient(kinetic_coefficient_),
          solid_density(solid_density_),
          linear_thermal_expansion_coefficient(
              linear_thermal_expansion_coefficient_),
          specific_heat_capacity(specific_heat_capacity_),
          thermal_conductivity(thermal_conductivity_),
          residual_thermal_conductivity(residual_thermal_conductivity_),
          reference_temperature(reference_temperature_),
          specific_body_force(specific_body_force_)
    {
    }

    ThermoMechanicalPhaseFieldProcessData(
        ThermoMechanicalPhaseFieldProcessData&& other) = default;

    //! Copies are forbidden.
    ThermoMechanicalPhaseFieldProcessData(
        ThermoMechanicalPhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ThermoMechanicalPhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ThermoMechanicalPhaseFieldProcessData&&) = delete;

    MeshLib::PropertyVector<int> const* const material_ids;

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
    double const reference_temperature;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt;
    double t;
};

}  // namespace ThermoMechanicalPhaseField
}  // namespace ProcessLib
