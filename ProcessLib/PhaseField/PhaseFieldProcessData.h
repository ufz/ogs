/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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

namespace PhaseField
{
template <int DisplacementDim>
struct PhaseFieldProcessData
{
    PhaseFieldProcessData(
        MeshLib::PropertyVector<int> const* const material_ids_,
        std::map<int,
                 std::unique_ptr<
                     MaterialLib::Solids::MechanicsBase<DisplacementDim>>>&&
            solid_materials_,
        Parameter<double> const& residual_stiffness_,
        Parameter<double> const& crack_resistance_,
        Parameter<double> const& crack_length_scale_,
        Parameter<double> const& kinetic_coefficient_,
        Parameter<double> const& solid_density_,
        Parameter<double>& history_field_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_,
        bool const propagating_crack_,
        bool const crack_pressure_)
        : material_ids(material_ids_),
          solid_materials{std::move(solid_materials_)},
          residual_stiffness(residual_stiffness_),
          crack_resistance(crack_resistance_),
          crack_length_scale(crack_length_scale_),
          kinetic_coefficient(kinetic_coefficient_),
          solid_density(solid_density_),
          history_field(history_field_),
          specific_body_force(specific_body_force_),
          propagating_crack(propagating_crack_),
          crack_pressure(crack_pressure_)
    {
    }

    PhaseFieldProcessData(PhaseFieldProcessData&& other) = default;

    //! Copies are forbidden.
    PhaseFieldProcessData(PhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(PhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(PhaseFieldProcessData&&) = delete;

    MeshLib::PropertyVector<int> const* const material_ids;

    std::map<int, std::unique_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    Parameter<double> const& residual_stiffness;
    Parameter<double> const& crack_resistance;
    Parameter<double> const& crack_length_scale;
    Parameter<double> const& kinetic_coefficient;
    Parameter<double> const& solid_density;
    Parameter<double>& history_field;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt = 0.0;
    double t = 0.0;
    double const unity_pressure = 1.0;
    double pressure = 0.0;
    double pressure_old = 0.0;
    double pressure_error = 0.0;
    double injected_volume = 0.0;
    double crack_volume = 0.0;
    double elastic_energy = 0.0;
    double surface_energy = 0.0;
    double pressure_work = 0.0;
    bool propagating_crack = false;
    bool crack_pressure = false;
};

}  // namespace PhaseField
}  // namespace ProcessLib
