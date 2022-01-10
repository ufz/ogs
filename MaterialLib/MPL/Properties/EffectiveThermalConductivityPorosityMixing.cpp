/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "EffectiveThermalConductivityPorosityMixing.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
//
// For 1D problems
//
template <>
EffectiveThermalConductivityPorosityMixing<1>::
    EffectiveThermalConductivityPorosityMixing(
        std::string name,
        ParameterLib::CoordinateSystem const* const local_coordinate_system)
    : local_coordinate_system_(local_coordinate_system)
{
    name_ = std::move(name);
}

template <>
void EffectiveThermalConductivityPorosityMixing<1>::checkScale() const
{
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'EffectiveThermalConductivityPorosityMixing' is "
            "implemented on the 'medium' scale only.");
    }
}

template <>
PropertyDataType EffectiveThermalConductivityPorosityMixing<1>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    auto const& medium = std::get<Medium*>(scale_);

    // Assuming there is either a gas phase or a liquid phase or both.
    auto const gas_phase =
        medium->hasPhase("Gas") ? &medium->phase("Gas") : nullptr;
    auto const liquid_phase = medium->hasPhase("AqueousLiquid")
                                  ? &medium->phase("AqueousLiquid")
                                  : nullptr;
    // Assuming there is always a solid phase.
    auto const& solid_phase = medium->phase("Solid");

    auto const gas_thermal_conductivity =
        gas_phase
            ? gas_phase
                  ->property(
                      MaterialPropertyLib::PropertyType::thermal_conductivity)
                  .template value<double>(variable_array, pos, t, dt)
            : 0.;

    auto const liquid_thermal_conductivity =
        liquid_phase
            ? liquid_phase
                  ->property(
                      MaterialPropertyLib::PropertyType::thermal_conductivity)
                  .template value<double>(variable_array, pos, t, dt)
            : 0.;

    auto const solid_thermal_conductivity =
        solid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template value<double>(variable_array, pos, t, dt);

    auto const porosity =
        std::get<double>(variable_array[static_cast<int>(Variable::porosity)]);
    auto const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);
    auto const S_G = 1. - S_L;

    auto const phi_G = porosity * S_G;
    auto const phi_L = porosity * S_L;
    auto const phi_S = 1. - porosity;

    double const effective_thermal_conductivity =
        phi_G * gas_thermal_conductivity + phi_L * liquid_thermal_conductivity +
        phi_S * solid_thermal_conductivity;

    return effective_thermal_conductivity;
}
template <>
PropertyDataType EffectiveThermalConductivityPorosityMixing<1>::dValue(
    VariableArray const&, Variable const, ParameterLib::SpatialPosition const&,
    double const, double const) const
{
    OGS_FATAL(
        "dValue is not implemented for "
        "EffectiveThermalConductivityPorosityMixing");
}
//
// For 2D and 3D problems
//
template <int GlobalDim>
EffectiveThermalConductivityPorosityMixing<GlobalDim>::
    EffectiveThermalConductivityPorosityMixing(
        std::string name,
        ParameterLib::CoordinateSystem const* const local_coordinate_system)
    : local_coordinate_system_(local_coordinate_system)
{
    name_ = std::move(name);
}

template <int GlobalDim>
void EffectiveThermalConductivityPorosityMixing<GlobalDim>::checkScale() const
{
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'EffectiveThermalConductivityPorosityMixing' is "
            "implemented on the 'medium' scale only.");
    }
}

template <int GlobalDim>
PropertyDataType EffectiveThermalConductivityPorosityMixing<GlobalDim>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    auto const& medium = std::get<Medium*>(scale_);
    // Assuming there is either a gas phase or a liquid phase or both.
    auto const gas_phase =
        medium->hasPhase("Gas") ? &medium->phase("Gas") : nullptr;
    auto const liquid_phase = medium->hasPhase("AqueousLiquid")
                                  ? &medium->phase("AqueousLiquid")
                                  : nullptr;
    // Assuming there is always a solid phase.
    auto const& solid_phase = medium->phase("Solid");

    auto const gas_thermal_conductivity =
        gas_phase
            ? gas_phase
                  ->property(
                      MaterialPropertyLib::PropertyType::thermal_conductivity)
                  .template value<double>(variable_array, pos, t, dt)
            : 0.;

    auto const liquid_thermal_conductivity =
        liquid_phase
            ? liquid_phase
                  ->property(
                      MaterialPropertyLib::PropertyType::thermal_conductivity)
                  .template value<double>(variable_array, pos, t, dt)
            : 0.;

    auto solid_thermal_conductivity = formEigenTensor<GlobalDim>(
        solid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, t, dt));

    auto const porosity =
        std::get<double>(variable_array[static_cast<int>(Variable::porosity)]);

    auto const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);
    auto const S_G = 1. - S_L;

    auto const phi_G = porosity * S_G;
    auto const phi_L = porosity * S_L;
    auto const phi_S = 1. - porosity;

    // Local coordinate transformation is only applied for the case that the
    // initial solid thermal conductivity is given with orthotropic assumption.
    if (local_coordinate_system_ &&
        (solid_thermal_conductivity.cols() == GlobalDim))
    {
        Eigen::Matrix<double, GlobalDim, GlobalDim> const e =
            local_coordinate_system_->transformation<GlobalDim>(pos);

        solid_thermal_conductivity =
            e.transpose() * solid_thermal_conductivity * e;
    }
    auto const I = Eigen::Matrix<double, GlobalDim, GlobalDim>::Identity();
    Eigen::Matrix<double, GlobalDim, GlobalDim> const
        effective_thermal_conductivity = (phi_G * gas_thermal_conductivity +
                                          phi_L * liquid_thermal_conductivity) *
                                             I +
                                         phi_S * solid_thermal_conductivity;
    return effective_thermal_conductivity;
}

template <int GlobalDim>
PropertyDataType EffectiveThermalConductivityPorosityMixing<GlobalDim>::dValue(
    VariableArray const&, Variable const, ParameterLib::SpatialPosition const&,
    double const, double const) const
{
    OGS_FATAL(
        "dValue is not implemented for "
        "EffectiveThermalConductivityPorosityMixing");
}
template class EffectiveThermalConductivityPorosityMixing<2>;
template class EffectiveThermalConductivityPorosityMixing<3>;
}  // namespace MaterialPropertyLib
