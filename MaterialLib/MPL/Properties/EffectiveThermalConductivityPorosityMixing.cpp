// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "EffectiveThermalConductivityPorosityMixing.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/PropertyType.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ParameterLib/CoordinateSystem.h"

namespace MaterialPropertyLib
{
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

    auto const porosity = variable_array.porosity;

    auto const S_L = variable_array.liquid_saturation;
    auto const S_G = 1. - S_L;

    auto const phi_G = porosity * S_G;
    auto const phi_L = porosity * S_L;
    auto const phi_S = 1. - porosity;

    if constexpr (GlobalDim == 1)
    {
        // For 1D, return scalar value
        double const effective_thermal_conductivity =
            phi_G * gas_thermal_conductivity +
            phi_L * liquid_thermal_conductivity +
            phi_S * solid_thermal_conductivity[0];

        return effective_thermal_conductivity;
    }
    else
    {
        // For 2D/3D, use tensor operations

        // Local coordinate transformation is only applied for the case that the
        // initial solid thermal conductivity is given with orthotropic
        // assumption.
        if (local_coordinate_system_)
        {
            solid_thermal_conductivity =
                local_coordinate_system_->rotateTensor<GlobalDim>(
                    solid_thermal_conductivity, pos);
        }
        auto const I = Eigen::Matrix<double, GlobalDim, GlobalDim>::Identity();
        Eigen::Matrix<double, GlobalDim, GlobalDim> const
            effective_thermal_conductivity =
                (phi_G * gas_thermal_conductivity +
                 phi_L * liquid_thermal_conductivity) *
                    I +
                phi_S * solid_thermal_conductivity;
        return effective_thermal_conductivity;
    }
}

template <int GlobalDim>
PropertyDataType EffectiveThermalConductivityPorosityMixing<GlobalDim>::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    if (variable != Variable::temperature)
    {
        OGS_FATAL(
            "The derivative of the "
            "EffectiveThermalConductivityPorosityMixing is implemented only "
            "w.r.t. temperature.");
    }

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

    auto const porosity = variable_array.porosity;

    auto const S_L = variable_array.liquid_saturation;
    auto const S_G = 1. - S_L;

    auto const phi_G = porosity * S_G;
    auto const phi_L = porosity * S_L;
    auto const phi_S = 1. - porosity;

    // Derivatives of thermal conductivities w.r.t. temperature
    auto const d_gas_thermal_conductivity_dT =
        gas_phase
            ? gas_phase
                  ->property(
                      MaterialPropertyLib::PropertyType::thermal_conductivity)
                  .template dValue<double>(variable_array,
                                           Variable::temperature, pos, t, dt)
            : 0.;

    auto const d_liquid_thermal_conductivity_dT =
        liquid_phase
            ? liquid_phase
                  ->property(
                      MaterialPropertyLib::PropertyType::thermal_conductivity)
                  .template dValue<double>(variable_array,
                                           Variable::temperature, pos, t, dt)
            : 0.;

    // For volume fractions, we need to consider derivatives w.r.t. porosity and
    // saturation. Assuming porosity and saturation may depend on temperature
    // d(phi_G)/dT =
    //  = d(porosity * S_G)/dT
    //  = d(porosity)/dT * S_G + porosity * d(S_G)/dT
    //  = d(porosity)/dT * S_G - porosity * d(S_L)/dT
    // d(phi_L)/dT =
    //  = d(porosity * S_L)/dT
    //  = d(porosity)/dT * S_L + porosity * d(S_L)/dT
    // d(phi_S)/dT =
    //  = d(1 - porosity)/dT
    //  = -d(porosity)/dT.
    double const d_porosity_dT =
        medium->property(MaterialPropertyLib::PropertyType::porosity)
            .template dValue<double>(variable_array, Variable::temperature, pos,
                                     t, dt);
    double const d_S_L_dT =
        // Some processes might not have saturation property (like THM) and
        // saturation passed in the variable_array is always 1.
        medium->hasProperty(MaterialPropertyLib::PropertyType::saturation)
            ? medium->property(MaterialPropertyLib::PropertyType::saturation)
                  .template dValue<double>(variable_array,
                                           Variable::temperature, pos, t, dt)
            : 0.;

    auto const d_phi_G_dT = d_porosity_dT * S_G - porosity * d_S_L_dT;
    auto const d_phi_L_dT = d_porosity_dT * S_L + porosity * d_S_L_dT;
    auto const d_phi_S_dT = -d_porosity_dT;

    auto solid_thermal_conductivity = formEigenTensor<GlobalDim>(
        solid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, t, dt));

    auto d_solid_thermal_conductivity_dT = formEigenTensor<GlobalDim>(
        solid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .dValue(variable_array, Variable::temperature, pos, t, dt));

    if constexpr (GlobalDim == 1)
    {
        // For 1D, return scalar derivative
        // Total derivative of effective thermal conductivity w.r.t. temperature
        double const d_effective_thermal_conductivity_dT =
            d_phi_G_dT * gas_thermal_conductivity +
            phi_G * d_gas_thermal_conductivity_dT +
            d_phi_L_dT * liquid_thermal_conductivity +
            phi_L * d_liquid_thermal_conductivity_dT +
            d_phi_S_dT * solid_thermal_conductivity[0] +
            phi_S * d_solid_thermal_conductivity_dT[0];

        return d_effective_thermal_conductivity_dT;
    }
    else
    {
        // For 2D/3D, use tensor operations

        // Local coordinate transformation is only applied for the case that the
        // initial solid thermal conductivity is given with orthotropic
        // assumption.
        if (local_coordinate_system_)
        {
            solid_thermal_conductivity =
                local_coordinate_system_->rotateTensor<GlobalDim>(
                    solid_thermal_conductivity, pos);
            d_solid_thermal_conductivity_dT =
                local_coordinate_system_->rotateTensor<GlobalDim>(
                    d_solid_thermal_conductivity_dT, pos);
        }

        auto const I = Eigen::Matrix<double, GlobalDim, GlobalDim>::Identity();
        Eigen::Matrix<double, GlobalDim, GlobalDim> const
            d_effective_thermal_conductivity_dT =
                (d_phi_G_dT * gas_thermal_conductivity +
                 phi_G * d_gas_thermal_conductivity_dT +
                 d_phi_L_dT * liquid_thermal_conductivity +
                 phi_L * d_liquid_thermal_conductivity_dT) *
                    I +
                d_phi_S_dT * solid_thermal_conductivity +
                phi_S * d_solid_thermal_conductivity_dT;

        return d_effective_thermal_conductivity_dT;
    }
}
template class EffectiveThermalConductivityPorosityMixing<1>;
template class EffectiveThermalConductivityPorosityMixing<2>;
template class EffectiveThermalConductivityPorosityMixing<3>;
}  // namespace MaterialPropertyLib
