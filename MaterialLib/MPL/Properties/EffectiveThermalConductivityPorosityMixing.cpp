/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
EffectiveThermalConductivityPorosityMixing<1>::EffectiveThermalConductivityPorosityMixing(
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
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    auto const porosity = medium->property(
                        MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(variable_array, pos, t, dt);
    auto const liquid_thermal_conductivity =
        liquid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template value<double>(variable_array, pos, t, dt);
    auto const solid_thermal_conductivity = solid_phase.property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template value<double>(variable_array, pos, t, dt);
    auto const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    double const effective_thermal_conductivity =
        (1.0 - porosity) * solid_thermal_conductivity +
        porosity * liquid_thermal_conductivity * S_L;
    return effective_thermal_conductivity;
}
template <>
PropertyDataType EffectiveThermalConductivityPorosityMixing<1>::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    auto const& medium = std::get<Medium*>(scale_);
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    auto const porosity = medium->property(
                        MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(variable_array, pos, t, dt);
   auto const liquid_thermal_conductivity =
        liquid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template value<double>(variable_array, pos, t, dt);
    auto const solid_thermal_conductivity = solid_phase.property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template value<double>(variable_array, pos, t, dt);
    auto const dporosity = medium->property(
                        MaterialPropertyLib::PropertyType::porosity)
                    .template dValue<double>(variable_array, variable, pos, t, dt);
    auto const dliquid_thermal_conductivity =
        liquid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template dValue<double>(variable_array, variable, pos, t, dt);
    auto const dsolid_thermal_conductivity = solid_phase.property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template dValue<double>(variable_array, variable, pos, t, dt);

    double const deffective_thermal_conductivity =
        -1.0 * dporosity * solid_thermal_conductivity +
        (1.0 - porosity) * dsolid_thermal_conductivity +
        dporosity * liquid_thermal_conductivity +
        porosity * dliquid_thermal_conductivity;
    return deffective_thermal_conductivity;

}
//
// For 2D and 3D problems
//
template <int GlobalDim>
EffectiveThermalConductivityPorosityMixing<GlobalDim>::EffectiveThermalConductivityPorosityMixing(
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
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    auto const porosity = medium->property(
                        MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(variable_array, pos, t, dt);
    auto const liquid_thermal_conductivity =
        liquid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .template value<double>(variable_array, pos, t, dt);
    auto solid_thermal_conductivity = formEigenTensor<GlobalDim>(
        solid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_conductivity)
            .value(variable_array, pos, t, dt));
    auto const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    // Local coordinate transformation is only applied for the case that the
    // initial solid thermal conductivity is given with orthotropic assumption.
    if (local_coordinate_system_ && (solid_thermal_conductivity.cols() == GlobalDim))
    {
        Eigen::Matrix<double, GlobalDim, GlobalDim> const e =
            local_coordinate_system_->transformation<GlobalDim>(pos);

        solid_thermal_conductivity = e.transpose() * solid_thermal_conductivity * e;
    }
    Eigen::Matrix<double, GlobalDim, GlobalDim> const
        effective_thermal_conductivity =
            (1.0 - porosity) * solid_thermal_conductivity +
            porosity * liquid_thermal_conductivity *
                Eigen::Matrix<double, GlobalDim, GlobalDim>::Identity() * S_L;
    return effective_thermal_conductivity;
}

template <int GlobalDim>
PropertyDataType EffectiveThermalConductivityPorosityMixing<GlobalDim>::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    auto const& medium = std::get<Medium*>(scale_);
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    auto const porosity = medium->property(
                        MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(variable_array, pos, t, dt);
    auto const liquid_thermal_conductivity =
            liquid_phase
                .property(
                    MaterialPropertyLib::PropertyType::thermal_conductivity)
                .template value<double>(variable_array, pos, t, dt);
    auto solid_thermal_conductivity =
        formEigenTensor<GlobalDim>(solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .value(variable_array, pos, t, dt));
    auto const dporosity = medium->property(
                        MaterialPropertyLib::PropertyType::porosity)
                    .template dValue<double>(variable_array, variable, pos, t, dt);
    auto const dliquid_thermal_conductivity =
            liquid_phase
                .property(
                    MaterialPropertyLib::PropertyType::thermal_conductivity)
                .template dValue<double>(variable_array, variable, pos, t, dt);
    auto dsolid_thermal_conductivity =
        formEigenTensor<GlobalDim>(solid_phase
                   .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .dValue(variable_array, variable, pos, t, dt));

    // Local coordinate transformation is only applied for the case that the
    // initial solid thermal conductivity is given with orthotropic assumption.
    if (local_coordinate_system_ && (solid_thermal_conductivity.cols() == GlobalDim))
    {
        Eigen::Matrix<double, GlobalDim, GlobalDim> const e =
            local_coordinate_system_->transformation<GlobalDim>(pos);

        solid_thermal_conductivity = e.transpose() * solid_thermal_conductivity * e;
        dsolid_thermal_conductivity =
            e.transpose() * dsolid_thermal_conductivity * e;
    }
    Eigen::Matrix<double, GlobalDim, GlobalDim> const
        deffective_thermal_conductivity =
            -1.0 * dporosity * solid_thermal_conductivity +
            (1.0 - porosity) * dsolid_thermal_conductivity +
            dporosity * liquid_thermal_conductivity *
                Eigen::Matrix<double, GlobalDim, GlobalDim>::Identity() +
            porosity * dliquid_thermal_conductivity *
                Eigen::Matrix<double, GlobalDim, GlobalDim>::Identity();
    return deffective_thermal_conductivity;
}
template class EffectiveThermalConductivityPorosityMixing<2>;
template class EffectiveThermalConductivityPorosityMixing<3>;
}  // namespace MaterialPropertyLib
