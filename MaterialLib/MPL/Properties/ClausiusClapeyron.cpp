/**
 * \file
 * \author Norbert Grunwald
 * \date   Dec 07, 2020
 * \brief
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/ClausiusClapeyron.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
ClausiusClapeyron::ClausiusClapeyron(std::string name,
                                     const double triple_temperature,
                                     const double triple_pressure,
                                     const double critical_temperature,
                                     const double critical_pressure,
                                     const double ref_temperature,
                                     const double ref_pressure)
    : T_triple_(triple_temperature),
      p_triple_(triple_pressure),
      T_critical_(critical_temperature),
      p_critical_(critical_pressure),
      T_ref_(ref_temperature),
      p_ref_(ref_pressure)
{
    name_ = std::move(name);
}

void ClausiusClapeyron::checkScale() const
{
    if (!(std::holds_alternative<Phase*>(scale_) ||
          std::holds_alternative<Component*>(scale_)))
    {
        OGS_FATAL(
            "The property 'ClausiusClapeyron' is implemented on 'phase' and "
            "'component' scales only.");
    }
}

double ClausiusClapeyron::molarMass(
    std::variant<Medium*, Phase*, Component*> const scale,
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    return std::visit(
        [&variable_array, &pos, t, dt](auto&& s) -> double {
            return s->property(PropertyType::molar_mass)
                .template value<double>(variable_array, pos, t, dt);
        },
        scale);
}

double ClausiusClapeyron::dMolarMass(
    std::variant<Medium*, Phase*, Component*> const scale,
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    return std::visit(
        [&variable_array, &primary_variable, &pos, t, dt](auto&& s) -> double {
            return s->property(PropertyType::molar_mass)
                .template dValue<double>(variable_array, primary_variable, pos,
                                         t, dt);
        },
        scale);
}

PropertyDataType ClausiusClapeyron::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);

    const double M = molarMass(scale_, variable_array, pos, t, dt);

    if (T >= T_critical_)
    {
        return p_critical_;
    }
    if (T <= T_triple_)
    {
        return p_triple_;
    }

    const double dh = std::get<double>(
        variable_array[static_cast<int>(Variable::enthalpy_of_evaporation)]);
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;

    return p_ref_ * std::exp((1. / T_ref_ - 1. / T) * M * dh / R);
}

PropertyDataType ClausiusClapeyron::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);

    const double M = molarMass(scale_, variable_array, pos, t, dt);
    const double dM =
        dMolarMass(scale_, variable_array, primary_variable, pos, t, dt);

    if (T > T_critical_)
    {
        return 0.;
    }
    if (T < T_triple_)
    {
        return 0.;
    }

    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double dh = std::get<double>(
        variable_array[static_cast<int>(Variable::enthalpy_of_evaporation)]);
    const double p_vap = std::get<double>(value(variable_array, pos, t, dt));

    if (primary_variable == Variable::temperature)
    {
        return p_vap * dh / R * ((1. / T_ref_ - 1. / T) * dM + M / (T * T));
    }
    if (primary_variable == Variable::phase_pressure)
    {
        return p_vap * dh / R * (1. / T_ref_ - 1. / T) * dM;
    }
    OGS_FATAL(
        "ClausiusClapeyron::dValue is implemented for derivatives with respect "
        "to phase pressure or temperature only.");
}

PropertyDataType ClausiusClapeyron::d2Value(
    VariableArray const& /*variable_array*/,
    Variable const /*primary_variable1*/, Variable const /*primary_variable2*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL("ClausiusClapeyron::d2Value is not implemented.");
}

}  // namespace MaterialPropertyLib
