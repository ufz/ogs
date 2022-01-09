/**
 * \file
 * \author Norbert Grunwald
 * \date   Dec 07, 2020
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * Vapour pressure as function of temperature based on Clausius-Clapeyron
 * equation.
 * This property must be either a phase or a component
 * property, it computes the saturation vapour pressure of a substance
 */
class ClausiusClapeyron final : public Property
{
public:
    explicit ClausiusClapeyron(std::string name,
                               const double triple_temperature,
                               const double triple_pressure,
                               const double critical_temperature,
                               const double critical_pressure,
                               const double ref_temperature,
                               const double ref_pressure);

    void checkScale() const override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const primary_variable1,
                             Variable const primary_variable2,
                             ParameterLib::SpatialPosition const& pos,
                             double const t,
                             double const dt) const override;

private:
    double T_triple_;
    double p_triple_;
    double T_critical_;
    double p_critical_;
    double T_ref_;
    double p_ref_;

    double molarMass(std::variant<Medium*, Phase*, Component*> const scale,
                     VariableArray const& variable_array,
                     ParameterLib::SpatialPosition const& pos, double const t,
                     double const dt) const;

    double dMolarMass(std::variant<Medium*, Phase*, Component*> const scale,
                      VariableArray const& variable_array,
                      Variable const primary_variable,
                      ParameterLib::SpatialPosition const& pos, double const t,
                      double const dt) const;
};

}  // namespace MaterialPropertyLib
