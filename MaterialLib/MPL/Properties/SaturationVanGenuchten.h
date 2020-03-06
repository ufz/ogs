/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <cmath>
#include <limits>

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

/// The van Genuchten soil characteristics function.
///
/// This property must be a medium property, it computes the saturation of the
/// wetting phase as function of capillary pressure.
class SaturationVanGenuchten final : public Property
{
private:
    Medium* _medium = nullptr;
    double const _S_L_res;
    double const _S_L_max;
    double const _m;
    double const _p_b;
    double const _pc_max; ///< Maximum capillary pressiure

public:
    SaturationVanGenuchten(double const residual_liquid_saturation,
                           double const residual_gas_saturation,
                           double const exponent,
                           double const entry_pressure,
                           double const max_capillary_pressure);

    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (!std::holds_alternative<Medium*>(scale_pointer))
        {
            OGS_FATAL(
                "The property 'SaturationVanGenuchten' is implemented on the "
                "'media' scale only.");
        }
        _medium = std::get<Medium*>(scale_pointer);
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property _values and _dValues.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/,
                           double const /*dt*/) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& /*pos*/,
                             double const /*t*/,
                             double const /*dt*/) const override;

    /// This method gets capillary pressure via saturation.
    PropertyDataType inverse_value(VariableArray const& variable_array,
                                   ParameterLib::SpatialPosition const& pos,
                                   double const t,
                                   double const dt) const override;
};
}  // namespace MaterialPropertyLib
