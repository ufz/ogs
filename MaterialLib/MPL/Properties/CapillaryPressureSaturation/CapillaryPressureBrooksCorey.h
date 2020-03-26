/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 25, 2020, 1:51 PM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

/// \copydoc MaterialPropertyLib::SaturationBrooksCorey
class CapillaryPressureBrooksCorey final : public Property
{
public:
    CapillaryPressureBrooksCorey(double const residual_liquid_saturation,
                                 double const maximum_liquid_saturation,
                                 double const exponent,
                                 double const entry_pressure,
                                 double const max_capillary_pressure)
        : _residual_liquid_saturation(residual_liquid_saturation),
          _max_liquid_saturation(maximum_liquid_saturation),
          _exponent(exponent),
          _entry_pressure(entry_pressure),
          _pc_max(max_capillary_pressure)
    {
    }
    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (!std::holds_alternative<Medium*>(scale_pointer))
        {
            OGS_FATAL(
                "The property 'CapillaryPressureBrooksCorey' is implemented on "
                "the "
                "'media' scale only.");
        }
        _medium = std::get<Medium*>(scale_pointer);
    }

    /// It returns \f$ p_c(S) \f$.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    /// It returns \f$ \frac{\partial p_c(S)}{\partial  S} \f$
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

    /// It returns \f$ \frac{\partial^2 p_c(S)}{\partial  S^2} \f$
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& pos,
                             double const t, double const dt) const override;

private:
    Medium* _medium = nullptr;

    /// Residual saturation of liquid phase.
    const double _residual_liquid_saturation;
    const double
        _max_liquid_saturation;    ///< Maximum saturation of liquid phase.
    const double _exponent;        ///< Exponent.
    const double _entry_pressure;  ///< Entry pressure.

    double const _pc_max;  ///< Maximum capillary pressiure
};
}  // namespace MaterialPropertyLib
