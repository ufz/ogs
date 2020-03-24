/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 20, 2020, 9:59 AM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

/// \copydoc MaterialPropertyLib::SaturationVanGenuchten
class CapillaryPressureVanGenuchten : public Property
{
public:
    CapillaryPressureVanGenuchten(double const residual_liquid_saturation,
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
                "The property 'CapillaryVanGenuchten' is implemented on the "
                "'media' scale only.");
        }
        this->_medium = std::get<Medium*>(scale_pointer);
    }

    /// gets \f$ p_c(S) \f$.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    /// gets \f$ \frac{\partial p_c(S)}{\partial  S} \f$
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

    /// gets \f$ \frac{\partial^2 p_c(S)}{\partial  S^2} \f$
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& pos,
                             double const t, double const dt) const override;

private:
    Medium* _medium = nullptr;
    /// Residual saturation of liquid phase.
    double const _residual_saturation;
    double const _maximuml_saturation;  ///< Maximum saturation of liquid phase.
    double const _m;                    ///< Exponent index.
    double const _p_b;                  ///< Entry pressure.
    double const _pc_max;               ///< Maximum capillary pressure.
};
}  // namespace MaterialPropertyLib
