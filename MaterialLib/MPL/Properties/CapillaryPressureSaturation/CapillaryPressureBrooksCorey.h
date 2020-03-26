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

#include "SaturationBrooksCorey.h"

namespace MaterialPropertyLib
{
/// \copydoc MaterialPropertyLib::SaturationBrooksCorey
class CapillaryPressureBrooksCorey final : public SaturationBrooksCorey
{
public:
    CapillaryPressureBrooksCorey(double const residual_liquid_saturation,
                                  double const residual_gas_saturation,
                                  double const exponent,
                                  double const entry_pressure,
                                  double const max_capillary_pressure)
        : SaturationBrooksCorey(residual_liquid_saturation,
                                 residual_gas_saturation, exponent,
                                 entry_pressure),
          _pc_max(max_capillary_pressure)
    {
    }

    void setScale(std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (!std::holds_alternative<Medium*>(scale_pointer))
        {
            OGS_FATAL(
                "The property 'CapillaryPressureBrooksCorey' is implemented on the "
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
    double const _pc_max;  ///< Maximum capillary pressiure
};
}  // namespace MaterialPropertyLib
