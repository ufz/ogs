/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 26, 2020, 4:15 PM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

/**
 * \copydoc MaterialPropertyLib::SaturationLiakopoulos
 *
 *  CapillaryPressureLiakopoulos handles the computations that are associated
 *  with
 *    \f[ p_\mathrm{c}(S)
 *       = \left(1.9722\cdot10^{11} (1-S)\right)^{1/2.4279}
 *    \f]
 */
class CapillaryPressureLiakopoulos final : public Property
{
public:
    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (!std::holds_alternative<Medium*>(scale_pointer))
        {
            OGS_FATAL(
                "The property 'CapillaryPressureLiakopoulos' is implemented on "
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
    const double _residual_liquid_saturation = 0.2;
    const double _parameter_a = 1.9722e-11;
    const double _parameter_b = 2.4279;
    const double _p_cap_max = std::pow(
        (1. - _residual_liquid_saturation) / _parameter_a, (1. / _parameter_b));
};
}  // namespace MaterialPropertyLib
