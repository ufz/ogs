/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <limits>
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \brief A well known soil characteristics function
 * \details This property must be a medium property, it
 * computes the saturation of the wetting phase as function
 * of capillary pressure.
 */
class SaturationBrooksCorey final : public Property
{
private:
    Medium* _medium = nullptr;
    const double _residual_liquid_saturation;
    const double _residual_gas_saturation;
    const double _exponent;
    const double _entry_pressure;

public:
    SaturationBrooksCorey(const double residual_liquid_saturation,
                          const double residual_gas_saturation,
                          const double exponent,
                          const double entry_pressure);

    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (!std::holds_alternative<Medium*>(scale_pointer))
        {
            OGS_FATAL(
                "The property 'SaturationBrooksCorey' is implemented on the "
                "'media' scale only.");
        }
        _medium = std::get<Medium*>(scale_pointer);
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property _values and _dValues.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/) const override;
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1,
                             Variable const variable2,
                             ParameterLib::SpatialPosition const& /*pos*/,
                             double const /*t*/) const override;
};
}  // namespace MaterialPropertyLib
