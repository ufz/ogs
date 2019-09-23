/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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
 * \class RelPermBrooksCorey
 * \brief Relative permeability function proposed by Brooks&Corey
 * \details This property must be a medium property, it
 * computes the permeability reduction due to saturation as function
 * of capillary pressure.
 */
class RelPermBrooksCorey final : public Property
{
private:
    Medium* _medium;
    const double _residual_liquid_saturation;
    const double _residual_gas_saturation;
    const double _min_relative_permeability_liquid;
    const double _min_relative_permeability_gas;
    const double _exponent;

public:
    RelPermBrooksCorey(const double /*residual_liquid_saturation*/,
                       const double /*residual_gas_saturation*/,
                       const double /*_min_relative_permeability_liquid*/,
                       const double /*_min_relative_permeability_gas*/,
                       const double /*exponent*/
    );
    /// This method assigns a pointer to the material object that is the owner
    /// of this property
    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (std::holds_alternative<Medium*>(scale_pointer))
        {
            _medium = std::get<Medium*>(scale_pointer);
        }
        else
        {
            OGS_FATAL(
                "The property 'RelPermBrooksCorey' is implemented on the "
                "'media' scale only.");
        }
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property _values and _dValues.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t) const override;
};

}  // namespace MaterialPropertyLib
