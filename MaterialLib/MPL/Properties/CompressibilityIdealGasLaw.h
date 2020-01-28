/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on January 28, 2020, 16:05 PM
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
 * \class CompressibilityIdealGasLaw
 * \brief Compressibility function for ideal gases
 * \details This property must be either a phase or a component property, it
 * computes the compressibility of an ideal gas as function of phase pressure
 */
class CompressibilityIdealGasLaw final : public Property
{
public:
    /// This method assigns a pointer to the material object that is the owner
    /// of this property
    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (std::holds_alternative<Phase*>(scale_pointer))
        {
            _phase = std::get<Phase*>(scale_pointer);
        }
        else if (std::holds_alternative<Component*>(scale_pointer))
        {
            _component = std::get<Component*>(scale_pointer);
        }
        else
        {
            OGS_FATAL(
                "The property 'CompressibilityIdealGasLaw' is implemented on the "
                "'phase' and 'component' scales only.");
        }
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
                             ParameterLib::SpatialPosition const& pos,
                             double const t) const override;

private:
    Phase* _phase = nullptr;
    Component* _component = nullptr;
};

}  // namespace MaterialPropertyLib
