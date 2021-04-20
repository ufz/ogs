/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
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
#pragma once

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \brief Density function for ideal gases.
 * \details This property must be either a phase or a component property, it
 * computes the density of an ideal gas as function of phase pressure and
 * temperature.
 */
class IdealGasLaw final : public Property
{
public:
    explicit IdealGasLaw(std::string name);

    void checkScale() const override
    {
        if (!(std::holds_alternative<Phase*>(scale_) ||
              std::holds_alternative<Component*>(scale_)))
        {
            OGS_FATAL(
                "The property 'IdealGasLaw' is implemented on the "
                "'phase' and 'component' scales only.");
        }
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property values_ and dValues_.
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
};

}  // namespace MaterialPropertyLib
