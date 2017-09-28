/**
 * \author Norbert Grunwald
 * \date   11.09.2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "../mpProperty.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * This property (usually a component property) computes a linear
 * density function of temperature based on a reference density,
 * a slope, and a reference temperature.
 */
class LinearTemperature final : public Property
{
private:
    /// This property is (currently) implemented to obtain component
    /// properties only.
    Component* _component;

public:
    /// Constructor that passes a pointer to the medium.
    LinearTemperature(Medium*);
    /// Constructor that passes a pointer to the phase.
    LinearTemperature(Phase*);
    /// Constructor that passes a pointer to the component.
    LinearTemperature(Component*);

    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const&) override;
};

}  // MaterialPropertyLib
