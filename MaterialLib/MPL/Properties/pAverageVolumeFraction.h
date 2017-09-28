/**
 * \author Norbert Grunwald
 * \date   13.09.2017
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
 * \class AverageVolumeFraction
 * \brief A function averaging a property by volume fraction
 * \details This property is usually a medium property, it
 * computes the average of individual phase properties
 * weighted by volume fraction.
 */
class AverageVolumeFraction final : public Property
{
private:
    /// A pointer to the medium object.
    Medium* _medium;

public:
    /// Constructor passing a pointer to the medium.
    AverageVolumeFraction(Medium*);
    /// Constructor passing a pointer to a phase.
    AverageVolumeFraction(Phase*);
    /// Constructor passing a pointer to a component.
    AverageVolumeFraction(Component*);

    /// This method overrides the base class implementation and
    /// actually computes and sets the property _value.
    PropertyDataType value(VariableArray const&) override;
};

}  // MaterialPropertyLib
