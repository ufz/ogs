/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
/// The constant property class. This property simply retrieves the stored
/// constant value. It accepts all datatypes defined in PropertyDataType
/// (currently: double, Vector, Tensor, std::string)
class Constant final : public Property
{
public:
    /// This constructor accepts single values of any data type defined in the
    /// PropertyDataType definition and sets the protected attribute _value of
    /// the base class Property to that value.
    explicit Constant(PropertyDataType const& v);
};
}  // namespace MaterialPropertyLib
