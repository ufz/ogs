/**
 * \file
 * \author Norbert Grunwald
 * \date   July 3rd, 2018
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "MaterialLib/MPL/mpProperty.h"

namespace MaterialPropertyLib
{
/// The constant property class. This property simply retrieves the stored
/// constant value. It accepts all datatypes defined in PropertyDataType.
class Undefined final : public Property
{
private:
    PropertyEnum thisPropertyEnum;
public:
    explicit Undefined(PropertyEnum const&);
    PropertyDataType value() const override;
    PropertyDataType value(VariableArray const&) override;
};
}  // MaterialPropertyLib
