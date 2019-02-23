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
 *
 */

#include "Property.h"

#include <string>

namespace MaterialPropertyLib
{

PropertyDataType Property::value() const
{
    return _value;
}
/// The default implementation of this method only returns the property value
/// without altering it.
PropertyDataType Property::value(VariableArray const& /*variable_array*/) const
{
    return _value;
}

/// The default implementation of this method only returns the
/// property value derivative without altering it.
PropertyDataType Property::dValue(VariableArray const& /*variable_array*/,
                                  Variable const /*variable*/) const
{
    return _dvalue;
}

/// Default implementation: 2nd derivative of any constant property is zero.
PropertyDataType Property::d2Value(VariableArray const& /*variable_array*/,
                                   Variable const /*variable*/,
                                   Variable const /*variable*/) const
{
    return 0.0;
}

void Property::notImplemented(const std::string& property,
                              const std::string& material) const
{
    OGS_FATAL("The property '%s' is not available on the '%s' scale",
              property.c_str(), material.c_str());
}
}  // namespace MaterialPropertyLib
