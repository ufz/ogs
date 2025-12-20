// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "PropertyType.h"

#include <boost/algorithm/string/predicate.hpp>

#include "BaseLib/Error.h"

namespace MaterialPropertyLib
{
PropertyType convertStringToProperty(std::string const& string)
{
    for (int i = 0; i < static_cast<int>(PropertyType::number_of_properties);
         ++i)
    {
        if (boost::iequals(string, property_enum_to_string[i]))
        {
            return static_cast<PropertyType>(i);
        }
    }

    OGS_FATAL(
        "The property name '{:s}' does not correspond to any known property",
        string);
}
}  // namespace MaterialPropertyLib
