// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

// This test verifies the correct conversion routines from enum to string and
// back. Errors in this conversions can occur easily when the corresponding
// string vectors and the property enum is not lined up correctly.
TEST(Material, convertProperties)
{
    for (int pIndex = 0; pIndex < MPL::PropertyType::number_of_properties;
         pIndex++)
    {
        auto const property_from_enum = static_cast<MPL::PropertyType>(pIndex);
        auto const property_from_string =
            MPL::convertStringToProperty(MPL::property_enum_to_string[pIndex]);

        // If this fails, the enum 'PropertyType' and the string array
        // 'property_enum_to_string' are not in sync!
        ASSERT_EQ(property_from_enum, property_from_string);
    }
}
