// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"

class MeshLibPropertyVector : public ::testing::Test
{
protected:
    MeshLib::PropertyVector<double>& createProperty(
        std::size_t const n_property_values)
    {
        auto* const property = properties.createNewPropertyVector<double>(
            "property", MeshLib::MeshItemType::Cell, n_property_values, 1);
        return *property;
    }

    MeshLib::Properties properties;
};

TEST_F(MeshLibPropertyVector, ExtendGrowsByCountAndPointsToTheFirstNewValue)
{
    auto& property = createProperty(2);
    property[0] = 1.0;
    property[1] = 2.0;

    auto const first_new = MeshLib::extend(property, 3);

    ASSERT_EQ(5, property.size());
    EXPECT_EQ(property.begin() + 2, first_new);
    EXPECT_EQ(1.0, property[0]);
    EXPECT_EQ(2.0, property[1]);
}

TEST_F(MeshLibPropertyVector, ExtendValueInitialisesTheNewValues)
{
    auto& property = createProperty(0);

    MeshLib::extend(property, 2);

    ASSERT_EQ(2, property.size());
    EXPECT_EQ(0.0, property[0]);
    EXPECT_EQ(0.0, property[1]);
}

TEST_F(MeshLibPropertyVector, ExtendByZeroLeavesThePropertyUnchanged)
{
    auto& property = createProperty(1);
    property[0] = 1.0;

    auto const first_new = MeshLib::extend(property, 0);

    ASSERT_EQ(1, property.size());
    EXPECT_EQ(property.end(), first_new);
    EXPECT_EQ(1.0, property[0]);
}

TEST_F(MeshLibPropertyVector, AppendValuesKeepsThePresentValues)
{
    auto& property = createProperty(1);
    property[0] = 1.0;

    MeshLib::appendValues(property, std::array{2.0, 3.0});

    ASSERT_EQ(3, property.size());
    EXPECT_EQ((std::vector<double>{1.0, 2.0, 3.0}),
              std::vector<double>(property.begin(), property.end()));
}

TEST_F(MeshLibPropertyVector, AppendValuesConvertsToThePropertyValueType)
{
    auto& property = createProperty(0);

    MeshLib::appendValues(property, std::vector<int>{1, 2});

    ASSERT_EQ(2, property.size());
    EXPECT_EQ((std::vector<double>{1.0, 2.0}),
              std::vector<double>(property.begin(), property.end()));
}

TEST_F(MeshLibPropertyVector, AppendValuesOfAnEmptyRangeIsANoOp)
{
    auto& property = createProperty(1);
    property[0] = 1.0;

    MeshLib::appendValues(property, std::vector<double>{});

    ASSERT_EQ(1, property.size());
    EXPECT_EQ(1.0, property[0]);
}
