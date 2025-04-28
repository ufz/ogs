/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#include "MeshLib/Utils/IntegrationPointMetaData.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "BaseLib/Error.h"

using namespace MeshLib;

TEST(MeshLibIntegrationPointMetaDataTest, ThrowsIfNoMatchingName)
{
    std::vector<IntegrationPointMetaData> data = {{"stress", 6, 2}};

    std::string const json_str = IntegrationPointMetaData::toJsonString(data);

    try
    {
        IntegrationPointMetaData::fromJsonString(json_str, "strain");
        FAIL() << "Expected std::runtime_error";
    }
    catch (const std::runtime_error& e)
    {
        EXPECT_THAT(
            e.what(),
            testing::HasSubstr(
                "No integration point meta data with name 'strain' found."));
    }
}

TEST(MeshLibIntegrationPointMetaDataTest, ThrowsIfMultipleMatchingNames)
{
    std::vector<IntegrationPointMetaData> data = {{"stress", 6, 2},
                                                  {"stress", 3, 1}};

    std::string const json_str = IntegrationPointMetaData::toJsonString(data);

    try
    {
        IntegrationPointMetaData::fromJsonString(json_str, "stress");
        FAIL() << "Expected std::runtime_error";
    }
    catch (const std::runtime_error& e)
    {
        EXPECT_THAT(
            e.what(),
            testing::HasSubstr("Expected exactly one integration point meta "
                               "data with name 'stress', found 2."));
    }
}

TEST(MeshLibIntegrationPointMetaDataTest, FindsAllUniqueQuantities)
{
    std::vector<IntegrationPointMetaData> data = {
        {"stress", 6, 2}, {"strain", 6, 2}, {"damage", 1, 1}};

    std::string const json_str = IntegrationPointMetaData::toJsonString(data);

    auto const stress =
        IntegrationPointMetaData::fromJsonString(json_str, "stress");
    EXPECT_EQ(stress.field_name, "stress");
    EXPECT_EQ(stress.n_components, 6);
    EXPECT_EQ(stress.integration_order, 2);

    auto const strain =
        IntegrationPointMetaData::fromJsonString(json_str, "strain");
    EXPECT_EQ(strain.field_name, "strain");
    EXPECT_EQ(strain.n_components, 6);
    EXPECT_EQ(strain.integration_order, 2);

    auto const damage =
        IntegrationPointMetaData::fromJsonString(json_str, "damage");
    EXPECT_EQ(damage.field_name, "damage");
    EXPECT_EQ(damage.n_components, 1);
    EXPECT_EQ(damage.integration_order, 1);
}
