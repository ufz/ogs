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

TEST(MeshLibIntegrationPointMetaDataSingleFieldTest, ThrowsIfNoMatchingName)
{
    std::vector<IntegrationPointMetaDataSingleField> data = {{"stress", 6, 2}};
    IntegrationPointMetaData const ip_meta_data{data};

    try
    {
        ip_meta_data["strain"];
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

TEST(MeshLibIntegrationPointMetaDataSingleFieldTest,
     ThrowsIfMultipleMatchingNames)
{
    std::vector<IntegrationPointMetaDataSingleField> data = {{"stress", 6, 2},
                                                             {"stress", 3, 1}};

    try
    {
        IntegrationPointMetaData const ip_meta_data{data};
        FAIL() << "Expected std::runtime_error";
    }
    catch (const std::runtime_error& e)
    {
        EXPECT_THAT(
            e.what(),
            testing::HasSubstr(
                "Duplicate integration point meta data names found: stress."));
    }
}

TEST(MeshLibIntegrationPointMetaDataSingleFieldTest, FindsAllUniqueQuantities)
{
    std::vector<IntegrationPointMetaDataSingleField> data = {
        {"stress", 6, 2}, {"strain", 6, 2}, {"damage", 1, 1}};

    IntegrationPointMetaData const ip_meta_data{data};

    auto const stress = ip_meta_data["stress"];
    EXPECT_EQ(stress.field_name, "stress");
    EXPECT_EQ(stress.n_components, 6);
    EXPECT_EQ(stress.integration_order, 2);

    auto const strain = ip_meta_data["strain"];
    EXPECT_EQ(strain.field_name, "strain");
    EXPECT_EQ(strain.n_components, 6);
    EXPECT_EQ(strain.integration_order, 2);

    auto const damage = ip_meta_data["damage"];
    EXPECT_EQ(damage.field_name, "damage");
    EXPECT_EQ(damage.n_components, 1);
    EXPECT_EQ(damage.integration_order, 1);
}

// TODO (naumov) Roundtrip tests for json conversion
