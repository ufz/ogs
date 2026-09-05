// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <vector>

#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"
#include "MeshLib/Utils/nextUnusedMaterialId.h"

class MeshLibNextUnusedMaterialId : public ::testing::Test
{
protected:
    MeshLib::PropertyVector<int>& createMaterialIds(
        std::vector<int> const& values)
    {
        auto* const material_ids = properties.createNewPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell, values.size(), 1);
        MeshLib::PropertyVector<int>& material_ids_ref = *material_ids;
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            material_ids_ref[i] = values[i];
        }
        return material_ids_ref;
    }

    MeshLib::Properties properties;
};

TEST_F(MeshLibNextUnusedMaterialId, EmptyMaterialIdsGiveZero)
{
    EXPECT_EQ(0, MeshLib::nextUnusedMaterialId(createMaterialIds({})));
}

TEST_F(MeshLibNextUnusedMaterialId, GivesOnePastTheLargestId)
{
    // The largest ID is neither the first nor the last entry, so an
    // implementation reading only one end of the vector would fail here.
    EXPECT_EQ(6, MeshLib::nextUnusedMaterialId(createMaterialIds({1, 5, 3})));
}

TEST_F(MeshLibNextUnusedMaterialId, IgnoresGapsAndRepetitions)
{
    EXPECT_EQ(5, MeshLib::nextUnusedMaterialId(createMaterialIds({0, 4, 4})));
}

TEST_F(MeshLibNextUnusedMaterialId, LeavesTheMaterialIdsUnchanged)
{
    auto const& material_ids = createMaterialIds({1, 2});

    MeshLib::nextUnusedMaterialId(material_ids);

    EXPECT_EQ((std::vector<int>{1, 2}),
              std::vector<int>(material_ids.begin(), material_ids.end()));
}
