/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <memory>
#include <numeric>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/getMeshElementsForMaterialIDs.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"

using namespace MeshLib;

TEST(MeshLib, getMeshElementsForMaterialIDs)
{
    auto input_mesh = std::unique_ptr<Mesh>(
        MeshToolsLib::MeshGenerator::generateLineMesh(1.0, 11));
    auto* material_ids = getOrCreateMeshProperty<int>(
        *input_mesh.get(), "MaterialIDs", MeshItemType::Cell, 1);
    std::iota(material_ids->begin(), material_ids->end(), 0);

    for (std::size_t i = 0; i < material_ids->size(); ++i)
    {
        auto elements = getMeshElementsForMaterialIDs(*input_mesh.get(),
                                                      {static_cast<int>(i)});
        ASSERT_EQ(1, elements.size());
    }

    std::fill(material_ids->begin(), material_ids->end(), 0);
    {
        auto elements = getMeshElementsForMaterialIDs(*input_mesh.get(), {0});
        ASSERT_EQ(material_ids->size(), elements.size());
    }
    {
        auto elements = getMeshElementsForMaterialIDs(*input_mesh.get(), {1});
        ASSERT_EQ(0, elements.size());
    }

    std::fill(material_ids->begin() + material_ids->size() / 2,
              material_ids->end(),
              1);
    {
        auto elements = getMeshElementsForMaterialIDs(*input_mesh.get(), {0});
        ASSERT_EQ(material_ids->size() / 2, elements.size());
    }
    {
        auto elements = getMeshElementsForMaterialIDs(*input_mesh.get(), {1});
        // if size of material ids is odd material_ids->size()/2 is equal to
        // (material_ids->size() - 1)/2
        ASSERT_EQ(material_ids->size() - (material_ids->size() / 2),
                  elements.size());
    }
    {
        auto elements =
            getMeshElementsForMaterialIDs(*input_mesh.get(), {0, 1});
        ASSERT_EQ(material_ids->size(), elements.size());
    }
    {
        auto elements =
            getMeshElementsForMaterialIDs(*input_mesh.get(), {2, 3});
        ASSERT_EQ(0, elements.size());
    }
}
