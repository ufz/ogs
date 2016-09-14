/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "BaseLib/BuildInfo.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Properties.h"

#include "Applications/FileIO/FEFLOW/FEFLOWMeshInterface.h"

TEST(FileIO, TestFEFLOWMeshInterface)
{
    std::string const file_name (BaseLib::BuildInfo::data_path + "/FileIO/small_cube_hex.fem");

    FileIO::FEFLOWMeshInterface feflowIO;
    std::unique_ptr<MeshLib::Mesh const> mesh(feflowIO.readFEFLOWFile(file_name));

    ASSERT_TRUE(static_cast<bool>(mesh));

    ASSERT_EQ(144u, mesh->getNumberOfNodes());
    ASSERT_EQ(75u, mesh->getNumberOfElements());
    auto opt_material_ids(mesh->getProperties().getPropertyVector<int>("MaterialIDs"));
    ASSERT_TRUE(static_cast<bool>(opt_material_ids));
    ASSERT_EQ(75u, opt_material_ids->size());

    MeshLib::Element const* e = mesh->getElement(0);
    ASSERT_EQ(MeshLib::CellType::HEX8, e->getCellType());
    ASSERT_EQ(1, (*opt_material_ids)[0]);
    ASSERT_EQ(0, (*opt_material_ids)[6]);
}
