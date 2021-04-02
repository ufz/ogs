/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <memory>

#include "Applications/FileIO/GMSInterface.h"
#include "InfoLib/TestInfo.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshInformation.h"
#include "gtest/gtest.h"

TEST(FileIO, TestGmsInterface)
{
    std::string const file_name(TestInfoLib::TestInfo::data_path +
                                "/FileIO/3DMeshData.3dm");
    std::unique_ptr<MeshLib::Mesh> mesh(
        FileIO::GMSInterface::readGMS3DMMesh(file_name));
    ASSERT_TRUE(mesh != nullptr);
    ASSERT_EQ(11795, mesh->getNumberOfNodes());
    ASSERT_EQ(19885, mesh->getNumberOfElements());

    auto const& types =
        MeshLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(1456, types.at(MeshLib::MeshElemType::TETRAHEDRON));  // tets
    ASSERT_EQ(1355, types.at(MeshLib::MeshElemType::PYRAMID));      // pyramids
    ASSERT_EQ(17074, types.at(MeshLib::MeshElemType::PRISM));       // prism
    auto const& bounds =
        MeshLib::MeshInformation::getValueBounds(*materialIDs(*mesh));
    ASSERT_TRUE(bounds.has_value());
    ASSERT_EQ(1, bounds->first);
    ASSERT_EQ(63, bounds->second);
}
