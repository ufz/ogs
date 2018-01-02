/**
 * \file   TestTetGenInterface.cpp
 * \author Karsten Rink
 * \date   2016-02-19
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <limits>
#include <vector>
#include <memory>

#include "gtest/gtest.h"

#include "Applications/FileIO/TetGenInterface.h"
#include "BaseLib/BuildInfo.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshInformation.h"

// read TetGen geometry
TEST(FileIO, TetGenSmeshReader)
{
    std::string const file_name (BaseLib::BuildInfo::data_path + "/FileIO/twolayermdl.smesh");
    GeoLib::GEOObjects geo_objects;
    FileIO::TetGenInterface tgi;
    bool const result (tgi.readTetGenGeometry(file_name, geo_objects));
    ASSERT_TRUE(result);

    std::vector<GeoLib::Point*> const& pnts (*geo_objects.getPointVec("twolayermdl"));
    ASSERT_EQ(744, pnts.size());
    std::vector<GeoLib::Surface*> const& sfcs (*geo_objects.getSurfaceVec("twolayermdl"));
    ASSERT_EQ(5, sfcs.size());
    ASSERT_EQ(468, sfcs[2]->getNumberOfTriangles());
    ASSERT_EQ(191, sfcs[3]->getNumberOfTriangles());
}

// existing mesh to TetGen geometry
#ifndef USE_PETSC
TEST(FileIO, TetGenSmeshInterface)
#else
TEST(FileIO, DISABLED_TetGenSmeshInterface)
#endif
{
    std::string const file_name (BaseLib::BuildInfo::data_path + "/FileIO/AmmerSubsurfaceCoarse.vtu");
    std::unique_ptr<MeshLib::Mesh const> const mesh (MeshLib::IO::readMeshFromFile(file_name));

    std::string const tg_new_name ("TestSmeshWriter");
    std::string const output_name(BaseLib::BuildInfo::tests_tmp_path + tg_new_name + ".smesh");
    std::vector<MeshLib::Node> attr_pnts;
    FileIO::TetGenInterface tgi;
    bool result (tgi.writeTetGenSmesh(output_name, *mesh, attr_pnts));
    ASSERT_TRUE(result);

    GeoLib::GEOObjects geo_objects;
    result = tgi.readTetGenGeometry(output_name, geo_objects);
    ASSERT_TRUE(result);
    std::string const ref_name(BaseLib::BuildInfo::data_path + "/FileIO/AmmerSubsurfaceCoarse.smesh");
    result = tgi.readTetGenGeometry(ref_name, geo_objects);
    ASSERT_TRUE(result);

    std::vector<GeoLib::Point*> const& ref_pnts (*geo_objects.getPointVec("AmmerSubsurfaceCoarse"));
    std::vector<GeoLib::Point*> const& new_pnts (*geo_objects.getPointVec(tg_new_name));
    ASSERT_EQ(ref_pnts.size(), new_pnts.size());

    std::vector<GeoLib::Surface*> const& ref_sfc (*geo_objects.getSurfaceVec("AmmerSubsurfaceCoarse"));
    std::vector<GeoLib::Surface*> const& new_sfc (*geo_objects.getSurfaceVec(tg_new_name));
    ASSERT_EQ(ref_sfc.size(), new_sfc.size());

    for (std::size_t i=0; i<ref_sfc.size(); ++i)
        ASSERT_EQ(ref_sfc[i]->getNumberOfTriangles(), new_sfc[i]->getNumberOfTriangles());

    std::remove(output_name.c_str());
}

// TetGen mesh with material array
TEST(FileIO, TetGenMeshReaderWithMaterials)
{
    std::string const node_name (BaseLib::BuildInfo::data_path + "/FileIO/twolayermdl.node");
    std::string const ele_name  (BaseLib::BuildInfo::data_path + "/FileIO/twolayermdl.ele");
    FileIO::TetGenInterface tgi;
    std::unique_ptr<MeshLib::Mesh> mesh (tgi.readTetGenMesh(node_name, ele_name));
    ASSERT_TRUE(mesh != nullptr);
    ASSERT_EQ(1378, mesh->getNumberOfNodes());
    ASSERT_EQ(5114, mesh->getNumberOfElements());

    std::pair<int, int> bounds (MeshLib::MeshInformation::getValueBounds<int>(*mesh, "MaterialIDs"));
    ASSERT_EQ(-20, bounds.first);
    ASSERT_EQ(-10, bounds.second);
}

// TetGen mesh without additional information
TEST(FileIO, TetGenMeshReaderWithoutMaterials)
{
    std::string const node_name (BaseLib::BuildInfo::data_path + "/FileIO/tetgen_example.node");
    std::string const ele_name  (BaseLib::BuildInfo::data_path + "/FileIO/tetgen_example.ele");
    FileIO::TetGenInterface tgi;
    std::unique_ptr<MeshLib::Mesh> mesh (tgi.readTetGenMesh(node_name, ele_name));
    ASSERT_TRUE(mesh != nullptr);
    ASSERT_EQ(202, mesh->getNumberOfNodes());
    ASSERT_EQ(650, mesh->getNumberOfElements());

    std::pair<int, int> bounds (MeshLib::MeshInformation::getValueBounds<int>(*mesh, "MaterialIDs"));
    ASSERT_EQ(std::numeric_limits<int>::max(), bounds.first);
    ASSERT_EQ(std::numeric_limits<int>::max(), bounds.second);
}
