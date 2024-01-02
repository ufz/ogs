/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdio>
#include <memory>

#include "GeoLib/IO/AsciiRasterInterface.h"
#include "GeoLib/Raster.h"
#include "InfoLib/TestInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshEditing/RasterDataToMesh.h"

class RasterDataToMeshTest : public ::testing::Test
{
public:
    RasterDataToMeshTest()
    {
        std::string const mesh_path =
            TestInfoLib::TestInfo::data_path + "/MeshLib/RasterDataToMesh.vtu";
        std::string const raster_path =
            TestInfoLib::TestInfo::data_path + "/MeshLib/RasterDataToMesh.asc";
        _mesh.reset(MeshLib::IO::VtuInterface::readVTUFile(mesh_path));
        _raster.reset(FileIO::AsciiRasterInterface::readRaster(raster_path));
    }

protected:
    std::unique_ptr<MeshLib::Mesh> _mesh;
    std::unique_ptr<GeoLib::Raster> _raster;
};

TEST_F(RasterDataToMeshTest, readRasterValuesToNodes)
{
    std::string const vec_name = "imgvalues";
    bool ret = MeshToolsLib::RasterDataToMesh::projectToNodes(*_mesh, *_raster,
                                                              0, vec_name);
    ASSERT_TRUE(ret);

    std::string const compare_path = TestInfoLib::TestInfo::data_path +
                                     "/MeshLib/RasterDataToMesh-Nodes.vtu";
    std::unique_ptr<MeshLib::Mesh> compare{
        MeshLib::IO::VtuInterface::readVTUFile(compare_path)};

    auto new_vec = _mesh->getProperties().getPropertyVector<double>(vec_name);
    auto org_vec =
        compare->getProperties().getPropertyVector<double>("RasterValues");

    bool is_equal =
        std::equal(org_vec->cbegin(), org_vec->cend(), new_vec->cbegin());

    ASSERT_TRUE(is_equal);
}

TEST_F(RasterDataToMeshTest, readRasterValuesToElements)
{
    std::string const vec_name = "imgvalues";
    bool ret = MeshToolsLib::RasterDataToMesh::projectToElements(
        *_mesh, *_raster, 0, vec_name);
    ASSERT_TRUE(ret);

    std::string const compare_path = TestInfoLib::TestInfo::data_path +
                                     "/MeshLib/RasterDataToMesh-Elements.vtu";
    std::unique_ptr<MeshLib::Mesh> compare{
        MeshLib::IO::VtuInterface::readVTUFile(compare_path)};

    auto new_vec = _mesh->getProperties().getPropertyVector<double>(vec_name);
    auto org_vec =
        compare->getProperties().getPropertyVector<double>("RasterValues");

    bool is_equal =
        std::equal(org_vec->cbegin(), org_vec->cend(), new_vec->cbegin());

    ASSERT_TRUE(is_equal);
}
