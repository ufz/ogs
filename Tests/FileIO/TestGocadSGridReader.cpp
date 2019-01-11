/**
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <memory>

#include "gtest/gtest.h"

#include "Applications/FileIO/GocadIO/GenerateFaceSetMeshes.h"
#include "Applications/FileIO/GocadIO/GocadSGridReader.h"
#include "BaseLib/BuildInfo.h"
#include "MeshLib/Mesh.h"

TEST(FileIOGocadSGridReader, flow_simulation_grid_klein)
{
    std::string const file_name(
        BaseLib::BuildInfo::data_path +
        "/FileIO/Gocad/IFG-Uni-Kiel_kleinesGrid/flow_simulation_grid_klein.sg");
    FileIO::Gocad::GocadSGridReader reader(file_name);
    std::unique_ptr<MeshLib::Mesh> mesh(reader.getMesh());

    ASSERT_EQ(220, mesh->getNumberOfNodes());
    ASSERT_EQ(120, mesh->getNumberOfElements());

    auto property_names(mesh->getProperties().getPropertyVectorNames());
    ASSERT_EQ(11, property_names.size());
    ASSERT_EQ(" heat conductivity1", property_names[0]);
    ASSERT_EQ(" specific heat capacity1", property_names[1]);
    ASSERT_EQ("CELL_VOLUME_CORRECTION_FACTOR", property_names[2]);
    ASSERT_EQ("GEOLOGIC_K", property_names[3]);
    ASSERT_EQ("RegionFlags", property_names[4]);
    ASSERT_EQ("UNIT", property_names[5]);
    ASSERT_EQ("VISIBLE_TETRA_RECORD_ID", property_names[6]);
    ASSERT_EQ("density1", property_names[7]);
    ASSERT_EQ("fault_blocks", property_names[8]);
    ASSERT_EQ("kf1", property_names[9]);
    ASSERT_EQ("porosity1", property_names[10]);
}

TEST(FileIOGocadSGridReader, flow_simulation_grid_klein_Rinne)
{
    std::string const file_name(
        BaseLib::BuildInfo::data_path +
        "/FileIO/Gocad/IFG-Uni-Kiel_kleinesGridRinne/flow_simulation_grid_klein_Rinne.sg");
    FileIO::Gocad::GocadSGridReader reader(file_name);
    std::unique_ptr<MeshLib::Mesh> mesh(reader.getMesh());

    ASSERT_EQ(406, mesh->getNumberOfNodes());
    ASSERT_EQ(180, mesh->getNumberOfElements());

    auto property_names(mesh->getProperties().getPropertyVectorNames());
    ASSERT_EQ(11, property_names.size());
    ASSERT_EQ("CELL_VOLUME_CORRECTION_FACTOR", property_names[0]);
    ASSERT_EQ("GEOLOGIC_K", property_names[1]);
    ASSERT_EQ("RegionFlags", property_names[2]);
    ASSERT_EQ("UNIT", property_names[3]);
    ASSERT_EQ("VISIBLE_TETRA_RECORD_ID", property_names[4]);
    ASSERT_EQ("density", property_names[5]);
    ASSERT_EQ("fault_blocks", property_names[6]);
    ASSERT_EQ("hc", property_names[7]);
    ASSERT_EQ("kf", property_names[8]);
    ASSERT_EQ("porosity", property_names[9]);
    ASSERT_EQ("shc", property_names[10]);
}

// this project uses split nodes
TEST(FileIOGocadSGridReader, IFG_Uni_Kiel_SGrid_Test)
{
    std::string const file_name(
        BaseLib::BuildInfo::data_path +
        "/FileIO/Gocad/IFG-Uni-Kiel_SGrid2OGS_Test/SGrid_Test.sg");
    FileIO::Gocad::GocadSGridReader reader(file_name);
    std::unique_ptr<MeshLib::Mesh> mesh(reader.getMesh());

    ASSERT_EQ(432141, mesh->getNumberOfNodes());
    ASSERT_EQ(388608, mesh->getNumberOfElements());

    auto property_names(mesh->getProperties().getPropertyVectorNames());
    ASSERT_EQ(6, property_names.size());
    ASSERT_EQ("CELL_VOLUME_CORRECTION_FACTOR", property_names[0]);
    ASSERT_EQ("GEOLOGIC_K", property_names[1]);
    ASSERT_EQ("RegionFlags", property_names[2]);
    ASSERT_EQ("UNIT", property_names[3]);
    ASSERT_EQ("VISIBLE_TETRA_RECORD_ID", property_names[4]);
    ASSERT_EQ("fault_blocks", property_names[5]);
}
