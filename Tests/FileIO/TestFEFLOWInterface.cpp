// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <memory>
#include <string>

#include "Applications/FileIO/FEFLOW/FEFLOWMeshInterface.h"
#include "InfoLib/TestInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Properties.h"

TEST(FileIO, TestFEFLOWMeshInterface)
{
    std::string const file_name(TestInfoLib::TestInfo::data_path +
                                "/FileIO/small_cube_hex.fem");

    FileIO::FEFLOWMeshInterface feflowIO;
    std::unique_ptr<MeshLib::Mesh const> mesh(
        feflowIO.readFEFLOWFile(file_name));

    ASSERT_TRUE(static_cast<bool>(mesh));

    ASSERT_EQ(144u, mesh->getNumberOfNodes());
    ASSERT_EQ(75u, mesh->getNumberOfElements());
    auto opt_material_ids(
        mesh->getProperties().getPropertyVector<int>("MaterialIDs"));
    ASSERT_TRUE(static_cast<bool>(opt_material_ids));
    ASSERT_EQ(75u, opt_material_ids->size());

    MeshLib::Element const* e = mesh->getElement(0);
    ASSERT_EQ(MeshLib::CellType::HEX8, e->getCellType());
    ASSERT_EQ(1, (*opt_material_ids)[0]);
    ASSERT_EQ(0, (*opt_material_ids)[6]);
}

TEST(FileIO, TestFEFLOWReadHexMesh)
{
    std::string const fname(TestInfoLib::TestInfo::data_path +
                            "/FileIO/FEFLOW/hex.fem");
    FileIO::FEFLOWMeshInterface feflowIO;
    std::unique_ptr<MeshLib::Mesh const> mesh(feflowIO.readFEFLOWFile(fname));

    EXPECT_EQ(12, mesh->getNumberOfNodes());
    EXPECT_EQ(2, mesh->getNumberOfElements());
    for (std::size_t k(0); k < mesh->getNumberOfElements(); ++k)
    {
        EXPECT_EQ(MeshLib::MeshElemType::HEXAHEDRON,
                  mesh->getElement(k)->getGeomType());
    }
}

TEST(FileIO, TestFEFLOWReadPrismMesh)
{
    std::string const fname(TestInfoLib::TestInfo::data_path +
                            "/FileIO/FEFLOW/prism.fem");
    FileIO::FEFLOWMeshInterface feflowIO;
    std::unique_ptr<MeshLib::Mesh const> mesh(feflowIO.readFEFLOWFile(fname));

    EXPECT_EQ(12, mesh->getNumberOfNodes());
    EXPECT_EQ(4, mesh->getNumberOfElements());
    for (std::size_t k(0); k < mesh->getNumberOfElements(); ++k)
    {
        EXPECT_EQ(MeshLib::MeshElemType::PRISM,
                  mesh->getElement(k)->getGeomType());
    }
}

TEST(FileIO, TestFEFLOWReadHexPrismMesh)
{
    std::string const fname(TestInfoLib::TestInfo::data_path +
                            "/FileIO/FEFLOW/hex_prism.fem");
    FileIO::FEFLOWMeshInterface feflowIO;
    std::unique_ptr<MeshLib::Mesh const> mesh(feflowIO.readFEFLOWFile(fname));

    EXPECT_EQ(12, mesh->getNumberOfNodes());
    EXPECT_EQ(3, mesh->getNumberOfElements());
    for (std::size_t k(0); k < 2; ++k)
    {
        EXPECT_EQ(MeshLib::MeshElemType::PRISM,
                  mesh->getElement(k)->getGeomType());
    }
    EXPECT_EQ(MeshLib::MeshElemType::HEXAHEDRON,
              mesh->getElement(2)->getGeomType());
}

TEST(FileIO, TestFEFLOWReadTetMesh)
{
    std::string const fname(TestInfoLib::TestInfo::data_path +
                            "/FileIO/FEFLOW/tet.fem");
    FileIO::FEFLOWMeshInterface feflowIO;
    std::unique_ptr<MeshLib::Mesh const> mesh(feflowIO.readFEFLOWFile(fname));

    EXPECT_EQ(13, mesh->getNumberOfNodes());
    EXPECT_EQ(22, mesh->getNumberOfElements());
    for (std::size_t k(0); k < mesh->getNumberOfElements(); ++k)
    {
        EXPECT_EQ(MeshLib::MeshElemType::TETRAHEDRON,
                  mesh->getElement(k)->getGeomType());
    }
}

TEST(FileIO, TestFEFLOWReadPrismTetMesh)
{
    std::string const fname(TestInfoLib::TestInfo::data_path +
                            "/FileIO/FEFLOW/prism_tet.fem");
    FileIO::FEFLOWMeshInterface feflowIO;
    std::unique_ptr<MeshLib::Mesh const> mesh(feflowIO.readFEFLOWFile(fname));

    EXPECT_EQ(16, mesh->getNumberOfNodes());
    EXPECT_EQ(20, mesh->getNumberOfElements());
    for (std::size_t k(0); k < 4; ++k)
    {
        EXPECT_EQ(MeshLib::MeshElemType::PRISM,
                  mesh->getElement(k)->getGeomType());
    }
    for (std::size_t k(4); k < mesh->getNumberOfElements(); ++k)
    {
        EXPECT_EQ(MeshLib::MeshElemType::TETRAHEDRON,
                  mesh->getElement(k)->getGeomType());
    }
}

TEST(FileIO, TestFEFLOWReadMeshWithMaterialProperties)
{
    std::string const fname(
        TestInfoLib::TestInfo::data_path +
        "/FileIO/FEFLOW/prism_with_material_properties.fem");
    FileIO::FEFLOWMeshInterface feflowIO;
    std::unique_ptr<MeshLib::Mesh const> mesh(feflowIO.readFEFLOWFile(fname));

    EXPECT_EQ(294u, mesh->getNumberOfNodes());
    EXPECT_EQ(380u, mesh->getNumberOfElements());
    EXPECT_EQ(MeshLib::MeshElemType::PRISM, mesh->getElement(0)->getGeomType());
    EXPECT_TRUE(
        mesh->getProperties().hasPropertyVector("Conductivity in x-direction"));
    EXPECT_TRUE(
        mesh->getProperties().hasPropertyVector("Conductivity in y-direction"));
    EXPECT_TRUE(
        mesh->getProperties().hasPropertyVector("Conductivity in z-direction"));
    EXPECT_TRUE(mesh->getProperties().hasPropertyVector(
        "Storativity (drain- or fillable) or density ratio"));
    EXPECT_TRUE(mesh->getProperties().hasPropertyVector(
        "Specific Storage due to compressibility effects"));

    auto opt_Kx(mesh->getProperties().getPropertyVector<double>(
        "Conductivity in x-direction"));
    ASSERT_EQ(380u, opt_Kx->size());
    ASSERT_EQ(1e-6, (*opt_Kx)[0]);
    ASSERT_EQ(1e-5, (*opt_Kx)[8]);
    ASSERT_EQ(1e-7, (*opt_Kx)[76]);

    auto opt_Ky(mesh->getProperties().getPropertyVector<double>(
        "Conductivity in y-direction"));
    ASSERT_EQ(380u, opt_Ky->size());
    ASSERT_EQ(2e-6, (*opt_Ky)[0]);
    ASSERT_EQ(2e-5, (*opt_Ky)[8]);
    ASSERT_EQ(2e-7, (*opt_Ky)[76]);
}
