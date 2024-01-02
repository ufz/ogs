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

#include <cstdio>
#include <memory>

#include "GeoLib/IO/AsciiRasterInterface.h"
#include "GeoLib/Raster.h"
#include "InfoLib/TestInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshGenerators/RasterToMesh.h"
#include "MeshToolsLib/MeshInformation.h"

#ifdef OGS_BUILD_GUI
#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include "Applications/DataExplorer/VtkVis/VtkGeoImageSource.h"
#include "Applications/DataExplorer/VtkVis/VtkRaster.h"
#endif

class RasterToMeshTest : public ::testing::Test
{
public:
    RasterToMeshTest()
        : _file_name(TestInfoLib::TestInfo::data_path +
                     "/MeshLib/testraster_selke.asc")
    {
        _raster.reset(FileIO::AsciiRasterInterface::readRaster(_file_name));
    }

protected:
    std::size_t const _n_pix = 542;
    std::size_t const _n_nodes = 626;
    double _spacing = 1000;
    std::string const _file_name;
    std::unique_ptr<GeoLib::Raster> _raster;
};

TEST_F(RasterToMeshTest, convertRasterToTriMeshElevation)
{
    auto const mesh = MeshToolsLib::RasterToMesh::convert(
        *_raster, MeshLib::MeshElemType::TRIANGLE,
        MeshLib::UseIntensityAs::ELEVATION, "test");
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(_n_nodes, mesh->getNodes().size());
    ASSERT_EQ(_n_nodes, mesh->getNumberOfNodes());

    ASSERT_EQ(0, mesh->getProperties().size());

    auto const& n_types =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(2 * _n_pix, n_types.at(MeshLib::MeshElemType::TRIANGLE));

    GeoLib::AABB const aabb =
        MeshToolsLib::MeshInformation::getBoundingBox(*mesh);
    ASSERT_NEAR(aabb.getMinPoint()[2], 0,
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(aabb.getMaxPoint()[2], 0.07,
                std::numeric_limits<double>::epsilon());
}

TEST_F(RasterToMeshTest, convertRasterToQuadMeshElevation)
{
    auto const mesh = MeshToolsLib::RasterToMesh::convert(
        *_raster, MeshLib::MeshElemType::QUAD,
        MeshLib::UseIntensityAs::ELEVATION, "test");
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(_n_nodes, mesh->computeNumberOfBaseNodes());

    ASSERT_EQ(0, mesh->getProperties().size());

    auto const& n_types =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(_n_pix, n_types.at(MeshLib::MeshElemType::QUAD));

    GeoLib::AABB const aabb =
        MeshToolsLib::MeshInformation::getBoundingBox(*mesh);
    ASSERT_NEAR(aabb.getMinPoint()[2], 0,
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(aabb.getMaxPoint()[2], 0.07,
                std::numeric_limits<double>::epsilon());
}

TEST_F(RasterToMeshTest, convertRasterTo3DMeshElevation)
{
    auto const mesh = MeshToolsLib::RasterToMesh::convert(
        *_raster, MeshLib::MeshElemType::PRISM,
        MeshLib::UseIntensityAs::ELEVATION, "test");
    ASSERT_TRUE(mesh == nullptr);

    auto const mesh2 = MeshToolsLib::RasterToMesh::convert(
        *_raster, MeshLib::MeshElemType::HEXAHEDRON,
        MeshLib::UseIntensityAs::ELEVATION, "test");
    ASSERT_TRUE(mesh2 == nullptr);
}

TEST_F(RasterToMeshTest, convertRasterToTriMeshValue)
{
    auto const mesh = MeshToolsLib::RasterToMesh::convert(
        *_raster, MeshLib::MeshElemType::TRIANGLE,
        MeshLib::UseIntensityAs::DATAVECTOR, "test");
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(_n_nodes, mesh->computeNumberOfBaseNodes());

    ASSERT_EQ(1, mesh->getProperties().size());

    MeshLib::PropertyVector<double>* const prop =
        mesh->getProperties().getPropertyVector<double>("test");
    ASSERT_TRUE(prop != nullptr);
    ASSERT_EQ(2 * _n_pix, prop->size());

    auto const& bounds = MeshToolsLib::MeshInformation::getValueBounds(*prop);
    ASSERT_TRUE(bounds.has_value());
    ASSERT_NEAR(0, bounds->first, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.07, bounds->second, std::numeric_limits<double>::epsilon());

    for (MeshLib::Node* n : mesh->getNodes())
    {
        ASSERT_NEAR(0, (*n)[2], std::numeric_limits<double>::epsilon());
    }

    auto const& n_types =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(2 * _n_pix, n_types.at(MeshLib::MeshElemType::TRIANGLE));
}

TEST_F(RasterToMeshTest, convertRasterToQuadMeshValue)
{
    auto const mesh = MeshToolsLib::RasterToMesh::convert(
        *_raster, MeshLib::MeshElemType::QUAD,
        MeshLib::UseIntensityAs::DATAVECTOR, "test");
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(_n_nodes, mesh->computeNumberOfBaseNodes());

    ASSERT_EQ(1, mesh->getProperties().size());

    MeshLib::PropertyVector<double>* const prop =
        mesh->getProperties().getPropertyVector<double>("test");
    ASSERT_TRUE(prop != nullptr);
    ASSERT_EQ(_n_pix, prop->size());

    auto const& bounds = MeshToolsLib::MeshInformation::getValueBounds(*prop);
    ASSERT_TRUE(bounds.has_value());
    ASSERT_NEAR(0, bounds->first, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.07, bounds->second, std::numeric_limits<double>::epsilon());

    for (MeshLib::Node* n : mesh->getNodes())
    {
        ASSERT_TRUE((*n)[2] == 0);
    }

    auto const& n_types =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(_n_pix, n_types.at(MeshLib::MeshElemType::QUAD));
}

TEST_F(RasterToMeshTest, convertRasterToPrismMeshValue)
{
    auto const mesh = MeshToolsLib::RasterToMesh::convert(
        *_raster, MeshLib::MeshElemType::PRISM,
        MeshLib::UseIntensityAs::DATAVECTOR, "test");
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(2 * _n_nodes, mesh->computeNumberOfBaseNodes());

    ASSERT_EQ(1, mesh->getProperties().size());

    MeshLib::PropertyVector<double>* const prop =
        mesh->getProperties().getPropertyVector<double>("test");
    ASSERT_TRUE(prop != nullptr);
    ASSERT_EQ(2 * _n_pix, prop->size());

    auto const& bounds = MeshToolsLib::MeshInformation::getValueBounds(*prop);
    ASSERT_TRUE(bounds.has_value());
    ASSERT_NEAR(0, bounds->first, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.07, bounds->second, std::numeric_limits<double>::epsilon());

    for (MeshLib::Node* n : mesh->getNodes())
    {
        ASSERT_TRUE(((*n)[2] == 0) || ((*n)[2] == _spacing));
    }

    auto const& n_types =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(2 * _n_pix, n_types.at(MeshLib::MeshElemType::PRISM));
}

TEST_F(RasterToMeshTest, convertRasterToHexMeshValue)
{
    auto const mesh = MeshToolsLib::RasterToMesh::convert(
        *_raster, MeshLib::MeshElemType::HEXAHEDRON,
        MeshLib::UseIntensityAs::MATERIALS, "MaterialIDs");
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(2 * _n_nodes, mesh->computeNumberOfBaseNodes());

    ASSERT_EQ(1, mesh->getProperties().size());

    MeshLib::PropertyVector<int>* const prop =
        mesh->getProperties().getPropertyVector<int>("MaterialIDs");
    ASSERT_TRUE(prop != nullptr);
    ASSERT_EQ(_n_pix, prop->size());

    auto const& bounds = MeshToolsLib::MeshInformation::getValueBounds(*prop);
    ASSERT_TRUE(bounds.has_value());
    ASSERT_NEAR(0, bounds->first, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0, bounds->second, std::numeric_limits<double>::epsilon());

    for (MeshLib::Node* n : mesh->getNodes())
    {
        ASSERT_TRUE(((*n)[2] == 0) || ((*n)[2] == _spacing));
    }

    auto const& n_types =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(_n_pix, n_types.at(MeshLib::MeshElemType::HEXAHEDRON));
}

TEST_F(RasterToMeshTest, convertRasterToQuadMeshNone)
{
    auto const mesh = MeshToolsLib::RasterToMesh::convert(
        *_raster, MeshLib::MeshElemType::QUAD, MeshLib::UseIntensityAs::NONE,
        "test");
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(_n_nodes, mesh->computeNumberOfBaseNodes());

    ASSERT_EQ(0, mesh->getProperties().size());

    for (MeshLib::Node* n : mesh->getNodes())
    {
        ASSERT_TRUE((*n)[2] == 0);
    }

    auto const& n_types =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(_n_pix, n_types.at(MeshLib::MeshElemType::QUAD));
}

#ifdef OGS_BUILD_GUI
TEST_F(RasterToMeshTest, vtkImage)
{
    vtkSmartPointer<vtkImageAlgorithm> const vtk_raster =
        VtkRaster::loadImage(_file_name);
    vtkSmartPointer<VtkGeoImageSource> geoImage = VtkGeoImageSource::New();
    geoImage->setImage(vtk_raster, "vtkImage");
    geoImage->Update();
    auto const ascii_raster = VtkGeoImageSource::convertToRaster(geoImage);

    ASSERT_TRUE(ascii_raster != std::nullopt);
    auto const mesh = MeshToolsLib::RasterToMesh::convert(
        *ascii_raster, MeshLib::MeshElemType::TRIANGLE,
        MeshLib::UseIntensityAs::DATAVECTOR, "test");
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(_n_nodes, mesh->computeNumberOfBaseNodes());

    ASSERT_EQ(1, mesh->getProperties().size());

    MeshLib::PropertyVector<double>* const prop =
        mesh->getProperties().getPropertyVector<double>("test");
    ASSERT_TRUE(prop != nullptr);
    ASSERT_EQ(2 * _n_pix, prop->size());

    auto const& bounds = MeshToolsLib::MeshInformation::getValueBounds(*prop);
    ASSERT_TRUE(bounds.has_value());
    ASSERT_NEAR(0, bounds->first, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.07, bounds->second, std::numeric_limits<float>::epsilon());

    for (MeshLib::Node* const n : mesh->getNodes())
    {
        ASSERT_TRUE((*n)[2] == 0);
    }

    auto const& n_types =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(2 * _n_pix, n_types.at(MeshLib::MeshElemType::TRIANGLE));
}
#endif
