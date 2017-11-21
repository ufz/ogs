/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <cstdio>
#include <memory>

#include "gtest/gtest.h"

#include "BaseLib/BuildInfo.h"

#include "Applications/FileIO/AsciiRasterInterface.h"
#include "GeoLib/Raster.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/RasterToMesh.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/Node.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"

#ifdef OGS_BUILD_GUI
#include "Applications/DataExplorer/VtkVis/VtkGeoImageSource.h"
#include "Applications/DataExplorer/VtkVis/VtkRaster.h"

#include <vtkImageData.h>
#endif

class RasterToMeshTest : public ::testing::Test
{
public:
    RasterToMeshTest()
        : _file_name(BaseLib::BuildInfo::data_path + "/MeshLib/testraster_selke.asc"),
          _mesh_name(BaseLib::BuildInfo::data_path + "/MeshLib/testraster_selke.vtu")
    {
        _raster.reset(FileIO::AsciiRasterInterface::readRaster(_file_name));
    }

    ~RasterToMeshTest() { std::remove(_mesh_name.c_str());  }

protected:
    std::size_t const _n_pix = 542;
    std::size_t const _n_nodes = 626;
    double _spacing = 1000;
    std::string const _file_name;
    std::string const _mesh_name;
    std::unique_ptr<GeoLib::Raster> _raster;
};

TEST_F(RasterToMeshTest, convertRasterToTriMeshElevation)
{
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::RasterToMesh::convert(
        *_raster, MeshLib::MeshElemType::TRIANGLE,
        MeshLib::UseIntensityAs::ELEVATION, "test"));
    ASSERT_TRUE(mesh != nullptr);

    MeshLib::IO::VtuInterface vtkio(mesh.get(), 0, false);
    std::string name(BaseLib::BuildInfo::data_path +
                     "/MeshLib/testraster_selke.vtu");
    vtkio.writeToFile(name);

    std::cout << mesh->getNodes().size() << ", " << mesh->getNumberOfNodes()
              << std::endl;
    ASSERT_EQ(_n_nodes, mesh->getNodes().size());
    ASSERT_EQ(_n_nodes, mesh->getNumberOfNodes());

    std::vector<std::string> names =
        mesh->getProperties().getPropertyVectorNames();
    ASSERT_TRUE(names.empty());

    std::array<unsigned, 7> n_types =
        MeshLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(2 * _n_pix, n_types[1]);

    GeoLib::AABB aabb = MeshLib::MeshInformation::getBoundingBox(*mesh);
    ASSERT_NEAR(aabb.getMinPoint()[2], 0,
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(aabb.getMaxPoint()[2], 0.07,
                std::numeric_limits<double>::epsilon());
}

TEST_F(RasterToMeshTest, convertRasterToQuadMeshElevation)
{
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::RasterToMesh::convert(
        *_raster, MeshLib::MeshElemType::QUAD,  MeshLib::UseIntensityAs::ELEVATION, "test"));
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(_n_nodes, mesh->getNumberOfBaseNodes());

    std::vector<std::string> names = mesh->getProperties().getPropertyVectorNames();
    ASSERT_TRUE(names.empty());

    std::array<unsigned, 7> n_types =
        MeshLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(_n_pix, n_types[2]);

    GeoLib::AABB aabb = MeshLib::MeshInformation::getBoundingBox(*mesh);
    ASSERT_NEAR(aabb.getMinPoint()[2], 0,
                std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(aabb.getMaxPoint()[2], 0.07,
                std::numeric_limits<double>::epsilon());
}

TEST_F(RasterToMeshTest, convertRasterTo3DMeshElevation)
{
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::RasterToMesh::convert(
        *_raster, MeshLib::MeshElemType::PRISM,
        MeshLib::UseIntensityAs::ELEVATION, "test"));
    ASSERT_TRUE(mesh == nullptr);

    std::unique_ptr<MeshLib::Mesh> mesh2(MeshLib::RasterToMesh::convert(
        *_raster.get(), MeshLib::MeshElemType::HEXAHEDRON,
        MeshLib::UseIntensityAs::ELEVATION, "test"));
    ASSERT_TRUE(mesh2 == nullptr);
}

TEST_F(RasterToMeshTest, convertRasterToTriMeshValue)
{
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::RasterToMesh::convert(
        *_raster.get(), MeshLib::MeshElemType::TRIANGLE,
        MeshLib::UseIntensityAs::DATAVECTOR, "test"));
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(_n_nodes, mesh->getNumberOfBaseNodes());

    std::vector<std::string> names =
        mesh->getProperties().getPropertyVectorNames();
    ASSERT_EQ(1, names.size());

    MeshLib::PropertyVector<double>* prop =
        mesh->getProperties().getPropertyVector<double>("test");
    ASSERT_EQ(2 * _n_pix, prop->size());

    std::pair<double, double> const& bounds =
        MeshLib::MeshInformation::getValueBounds<double>(*mesh, "test");
    ASSERT_NEAR(0, bounds.first, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.07, bounds.second, std::numeric_limits<double>::epsilon());

    std::vector<MeshLib::Node*> const& nodes = mesh->getNodes();
    for (MeshLib::Node* n : nodes)
        ASSERT_NEAR(0, (*n)[2], std::numeric_limits<double>::epsilon());

    std::array<unsigned, 7> n_types =
        MeshLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(2 * _n_pix, n_types[1]);
}

TEST_F(RasterToMeshTest, convertRasterToQuadMeshValue)
{
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::RasterToMesh::convert(
        *_raster.get(), MeshLib::MeshElemType::QUAD,
        MeshLib::UseIntensityAs::DATAVECTOR, "test"));
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(_n_nodes, mesh->getNumberOfBaseNodes());

    std::vector<std::string> names =
        mesh->getProperties().getPropertyVectorNames();
    ASSERT_EQ(1, names.size());

    MeshLib::PropertyVector<double>* prop =
        mesh->getProperties().getPropertyVector<double>("test");
    ASSERT_EQ(_n_pix, prop->size());

    std::pair<double, double> const& bounds =
        MeshLib::MeshInformation::getValueBounds<double>(*mesh, "test");
    ASSERT_NEAR(0, bounds.first, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.07, bounds.second, std::numeric_limits<double>::epsilon());

    std::vector<MeshLib::Node*> const& nodes = mesh->getNodes();
    for (MeshLib::Node* n : nodes)
        ASSERT_TRUE((*n)[2] == 0);

    std::array<unsigned, 7> n_types =
        MeshLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(_n_pix, n_types[2]);
}

TEST_F(RasterToMeshTest, convertRasterToPrismMeshValue)
{
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::RasterToMesh::convert(
        *_raster.get(), MeshLib::MeshElemType::PRISM,
        MeshLib::UseIntensityAs::DATAVECTOR, "test"));
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(2 * _n_nodes, mesh->getNumberOfBaseNodes());

    std::vector<std::string> names =
        mesh->getProperties().getPropertyVectorNames();
    ASSERT_EQ(1, names.size());

    MeshLib::PropertyVector<double>* prop =
        mesh->getProperties().getPropertyVector<double>("test");
    ASSERT_EQ(2 * _n_pix, prop->size());

    std::pair<double, double> const& bounds =
        MeshLib::MeshInformation::getValueBounds<double>(*mesh, "test");
    ASSERT_NEAR(0, bounds.first, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.07, bounds.second, std::numeric_limits<double>::epsilon());

    std::vector<MeshLib::Node*> const& nodes = mesh->getNodes();
    for (MeshLib::Node* n : nodes)
        ASSERT_TRUE(((*n)[2] == 0) || ((*n)[2] == _spacing));

    std::array<unsigned, 7> n_types =
        MeshLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(2 * _n_pix, n_types[6]);
}

TEST_F(RasterToMeshTest, convertRasterToHexMeshValue)
{
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::RasterToMesh::convert(
        *_raster.get(), MeshLib::MeshElemType::HEXAHEDRON,
        MeshLib::UseIntensityAs::MATERIALS, "MaterialIDs"));
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(2 * _n_nodes, mesh->getNumberOfBaseNodes());

    std::vector<std::string> names =
        mesh->getProperties().getPropertyVectorNames();
    ASSERT_EQ(1, names.size());

    MeshLib::PropertyVector<int>* prop =
        mesh->getProperties().getPropertyVector<int>("MaterialIDs");
    ASSERT_EQ(_n_pix, prop->size());

    std::pair<int, int> const& bounds =
        MeshLib::MeshInformation::getValueBounds<int>(*mesh, "MaterialIDs");
    ASSERT_NEAR(0, bounds.first, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0, bounds.second, std::numeric_limits<double>::epsilon());

    std::vector<MeshLib::Node*> const& nodes = mesh->getNodes();
    for (MeshLib::Node* n : nodes)
        ASSERT_TRUE(((*n)[2] == 0) || ((*n)[2] == _spacing));

    std::array<unsigned, 7> n_types =
        MeshLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(_n_pix, n_types[4]);
}

TEST_F(RasterToMeshTest, convertRasterToQuadMeshNone)
{
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::RasterToMesh::convert(
        *_raster.get(), MeshLib::MeshElemType::QUAD,
        MeshLib::UseIntensityAs::NONE, "test"));
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(_n_nodes, mesh->getNumberOfBaseNodes());

    std::vector<std::string> names =
        mesh->getProperties().getPropertyVectorNames();
    ASSERT_TRUE(names.empty());

    std::vector<MeshLib::Node*> const& nodes = mesh->getNodes();
    for (MeshLib::Node* n : nodes)
        ASSERT_TRUE((*n)[2] == 0);

    std::array<unsigned, 7> n_types =
        MeshLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(_n_pix, n_types[2]);
}

#ifdef OGS_BUILD_GUI
TEST_F(RasterToMeshTest, vtkImage)
{
    double x0, y0, spacing;
    vtkImageAlgorithm* raster =
        VtkRaster::loadImage(_file_name, x0, y0, spacing);
    double origin[3];
    raster->GetOutput()->GetOrigin(origin);

    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::RasterToMesh::convert(
        raster->GetOutput(), origin, spacing, MeshLib::MeshElemType::TRIANGLE,
        MeshLib::UseIntensityAs::DATAVECTOR, "test"));
    ASSERT_TRUE(mesh != nullptr);

    ASSERT_EQ(_n_nodes, mesh->getNumberOfBaseNodes());

    std::vector<std::string> names =
        mesh->getProperties().getPropertyVectorNames();
    ASSERT_EQ(1, names.size());

    MeshLib::PropertyVector<double>* prop =
        mesh->getProperties().getPropertyVector<double>("test");
    ASSERT_EQ(2 * _n_pix, prop->size());

    std::pair<double, double> const& bounds =
        MeshLib::MeshInformation::getValueBounds<double>(*mesh, "test");
    ASSERT_NEAR(0, bounds.first, std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.07, bounds.second, std::numeric_limits<float>::epsilon());

    std::vector<MeshLib::Node*> const& nodes = mesh->getNodes();
    for (MeshLib::Node* n : nodes)
        ASSERT_TRUE((*n)[2] == 0);

    std::array<unsigned, 7> n_types =
        MeshLib::MeshInformation::getNumberOfElementTypes(*mesh);
    ASSERT_EQ(2 * _n_pix, n_types[1]);
}
#endif
