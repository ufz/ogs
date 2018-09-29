/**
 * @file TestMeshValidation.cpp
 * @author Karsten Rink
 * @date 2015-01-29
 * @brief Tests for MeshRevision class
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include <memory>

#include "MathLib/Point3d.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshQuality/MeshValidation.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshGenerators/RasterToMesh.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/Elements/Element.h"

#include "GeoLib/Raster.h"

void
detectHoles(MeshLib::Mesh const& mesh,
    std::vector<std::size_t> erase_elems,
    std::size_t const expected_n_holes)
{
    std::vector<MeshLib::Node*> nodes = MeshLib::copyNodeVector(mesh.getNodes());
    std::vector<MeshLib::Element*> elems = MeshLib::copyElementVector(mesh.getElements(),nodes);
    for (auto pos : erase_elems)
    {
        delete elems[pos];
        elems.erase(elems.begin()+pos);
    }
    MeshLib::Mesh mesh2("mesh2", nodes, elems);
    ASSERT_EQ(expected_n_holes, MeshLib::MeshValidation::detectHoles(mesh2));
};

TEST(MeshValidation, DetectHolesTri)
{
    std::array<double, 12> pix = {{0,0.1,0.2,0.1,0,0,0.1,0,0,0,-0.1,0}};
    GeoLib::RasterHeader const header = {
        4, 3, 1, MathLib::Point3d(std::array<double, 3>{{0, 0, 0}}), 1, -9999};
    GeoLib::Raster const raster(header ,pix.begin(), pix.end());
    std::unique_ptr<MeshLib::Mesh> mesh (MeshLib::RasterToMesh::convert(
        raster, MeshLib::MeshElemType::TRIANGLE, MeshLib::UseIntensityAs::ELEVATION));
    ASSERT_EQ(0, MeshLib::MeshValidation::detectHoles(*mesh));

    detectHoles(*mesh, {12}, 1);
    detectHoles(*mesh, {11, 11}, 1);
    detectHoles(*mesh, {10, 12}, 2);
}

TEST(MeshValidation, DetectHolesHex)
{
    auto mesh = std::unique_ptr<MeshLib::Mesh>{
        MeshLib::MeshGenerator::generateRegularHexMesh(
            5, 4, 4, 1.0, 1.0, 1.0, MathLib::ORIGIN, "mesh")};
    ASSERT_EQ(0, MeshLib::MeshValidation::detectHoles(*mesh));

    detectHoles(*mesh, {27}, 1);
    detectHoles(*mesh, {28, 27}, 1);
    detectHoles(*mesh, {29, 27}, 1);
}




