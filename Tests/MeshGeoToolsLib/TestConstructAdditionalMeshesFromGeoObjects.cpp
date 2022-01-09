/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <ctime>
#include <memory>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "MeshGeoToolsLib/ConstructMeshesFromGeometries.h"
#include "MeshGeoToolsLib/SearchLength.h"
#include "MeshLib/Mesh.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"
#endif

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "Tests/GeoLib/CreateTestPoints.h"

TEST(ConstructAdditionalMeshesFromGeoObjects, PointMesh)
{
    // create 10x10x10 mesh using 1000 hexahedra
    const double length = 10.0;
    const std::size_t n_subdivisions = 10;

#ifdef USE_PETSC
    std::unique_ptr<MeshLib::NodePartitionedMesh> mesh =
        std::make_unique<MeshLib::NodePartitionedMesh>(
            *MeshLib::MeshGenerator::generateRegularHexMesh(length,
                                                            n_subdivisions));
#else
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::MeshGenerator::generateRegularHexMesh(length, n_subdivisions));
#endif

    // create geometry: for every mesh node exactly one point
    GeoLib::GEOObjects geometries;
    GeoLib::Point shift(0.0, 0.0, 0.0);
    std::string geometry_name("AllMeshNodes");
    // all points have a name
    int const points_per_edge = 10;
    createSetOfTestPointsAndAssociatedNames(geometries, geometry_name,
                                            points_per_edge, shift);

    // construct meshes from the points
    double const search_length = std::numeric_limits<double>::epsilon();
    bool const multiple_nodes_allowed = false;
    auto const meshes =
        MeshGeoToolsLib::constructAdditionalMeshesFromGeoObjects(
            geometries, *mesh,
            std::make_unique<MeshGeoToolsLib::SearchLength>(search_length),
            multiple_nodes_allowed);

    // the number of constructed meshes have to equal the number of points
    ASSERT_EQ(geometries.getPointVec(geometry_name)->size(), meshes.size());
}

TEST(ConstructAdditionalMeshesFromGeoObjects, PointMeshLargeSearchRadius)
{
    // create 10x10x10 mesh using 1000 hexahedra
    const double length = 10.0;
    const std::size_t n_subdivisions = 10;
#ifdef USE_PETSC

    std::unique_ptr<MeshLib::NodePartitionedMesh> mesh =
        std::make_unique<MeshLib::NodePartitionedMesh>(
            *MeshLib::MeshGenerator::generateRegularHexMesh(length,
                                                            n_subdivisions));
#else
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::MeshGenerator::generateRegularHexMesh(length, n_subdivisions));
#endif

    // create geometry: for every mesh node exactly one point
    GeoLib::GEOObjects geometries;
    GeoLib::Point shift(0.0, 0.0, 0.0);
    std::string geometry_name("AllMeshNodes");
    // all points have a name
    int const points_per_edge = 2;
    createSetOfTestPointsAndAssociatedNames(geometries, geometry_name,
                                            points_per_edge, shift);

    // construct meshes from the points
    double const search_length = 1 + std::numeric_limits<double>::epsilon();
    bool multiple_nodes_allowed = false;

    EXPECT_ANY_THROW(
        auto const meshes =
            MeshGeoToolsLib::constructAdditionalMeshesFromGeoObjects(
                geometries, *mesh,
                std::make_unique<MeshGeoToolsLib::SearchLength>(search_length),
                multiple_nodes_allowed));

    multiple_nodes_allowed = true;
    auto const meshes_from_multiple_nodes =
        MeshGeoToolsLib::constructAdditionalMeshesFromGeoObjects(
            geometries, *mesh,
            std::make_unique<MeshGeoToolsLib::SearchLength>(search_length),
            multiple_nodes_allowed);

    // the number of constructed meshes have to equal the number of points
    ASSERT_EQ(geometries.getPointVec(geometry_name)->size(),
              meshes_from_multiple_nodes.size());
}
