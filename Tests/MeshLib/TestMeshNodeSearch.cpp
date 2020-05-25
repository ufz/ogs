/**
 * \file
 * \date Oct 28, 2013
 * \brief Test the implementation of class MeshNodeSearch.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "gtest/gtest.h"

#include <memory>

#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/HeuristicSearchLength.h"

using namespace MeshLib;

class MeshLibMeshNodeSearchInSimpleQuadMesh : public testing::Test
{
public:
    MeshLibMeshNodeSearchInSimpleQuadMesh()
        : quad_mesh_(MeshGenerator::generateRegularQuadMesh(
              geometric_size_, number_of_subdivisions_per_direction_))
    {}

    ~MeshLibMeshNodeSearchInSimpleQuadMesh() override { delete quad_mesh_; }

protected:
    const double geometric_size_{10.0};
    const std::size_t number_of_subdivisions_per_direction_{99};
    Mesh* quad_mesh_;
};

class MeshLibMeshNodeSearchInSimpleHexMesh : public testing::Test
{
public:
    MeshLibMeshNodeSearchInSimpleHexMesh()
        : hex_mesh_(MeshGenerator::generateRegularHexMesh(
              geometric_size_, number_of_subdivisions_per_direction_))
    {}

    ~MeshLibMeshNodeSearchInSimpleHexMesh() override { delete hex_mesh_; }

protected:
    const double geometric_size_{10.0};
    const std::size_t number_of_subdivisions_per_direction_{10};
    Mesh* hex_mesh_;
};

TEST_F(MeshLibMeshNodeSearchInSimpleQuadMesh, PointSearchEpsHalfEdge)
{
    double dx_ = geometric_size_ / number_of_subdivisions_per_direction_;
    double dx_half = dx_*0.5;

    ASSERT_TRUE(quad_mesh_ != nullptr);

    // 2 perform search and compare results with expected vals
    auto search_length =
        std::make_unique<MeshGeoToolsLib::SearchLength>(dx_half);
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(*quad_mesh_,
        std::move(search_length), MeshGeoToolsLib::SearchAllNodes::Yes);

    GeoLib::Point p1(0.0, 0.0, 0.0);
    EXPECT_EQ(1u, mesh_node_searcher.getMeshNodeIDsForPoint(p1).size());
    EXPECT_EQ(0u, mesh_node_searcher.getMeshNodeIDsForPoint(p1)[0]);

    GeoLib::Point p2(dx_half*0.99, 0.0, 0.0);
    EXPECT_EQ(1u, mesh_node_searcher.getMeshNodeIDsForPoint(p2).size());
    EXPECT_EQ(0u, mesh_node_searcher.getMeshNodeIDsForPoint(p2)[0]);

    GeoLib::Point p3(dx_half, 0.0, 0.0);
    EXPECT_EQ(0u, mesh_node_searcher.getMeshNodeIDsForPoint(p3).size());

    GeoLib::Point p4(dx_half*1.01, 0.0, 0.0);
    ASSERT_EQ(1u, mesh_node_searcher.getMeshNodeIDsForPoint(p4).size());
    ASSERT_EQ(1u, mesh_node_searcher.getMeshNodeIDsForPoint(p4)[0]);
}

TEST_F(MeshLibMeshNodeSearchInSimpleQuadMesh, PointSearchZeroEps)
{
    ASSERT_TRUE(quad_mesh_ != nullptr);
    // 1 create a geometry

    // 2 perform search and compare results with expected vals
    auto search_length = std::make_unique<MeshGeoToolsLib::SearchLength>();
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(*quad_mesh_,
        std::move(search_length), MeshGeoToolsLib::SearchAllNodes::Yes);

    // find ORIGIN
    GeoLib::Point pnt1(0.0, 0.0, 0.0);
    ASSERT_EQ(0u, mesh_node_searcher.getMeshNodeIDsForPoint(pnt1)[0]);

    GeoLib::Point pnt2(0.049, 0.049, 0.0);
    ASSERT_EQ(0u, mesh_node_searcher.getMeshNodeIDsForPoint(pnt2).size());

    GeoLib::Point pnt3(0.051, 0.049, 0.0);
    ASSERT_EQ(0u, mesh_node_searcher.getMeshNodeIDsForPoint(pnt3).size());

    GeoLib::Point pnt4(0.049, 0.051, 0.0);
    ASSERT_EQ(0u, mesh_node_searcher.getMeshNodeIDsForPoint(pnt4).size());

    GeoLib::Point pnt5(0.051, 0.051, 0.0);
    ASSERT_EQ(0u, mesh_node_searcher.getMeshNodeIDsForPoint(pnt5).size());

    GeoLib::Point pnt6(10, 10, 0.0);
    ASSERT_EQ(1u, mesh_node_searcher.getMeshNodeIDsForPoint(pnt6).size());
    EXPECT_EQ(
        pow(number_of_subdivisions_per_direction_+1,2)-1,
        mesh_node_searcher.getMeshNodeIDsForPoint(pnt6)[0]);
}

TEST_F(MeshLibMeshNodeSearchInSimpleQuadMesh, PolylineSearch)
{
    ASSERT_TRUE(quad_mesh_ != nullptr);
    // create geometry
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(geometric_size_, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, 0.049, 0.0));
    pnts.push_back(new GeoLib::Point(geometric_size_, 0.049, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, geometric_size_, 0.0));
    pnts.push_back(new GeoLib::Point(geometric_size_, geometric_size_, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, geometric_size_ - 0.049, 0.0));
    pnts.push_back(new GeoLib::Point(geometric_size_, geometric_size_ - 0.049, 0.0));

    GeoLib::Polyline ply0(pnts);
    ply0.addPoint(0);
    ply0.addPoint(1);

    // perform search and compare results with expected vals
    auto search_length =
        std::make_unique<MeshGeoToolsLib::HeuristicSearchLength>(*quad_mesh_);
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(*quad_mesh_,
        std::move(search_length), MeshGeoToolsLib::SearchAllNodes::Yes);
    std::vector<std::size_t> const& found_ids_ply0(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply0));

    ASSERT_EQ(100u, found_ids_ply0.size());
    for (std::size_t k(0); k < found_ids_ply0.size(); k++)
    {
        ASSERT_EQ(k, found_ids_ply0[k]);
    }

    GeoLib::Polyline ply1(pnts);
    ply1.addPoint(2);
    ply1.addPoint(3);
    std::vector<std::size_t> const& found_ids_ply1(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply1));

    ASSERT_EQ(100u, found_ids_ply1.size());
    for (std::size_t k(0); k < found_ids_ply1.size(); k++)
    {
        ASSERT_EQ(k, found_ids_ply1[k]);
    }

    GeoLib::Polyline ply2(pnts);
    ply2.addPoint(4);
    ply2.addPoint(5);
    std::vector<std::size_t> const& found_ids_ply2(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply2));

    std::size_t offset((number_of_subdivisions_per_direction_+1)*number_of_subdivisions_per_direction_);
    ASSERT_EQ(100u, found_ids_ply2.size());
    for (std::size_t k(0); k < found_ids_ply2.size(); k++)
    {
        ASSERT_EQ(offset + k, found_ids_ply2[k]);
    }

    GeoLib::Polyline ply3(pnts);
    ply3.addPoint(6);
    ply3.addPoint(7);
    std::vector<std::size_t> const& found_ids_ply3(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply3));

    ASSERT_EQ(100u, found_ids_ply3.size());
    for (std::size_t k(0); k < found_ids_ply3.size(); k++)
    {
        ASSERT_EQ(offset + k, found_ids_ply3[k]);
    }

    // left border
    GeoLib::Polyline ply4(pnts);
    ply4.addPoint(0);
    ply4.addPoint(6);
    std::vector<std::size_t> const& found_ids_ply4(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply4));

    ASSERT_EQ(100u, found_ids_ply4.size());
    for (std::size_t k(0); k < found_ids_ply4.size(); k++)
    {
        ASSERT_EQ(k * (number_of_subdivisions_per_direction_ + 1),
                  found_ids_ply4[k]);
    }

    // right border
    GeoLib::Polyline ply5(pnts);
    ply5.addPoint(1);
    ply5.addPoint(7);
    std::vector<std::size_t> const& found_ids_ply5(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply5));

    ASSERT_EQ(100u, found_ids_ply5.size());
    for (std::size_t k(0); k < found_ids_ply5.size(); k++)
    {
        ASSERT_EQ(k * (number_of_subdivisions_per_direction_ + 1) +
                      number_of_subdivisions_per_direction_,
                  found_ids_ply5[k]);
    }

    // diagonal
    GeoLib::Polyline ply6(pnts);
    ply6.addPoint(0);
    ply6.addPoint(5);
    std::vector<std::size_t> const& found_ids_ply6(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply6));

    ASSERT_EQ(100u, found_ids_ply6.size());
    for (std::size_t k(0); k < found_ids_ply6.size(); k++)
    {
        ASSERT_EQ(k * (number_of_subdivisions_per_direction_ + 1) + k,
                  found_ids_ply6[k]);
    }

    std::for_each(pnts.begin(), pnts.end(), [](GeoLib::Point* pnt) { delete pnt; });
}

TEST_F(MeshLibMeshNodeSearchInSimpleQuadMesh, SurfaceSearch)
{
    ASSERT_TRUE(quad_mesh_ != nullptr);
    // create geometry
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(geometric_size_, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(geometric_size_, geometric_size_, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, geometric_size_, 0.0));
    pnts.push_back(new GeoLib::Point(geometric_size_, 0.5*geometric_size_, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, 0.5*geometric_size_, 0.0));

    auto search_length = std::make_unique<MeshGeoToolsLib::SearchLength>();
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(*quad_mesh_,
        std::move(search_length), MeshGeoToolsLib::SearchAllNodes::Yes);

    // entire domain
    GeoLib::Surface sfc0(pnts);
    sfc0.addTriangle(0, 1, 2);
    sfc0.addTriangle(0, 2, 3);

    std::vector<std::size_t> const& found_ids_sfc0(
        mesh_node_searcher.getMeshNodeIDsAlongSurface(sfc0));

    ASSERT_EQ(quad_mesh_->getNumberOfNodes(), found_ids_sfc0.size());
    for (std::size_t k(0); k < found_ids_sfc0.size(); k++)
    {
        ASSERT_EQ(k, found_ids_sfc0[k]);
    }

    // bottom half domain
    GeoLib::Surface sfc1(pnts);
    sfc1.addTriangle(0, 1, 4);
    sfc1.addTriangle(0, 4, 5);

    std::vector<std::size_t> const& found_ids_sfc1(
        mesh_node_searcher.getMeshNodeIDsAlongSurface(sfc1));

    ASSERT_EQ(quad_mesh_->getNumberOfNodes()/2, found_ids_sfc1.size());
    for (std::size_t k(0); k < found_ids_sfc1.size(); k++)
    {
        ASSERT_EQ(k, found_ids_sfc1[k]);
    }

    std::for_each(pnts.begin(), pnts.end(), [](GeoLib::Point* pnt) { delete pnt; });
}

TEST_F(MeshLibMeshNodeSearchInSimpleHexMesh, SurfaceSearch)
{
    ASSERT_TRUE(hex_mesh_ != nullptr);
    // create geometry
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(geometric_size_, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(geometric_size_, geometric_size_, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, geometric_size_, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, 0.0, geometric_size_));
    pnts.push_back(new GeoLib::Point(geometric_size_, 0.0, geometric_size_));
    pnts.push_back(new GeoLib::Point(geometric_size_, geometric_size_, geometric_size_));
    pnts.push_back(new GeoLib::Point(0.0, geometric_size_, geometric_size_));

    auto search_length = std::make_unique<MeshGeoToolsLib::SearchLength>();
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(*hex_mesh_,
        std::move(search_length), MeshGeoToolsLib::SearchAllNodes::Yes);

    const std::size_t n_nodes_1d = number_of_subdivisions_per_direction_ + 1;
    const std::size_t n_nodes_2d = n_nodes_1d * n_nodes_1d;
    const std::size_t n_nodes_3d = n_nodes_1d * n_nodes_1d * n_nodes_1d;

    // bottom surface
    GeoLib::Surface sfc_bottom(pnts);
    sfc_bottom.addTriangle(0, 1, 2);
    sfc_bottom.addTriangle(0, 2, 3);

    std::vector<std::size_t> const& found_ids_sfc_b(
        mesh_node_searcher.getMeshNodeIDsAlongSurface(sfc_bottom));
    ASSERT_EQ(n_nodes_2d, found_ids_sfc_b.size());
    for (std::size_t k(0); k < found_ids_sfc_b.size(); k++)
    {
        ASSERT_EQ(k, found_ids_sfc_b[k]);
    }

    // top surface
    GeoLib::Surface sfc_top(pnts);
    sfc_top.addTriangle(4, 5, 6);
    sfc_top.addTriangle(4, 6, 7);

    std::vector<std::size_t> const& found_ids_sfc_t(
        mesh_node_searcher.getMeshNodeIDsAlongSurface(sfc_top));
    ASSERT_EQ(n_nodes_2d, found_ids_sfc_t.size());
    const std::size_t offset_t = n_nodes_3d - n_nodes_2d;
    for (std::size_t k(0); k < found_ids_sfc_t.size(); k++)
    {
        ASSERT_EQ(offset_t + k, found_ids_sfc_t[k]);
    }

    // front
    GeoLib::Surface sfc_front(pnts);
    sfc_front.addTriangle(0, 1, 5);
    sfc_front.addTriangle(0, 5, 4);

    std::vector<std::size_t> const& found_ids_sfc_f(
        mesh_node_searcher.getMeshNodeIDsAlongSurface(sfc_front));
    ASSERT_EQ(n_nodes_2d, found_ids_sfc_f.size());
    std::size_t cnt=0;
    for (std::size_t k(0); k<n_nodes_1d; k++) {
        for (std::size_t i(0); i < n_nodes_1d; i++)
        {
            ASSERT_EQ(k * n_nodes_2d + i, found_ids_sfc_f[cnt++]);
        }
    }

    // back
    GeoLib::Surface sfc_back(pnts);
    sfc_back.addTriangle(3, 2, 6);
    sfc_back.addTriangle(3, 6, 7);

    std::vector<std::size_t> const& found_ids_sfc_back(
        mesh_node_searcher.getMeshNodeIDsAlongSurface(sfc_back));
    ASSERT_EQ(n_nodes_2d, found_ids_sfc_back.size());
    cnt = 0;
    const std::size_t y_offset = n_nodes_1d*(n_nodes_1d-1);
    for (std::size_t k(0); k<n_nodes_1d; k++) {
        const std::size_t z_offset = n_nodes_2d*k;
        for (std::size_t i(0); i < n_nodes_1d; i++)
        {
            ASSERT_EQ(z_offset + y_offset + i, found_ids_sfc_back[cnt++]);
        }
    }

    // left
    GeoLib::Surface sfc_left(pnts);
    sfc_left.addTriangle(0, 3, 7);
    sfc_left.addTriangle(0, 7, 4);

    std::vector<std::size_t> const& found_ids_sfc_left(
        mesh_node_searcher.getMeshNodeIDsAlongSurface(sfc_left));
    ASSERT_EQ(n_nodes_2d, found_ids_sfc_left.size());
    cnt = 0;
    for (std::size_t k(0); k<n_nodes_1d; k++) {
        const std::size_t z_offset = n_nodes_2d*k;
        for (std::size_t j(0); j < n_nodes_1d; j++)
        {
            ASSERT_EQ(z_offset + j * n_nodes_1d, found_ids_sfc_left[cnt++]);
        }
    }

    // right
    GeoLib::Surface sfc_right(pnts);
    sfc_right.addTriangle(1, 2, 6);
    sfc_right.addTriangle(1, 6, 5);

    std::vector<std::size_t> const& found_ids_sfc_right(
        mesh_node_searcher.getMeshNodeIDsAlongSurface(sfc_right));
    ASSERT_EQ(n_nodes_2d, found_ids_sfc_right.size());
    cnt = 0;
    for (std::size_t k(0); k<n_nodes_1d; k++) {
        const std::size_t z_offset = n_nodes_2d*k;
        for (std::size_t j(0); j < n_nodes_1d; j++)
        {
            ASSERT_EQ(z_offset + (j + 1) * n_nodes_1d - 1,
                      found_ids_sfc_right[cnt++]);
        }
    }


    std::for_each(pnts.begin(), pnts.end(), [](GeoLib::Point* pnt) { delete pnt; });
}

