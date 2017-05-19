/**
 * @file TestMeshNodeSearch.cpp
 * @author git blame TestMeshNodeSearch.cpp
 * @date Oct 28, 2013
 * @brief Test the implementation of class MeshNodeSearch.
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
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
    MeshLibMeshNodeSearchInSimpleQuadMesh() :
        _geometric_size(10.0), _number_of_subdivisions_per_direction(99),
        _quad_mesh(MeshGenerator::generateRegularQuadMesh(_geometric_size, _number_of_subdivisions_per_direction))
    {}

    ~MeshLibMeshNodeSearchInSimpleQuadMesh() override { delete _quad_mesh; }

protected:
    const double _geometric_size;
    const std::size_t _number_of_subdivisions_per_direction;
    Mesh* _quad_mesh;
};

class MeshLibMeshNodeSearchInSimpleHexMesh : public testing::Test
{
public:
    MeshLibMeshNodeSearchInSimpleHexMesh() :
        _geometric_size(10.0), _number_of_subdivisions_per_direction(10),
        _hex_mesh(MeshGenerator::generateRegularHexMesh(_geometric_size, _number_of_subdivisions_per_direction))
    {}

    ~MeshLibMeshNodeSearchInSimpleHexMesh() override { delete _hex_mesh; }

protected:
    const double _geometric_size;
    const std::size_t _number_of_subdivisions_per_direction;
    Mesh* _hex_mesh;
};

TEST_F(MeshLibMeshNodeSearchInSimpleQuadMesh, PointSearchEpsHalfEdge)
{
    double _dx = _geometric_size / _number_of_subdivisions_per_direction;
    double dx_half = _dx*0.5;

    ASSERT_TRUE(_quad_mesh != nullptr);

    // 2 perform search and compare results with expected vals
    std::unique_ptr<MeshGeoToolsLib::SearchLength> search_length{
        new MeshGeoToolsLib::SearchLength(dx_half)};
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(*_quad_mesh,
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
    ASSERT_TRUE(_quad_mesh != nullptr);
    // 1 create a geometry

    // 2 perform search and compare results with expected vals
    std::unique_ptr<MeshGeoToolsLib::SearchLength> search_length{
        new MeshGeoToolsLib::SearchLength};
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(*_quad_mesh,
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
        pow(_number_of_subdivisions_per_direction+1,2)-1,
        mesh_node_searcher.getMeshNodeIDsForPoint(pnt6)[0]);
}

TEST_F(MeshLibMeshNodeSearchInSimpleQuadMesh, PolylineSearch)
{
    ASSERT_TRUE(_quad_mesh != nullptr);
    // create geometry
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, 0.049, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, 0.049, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, _geometric_size, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, _geometric_size, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, _geometric_size - 0.049, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, _geometric_size - 0.049, 0.0));

    GeoLib::Polyline ply0(pnts);
    ply0.addPoint(0);
    ply0.addPoint(1);

    // perform search and compare results with expected vals
    std::unique_ptr<MeshGeoToolsLib::HeuristicSearchLength> search_length{
        new MeshGeoToolsLib::HeuristicSearchLength(*_quad_mesh)};
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(*_quad_mesh,
        std::move(search_length), MeshGeoToolsLib::SearchAllNodes::Yes);
    std::vector<std::size_t> const& found_ids_ply0(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply0));

    ASSERT_EQ(100u, found_ids_ply0.size());
    for (std::size_t k(0); k<found_ids_ply0.size(); k++)
        ASSERT_EQ(k, found_ids_ply0[k]);

    GeoLib::Polyline ply1(pnts);
    ply1.addPoint(2);
    ply1.addPoint(3);
    std::vector<std::size_t> const& found_ids_ply1(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply1));

    ASSERT_EQ(100u, found_ids_ply1.size());
    for (std::size_t k(0); k<found_ids_ply1.size(); k++)
        ASSERT_EQ(k, found_ids_ply1[k]);

    GeoLib::Polyline ply2(pnts);
    ply2.addPoint(4);
    ply2.addPoint(5);
    std::vector<std::size_t> const& found_ids_ply2(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply2));

    std::size_t offset((_number_of_subdivisions_per_direction+1)*_number_of_subdivisions_per_direction);
    ASSERT_EQ(100u, found_ids_ply2.size());
    for (std::size_t k(0); k<found_ids_ply2.size(); k++)
        ASSERT_EQ(offset + k, found_ids_ply2[k]);

    GeoLib::Polyline ply3(pnts);
    ply3.addPoint(6);
    ply3.addPoint(7);
    std::vector<std::size_t> const& found_ids_ply3(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply3));

    ASSERT_EQ(100u, found_ids_ply3.size());
    for (std::size_t k(0); k<found_ids_ply3.size(); k++)
        ASSERT_EQ(offset + k, found_ids_ply3[k]);

    // left border
    GeoLib::Polyline ply4(pnts);
    ply4.addPoint(0);
    ply4.addPoint(6);
    std::vector<std::size_t> const& found_ids_ply4(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply4));

    ASSERT_EQ(100u, found_ids_ply4.size());
    for (std::size_t k(0); k<found_ids_ply4.size(); k++)
        ASSERT_EQ(k*(_number_of_subdivisions_per_direction+1), found_ids_ply4[k]);

    // right border
    GeoLib::Polyline ply5(pnts);
    ply5.addPoint(1);
    ply5.addPoint(7);
    std::vector<std::size_t> const& found_ids_ply5(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply5));

    ASSERT_EQ(100u, found_ids_ply5.size());
    for (std::size_t k(0); k<found_ids_ply5.size(); k++)
        ASSERT_EQ(k*(_number_of_subdivisions_per_direction+1)+_number_of_subdivisions_per_direction, found_ids_ply5[k]);

    // diagonal
    GeoLib::Polyline ply6(pnts);
    ply6.addPoint(0);
    ply6.addPoint(5);
    std::vector<std::size_t> const& found_ids_ply6(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply6));

    ASSERT_EQ(100u, found_ids_ply6.size());
    for (std::size_t k(0); k<found_ids_ply6.size(); k++)
        ASSERT_EQ(k*(_number_of_subdivisions_per_direction+1)+k, found_ids_ply6[k]);

    std::for_each(pnts.begin(), pnts.end(), [](GeoLib::Point* pnt) { delete pnt; });
}

TEST_F(MeshLibMeshNodeSearchInSimpleQuadMesh, SurfaceSearch)
{
    ASSERT_TRUE(_quad_mesh != nullptr);
    // create geometry
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, _geometric_size, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, _geometric_size, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, 0.5*_geometric_size, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, 0.5*_geometric_size, 0.0));

    std::unique_ptr<MeshGeoToolsLib::SearchLength> search_length{
        new MeshGeoToolsLib::SearchLength};
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(*_quad_mesh,
        std::move(search_length), MeshGeoToolsLib::SearchAllNodes::Yes);

    // entire domain
    GeoLib::Polyline ply0(pnts);
    ply0.addPoint(0);
    ply0.addPoint(1);
    ply0.addPoint(2);
    ply0.addPoint(3);
    ply0.addPoint(0);
    std::unique_ptr<GeoLib::Surface> sfc0(GeoLib::Surface::createSurface(ply0));

    std::vector<std::size_t> const& found_ids_sfc0(mesh_node_searcher.getMeshNodeIDsAlongSurface(*sfc0));

    ASSERT_EQ(_quad_mesh->getNumberOfNodes(), found_ids_sfc0.size());
    for (std::size_t k(0); k<found_ids_sfc0.size(); k++)
        ASSERT_EQ(k, found_ids_sfc0[k]);

    // bottom half domain
    GeoLib::Polyline ply1(pnts);
    ply1.addPoint(0);
    ply1.addPoint(1);
    ply1.addPoint(4);
    ply1.addPoint(5);
    ply1.addPoint(0);
    std::unique_ptr<GeoLib::Surface> sfc1(GeoLib::Surface::createSurface(ply1));

    std::vector<std::size_t> const& found_ids_sfc1(mesh_node_searcher.getMeshNodeIDsAlongSurface(*sfc1));

    ASSERT_EQ(_quad_mesh->getNumberOfNodes()/2, found_ids_sfc1.size());
    for (std::size_t k(0); k<found_ids_sfc1.size(); k++)
        ASSERT_EQ(k, found_ids_sfc1[k]);

    std::for_each(pnts.begin(), pnts.end(), [](GeoLib::Point* pnt) { delete pnt; });
}

TEST_F(MeshLibMeshNodeSearchInSimpleHexMesh, SurfaceSearch)
{
    ASSERT_TRUE(_hex_mesh != nullptr);
    // create geometry
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, _geometric_size, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, _geometric_size, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, 0.0, _geometric_size));
    pnts.push_back(new GeoLib::Point(_geometric_size, 0.0, _geometric_size));
    pnts.push_back(new GeoLib::Point(_geometric_size, _geometric_size, _geometric_size));
    pnts.push_back(new GeoLib::Point(0.0, _geometric_size, _geometric_size));

    std::unique_ptr<MeshGeoToolsLib::SearchLength> search_length{
        new MeshGeoToolsLib::SearchLength};
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(*_hex_mesh,
        std::move(search_length), MeshGeoToolsLib::SearchAllNodes::Yes);

    const std::size_t n_nodes_1d = _number_of_subdivisions_per_direction + 1;
    const std::size_t n_nodes_2d = n_nodes_1d * n_nodes_1d;
    const std::size_t n_nodes_3d = n_nodes_1d * n_nodes_1d * n_nodes_1d;

    // bottom surface
    GeoLib::Polyline ply_bottom(pnts);
    ply_bottom.addPoint(0);
    ply_bottom.addPoint(1);
    ply_bottom.addPoint(2);
    ply_bottom.addPoint(3);
    ply_bottom.addPoint(0);
    std::unique_ptr<GeoLib::Surface> sfc_bottom(GeoLib::Surface::createSurface(ply_bottom));

    std::vector<std::size_t> const& found_ids_sfc_b(mesh_node_searcher.getMeshNodeIDsAlongSurface(*sfc_bottom));
    ASSERT_EQ(n_nodes_2d, found_ids_sfc_b.size());
    for (std::size_t k(0); k<found_ids_sfc_b.size(); k++)
        ASSERT_EQ(k, found_ids_sfc_b[k]);

    // top surface
    GeoLib::Polyline ply_top(pnts);
    ply_top.addPoint(4);
    ply_top.addPoint(5);
    ply_top.addPoint(6);
    ply_top.addPoint(7);
    ply_top.addPoint(4);
    std::unique_ptr<GeoLib::Surface> sfc_top(GeoLib::Surface::createSurface(ply_top));

    std::vector<std::size_t> const& found_ids_sfc_t(mesh_node_searcher.getMeshNodeIDsAlongSurface(*sfc_top));
    ASSERT_EQ(n_nodes_2d, found_ids_sfc_t.size());
    const std::size_t offset_t = n_nodes_3d - n_nodes_2d;
    for (std::size_t k(0); k<found_ids_sfc_t.size(); k++)
        ASSERT_EQ(offset_t + k, found_ids_sfc_t[k]);

    // front
    GeoLib::Polyline ply_front(pnts);
    ply_front.addPoint(0);
    ply_front.addPoint(1);
    ply_front.addPoint(5);
    ply_front.addPoint(4);
    ply_front.addPoint(0);
    std::unique_ptr<GeoLib::Surface> sfc_front(GeoLib::Surface::createSurface(ply_front));

    std::vector<std::size_t> const& found_ids_sfc_f(mesh_node_searcher.getMeshNodeIDsAlongSurface(*sfc_front));
    ASSERT_EQ(n_nodes_2d, found_ids_sfc_f.size());
    std::size_t cnt=0;
    for (std::size_t k(0); k<n_nodes_1d; k++) {
        for (std::size_t i(0); i<n_nodes_1d; i++)
            ASSERT_EQ(k*n_nodes_2d+i, found_ids_sfc_f[cnt++]);
    }

    // back
    GeoLib::Polyline ply_back(pnts);
    ply_back.addPoint(3);
    ply_back.addPoint(2);
    ply_back.addPoint(6);
    ply_back.addPoint(7);
    ply_back.addPoint(3);
    std::unique_ptr<GeoLib::Surface> sfc_back(GeoLib::Surface::createSurface(ply_back));
    std::vector<std::size_t> const& found_ids_sfc_back(mesh_node_searcher.getMeshNodeIDsAlongSurface(*sfc_back));
    ASSERT_EQ(n_nodes_2d, found_ids_sfc_back.size());
    cnt = 0;
    const std::size_t y_offset = n_nodes_1d*(n_nodes_1d-1);
    for (std::size_t k(0); k<n_nodes_1d; k++) {
        const std::size_t z_offset = n_nodes_2d*k;
        for (std::size_t i(0); i<n_nodes_1d; i++)
            ASSERT_EQ(z_offset + y_offset + i, found_ids_sfc_back[cnt++]);
    }

    // left
    GeoLib::Polyline ply_left(pnts);
    ply_left.addPoint(0);
    ply_left.addPoint(3);
    ply_left.addPoint(7);
    ply_left.addPoint(4);
    ply_left.addPoint(0);
    std::unique_ptr<GeoLib::Surface> sfc_left(GeoLib::Surface::createSurface(ply_left));
    std::vector<std::size_t> const& found_ids_sfc_left(mesh_node_searcher.getMeshNodeIDsAlongSurface(*sfc_left));
    ASSERT_EQ(n_nodes_2d, found_ids_sfc_left.size());
    cnt = 0;
    for (std::size_t k(0); k<n_nodes_1d; k++) {
        const std::size_t z_offset = n_nodes_2d*k;
        for (std::size_t j(0); j<n_nodes_1d; j++)
            ASSERT_EQ(z_offset + j*n_nodes_1d, found_ids_sfc_left[cnt++]);
    }

    // right
    GeoLib::Polyline ply_right(pnts);
    ply_right.addPoint(1);
    ply_right.addPoint(2);
    ply_right.addPoint(6);
    ply_right.addPoint(5);
    ply_right.addPoint(1);
    std::unique_ptr<GeoLib::Surface> sfc_right(GeoLib::Surface::createSurface(ply_right));
    std::vector<std::size_t> const& found_ids_sfc_right(mesh_node_searcher.getMeshNodeIDsAlongSurface(*sfc_right));
    ASSERT_EQ(n_nodes_2d, found_ids_sfc_right.size());
    cnt = 0;
    for (std::size_t k(0); k<n_nodes_1d; k++) {
        const std::size_t z_offset = n_nodes_2d*k;
        for (std::size_t j(0); j<n_nodes_1d; j++)
            ASSERT_EQ(z_offset + (j+1)*n_nodes_1d-1, found_ids_sfc_right[cnt++]);
    }


    std::for_each(pnts.begin(), pnts.end(), [](GeoLib::Point* pnt) { delete pnt; });
}

