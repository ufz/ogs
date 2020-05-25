/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "gtest/gtest.h"

#include <memory>
#include <numeric>

#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"

#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"
#include "MeshGeoToolsLib/HeuristicSearchLength.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshGenerators/QuadraticMeshGenerator.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/Node.h"

using namespace MeshLib;

class MeshLibBoundaryElementSearchInSimpleQuadMesh : public testing::Test
{
public:
    MeshLibBoundaryElementSearchInSimpleQuadMesh()
        : quad_mesh_(MeshGenerator::generateRegularQuadMesh(
              geometric_size_, number_of_subdivisions_per_direction_))
    {}

    ~MeshLibBoundaryElementSearchInSimpleQuadMesh() override
    {
        for (auto p : pnts_)
        {
            delete p;
        }
    }

protected:
    const double geometric_size_{10.0};
    const std::size_t number_of_subdivisions_per_direction_{10};
    std::unique_ptr<Mesh> quad_mesh_;
    std::vector<GeoLib::Point*> pnts_;
};

class MeshLibBoundaryElementSearchInSimpleHexMesh : public testing::Test
{
public:
    MeshLibBoundaryElementSearchInSimpleHexMesh()
        : hex_mesh_(MeshGenerator::generateRegularHexMesh(
              geometric_size_, number_of_subdivisions_per_direction_))
    {
        // Points for the surfaces. Corners of a cube.
        pnts_.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
        pnts_.push_back(new GeoLib::Point(geometric_size_, 0.0, 0.0));
        pnts_.push_back(
            new GeoLib::Point(geometric_size_, geometric_size_, 0.0));
        pnts_.push_back(new GeoLib::Point(0.0, geometric_size_, 0.0));
        pnts_.push_back(
            new GeoLib::Point(geometric_size_, 0.0, geometric_size_));
        pnts_.push_back(new GeoLib::Point(0.0, 0.0, geometric_size_));
    }

    ~MeshLibBoundaryElementSearchInSimpleHexMesh() override
    {
        for (auto p : pnts_)
        {
            delete p;
        }
    }

    void TestHexSurfacesBoundaryElementSearcher(
        Mesh const& mesh,
        std::size_t const n_nodes_2d,
        std::size_t const n_eles_2d) const;

protected:
    const double geometric_size_{10.0};
    const std::size_t number_of_subdivisions_per_direction_{10};
    std::unique_ptr<Mesh> hex_mesh_;
    std::vector<GeoLib::Point*> pnts_;
};

TEST_F(MeshLibBoundaryElementSearchInSimpleQuadMesh, PolylineSearch)
{
    ASSERT_TRUE(quad_mesh_ != nullptr);

    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(
        *quad_mesh_,
        std::make_unique<MeshGeoToolsLib::HeuristicSearchLength>(*quad_mesh_),
        MeshGeoToolsLib::SearchAllNodes::Yes);
    MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(
        *quad_mesh_, mesh_node_searcher);
    bool const multiple_nodes_allowed = false;

    const unsigned n_eles_per_dir = number_of_subdivisions_per_direction_;
    const unsigned n_nodes_per_dir = number_of_subdivisions_per_direction_ + 1;

    // points for the polylines
    pnts_.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
    pnts_.push_back(new GeoLib::Point(geometric_size_, 0.0, 0.0));
    pnts_.push_back(new GeoLib::Point(geometric_size_, geometric_size_, 0.0));
    pnts_.push_back(new GeoLib::Point(0.0, geometric_size_, 0.0));

    GeoLib::Polyline ply_bottom(pnts_);
    ply_bottom.addPoint(0);
    ply_bottom.addPoint(1);
    {
        std::vector<MeshLib::Element*> const& found_edges_ply_bottom(
            boundary_element_searcher.getBoundaryElements(
                ply_bottom, multiple_nodes_allowed));

        // check the total number of found edges
        ASSERT_EQ(n_eles_per_dir, found_edges_ply_bottom.size());
        for (unsigned i = 0; i < n_eles_per_dir; i++)
        {
            // edges found on a bottom line
            auto* edge = found_edges_ply_bottom[i];
            ASSERT_EQ(2u, edge->getNumberOfBaseNodes());
            ASSERT_EQ(i, edge->getNodeIndex(0));
            ASSERT_EQ(i + 1, edge->getNodeIndex(1));
        }
    }

    GeoLib::Polyline ply_right(pnts_);
    ply_right.addPoint(1);
    ply_right.addPoint(2);
    {
        auto const& found_edges_ply_right(
            boundary_element_searcher.getBoundaryElements(
                ply_right, multiple_nodes_allowed));

        // check the total number of found edges
        ASSERT_EQ(n_eles_per_dir, found_edges_ply_right.size());
        for (unsigned i = 0; i < n_eles_per_dir; i++)
        {
            // edges found on a right line
            auto* edge = found_edges_ply_right[i];
            ASSERT_EQ(2u, edge->getNumberOfBaseNodes());
            ASSERT_EQ((i + 1) * n_nodes_per_dir - 1, edge->getNodeIndex(0));
            ASSERT_EQ((i + 2) * n_nodes_per_dir - 1, edge->getNodeIndex(1));
        }
    }

    GeoLib::Polyline ply_top(pnts_);
    ply_top.addPoint(2);
    ply_top.addPoint(3);
    {
        auto const& found_edges_ply_top(
            boundary_element_searcher.getBoundaryElements(
                ply_top, multiple_nodes_allowed));

        // check the total number of found edges
        ASSERT_EQ(n_eles_per_dir, found_edges_ply_top.size());
        for (unsigned i = 0; i < n_eles_per_dir; i++)
        {
            // edges found on a top line
            auto* edge = found_edges_ply_top[i];
            ASSERT_EQ(2u, edge->getNumberOfBaseNodes());
            ASSERT_EQ(n_nodes_per_dir * n_nodes_per_dir - 1 - i,
                      edge->getNodeIndex(0));
            ASSERT_EQ(n_nodes_per_dir * n_nodes_per_dir - 2 - i,
                      edge->getNodeIndex(1));
        }
    }

    GeoLib::Polyline ply_left(pnts_);
    ply_left.addPoint(3);
    ply_left.addPoint(0);
    {
        auto const& found_edges_ply_left(
            boundary_element_searcher.getBoundaryElements(
                ply_left, multiple_nodes_allowed));

        // check the total number of found edges
        ASSERT_EQ(n_eles_per_dir, found_edges_ply_left.size());
        for (unsigned i = 0; i < n_eles_per_dir; i++)
        {
            // edges found on a left line
            auto* edge = found_edges_ply_left[i];
            ASSERT_EQ(2u, edge->getNumberOfBaseNodes());
            ASSERT_EQ(n_nodes_per_dir * (n_nodes_per_dir - 1 - i),
                      edge->getNodeIndex(0));
            ASSERT_EQ(n_nodes_per_dir * (n_nodes_per_dir - 2 - i),
                      edge->getNodeIndex(1));
        }
    }
}

template <typename ElementIterator>
double computeAreaOfFaceElements(ElementIterator first, ElementIterator last)
{
    return std::accumulate(first, last, 0.0, [](double v, MeshLib::Element* e) {
        return v + e->getContent();
    });
}

void MeshLibBoundaryElementSearchInSimpleHexMesh::
    TestHexSurfacesBoundaryElementSearcher(Mesh const& mesh,
                                           std::size_t const n_nodes_2d,
                                           std::size_t const n_eles_2d) const
{
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(
        mesh,
        std::make_unique<MeshGeoToolsLib::SearchLength>(),
        MeshGeoToolsLib::SearchAllNodes::Yes);
    MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(
        mesh, mesh_node_searcher);

    // perform search on the bottom surface
    GeoLib::Surface sfc_bottom(pnts_);
    sfc_bottom.addTriangle(0, 1, 2);
    sfc_bottom.addTriangle(0, 2, 3);

    bool const multiple_nodes_allowed = false;
    std::vector<MeshLib::Element*> const& found_faces_sfc_b(
        boundary_element_searcher.getBoundaryElements(sfc_bottom,
                                                      multiple_nodes_allowed));
    ASSERT_EQ(n_eles_2d, found_faces_sfc_b.size());
    ASSERT_EQ(geometric_size_ * geometric_size_,
              computeAreaOfFaceElements(found_faces_sfc_b.begin(),
                                        found_faces_sfc_b.end()));
    auto connected_nodes_b = MeshLib::getUniqueNodes(found_faces_sfc_b);
    ASSERT_EQ(n_nodes_2d, connected_nodes_b.size());
    for (auto node : connected_nodes_b)
    {
        ASSERT_EQ(0.0, (*node)[2]);  // check z coordinates
    }

    // perform search on the front surface
    GeoLib::Surface sfc_front(pnts_);
    sfc_front.addTriangle(0, 1, 4);
    sfc_front.addTriangle(0, 4, 5);

    std::vector<MeshLib::Element*> const& found_faces_sfc_f(
        boundary_element_searcher.getBoundaryElements(sfc_front,
                                                      multiple_nodes_allowed));
    ASSERT_EQ(n_eles_2d, found_faces_sfc_f.size());
    ASSERT_EQ(geometric_size_ * geometric_size_,
              computeAreaOfFaceElements(found_faces_sfc_f.begin(),
                                        found_faces_sfc_f.end()));
    auto connected_nodes_f = MeshLib::getUniqueNodes(found_faces_sfc_f);
    ASSERT_EQ(n_nodes_2d, connected_nodes_f.size());
    for (auto node : connected_nodes_f)
    {
        ASSERT_EQ(0.0, (*node)[1]);  // check y coordinates
    }
}

TEST_F(MeshLibBoundaryElementSearchInSimpleHexMesh, SurfaceSearch)
{
    ASSERT_TRUE(hex_mesh_ != nullptr);
    const std::size_t& s = number_of_subdivisions_per_direction_;
    const std::size_t n_nodes_2d = (s + 1) * (s + 1);
    const std::size_t n_eles_2d = s * s;

    TestHexSurfacesBoundaryElementSearcher(*hex_mesh_, n_nodes_2d, n_eles_2d);
}

// This is identical to the above
// MeshLibBoundaryElementSearchInSimpleHexMesh.SurfaceSearch test but
// creates a quadratic mesh from the original hex mesh.
TEST_F(MeshLibBoundaryElementSearchInSimpleHexMesh, QuadElementsSurfaceSearch)
{
    ASSERT_TRUE(hex_mesh_ != nullptr);
    auto mesh = MeshLib::createQuadraticOrderMesh(*hex_mesh_);

    const std::size_t& s = number_of_subdivisions_per_direction_;
    const std::size_t n_nodes_2d = (s + 1) * (3 * s + 1);
    const std::size_t n_eles_2d = s * s;

    TestHexSurfacesBoundaryElementSearcher(*mesh, n_nodes_2d, n_eles_2d);
}
