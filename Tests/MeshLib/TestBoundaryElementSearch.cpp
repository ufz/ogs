/**
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
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
    MeshLibBoundaryElementSearchInSimpleQuadMesh() :
        _geometric_size(10.0), _number_of_subdivisions_per_direction(10),
        _quad_mesh(MeshGenerator::generateRegularQuadMesh(_geometric_size, _number_of_subdivisions_per_direction))
    {}

    ~MeshLibBoundaryElementSearchInSimpleQuadMesh() override
    {
        for (auto p : _pnts)
            delete p;
    }

protected:
    const double _geometric_size;
    const std::size_t _number_of_subdivisions_per_direction;
    std::unique_ptr<Mesh> _quad_mesh;
    std::vector<GeoLib::Point*> _pnts;
};

class MeshLibBoundaryElementSearchInSimpleHexMesh : public testing::Test
{
public:
    MeshLibBoundaryElementSearchInSimpleHexMesh() :
        _geometric_size(10.0), _number_of_subdivisions_per_direction(10),
        _hex_mesh(MeshGenerator::generateRegularHexMesh(_geometric_size, _number_of_subdivisions_per_direction))
    {
        // Points for the surfaces. Corners of a cube.
        _pnts.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
        _pnts.push_back(new GeoLib::Point(_geometric_size, 0.0, 0.0));
        _pnts.push_back(
            new GeoLib::Point(_geometric_size, _geometric_size, 0.0));
        _pnts.push_back(new GeoLib::Point(0.0, _geometric_size, 0.0));
        _pnts.push_back(
            new GeoLib::Point(_geometric_size, 0.0, _geometric_size));
        _pnts.push_back(new GeoLib::Point(0.0, 0.0, _geometric_size));
    }

    ~MeshLibBoundaryElementSearchInSimpleHexMesh() override
    {
        for (auto p : _pnts)
            delete p;
    }

    void TestHexSurfacesBoundaryElementSearcher(
        Mesh const& mesh,
        std::size_t const n_nodes_2d,
        std::size_t const n_eles_2d) const;

protected:
    const double _geometric_size;
    const std::size_t _number_of_subdivisions_per_direction;
    std::unique_ptr<Mesh> _hex_mesh;
    std::vector<GeoLib::Point*> _pnts;
};

TEST_F(MeshLibBoundaryElementSearchInSimpleQuadMesh, PolylineSearch)
{
    ASSERT_TRUE(_quad_mesh != nullptr);
    const unsigned n_eles_per_dir = _number_of_subdivisions_per_direction;
    const unsigned n_nodes_per_dir = _number_of_subdivisions_per_direction + 1;

    // create a polyline (closed)
    _pnts.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
    _pnts.push_back(new GeoLib::Point(_geometric_size, 0.0, 0.0));
    _pnts.push_back(new GeoLib::Point(_geometric_size, _geometric_size, 0.0));
    _pnts.push_back(new GeoLib::Point(0.0, _geometric_size, 0.0));
    GeoLib::Polyline ply0(_pnts);
    ply0.addPoint(0);
    ply0.addPoint(1);
    ply0.addPoint(2);
    ply0.addPoint(3);
    ply0.addPoint(0);

    // perform search on the polyline
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(
        *_quad_mesh,
        std::make_unique<MeshGeoToolsLib::HeuristicSearchLength>(*_quad_mesh),
        MeshGeoToolsLib::SearchAllNodes::Yes);
    MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(*_quad_mesh, mesh_node_searcher);
    std::vector<MeshLib::Element*> const& found_edges_ply0(boundary_element_searcher.getBoundaryElements(ply0));

    // check the total number of found edges
    ASSERT_EQ(n_eles_per_dir*4u, found_edges_ply0.size());
    // check node IDs of edges
    for (unsigned i=0; i<n_eles_per_dir; i++) {
        // edge found on a bottom line
        auto* edge0 = found_edges_ply0[i];
        ASSERT_EQ(2u, edge0->getNumberOfBaseNodes());
        ASSERT_EQ(i, edge0->getNodeIndex(0));
        ASSERT_EQ(i+1, edge0->getNodeIndex(1));
        // edge found on a right line
        auto* edge1 = found_edges_ply0[i+n_eles_per_dir];
        ASSERT_EQ(2u, edge1->getNumberOfBaseNodes());
        ASSERT_EQ(n_nodes_per_dir*i+n_nodes_per_dir-1, edge1->getNodeIndex(0));
        ASSERT_EQ(n_nodes_per_dir*(i+1)+n_nodes_per_dir-1, edge1->getNodeIndex(1));
        // edge found on a top line
        auto* edge2 = found_edges_ply0[i+n_eles_per_dir*2];
        ASSERT_EQ(2u, edge2->getNumberOfBaseNodes());
        ASSERT_EQ(n_nodes_per_dir*n_nodes_per_dir-1-i, edge2->getNodeIndex(0));
        ASSERT_EQ(n_nodes_per_dir*n_nodes_per_dir-1-(i+1), edge2->getNodeIndex(1));
        // edge found on a left line
        auto* edge3 = found_edges_ply0[i+n_eles_per_dir*3];
        ASSERT_EQ(2u, edge3->getNumberOfBaseNodes());
        ASSERT_EQ(n_nodes_per_dir*(n_nodes_per_dir-1)-n_nodes_per_dir*i, edge3->getNodeIndex(0));
        ASSERT_EQ(n_nodes_per_dir*(n_nodes_per_dir-1)-n_nodes_per_dir*(i+1), edge3->getNodeIndex(1));
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
    GeoLib::Surface sfc_bottom(_pnts);
    sfc_bottom.addTriangle(0, 1, 2);
    sfc_bottom.addTriangle(0, 2, 3);

    std::vector<MeshLib::Element*> const& found_faces_sfc_b(
        boundary_element_searcher.getBoundaryElements(sfc_bottom));
    ASSERT_EQ(n_eles_2d, found_faces_sfc_b.size());
    ASSERT_EQ(_geometric_size * _geometric_size,
              computeAreaOfFaceElements(found_faces_sfc_b.begin(),
                                        found_faces_sfc_b.end()));
    auto connected_nodes_b = MeshLib::getUniqueNodes(found_faces_sfc_b);
    ASSERT_EQ(n_nodes_2d, connected_nodes_b.size());
    for (auto node : connected_nodes_b)
        ASSERT_EQ(0.0, (*node)[2]);  // check z coordinates

    // perform search on the front surface
    GeoLib::Surface sfc_front(_pnts);
    sfc_front.addTriangle(0, 1, 4);
    sfc_front.addTriangle(0, 4, 5);

    std::vector<MeshLib::Element*> const& found_faces_sfc_f(
        boundary_element_searcher.getBoundaryElements(sfc_front));
    ASSERT_EQ(n_eles_2d, found_faces_sfc_f.size());
    ASSERT_EQ(_geometric_size * _geometric_size,
              computeAreaOfFaceElements(found_faces_sfc_f.begin(),
                                        found_faces_sfc_f.end()));
    auto connected_nodes_f = MeshLib::getUniqueNodes(found_faces_sfc_f);
    ASSERT_EQ(n_nodes_2d, connected_nodes_f.size());
    for (auto node : connected_nodes_f)
        ASSERT_EQ(0.0, (*node)[1]);  // check y coordinates
}

TEST_F(MeshLibBoundaryElementSearchInSimpleHexMesh, SurfaceSearch)
{
    ASSERT_TRUE(_hex_mesh != nullptr);
    const std::size_t& s = _number_of_subdivisions_per_direction;
    const std::size_t n_nodes_2d = (s + 1) * (s + 1);
    const std::size_t n_eles_2d = s * s;

    TestHexSurfacesBoundaryElementSearcher(*_hex_mesh, n_nodes_2d, n_eles_2d);
}

// This is identical to the above
// MeshLibBoundaryElementSearchInSimpleHexMesh.SurfaceSearch test but
// creates a quadratic mesh from the original hex mesh.
TEST_F(MeshLibBoundaryElementSearchInSimpleHexMesh, QuadElementsSurfaceSearch)
{
    ASSERT_TRUE(_hex_mesh != nullptr);
    auto mesh = MeshLib::createQuadraticOrderMesh(*_hex_mesh);

    const std::size_t& s = _number_of_subdivisions_per_direction;
    const std::size_t n_nodes_2d = (s + 1) * (3 * s + 1);
    const std::size_t n_eles_2d = s * s;

    TestHexSurfacesBoundaryElementSearcher(*mesh, n_nodes_2d, n_eles_2d);
}
