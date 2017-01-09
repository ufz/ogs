/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include <memory>
#include <numeric>

#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/HeuristicSearchLength.h"
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"

using namespace MeshLib;

class MeshLibBoundaryElementSearchInSimpleQuadMesh : public testing::Test
{
public:
    MeshLibBoundaryElementSearchInSimpleQuadMesh() :
        _geometric_size(10.0), _number_of_subdivisions_per_direction(10),
        _quad_mesh(MeshGenerator::generateRegularQuadMesh(_geometric_size, _number_of_subdivisions_per_direction))
    {}

protected:
    const double _geometric_size;
    const std::size_t _number_of_subdivisions_per_direction;
    std::unique_ptr<Mesh> _quad_mesh;
};

class MeshLibBoundaryElementSearchInSimpleHexMesh : public testing::Test
{
public:
    MeshLibBoundaryElementSearchInSimpleHexMesh() :
        _geometric_size(10.0), _number_of_subdivisions_per_direction(10),
        _hex_mesh(MeshGenerator::generateRegularHexMesh(_geometric_size, _number_of_subdivisions_per_direction))
    {}

    ~MeshLibBoundaryElementSearchInSimpleHexMesh()
    {
        delete _hex_mesh;
    }

protected:
    const double _geometric_size;
    const std::size_t _number_of_subdivisions_per_direction;
    Mesh* _hex_mesh;
};

TEST_F(MeshLibBoundaryElementSearchInSimpleQuadMesh, PolylineSearch)
{
    ASSERT_TRUE(_quad_mesh != nullptr);
    const unsigned n_eles_per_dir = _number_of_subdivisions_per_direction;
    const unsigned n_nodes_per_dir = _number_of_subdivisions_per_direction + 1;

    // create a polyline (closed)
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, _geometric_size, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, _geometric_size, 0.0));
    GeoLib::Polyline ply0(pnts);
    ply0.addPoint(0);
    ply0.addPoint(1);
    ply0.addPoint(2);
    ply0.addPoint(3);
    ply0.addPoint(0);

    // perform search on the polyline
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(*_quad_mesh,
        MeshGeoToolsLib::HeuristicSearchLength(*_quad_mesh));
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

    for (auto p : pnts)
        delete p;
}

TEST_F(MeshLibBoundaryElementSearchInSimpleHexMesh, SurfaceSearch)
{
    ASSERT_TRUE(_hex_mesh != nullptr);
    const std::size_t n_nodes_1d = _number_of_subdivisions_per_direction + 1;
    const std::size_t n_nodes_2d = n_nodes_1d * n_nodes_1d;
    const std::size_t n_eles_2d = (n_nodes_1d-1) * (n_nodes_1d-1);

    // create bottom and front surfaces of a cubic
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, 0.0, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, _geometric_size, 0.0));
    pnts.push_back(new GeoLib::Point(0.0, _geometric_size, 0.0));
    pnts.push_back(new GeoLib::Point(_geometric_size, 0.0, _geometric_size));
    pnts.push_back(new GeoLib::Point(0.0, 0.0, _geometric_size));

    GeoLib::Polyline ply_bottom(pnts);
    ply_bottom.addPoint(0);
    ply_bottom.addPoint(1);
    ply_bottom.addPoint(2);
    ply_bottom.addPoint(3);
    ply_bottom.addPoint(0);
    std::unique_ptr<GeoLib::Surface> sfc_bottom(GeoLib::Surface::createSurface(ply_bottom));

    GeoLib::Polyline ply_front(pnts);
    ply_front.addPoint(0);
    ply_front.addPoint(1);
    ply_front.addPoint(4);
    ply_front.addPoint(5);
    ply_front.addPoint(0);
    std::unique_ptr<GeoLib::Surface> sfc_front(GeoLib::Surface::createSurface(ply_front));

    // perform search on the bottom surface
    MeshGeoToolsLib::MeshNodeSearcher mesh_node_searcher(*_hex_mesh,
        MeshGeoToolsLib::HeuristicSearchLength(*_hex_mesh));
    MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(*_hex_mesh, mesh_node_searcher);
    std::vector<MeshLib::Element*> const& found_faces_sfc_b(boundary_element_searcher.getBoundaryElements(*sfc_bottom));
    ASSERT_EQ(n_eles_2d, found_faces_sfc_b.size());
    double sum_area_b = std::accumulate(found_faces_sfc_b.begin(), found_faces_sfc_b.end(), 0.0,
                [](double v, MeshLib::Element*e){return v+e->getContent();});
    ASSERT_EQ(_geometric_size*_geometric_size, sum_area_b);
    auto connected_nodes_b = MeshLib::getUniqueNodes(found_faces_sfc_b);
    ASSERT_EQ(n_nodes_2d, connected_nodes_b.size());
    for (auto node : connected_nodes_b)
        ASSERT_EQ(0.0, (*node)[2]); // check z coordinates

    // perform search on the front surface
    std::vector<MeshLib::Element*> const& found_faces_sfc_f(boundary_element_searcher.getBoundaryElements(*sfc_front));
    ASSERT_EQ(n_eles_2d, found_faces_sfc_f.size());
    double sum_area_f = std::accumulate(found_faces_sfc_f.begin(), found_faces_sfc_f.end(), 0.0,
                [](double v, MeshLib::Element*e){return v+e->getContent();});
    ASSERT_EQ(_geometric_size*_geometric_size, sum_area_f);
    auto connected_nodes_f = MeshLib::getUniqueNodes(found_faces_sfc_f);
    ASSERT_EQ(n_nodes_2d, connected_nodes_f.size());
    for (auto node : connected_nodes_f)
        ASSERT_EQ(0.0, (*node)[1]); // check y coordinates

    for (auto p : pnts)
        delete p;
}

