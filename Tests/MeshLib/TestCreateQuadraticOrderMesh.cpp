/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 29, 2022, 2:05 PM
 */

#include <gtest/gtest.h>

#include <fstream>
#include <memory>

#include "InfoLib/TestInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/HexRule20.h"
#include "MeshLib/Elements/LineRule3.h"
#include "MeshLib/Elements/PrismRule15.h"
#include "MeshLib/Elements/PyramidRule13.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/QuadRule8.h"
#include "MeshLib/Elements/TetRule10.h"
#include "MeshLib/Elements/TriRule6.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/MeshGenerators/QuadraticMeshGenerator.h"
#include "MeshLib/Node.h"

const unsigned* getEdgeNodeLocalIDs(MeshLib::Element const& element,
                                    unsigned const edge_id)
{
    switch (element.getGeomType())
    {
        case MeshLib::MeshElemType::LINE:
            return MeshLib::LineRule3::edge_nodes[edge_id];
        case MeshLib::MeshElemType::TRIANGLE:
            return MeshLib::TriRule6::edge_nodes[edge_id];
        case MeshLib::MeshElemType::QUAD:
            return MeshLib::QuadRule8::edge_nodes[edge_id];
        case MeshLib::MeshElemType::TETRAHEDRON:
            return MeshLib::TetRule10::edge_nodes[edge_id];
        case MeshLib::MeshElemType::HEXAHEDRON:
            return MeshLib::HexRule20::edge_nodes[edge_id];
        case MeshLib::MeshElemType::PRISM:
            return MeshLib::PrismRule15::edge_nodes[edge_id];
        case MeshLib::MeshElemType::PYRAMID:
            return MeshLib::PyramidRule13::edge_nodes[edge_id];
        default:
            OGS_FATAL(
                "Element with the geometry type of {:s} is not supported "
                "in createQuadraticOrderMesh.",
                MeshLib::MeshElemType2String(element.getGeomType()));
    }
}

void runCreateQuadraticOrderMeshTest(MeshLib::Mesh const& linear_mesh,
                                     MeshLib::Mesh const& quadratic_mesh)
{
    for (std::size_t e_id = 0; e_id < quadratic_mesh.getNumberOfElements();
         e_id++)
    {
        auto const linear_element = linear_mesh.getElement(e_id);
        auto const quadratic_element = quadratic_mesh.getElement(e_id);

        // Compare vertex nodes
        for (std::size_t i = 0; i < linear_element->getNumberOfNodes(); i++)
        {
            auto const x = linear_element->getNode(i)->asEigenVector3d();
            auto const x_q = quadratic_element->getNode(i)->asEigenVector3d();

            ASSERT_EQ((x - x_q).norm(), 0.0);
        }

        // Compare edge nodes:
        auto const nodes = linear_element->getNodes();
        auto const nodes_q = quadratic_element->getNodes();
        for (unsigned i = 0; i < quadratic_element->getNumberOfEdges(); i++)
        {
            auto const edge_node_ids =
                getEdgeNodeLocalIDs(*quadratic_element, i);
            auto const x_vertex_a = nodes[edge_node_ids[0]]->asEigenVector3d();
            auto const x_vertex_b = nodes[edge_node_ids[1]]->asEigenVector3d();
            auto const x_edge = nodes_q[edge_node_ids[2]]->asEigenVector3d();
            ASSERT_EQ((0.5 * (x_vertex_a + x_vertex_b) - x_edge).norm(), 0.0);
        }
    }
}

TEST(MeshLib, createQuadraticOrderMesh)
{
    // Use mesh linear_mesh.vtu that contains 7 elements with types of: line,
    // triangle, quadrilateral, tetrahedron, hexahedron, prism and pyramid:
    std::string const file_name =
        TestInfoLib::TestInfo::data_path + "/Utils/GMSH2OGS/linear_mesh.vtu";
    auto mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshLib::IO::readMeshFromFile(file_name));
    if (!mesh)
    {
        OGS_FATAL("Could not read mesh from '{:s}' file. No mesh created.",
                  file_name);
    }

    auto const quadratic_mesh =
        MeshLib::createQuadraticOrderMesh(*mesh, false /* add centre node*/);

    runCreateQuadraticOrderMeshTest(*mesh, *quadratic_mesh);
}

// Special case: quad9
TEST(MeshLib, createQuadraticOrderMesh_Quad9)
{
    // Note the memory allocated for nodes and elements are released in the
    // destructor of Mesh.
    std::array<MeshLib::Node*, 4> nodes{
        // nodes on base
        new MeshLib::Node(0.0, 0.0, 0.0), new MeshLib::Node(1.0, 0.0, 0.0),
        new MeshLib::Node(1.0, 0.0, 0.0), new MeshLib::Node(0.0, 1.0, 0.0)};

    std::vector<MeshLib::Element*> elements{new MeshLib::Quad(nodes)};
    MeshLib::Mesh mesh(
        "test_mesh", std::vector(begin(nodes), end(nodes)), elements);

    auto const quadratic_mesh =
        MeshLib::createQuadraticOrderMesh(mesh, true /* add centre node*/);

    runCreateQuadraticOrderMeshTest(mesh, *quadratic_mesh);

    // Test centre node:
    auto const node_c = quadratic_mesh->getElement(0)->getNode(8);
    auto x_centre_diff = node_c->asEigenVector3d();

    // Loop over four vertex nodes:
    for (auto const& node : nodes)
    {
        auto const x = node->asEigenVector3d();
        x_centre_diff -= 0.25 * x;
    }

    ASSERT_EQ(x_centre_diff.norm(), 0.0);
}
