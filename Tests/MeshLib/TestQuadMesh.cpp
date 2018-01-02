/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include <iterator>
#include "gtest/gtest.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Quad.h"

class MeshLibQuadMesh : public ::testing::Test
{
    public:
    MeshLibQuadMesh()
        : mesh(nullptr)
    {
        mesh = MeshLib::MeshGenerator::generateRegularQuadMesh(1.0, n_elements);
    }

    ~MeshLibQuadMesh() override { delete mesh; }
    static std::size_t const n_elements = 4;
    static std::size_t const n_nodes = n_elements + 1;
    static std::size_t const elements_stride = n_elements - 1;
    static std::size_t const nodes_stride = n_nodes - 1;

    MeshLib::Mesh const* mesh;

    public:
    // Helper functions to access elements and nodes.

    MeshLib::Element* getElement(std::size_t const i, std::size_t const j)
    {
        return mesh->getElements()[i * n_elements + j];
    }

    MeshLib::Node* getNode(std::size_t const i, std::size_t const j)
    {
        return mesh->getNodes()[i * n_nodes + j];
    }

    using Indices = std::list<std::size_t>;
    Indices getNeighbor(std::size_t const i) const
    {
        std::list<std::size_t> result;
        switch (i) {
        case 0:
            result.push_back(i + 1);
            break;
        case elements_stride:
            result.push_back(i - 1);
            break;
        default:
            result.push_back(i - 1);
            result.push_back(i + 1);
            break;
        }
        return result;
    }

    template <typename F>
    void
    testCornerElements(F const&& f)
    {
        f(getElement(0, 0), 0, 0);                            // BL
        f(getElement(0, elements_stride), 0, elements_stride);              // BR
        f(getElement(elements_stride, 0), elements_stride, 0);              // BL
        f(getElement(elements_stride, elements_stride), elements_stride, elements_stride);// BR
    }

    template <typename F>
    void
    testBoundaryElements(F const&& f)
    {
        // left
        for (std::size_t j = 1; j < elements_stride; ++j)
            f(getElement(0, j), 0, j);
        // right
        for (std::size_t j = 1; j < elements_stride; ++j)
            f(getElement(elements_stride, j), elements_stride, j);
        // bottom
        for (std::size_t i = 1; i < elements_stride; ++i)
            f(getElement(i, 0), i, 0);
        // top
        for (std::size_t i = 1; i < elements_stride; ++i)
            f(getElement(i, elements_stride), i, elements_stride);
    }

    template <typename F>
    void
    testInsideElements(F const&& f)
    {
        for (std::size_t i = 1; i < elements_stride; ++i)
            for (std::size_t j = 1; j < elements_stride; ++j)
                f(getElement(i, j), i, j);
    }

    template <typename F>
    void
    testAllElements(F const&& f)
    {
        for (std::size_t i = 0; i < elements_stride; ++i)
            for (std::size_t j = 0; j < elements_stride; ++j)
                f(getElement(i, j), i, j);
    }

    template <typename F>
    void
    testCornerNodes(F const&& f)
    {
        f(getNode(0, 0), 0, 0);
        f(getNode(0, nodes_stride), 0, nodes_stride);
        f(getNode(nodes_stride, 0), nodes_stride, 0);
        f(getNode(nodes_stride, nodes_stride),
                nodes_stride, nodes_stride);
    }

    template <typename F>
    void
    testBoundaryNodes(F const&& f)
    {
        // left
        for (std::size_t j = 1; j < nodes_stride; ++j)
            f(getNode(0, j), 0, j);
        // right
        for (std::size_t j = 1; j < nodes_stride; ++j)
            f(getNode(nodes_stride, j), nodes_stride, j);
        // bottom
        for (std::size_t i = 1; i < nodes_stride; ++i)
            f(getNode(i, 0), i, 0);
        // top
        for (std::size_t i = 1; i < nodes_stride; ++i)
            f(getNode(i, nodes_stride), i, nodes_stride);
    }

    template <typename F>
    void
    testInsideNodes(F const&& f)
    {
        for (std::size_t i = 1; i < nodes_stride; ++i)
            for (std::size_t j = 1; j < nodes_stride; ++j)
                f(getNode(i, j), i, j);
    }
};

TEST_F(MeshLibQuadMesh, Construction)
{
    ASSERT_TRUE(mesh != nullptr);

    // There are n_elements^2 elements in the mesh.
    ASSERT_EQ(n_elements * n_elements, mesh->getNumberOfElements());

    // There are n_nodes^2 nodes in the mesh.
    ASSERT_EQ(n_nodes*n_nodes, mesh->getNumberOfNodes());

    // All elements have maximum four neighbors.
    testAllElements([](MeshLib::Element const* const e, ...)
        {
            ASSERT_EQ(4u, e->getNumberOfNeighbors());
        });
}

TEST_F(MeshLibQuadMesh, ElementNeighbors)
{
    auto count_neighbors = [](MeshLib::Element const* const e)
        {
            unsigned count = 0;
            for (int i = 0; i < 4; i++)
                if (e->getNeighbor(i) != nullptr)
                    count++;
            return count;
        };

    auto getNeighborIndices = [this](std::size_t const i, std::size_t const j)
    {
        return std::make_pair(getNeighbor(i), getNeighbor(j));
    };

    auto testNeighbors = [this](
        MeshLib::Element const* const e,
        std::size_t const i,
        std::size_t const j,
        std::pair<Indices, Indices> const& neighbors)
    {
        for (auto i_neighbor : neighbors.first)
            ASSERT_TRUE(e->hasNeighbor(getElement(i_neighbor, j)));

        for (auto j_neighbor : neighbors.second)
            ASSERT_TRUE(e->hasNeighbor(getElement(i, j_neighbor)));
    };

    // Two neighbors for corner elements.
    testCornerElements([&](MeshLib::Element const* const e,
                       std::size_t const i, std::size_t const j)
        {
            EXPECT_EQ(2u, count_neighbors(e));

            std::pair<Indices, Indices> const ij_neighbors = getNeighborIndices(i, j);
            // Test the test
            EXPECT_EQ(1u, ij_neighbors.first.size());
            EXPECT_EQ(1u, ij_neighbors.second.size());
            ASSERT_TRUE(e->isBoundaryElement());

            testNeighbors(e, i, j, ij_neighbors);
        });

    // Three neighbors for boundary elements.
    testBoundaryElements([&](MeshLib::Element const* const e,
                         std::size_t const i, std::size_t const j)
        {
            EXPECT_EQ(3u, count_neighbors(e));

            std::pair<Indices, Indices> const ij_neighbors = getNeighborIndices(i, j);
            // Test the test
            EXPECT_EQ(3u, ij_neighbors.first.size() + ij_neighbors.second.size());
            ASSERT_TRUE(e->isBoundaryElement());

            testNeighbors(e, i, j, ij_neighbors);
        });

    // Four neighbors inside mesh.
    testInsideElements([&](MeshLib::Element const* const e,
                       std::size_t const i, std::size_t const j)
        {
            EXPECT_EQ(4u, count_neighbors(e));

            std::pair<Indices, Indices> const ij_neighbors = getNeighborIndices(i, j);
            // Test the test
            EXPECT_EQ(2u, ij_neighbors.first.size());
            EXPECT_EQ(2u, ij_neighbors.second.size());
            ASSERT_FALSE(e->isBoundaryElement());

            testNeighbors(e, i, j, ij_neighbors);
        });
}

TEST_F(MeshLibQuadMesh, ElementToNodeConnectivity)
{
    // An element (i,j) consists of four nodes (i,j), (i+1,j), (i+1, j+1),
    // and (i, j+1).
    testAllElements([this](
        MeshLib::Element const* const e,
        std::size_t const i,
        std::size_t const j)
        {
            EXPECT_EQ(4u, e->getNumberOfBaseNodes());
            EXPECT_EQ(getNode(i,   j),   e->getNode(0));
            EXPECT_EQ(getNode(i,   j+1), e->getNode(1));
            EXPECT_EQ(getNode(i+1, j+1), e->getNode(2));
            EXPECT_EQ(getNode(i+1, j),   e->getNode(3));
        });
}

// A node is connected to four elements inside the mesh, two on the boundary,
// and one in the corner.
TEST_F(MeshLibQuadMesh, NodeToElementConnectivity)
{
    testCornerNodes(
        [this](MeshLib::Node const* const node,
            std::size_t i,
            std::size_t j)
        {
            EXPECT_EQ(1u, node->getNumberOfElements());

            if (i == nodes_stride)
                i--;
            if (j == nodes_stride)
                j--;

            EXPECT_EQ(getElement(i, j), node->getElement(0));
        });

    testBoundaryNodes([this](
        MeshLib::Node const* const node,
        std::size_t i,
        std::size_t j)
        {
            EXPECT_EQ(2u, node->getNumberOfElements());

            if (i == 0)
            {
                EXPECT_EQ(getElement(i, j-1), node->getElement(0));
                EXPECT_EQ(getElement(i, j), node->getElement(1));
            }
            if (i == nodes_stride)
            {
                EXPECT_EQ(getElement(elements_stride, j-1), node->getElement(0));
                EXPECT_EQ(getElement(elements_stride, j), node->getElement(1));
            }
            if (j == 0)
            {
                EXPECT_EQ(getElement(i-1, j), node->getElement(0));
                EXPECT_EQ(getElement(i, j), node->getElement(1));
            }
            if (j == nodes_stride)
            {
                j--;
                EXPECT_EQ(getElement(i-1, j), node->getElement(0));
                EXPECT_EQ(getElement(i, j), node->getElement(1));
            }
        });

    testInsideNodes([this](
        MeshLib::Node const* const node,
        std::size_t const i,
        std::size_t const j)
        {
            EXPECT_EQ(4u, node->getNumberOfElements());

            EXPECT_EQ(getElement(i-1, j-1), node->getElement(0));
            EXPECT_EQ(getElement(i-1, j  ), node->getElement(1));
            EXPECT_EQ(getElement(i,   j-1), node->getElement(2));
            EXPECT_EQ(getElement(i,   j  ), node->getElement(3));
        });
}
