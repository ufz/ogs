/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include <random>
#include <variant>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/ConvertToLinearMesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshGenerators/QuadraticMeshGenerator.h"
#include "MeshLib/Node.h"

using namespace MeshLib;
using namespace std::string_literals;

class ConvertToLinearMesh : public ::testing::Test
{
public:
    ConvertToLinearMesh()
    {
        linear_mesh.reset(
            MeshGenerator::generateRegularHexMesh(1.0, mesh_size));
        quadratic_mesh =
            createQuadraticOrderMesh(*linear_mesh, false /* add centre node*/);
    }

    static constexpr int mesh_size = 5;

    std::unique_ptr<Mesh> linear_mesh;
    std::unique_ptr<Mesh> quadratic_mesh;
};

// Returns the (inverse) permutation vector of b_nodes to a_nodes, or an error
// message.
std::variant<std::vector<std::size_t>, std::string> compareNodes(
    std::vector<Node*> const& a_nodes, std::vector<Node*> const& b_nodes)
{
    // For each 'a' mesh node find a corresponding 'b' mesh node, not checking
    // the ids, but store the id's in a map.
    std::vector<std::size_t> b_node_id(a_nodes.size());

    std::size_t const nnodes = a_nodes.size();
    for (std::size_t i = 0; i < nnodes; ++i)
    {
        Node* a_node = a_nodes[i];
        auto const b_it =
            std::find_if(begin(b_nodes), end(b_nodes), [&](Node* const b_node) {
                return *a_node ==
                       *b_node;  // coordinate comparison up to epsilon
            });
        if (b_it == end(b_nodes))
        {
            return "Could not find a matching node for node " +
                   std::to_string(a_node->getID()) + ".";
        }

        b_node_id[i] = distance(begin(b_nodes), b_it);
    }

    return b_node_id;
}

// Two elements are equal if their corresponding node coordinates are equal up
// to a (shift) permutation.
// Returns a string in case of error containing a descriptive message.
std::optional<std::string> equal(Element const& a, Element const& b)
{
    if (typeid(a) != typeid(b))  // => (a.nodes.size == b.nodes.size)
    {
        return "Elements have different types."s;
    }

    auto const a_nodes = a.getNodes();
    auto const b_nodes = b.getNodes();
    auto const nnodes = a.getNumberOfNodes();

    std::size_t first_b_node = 0;
    while (first_b_node < nnodes)
    {
        if (*a_nodes[0] == *b_nodes[first_b_node])
        {
            break;
        }
        ++first_b_node;
    }
    if (first_b_node == nnodes)
    {
        return "Did not find first node of the element " +
               std::to_string(a.getID()) + " in the nodes of the element " +
               std::to_string(b.getID()) + ".";
    }

    for (std::size_t i = 1; i < nnodes; ++i)
    {
        if (!(*a_nodes[i] == *b_nodes[(first_b_node + i) % nnodes]))
        {
            return "node " + std::to_string(i) + " is not equal to the node " +
                   std::to_string((first_b_node + i) % nnodes) + ".";
        }
    }

    return {};  // All nodes are equal.
}

// Returns a vector filled with numbers [0, size) randomly permuted.
std::vector<std::size_t> generateRandomPermutationVector(std::size_t const size)
{
    std::random_device rd;
    std::mt19937 random_generator(rd());

    std::vector<std::size_t> permutation(size);
    std::iota(begin(permutation), end(permutation), 0);
    std::shuffle(begin(permutation), end(permutation), random_generator);
    return permutation;
}

// Test the inverse permutation for being the inverse: p^-1(p(.)) = Id.
// In case of failure an error message is returned.
std::optional<std::string> inversePermutationIdentityTest(
    std::vector<std::size_t> const& permutation,
    std::vector<std::size_t> const& inverse_permutation)
{
    if (permutation.size() != inverse_permutation.size())
    {
        return "Inverse permutation size differs from the permutation size: " +
               std::to_string(inverse_permutation.size()) + " vs. " +
               std::to_string(permutation.size()) + ".";
    }

    for (std::size_t i = 0; i < permutation.size(); ++i)
    {
        auto const p = permutation[i];
        if (inverse_permutation[p] == i)
        {
            continue;
        }
        return "Inverse permutation element " + std::to_string(p) +
               " does not point to " + std::to_string(i) +
               "-th element but to " + std::to_string(inverse_permutation[p]) +
               ".";
    }
    return {};
}

// Create new nodes in permuted order with new ids.
std::vector<Node*> permuteNodes(std::vector<std::size_t> const& permutation,
                                std::vector<Node*> const& nodes)
{
    std::vector<Node*> permuted_nodes;
    for (std::size_t i = 0; i < permutation.size(); ++i)
    {
        auto const node_p_i = nodes[permutation[i]];
        permuted_nodes.push_back(new Node{node_p_i->getCoords(), i});
    }
    return permuted_nodes;
}

// Tests the nodes of two meshes and returns a string in case of error with
// message, otherwise a permutation vector as in compareNodes().
std::variant<std::vector<std::size_t>, std::string> compareNodeVectors(
    std::vector<Node*> const& a, std::vector<Node*> const& b)
{
    if (a.size() != b.size())
    {
        return "There are different number of nodes: " +
               std::to_string(a.size()) + " and " + std::to_string(b.size()) +
               ".";
    }

    auto const compare_nodes_result = compareNodes(a, b);
    if (std::holds_alternative<std::string>(compare_nodes_result))
    {
        return "Node comparison failed: " +
               std::get<std::string>(compare_nodes_result);
    }

    return std::get<std::vector<std::size_t>>(compare_nodes_result);
}

TEST_F(ConvertToLinearMesh, GeneratedHexMeshRandomizedNodes)
{
    ASSERT_TRUE(linear_mesh != nullptr);
    ASSERT_TRUE(quadratic_mesh != nullptr);

    auto const permutation =
        generateRandomPermutationVector(quadratic_mesh->getNumberOfNodes());

    // Create new nodes in permuted order with new ids.
    auto const permuted_quadratic_nodes =
        permuteNodes(permutation, quadratic_mesh->getNodes());

    //
    // Test the permuted nodes and compute the inverse permutation map.
    //
    auto const compare_nodes_result = compareNodeVectors(
        quadratic_mesh->getNodes(), permuted_quadratic_nodes);
    ASSERT_FALSE(std::holds_alternative<std::string>(compare_nodes_result))
        << "Quadratic mesh nodes comparison with permuted quadratic nodes "
           "failed.\n"
        << std::get<std::string>(compare_nodes_result);

    auto const& inverse_permutation =
        std::get<std::vector<std::size_t>>(compare_nodes_result);
    {
        auto const result =
            inversePermutationIdentityTest(permutation, inverse_permutation);
        ASSERT_TRUE(result == boost::none)
            << "Quadratic mesh nodes permutation test failed: " << *result;
    }

    auto clone_element_using_permuted_nodes = [&](Element* const e) {
        Node** nodes = new Node*[e->getNumberOfNodes()];
        for (std::size_t i = 0; i < e->getNumberOfNodes(); ++i)
        {
            auto const q_i = e->getNode(i)->getID();
            auto const r = inverse_permutation[q_i];
            nodes[i] = permuted_quadratic_nodes[r];
        }
        return e->clone(nodes, e->getID());
    };

    // Create new elements in the same order, but using new nodes.
    std::vector<Element*> new_quadratic_elements;
    std::transform(begin(quadratic_mesh->getElements()),
                   end(quadratic_mesh->getElements()),
                   back_inserter(new_quadratic_elements),
                   clone_element_using_permuted_nodes);

    Mesh permuted_nodes_quadratic_mesh{"permuted_quadratic_mesh",
                                       permuted_quadratic_nodes,
                                       new_quadratic_elements};
    //
    // Test the quadratic meshes first.
    //
    {
        ASSERT_EQ(quadratic_mesh->getNumberOfElements(),
                  permuted_nodes_quadratic_mesh.getNumberOfElements());

        auto const nelements = quadratic_mesh->getNumberOfElements();
        for (std::size_t i = 0; i < nelements; ++i)
        {
            auto const& element_a = *quadratic_mesh->getElement(i);
            auto const& element_b =
                *permuted_nodes_quadratic_mesh.getElement(i);
            auto const elements_are_equal = equal(element_a, element_b);
            ASSERT_TRUE(elements_are_equal == boost::none)
                << *elements_are_equal << " For the element " << i << ".";
        }
    }

    {
        auto const converted_mesh = convertToLinearMesh(
            permuted_nodes_quadratic_mesh, "back_to_linear_mesh");

        //
        // Two meshes are equal if all of their nodes have exactly same
        // coordinates.
        //
        auto const compare_nodes_result = compareNodeVectors(
            linear_mesh->getNodes(), converted_mesh->getNodes());
        ASSERT_FALSE(std::holds_alternative<std::string>(compare_nodes_result))
            << "Linear mesh nodes comparison with converted linear mesh nodes "
               "failed.\n"
            << std::get<std::string>(compare_nodes_result);
        // TODO (naumov) check inverse permutation, but it requires the reduced
        // to base nodes forward permutation vector first.

        // Two meshes are equal if their elements share same nodes.
        ASSERT_EQ(linear_mesh->getNumberOfElements(),
                  converted_mesh->getNumberOfElements());

        auto const linear_mesh_nelements = linear_mesh->getNumberOfElements();
        for (std::size_t i = 0; i < linear_mesh_nelements; ++i)
        {
            auto const& element_a = *linear_mesh->getElement(i);
            auto const& element_b = *converted_mesh->getElement(i);
            auto const elements_are_equal = equal(element_a, element_b);
            ASSERT_TRUE(elements_are_equal == boost::none)
                << *elements_are_equal << " For the element " << i << ".";
        }
    }
}

TEST_F(ConvertToLinearMesh, GeneratedHexMeshBackToLinear)
{
    ASSERT_TRUE(linear_mesh != nullptr);
    ASSERT_TRUE(quadratic_mesh != nullptr);

    auto const converted_mesh =
        convertToLinearMesh(*quadratic_mesh, "back_to_linear_mesh");

    //
    // Two meshes are equal if all of their nodes have exactly same coordinates.
    //
    auto const compare_nodes_result =
        compareNodeVectors(linear_mesh->getNodes(), converted_mesh->getNodes());
    ASSERT_FALSE(std::holds_alternative<std::string>(compare_nodes_result))
        << "Linear mesh nodes comparison with converted linear mesh nodes "
           "failed.\n"
        << std::get<std::string>(compare_nodes_result);

    // Two meshes are equal if their elements share same nodes.
    ASSERT_EQ(linear_mesh->getNumberOfElements(),
              converted_mesh->getNumberOfElements());

    auto const linear_mesh_nelements = linear_mesh->getNumberOfElements();
    for (std::size_t i = 0; i < linear_mesh_nelements; ++i)
    {
        auto const& linear_mesh_element = *linear_mesh->getElement(i);
        auto const& converted_mesh_element = *converted_mesh->getElement(i);
        auto const elements_are_equal =
            equal(linear_mesh_element, converted_mesh_element);
        ASSERT_TRUE(elements_are_equal == boost::none)
            << *elements_are_equal << " For the element " << i << ".";
    }
}
