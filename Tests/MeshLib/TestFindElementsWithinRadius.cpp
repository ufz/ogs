/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>
#include <autocheck/autocheck.hpp>

#include "BaseLib/makeVectorUnique.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Node.h"
#include "MeshLib/findElementsWithinRadius.h"

#include "Tests/AutoCheckTools.h"

using namespace MeshLib;
namespace ac = autocheck;

struct MeshLibFindElementWithinRadius : public ::testing::Test
{
    ac::gtest_reporter gtest_reporter;
};

// For zero radius only the starting element is alway returned.
TEST_F(MeshLibFindElementWithinRadius, ZeroRadius)
{
    auto mesh = MeshGenerator::generateRegularQuadMesh(10., 10);

    auto same_element_returned = [&mesh](std::size_t& element_id) -> bool {
        auto result =
            findElementsWithinRadius(*mesh->getElement(element_id), 0.);
        return (result.size() == 1) && (result[0] == element_id);
    };

    ac::check<std::size_t>(same_element_returned, 100,
                           ac::make_arbitrary<std::size_t>(), gtest_reporter);
}

std::vector<std::size_t> findElementIdsConnectedToNode(Node const& node)
{
    std::vector<std::size_t> connected_elements;
    auto const& elements_of_node = node.getElements();
    std::transform(std::begin(elements_of_node),
                   std::end(elements_of_node),
                   std::back_inserter(connected_elements),
                   [](auto const& e) { return e->getID(); });
    return connected_elements;
}

// Find all elements connected through all nodes of the element
std::vector<std::size_t> findByNodeConnectedElements(Element const& element)
{
    std::vector<std::size_t> connected_elements;
    for (unsigned n = 0; n < element.getNumberOfBaseNodes(); ++n)
    {
        auto const& node = *element.getNode(n);
        auto const& elements = findElementIdsConnectedToNode(node);
        std::copy(std::begin(elements), std::end(elements),
                  std::back_inserter(connected_elements));
    }

    BaseLib::makeVectorUnique(connected_elements);
    return connected_elements;
}

// Find nodes in radius of a given node.
std::vector<Node const*> findNodesWithinRadius(Mesh const& mesh,
                                               Node const& node,
                                               double const radius)
{
    std::vector<Node const*> nodes;

    for (std::size_t i = 0; i < mesh.getNumberOfNodes(); ++i)
    {
        auto const* n = mesh.getNode(i);
        auto const distance_squared = MathLib::sqrDist(node, *n);
        if (distance_squared > radius * radius)
            continue;

        nodes.push_back(n);
    }
    return nodes;
}

// Find nodes in radius of all nodes of a given element.
std::vector<Node const*> findNodesWithinRadius(Mesh const& mesh,
                                               Element const& element,
                                               double const radius)
{
    std::vector<Node const*> nodes_in_radius;

    for (unsigned n = 0; n < element.getNumberOfBaseNodes(); ++n)
    {
        auto const& node = *element.getNode(n);
        auto const& nodes = findNodesWithinRadius(mesh, node, radius);
        std::copy(std::begin(nodes), std::end(nodes),
                  std::back_inserter(nodes_in_radius));
    }

    BaseLib::makeVectorUnique(nodes_in_radius);
    return nodes_in_radius;
}

// Brute force search for checking the actual algorithm.
std::vector<std::size_t> bruteForceFindElementIdsWithinRadius(
    Mesh const& mesh, Element const& element, double const radius)
{
    std::vector<std::size_t> connected_elements;
    auto const nodes_in_radius = findNodesWithinRadius(mesh, element, radius);
    for (auto n : nodes_in_radius)
    {
        auto const& elements = findElementIdsConnectedToNode(*n);
        std::copy(std::begin(elements), std::end(elements),
                  std::back_inserter(connected_elements));
    }

    BaseLib::makeVectorUnique(connected_elements);
    return connected_elements;
}

// For a small radius only the element and its neighbors (through all nodes) are
// expected.
TEST_F(MeshLibFindElementWithinRadius, VerySmallRadius)
{
    auto mesh = MeshGenerator::generateRegularQuadMesh(10., 10);

    auto neighboring_elements_returned =
        [&mesh](std::size_t& element_id) -> bool {
        auto const& element = *mesh->getElement(element_id);
        auto result = findElementsWithinRadius(element, 1e-5);

        std::sort(std::begin(result), std::end(result));
        auto const expected_elements = findByNodeConnectedElements(element);

        return result.size() == expected_elements.size() &&
               std::includes(std::begin(result), std::end(result),
                             std::begin(expected_elements),
                             std::end(expected_elements));
    };

    ac::check<std::size_t>(neighboring_elements_returned, 100,
                           ac::make_arbitrary<std::size_t>(), gtest_reporter);
}

// For radii large enough to cover all of the mesh all of the elements are
// expected to be found.
TEST_F(MeshLibFindElementWithinRadius, VeryLargeRadius)
{
    auto mesh = *MeshGenerator::generateRegularQuadMesh(10., 10);

    auto all_elements_returned = [&mesh](std::size_t& element_id) -> bool {
        auto result =
            findElementsWithinRadius(*mesh.getElement(element_id), 1e5);
        BaseLib::makeVectorUnique(result);

        return result.size() == mesh.getNumberOfElements();
    };

    ac::check<std::size_t>(all_elements_returned, 100,
                           ac::make_arbitrary<std::size_t>(), gtest_reporter);
}

// Random test (element_id and radius > 0); comparison with brute-force search
// algorithm.
TEST_F(MeshLibFindElementWithinRadius, RandomPositiveRadius2d)
{
    auto mesh = *MeshGenerator::generateRegularQuadMesh(10., 10);

    auto compare_to_brute_force_search = [&mesh](std::size_t& element_id,
                                                 double const radius) -> bool {
        auto const& element = *mesh.getElement(element_id);

        auto result = findElementsWithinRadius(element, radius);
        std::sort(std::begin(result), std::end(result));

        auto const expected_elements =
            bruteForceFindElementIdsWithinRadius(mesh, element, radius);

        return result.size() == expected_elements.size() &&
               std::includes(std::begin(result), std::end(result),
                             std::begin(expected_elements),
                             std::end(expected_elements));
    };

    ac::check<std::size_t, double>(
        compare_to_brute_force_search, 100,
        ac::make_arbitrary(ac::generator<std::size_t>(),
                           ac::map(&ac::absoluteValue, ac::generator<double>()))
            .discard_if(
                [](std::size_t, double const r) { return (r < 1e-16); }),
        gtest_reporter);
}

// Random test (element_id and radius > 0); comparison with brute-force search
// algorithm.
TEST_F(MeshLibFindElementWithinRadius, RandomPositiveRadius3d)
{
    auto mesh = *MeshGenerator::generateRegularHexMesh(10., 10);

    auto compare_to_brute_force_search = [&mesh](std::size_t& element_id,
                                                 double const radius) -> bool {
        auto const& element = *mesh.getElement(element_id);

        auto result = findElementsWithinRadius(element, radius);
        std::sort(std::begin(result), std::end(result));

        auto const expected_elements =
            bruteForceFindElementIdsWithinRadius(mesh, element, radius);

        return result.size() == expected_elements.size() &&
               std::includes(std::begin(result), std::end(result),
                             std::begin(expected_elements),
                             std::end(expected_elements));
    };

    ac::check<std::size_t, double>(
        compare_to_brute_force_search, 100,
        ac::make_arbitrary(ac::generator<std::size_t>(),
                           ac::map(&ac::absoluteValue, ac::generator<double>()))
            .discard_if(
                [](std::size_t, double const r) { return (r < 1e-16); }),
        gtest_reporter);
}
